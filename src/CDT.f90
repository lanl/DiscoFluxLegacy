module CDT_routines
private
public :: CDT_Combine_State, CDT_State_Init, CDT_Update_State
contains
SUBROUTINE CDT_State_Init(cell_sizes, cdt_state)
  USE GlobalParams
  USE CDT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nel), INTENT(IN)  :: cell_sizes
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_CDT_STATE),   INTENT(INOUT) :: cdt_state
!-----------------------------------------------------------------------
!~   call cdt_state%allocatememory()
  cdt_state%cell_sizes(:) = cell_sizes

  RETURN
END SUBROUTINE CDT_State_Init

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE CDT_Update_State(cdt_state_old, ddc_state_new, dmb_state_new,  &
           all_cdt_props, dt, flux_bc, cdt_state_new)
  USE GlobalParams
  USE CDT
  USE DDC
  USE DMB
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_CDT_State),         INTENT(IN)  :: cdt_state_old
  TYPE(Type_DDC_State),         INTENT(IN)  :: ddc_state_new
  TYPE(Type_DMB_State),         INTENT(IN)  :: dmb_state_new
  TYPE(Type_CDT_PROPS),         INTENT(IN)  :: all_cdt_props(Nregion)
  REAL(KIND=8),                       INTENT(IN)  :: dt
  REAL(KIND=8), DIMENSION(Nregion+1), INTENT(IN)  :: flux_bc
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_CDT_STATE),         INTENT(INOUT) :: cdt_state_new
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: varrho_src
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: velocity
  REAL(KIND=8), DIMENSION(Nel)       :: cell_sizes
  TYPE(Type_CDT_PROPS)               :: cdt_props
  INTEGER :: islip, iel, ich
!-----------------------------------------------------------------------
  allocate(varrho_src(Nel,Nslip,Nchar), velocity(Nel,Nslip,Nchar))
!~   call cdt_state_new%allocatememory()
  cdt_props = all_cdt_props(1)
  
  cell_sizes(:) = cdt_state_old%cell_sizes(:)

  call CDT_State_Init(cell_sizes,cdt_state_new)

  call CDT_Src_Calc(cdt_state_old,ddc_state_new,dmb_state_new,all_cdt_props, &
       varrho_src, velocity)

  if (cdt_props%no_flux) velocity = 0.d0

  if (cdt_props%update_method == 'FVM') then

    do ich=1,Nchar
      !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(islip)
      do islip=1,Nslip
        call FV_Update(cdt_state_old%varrho_field(:,islip,ich,1), velocity(:,islip,ich),        &
                     flux_bc, varrho_src(:,islip,ich), cell_sizes, dt, cdt_props%scheme,  &
                     cdt_props%limiter, cdt_state_new%varrho_field(:,islip,ich,1))

        call FV_Update(cdt_state_old%varrho_field(:,islip,ich,2), -velocity(:,islip,ich),       &
                     flux_bc, varrho_src(:,islip,ich), cell_sizes, dt, cdt_props%scheme,  &
                     cdt_props%limiter, cdt_state_new%varrho_field(:,islip,ich,2))
      end do ! islip
      !$OMP END PARALLEL DO
    end do  !ich

  elseif (cdt_props%update_method == 'FEM') then

    do ich=1,Nchar
      !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(islip)
      do islip=1,Nslip
        call FEM_Update(cdt_state_old%varrho_node(:,islip,ich,1), velocity(:,islip,ich),         &
                      flux_bc, varrho_src(:,islip,ich), cell_sizes, dt, cdt_props%scheme,  &
                      cdt_state_new%varrho_node(:,islip,ich,1))
          
        call FEM_Update(cdt_state_old%varrho_node(:,islip,ich,2), -velocity(:,islip,ich),        &
                      flux_bc, varrho_src(:,islip,ich), cell_sizes, dt, cdt_props%scheme,  &
                      cdt_state_new%varrho_node(:,islip,ich,2))
      end do ! islip
      !$OMP END PARALLEL DO
    end do  !ich

  else

    write(*,100) ">>>>>>>>>> Unknown numerical method '",cdt_props%update_method,"' <<<<<<<<< <"
    call Fatal('Check setup.')

  end if

100 format(3a)
  
! Dislocation densities evaluated at element centers are stored in varrho 
  if (cdt_props%update_method == 'FEM') then

    do iel=1,Nel
      cdt_state_new%varrho_field(iel,:,:,:) = 0.5d0*(cdt_state_new%varrho_node(iel,:,:,:) + cdt_state_new%varrho_node(iel+1,:,:,:))
    end do  !iel

  end if

  RETURN
END SUBROUTINE CDT_Update_State

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE FV_Update(concentration, velocity, flux_bc, src, Lc, dt, nsc, &
           limiter, c_new)
  USE GlobalParams
  USE CDT
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: concentration
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: velocity
  REAL(KIND=8), DIMENSION(Nregion+1), INTENT(IN)  :: flux_bc
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: src
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: Lc
  REAL(KIND=8),                       INTENT(IN)  :: dt
  CHARACTER(32),                INTENT(IN)  :: nsc
  CHARACTER(32),                INTENT(IN)  :: limiter
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nel),       INTENT(OUT) :: c_new
!-----------------------------------------------------------------------
!  Locals:
  INTEGER                  :: iEdge, iNode, iRegion
  REAL(KIND=8), DIMENSION(Nnode) :: flux
  REAL(KIND=8), DIMENSION(Nel+2) :: cfl_cell, slopes
  REAL(KIND=8), DIMENSION(Nel)   :: c_dot
  REAL(KIND=8), DIMENSION(2)     :: c,v
  REAL(KIND=8)                   :: dX
!-----------------------------------------------------------------------
!~ !  Function interfaces:
!~   INTERFACE
!~     FUNCTION numFlux(c,v,dX,dt,nsc)
!~       REAL(KIND=8)                                    :: numFlux
!~       REAL(KIND=8), DIMENSION(2),         INTENT(IN)  :: c,v
!~       REAL(KIND=8),                       INTENT(IN)  :: dX,dt
!~       CHARACTER(32),                INTENT(IN)  :: nsc
!~     END FUNCTION numFlux
!~   END INTERFACE
!~ !-----------------------------------------------------------------------

! CFL condition
  cfl_cell(1)       = velocity(1)*dt/Lc(1)
  cfl_cell(2:Nel+1) = velocity(1:Nel)*dt/Lc(1:Nel)
  cfl_cell(Nel+2)   = velocity(Nel)*dt/Lc(Nel)

  if ( (ANY(cfl_cell > 0.5d0)) .or. (ANY(cfl_cell < -0.5d0)) ) then
    write(*,*) ">>>>>>>>>> FV_UPDATE: CFL condition violated <<<<<<<<<<"
    call Fatal('Bad CFL.')
  end if

! Reconstruction
  call reconstruct(concentration,Lc,limiter, slopes)

  dX   = Lc(1)
  v(1) = velocity(1)
  v(2) = velocity(1)
  c(1) = concentration(1) + slopes(1) * 0.5d0*dX
  c(2) = concentration(1) - slopes(2) * 0.5d0*dX
  flux(1) = numFlux(c,v,dX,dt,nsc) ! Numerical flux
  
  do iEdge=2,Nel
    dX   = 0.5d0*(Lc(iEdge)+Lc(iEdge-1))
    v(1) = velocity(iEdge-1)
    v(2) = velocity(iEdge)
    c(1) = concentration(iEdge-1) + slopes(iEdge  ) * 0.5d0*Lc(iEdge-1)
    c(2) = concentration(iEdge)   - slopes(iEdge+1) * 0.5d0*Lc(iEdge)

!   Numerical flux
    flux(iEdge) = numFlux(c,v,dX,dt,nsc)
  end do  !iEdge

  dX   = Lc(Nel)
  v(1) = velocity(Nel)
  v(2) = velocity(Nel)
  c(1) = concentration(Nel) + slopes(Nel+1) * 0.5d0*dX
  c(2) = concentration(Nel) - slopes(Nel+2) * 0.5d0*dX
  flux(Nel+1) = numFlux(c,v,dX,dt,nsc) ! Numerical flux
 
! Flux BCs
  if(abs(rnul-flux_bc(1))>rzero) flux(1)     = flux_bc(1)
  do iregion=1,Nregion
    iNode = region_mask(iregion,2) + 1
    if(abs(rnul-flux_bc(iregion+1))>rzero) flux(iNode) = flux_bc(iregion+1)
  end do
  
! Evolution
  c_dot(1:Nel) = src(1:Nel) - (flux(2:Nel+1)-flux(1:Nel)) / Lc(1:Nel)
  c_new(:) = concentration(:) + c_dot(:) * dt
  
  RETURN
END SUBROUTINE FV_Update

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

SUBROUTINE reconstruct(concentration,Lc,limiter, slopes)
  USE GlobalParams
  USE CDT
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: concentration
!~   REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: velocity
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: Lc
  CHARACTER(32),                INTENT(IN)  :: limiter
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nel+2),     INTENT(OUT) :: slopes
!-----------------------------------------------------------------------
!  Locals:
  INTEGER                                   :: iel
  REAL(KIND=8), DIMENSION(Nel+4)                  :: Q,h
  REAL(KIND=8)                                    :: sed,awk
!-----------------------------------------------------------------------
!~ !  Function interfaces:
!~   INTERFACE
!~     FUNCTION minmod(a,b)
!~       REAL(KIND=8)                                :: minmod
!~       REAL(KIND=8),                   INTENT(IN)  :: a,b
!~     END FUNCTION minmod
!~     FUNCTION maxmod(a,b)
!~       REAL(KIND=8)                                :: maxmod
!~       REAL(KIND=8),                   INTENT(IN)  :: a,b
!~     END FUNCTION maxmod
!~   END INTERFACE
!~ !-----------------------------------------------------------------------

  slopes  = 0.d0

  if ( limiter == 'zero' ) then        ! 1st order (Godunov)
    return ! do nothing
  else
    Q(1:2)         = concentration(1)
    Q(3:Nel+2)     = concentration(1:Nel)
    Q(Nel+3:Nel+4) = concentration(Nel)
    h(1:2)         = Lc(1)
    h(3:Nel+2)     = Lc(1:Nel)
    h(Nel+3:Nel+4) = Lc(Nel)
  end if

  if ( limiter == 'Fromm' ) then   ! centered slope (Fromm)
    slopes(1:Nel+2) = ( Q(3:Nel+4)-Q(1:Nel+2) ) / ( h(2:Nel+3) + 0.5d0*(h(1:Nel+2)+h(3:Nel+4)) )
  elseif ( limiter == 'BW' ) then      ! upwind slope (Beam-Warming)
    slopes(1:Nel+2) = ( Q(2:Nel+3)-Q(1:Nel+2) ) / ( 0.5d0*(h(1:Nel+2)+h(2:Nel+3)) )
  elseif ( limiter == 'LxW' ) then     ! downwind slope (Lax-Wendroff)
    slopes(1:Nel+2) = ( Q(3:Nel+4)-Q(2:Nel+3) ) / ( 0.5d0*(h(2:Nel+3)+h(3:Nel+4)) )
  elseif ( limiter == 'MMod' ) then    ! Minmod limiter
    do iel = 1, Nel+2
        slopes(iel) = minmod( ( Q(iel+1)-Q(iel  ) ) / ( 0.5d0*(h(iel+1)+h(iel  )) ), &
                              ( Q(iel+2)-Q(iel+1) ) / ( 0.5d0*(h(iel+2)+h(iel+1)) ) )
    end do ! iel
  elseif ( limiter == 'SBee' ) then    ! Superbee limiter
    do iel = 1, Nel+2
        sed         = minmod(      ( Q(iel+2)-Q(iel+1) ) / ( 0.5d0*(h(iel+2)+h(iel+1)) ) , &
                              2.d0*( Q(iel+1)-Q(iel  ) ) / ( 0.5d0*(h(iel+1)+h(iel  )) ) )
        awk         = minmod( 2.d0*( Q(iel+2)-Q(iel+1) ) / ( 0.5d0*(h(iel+2)+h(iel+1)) ) , &
                                   ( Q(iel+1)-Q(iel  ) ) / ( 0.5d0*(h(iel+1)+h(iel  )) ) )
        slopes(iel) = maxmod( sed, awk )
    end do ! iel
  elseif ( limiter == 'MC' ) then      ! Monotonized central-diff. limiter
    do iel = 1, Nel+2
        sed         = minmod( 2.d0*( Q(iel+2)-Q(iel+1) ) / ( 0.5d0*(h(iel+2)+h(iel+1)) ) , &
                              2.d0*( Q(iel+1)-Q(iel  ) ) / ( 0.5d0*(h(iel+1)+h(iel  )) ) )
        slopes(iel) = minmod( sed, ( Q(iel+2)-Q(iel  ) ) / ( 0.5d0*(h(iel+2)+h(iel  ))+h(iel+1) ) )
    end do ! iel
  else
    write(*,100) ">>>>>>>>>> FV_UPDATE: Unknown limiter '",limiter,"' <<<<<<<<<<"
    call Fatal('Bad FV scheme.')
  end if


100 format(3a)

END SUBROUTINE reconstruct

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

FUNCTION minmod(a,b)
  REAL(KIND=8)                                :: minmod
  REAL(KIND=8),                   INTENT(IN)  :: a,b

  if ( a*b <= 0.d0 ) then
    minmod = 0.d0
  else
    if ( abs(a) < abs(b) ) then
      minmod = a
    else
      minmod = b
    end if
  end if

END FUNCTION minmod

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

FUNCTION maxmod(a,b)
  REAL(KIND=8)                                :: maxmod
  REAL(KIND=8),                   INTENT(IN)  :: a,b

  if ( a*b <= 0.d0 ) then
    maxmod = 0.d0
  else
    if ( abs(a) > abs(b) ) then
      maxmod = a
    else
      maxmod = b
    end if
  end if

END FUNCTION maxmod

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

FUNCTION numFlux(c,v,dX,dt,nsc)
  use utilities, only: Fatal
  IMPLICIT NONE
  REAL(KIND=8)                                    :: numFlux
!-----------------------------------------------------------------------
!  Inputs: 
  REAL(KIND=8), DIMENSION(2),         INTENT(IN)  :: c,v
  REAL(KIND=8),                       INTENT(IN)  :: dX,dt
  CHARACTER(32),                INTENT(IN)  :: nsc

!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(2)                      :: f
  REAL(KIND=8)                                    :: vlf

  f(1) = v(1)*c(1)
  f(2) = v(2)*c(2)
  
! Central difference flux (Central)
  if (nsc == "Central") then
    numFlux = 0.5d0*( f(1)+f(2) )

! Godunov flux (Godunov)
  elseif (nsc == "Godunov") then
    if ( c(1) < c(2) ) then
      numFlux = min(f(1),f(2))
    else
      numFlux = max(f(1),f(2))
    end if

! Lax-Friedrichs flux (LxF)
  elseif (nsc == "LxF") then

    vlf = dX/dt           ! Add an excessive amount of artificial diffusion
!   vlf = dX/dt*0.35      ! or add just a fraction (say 35%) of that
    numFlux = 0.5d0*( f(1)+f(2) ) - 0.5d0*vlf*( c(2)-c(1) )

! Local Lax-Friedrichs flux (LLF)
  elseif (nsc == "LLF") then
    vlf = max( abs(v(1)), abs(v(2)) )
    numFlux = 0.5d0*( f(1)+f(2) ) - 0.5d0*vlf*( c(2)-c(1) )

! Otherwise, complain bitterly
  else
    write(*,100) ">>>>>>>>>> FV_UPDATE: Unknown numerical flux '",nsc,"' <<<<<<<<<<"
    call Fatal('Bad FV scheme.')
  end if

100 format(3a)

END FUNCTION numFlux

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE FEM_Update(concentration, velocity, flux_bc, src, Lc, dt, nsc, &
           c_new)
  USE GlobalParams
  USE CDT
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nnode),     INTENT(IN)  :: concentration
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: velocity
  REAL(KIND=8), DIMENSION(Nregion+1), INTENT(IN)  :: flux_bc
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: src
  REAL(KIND=8), DIMENSION(Nel),       INTENT(IN)  :: Lc
  REAL(KIND=8),                       INTENT(IN)  :: dt
  CHARACTER(32),                INTENT(IN)  :: nsc
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nnode),     INTENT(OUT) :: c_new
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(Nel)   :: cfl_elem
  REAL(KIND=8), DIMENSION(Nnode) :: flux, source, M_lumped, rhs, c_dot
  REAL(KIND=8), DIMENSION(Nnode) :: flux_supg, source_supg, M_supg_lumped
!~   REAL(KIND=8), DIMENSION(Nnode,Nnode) :: M_consis, M_supg
  ! replace M_consis, M_supg with sparse versions (don't waste memory/performance), i.e.
  ! M_consis_sprs(e,1)=M_consis(e,e), M_consis_sprs(e,2)=M_consis(e+1,e), M_consis_sprs(e,3)=M_consis(e,e+1)
  REAL(KIND=8), DIMENSION(Nnode,3) :: M_consis_sprs, M_supg_sprs
!~   REAL(KIND=8), DIMENSION(2,2)   :: matrix1
  REAL(KIND=8), DIMENSION(Nel)   :: conc_ip, f_ip
  REAL(KIND=8), DIMENSION(Nnode) :: amid
  REAL(KIND=8), DIMENSION(Nel)   :: alow, aup
  INTEGER :: e, INFO
!-----------------------------------------------------------------------

! Check CFL condition
  cfl_elem(:)  = velocity(:)*dt/Lc(:)
  if ( (ANY(cfl_elem > 0.5d0)) .or. (ANY(cfl_elem < -0.5d0)) ) then
    write(*,*) ">>>>>>>>>> FEM_UPDATE: CFL condition violated <<<<<<<<<<"
    call Fatal('Bad CFL.')
  end if

  flux         = 0.d0
  source       = 0.d0
  M_lumped     = 0.d0
!~   M_consis     = 0.d0
  M_consis_sprs = 0.d0
!~   matrix1(1,:) = (/ 1.d0/3.d0, 1.d0/6.d0 /)
!~   matrix1(2,:) = (/ 1.d0/6.d0, 1.d0/3.d0 /)

! Concentration and flux at element centers
  conc_ip(1:Nel)  = 0.5d0*(concentration(2:Nel+1)+concentration(1:Nel))

  f_ip(:) = velocity(:) * conc_ip(:)

! Loop over elements
!~   !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(e)
  do e=1,Nel
!~     f_ip(e) = velocity(e) * conc_ip(e)
    flux(e)         = flux(e)   - f_ip(e)
    flux(e+1)       = flux(e+1) + f_ip(e)
    source(e:e+1)   = source(e:e+1) + 0.5d0*Lc(e)*src(e)
    M_lumped(e:e+1) = M_lumped(e:e+1) + 0.5d0*Lc(e)
!~     M_consis(e:e+1,e:e+1) = M_consis(e:e+1,e:e+1) + Lc(e)*matrix1(1:2,1:2)
    M_consis_sprs(e,2:3)   = M_consis_sprs(e,2:3) + Lc(e)/6.d0
    M_consis_sprs(e:e+1,1) = M_consis_sprs(e:e+1,1) + Lc(e)/3.d0
  end do  !e
!~   !$OMP END PARALLEL DO

! Standard Galerkin method (Bubnov)
  if (nsc == 'Bubnov') then

    continue               ! We don't need anything extra

! Upwind Petrov Galerkin method (SUPG)
  elseif (nsc == 'SUPG') then

  ! Compute SUPG arrays
!~     call SUPG_Arrays(velocity, concentration, src, Lc, &
!~                      M_supg, M_supg_lumped, flux_supg, source_supg)
    call SUPG_Arrays(velocity, concentration, src, Lc, &
                     M_supg_sprs, M_supg_lumped, flux_supg, source_supg)

  ! Add SUPG correction to FEM arrays
!~     M_consis = M_consis + M_supg
    M_consis_sprs = M_consis_sprs + M_supg_sprs
    M_lumped = M_lumped + M_supg_lumped
    flux     = flux     - flux_supg
    source   = source   + source_supg

! Otherwise, complain
  else
    write(*,100) ">>>>>>>>>> FEM_UPDATE: Unknown numerical scheme '",nsc,"' <<<<<<<<<<"
    call Fatal('Bad FEM scheme.')
  end if
100 format(3a)

! Flux BCs (convention: +ive flux BC indicates disloc. flow from left to right)
  if(abs(rnul-flux_bc(1)) < rzero) then
    flux(1) = 0.d0
  else
    flux(1) = flux(1) + flux_bc(1)
  end if
  if(abs(rnul-flux_bc(Nregion+1)) < rzero) then
    flux(Nnode) = 0.d0
  else
    flux(Nnode) = flux(Nnode) - flux_bc(Nregion+1)
  end if

  rhs(:)   = flux(:) + source(:)

!~   !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(e)
  do e=1,Nel
!~     aup(e)  = M_consis(e  ,e+1)
!~     amid(e) = M_consis(e  ,e  )
!~     alow(e) = M_consis(e+1,e  )
    aup(e)  = M_consis_sprs(e,3)
    amid(e) = M_consis_sprs(e,1)
    alow(e) = M_consis_sprs(e,2)
  end do
!~   !$OMP END PARALLEL DO
!~   amid(Nel+1) = M_consis(Nel+1,Nel+1)
  amid(Nel+1) = M_consis_sprs(Nel+1,1)
  call DGTSV(Nnode,1,alow,amid,aup,rhs,Nnode,INFO)
  if (INFO/=0) call Fatal('FEM_Update: Problem w/ linear solver.')
  c_dot(:) = rhs(:)

! alow(:) = 1.d0/6.d0 * Lc(:)
! aup(:) = alow(:)
! amid(1) = 1.d0/3.d0 * Lc(1)
! amid(2:Nel) = 1.d0/3.d0 * (Lc(2:Nel) + Lc(1:Nel-1))
! amid(Nnode) = 1.d0/3.d0 * Lc(Nel)
! call DGTSV(Nnode,1,alow,amid,aup,rhs,Nnode,INFO)
! c_dot(:) = rhs(:)
! c_dot(:) = rhs(:) / M_lumped(:)   !  Appears more diffusive (?)

  c_new(:) = concentration(:) + c_dot(:) * dt

  RETURN
END SUBROUTINE FEM_Update

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

SUBROUTINE SUPG_Arrays(v, c, src, Lc, &
                       M_supg, M_supg_lumped, flux_supg, source_supg)

  USE GlobalParams
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nel),         INTENT(IN)  :: v
  REAL(KIND=8), DIMENSION(Nnode),       INTENT(IN)  :: c
  REAL(KIND=8), DIMENSION(Nel),         INTENT(IN)  :: src
  REAL(KIND=8), DIMENSION(Nel),         INTENT(IN)  :: Lc
!-----------------------------------------------------------------------
!  Outputs:
!~   REAL(KIND=8), DIMENSION(Nnode,Nnode), INTENT(OUT) :: M_supg
  ! replace M_supg with sparse versions (don't waste memory/performance)
  REAL(KIND=8), DIMENSION(Nnode,3),     INTENT(OUT) :: M_supg
  REAL(KIND=8), DIMENSION(Nnode),       INTENT(OUT) :: M_supg_lumped, flux_supg, source_supg
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: e
  REAL(KIND=8)  :: kappa,alpha,xi,tauv
!-----------------------------------------------------------------------
  
  M_supg        = 0.d0
  M_supg_lumped = 0.d0
  flux_supg     = 0.d0
  source_supg   = 0.d0

!~   !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(e,kappa,alpha,xi,tauv)
  do e = 1, Nel

    ! Assume a really small (physical) diffusivity (instead of zero)
    kappa = 1.0d-14

    ! Evaluate the element Peclet number
    alpha = 0.5d0*abs(v(e))*Lc(e)/kappa

    ! Optimal numerical diffusivity (Hughes, Mallet and Mizukami; 1986)
    if (abs(alpha)<rzero) then
      xi = 0.d0
    else
      xi = 1.d0/tanh(alpha) - 1.d0/alpha
    end if

    ! evaluate SUPG stabilization parameter, tau=L/(2*|v|), multiplied by the velocity (v)
    tauv = sign(0.5d0*Lc(e)*xi,v(e))

    ! evaluate and assemble SUPG mass matrix
!~     M_supg(e,e)        = M_supg(e,e)     - 0.5d0*tauv
!~     M_supg(e,e+1)      = M_supg(e,e+1)   - 0.5d0*tauv
!~     M_supg(e+1,e)      = M_supg(e+1,e)   + 0.5d0*tauv
!~     M_supg(e+1,e+1)    = M_supg(e+1,e+1) + 0.5d0*tauv
    M_supg(e,1)        = M_supg(e,1)   - 0.5d0*tauv
    M_supg(e,3)        = M_supg(e,3)   - 0.5d0*tauv
    M_supg(e,2)        = M_supg(e,2)   + 0.5d0*tauv
    M_supg(e+1,1)      = M_supg(e+1,1) + 0.5d0*tauv

    ! Mass lumping via row summing
    M_supg_lumped(e)   = M_supg_lumped(e)   - tauv
    M_supg_lumped(e+1) = M_supg_lumped(e+1) + tauv

    ! evaluate and assemble SUPG flux correction
    flux_supg(e)       = flux_supg(e)   - tauv*v(e)/Lc(e) * (c(e+1)-c(e))
    flux_supg(e+1)     = flux_supg(e+1) + tauv*v(e)/Lc(e) * (c(e+1)-c(e))

    ! evaluate and assemble SUPG source correction
    source_supg(e)     = source_supg(e)   - tauv * src(e)
    source_supg(e+1)   = source_supg(e+1) + tauv * src(e)

  end do  !e
!~   !$OMP END PARALLEL DO

  RETURN
END SUBROUTINE SUPG_Arrays

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE CDT_Combine_State(state_A, state_B, mult_A, mult_B, state_new)
  USE GlobalParams
  USE CDT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_CDT_STATE),   INTENT(IN)  :: state_A, state_B
  REAL(KIND=8), INTENT(IN)                  :: mult_A, mult_B
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_CDT_STATE),   INTENT(OUT) :: state_new
!-----------------------------------------------------------------------

  call state_new%allocatememory()
  state_new%cell_sizes   =        state_A%cell_sizes
  state_new%varrho_field = mult_A*state_A%varrho_field + mult_B*state_B%varrho_field
  state_new%varrho_node  = mult_A*state_A%varrho_node  + mult_B*state_B%varrho_node

  RETURN
END SUBROUTINE CDT_Combine_State

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!> Compute the dislocation source terms.
SUBROUTINE CDT_Src_Calc(cdt_state_old, ddc_state_new, dmb_state_new, all_cdt_props, &
           varrho_src, velocity)
  USE GlobalParams
  USE CDT
  USE DDC
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_CDT_State),         INTENT(IN)  :: cdt_state_old
  TYPE(Type_DDC_State),         INTENT(IN)  :: ddc_state_new
  TYPE(Type_DMB_State),         INTENT(IN)  :: dmb_state_new
  TYPE(Type_CDT_PROPS),         INTENT(IN)  :: all_cdt_props(Nregion)
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nel,Nslip,Nchar),  INTENT(OUT) :: varrho_src
  REAL(KIND=8), DIMENSION(Nel,Nslip,Nchar),  INTENT(OUT) :: velocity
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: iel,islip,jslip,ich,jch, iregion
  REAL(KIND=8) :: M,Aint,A,Y_e,G_n0,tau_n0,rhoDot_n0,pn,qn,boltz,arg1
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: s_mult, s_ann, s_nuc_hom
  REAL(KIND=8), DIMENSION(:,:), allocatable   :: tau, tau_eff, varrho_obs
  REAL(KIND=8), DIMENSION(Nel,1)        :: tmp
  REAL(KIND=8), DIMENSION(:,:,:), allocatable  :: varrho_pos, varrho_neg, varrho_tot
  REAL(KIND=8), DIMENSION(Nslip, Nslip, Nchar) :: Amat
!-----------------------------------------------------------------------
  allocate(s_mult(Nel,Nslip,Nchar), s_ann(Nel,Nslip,Nchar), s_nuc_hom(Nel,Nslip,Nchar))
  allocate(tau(Nel,Nslip), tau_eff(Nel,Nslip), varrho_obs(Nel,Nslip))
  allocate(varrho_pos(Nel, Nslip, Nchar), varrho_neg(Nel, Nslip, Nchar), varrho_tot(Nel, Nslip, Nchar))

!  Initialize source/sink arrays.
  s_mult    = 0.d0
  s_ann     = 0.d0
  s_nuc_hom = 0.d0

!  Unpack state variables needed for source/sink calculations    
  velocity   = dmb_state_new%vel_x
  tau        = dmb_state_new%tau
  tau_eff    = tau - ddc_state_new%tau_b
  tmp        = dmb_state_new%T
  varrho_pos = cdt_state_old%varrho_field(:,:,:,1)
  varrho_neg = cdt_state_old%varrho_field(:,:,:,2)
  varrho_tot = varrho_pos + varrho_neg

  !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(iregion, M, Aint, A, Y_e, G_n0, tau_n0, rhoDot_n0, pn, qn, &
  !$OMP                              boltz, Amat, iel, islip, jch, jslip, ich, arg1)
  do iregion=1,nregion

!  Unpack CDT props
  M         = 1.d0/all_cdt_props(iregion)%K
  Aint      = all_cdt_props(iregion)%Aint
  A         = all_cdt_props(iregion)%A
  Y_e       = all_cdt_props(iregion)%Y_e
  G_n0      = all_cdt_props(iregion)%G_n0
  tau_n0    = all_cdt_props(iregion)%tau_n0
  rhoDot_n0 = all_cdt_props(iregion)%rhoDot_n0
  pn        = all_cdt_props(iregion)%pn
  qn        = all_cdt_props(iregion)%qn
  boltz     = all_cdt_props(iregion)%boltz
!  Amat = Aint*nmp.ones((Nslip,Nslip))
  if (Aint < 0.d0) then
    Amat = all_cdt_props(iregion)%Aforest !! use generalization of Aint if user provided a negative value
  else
    Amat(:,:,:)  = Aint
  end if
!  varrho_obs = nmp.transpose(nmp.dot(Amat,nmp.transpose(varrho_tot)))
  do iel= region_mask(iregion,1), region_mask(iregion,2)
      do islip=1,Nslip
        varrho_obs(iel,islip) = 0.d0
        do jch=1,Nchar
          do jslip=1,Nslip
            varrho_obs(iel,islip) = varrho_obs(iel,islip) + Amat(islip,jslip,jch) * varrho_tot(iel,jslip,jch)
          end do !jslip
        end do  !jch
      end do !islip
 
!  Compute source/sink terms       
    do ich=1,Nchar
      do islip=1,Nslip
        s_mult(iel,islip,ich) =  M*sqrt(varrho_obs(iel,islip))*varrho_tot(iel,islip,ich)*abs(velocity(iel,islip,ich))
        s_ann(iel,islip,ich)  = -A*Y_e*varrho_pos(iel,islip,ich)*varrho_neg(iel,islip,ich)*abs(velocity(iel,islip,ich))
        arg1 = (abs(tau_eff(iel,islip))/tau_n0)**pn
        if (arg1 < 1.d0) then
          s_nuc_hom(iel,islip,ich) = rhoDot_n0*exp(-G_n0*(1.d0 - arg1)**qn/boltz/tmp(iel,1))
        else
          s_nuc_hom(iel,islip,ich) = rhoDot_n0
        end if
      end do  !islip
    end do  !ich
  end do  !iel
  end do !iregion
  !$OMP END PARALLEL DO

  varrho_src = s_mult + s_ann + s_nuc_hom

  RETURN
END SUBROUTINE CDT_Src_Calc

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------
end module CDT_routines
