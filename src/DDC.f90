module DDC_routines
public
contains
SUBROUTINE DDC_State_Init(cell_sizes, ddc_state)
  USE GlobalParams
  USE DDC
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nel), INTENT(IN)  :: cell_sizes
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_DDC_STATE),   INTENT(INOUT) :: ddc_state
!-----------------------------------------------------------------------
!~   call ddc_state%allocatememory()
  ddc_state%cell_sizes(:) = cell_sizes

  RETURN
END SUBROUTINE DDC_State_Init

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DDC_Update(ddc_state_old, dmb_state, cdt_state,  &
           all_props_ddc, all_props_dmb, &
           ddc_state_new)
  USE GlobalParams
  USE DDC
  USE CDT
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DDC_STATE),   INTENT(IN) :: ddc_state_old
  TYPE(Type_DMB_STATE),   INTENT(IN) :: dmb_state
  TYPE(Type_CDT_STATE),   INTENT(IN) :: cdt_state
  TYPE(Type_DDC_PROPS),   INTENT(IN) :: all_props_ddc(Nregion)
  TYPE(Type_DMB_PROPS),   INTENT(INOUT) :: all_props_dmb(Nregion)
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_DDC_STATE),   INTENT(INOUT) :: ddc_state_new
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(Nnode)     :: X_ref
  REAL(KIND=8), DIMENSION(3,3,3,3)   :: Cinv
  REAL(KIND=8), DIMENSION(Nslip,3,3) :: schmid0
!~   REAL(KIND=8), DIMENSION(Nel,Nslip,Nchar) :: tau_internal
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: pk2_i
!  REAL(KIND=8), DIMENSION(3,3)       :: Csum
  INTEGER :: iip,islip,i,j,ich,iregion,e1,e2,n
  REAL(KIND=8) :: sum
!-----------------------------------------------------------------------
  allocate(pk2_i(Nel,3,3))
!~   call ddc_state_new%allocatememory()
  
  X_ref   = dmb_state%X_ref

  ! initialize state:
  call DDC_State_Init(ddc_state_old%cell_sizes,ddc_state_new)

  ! compute slip-system specific internal stress
!~   if (all_props_ddc(1)%backstress_model.eq.'zaiser') then
!~     print*,dmb_state%vel(1,1,1),dmb_state%vel(Nel/4,1,1)
!~   end if
  
  !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(iregion,e1,e2,n)
  do iregion=1,Nregion
    e1 = region_mask(iregion,1)      ! region-left-most element ID
    e2 = region_mask(iregion,2)      ! region-right-most element ID
    n = e2-e1 + 1                    ! number of elements in region
    
    call DDC_Compute_Backstress(n,cdt_state%varrho_field(e1:e2,:,:,:), &
                             all_props_ddc(iregion), all_props_dmb(iregion), dmb_state%vel(e1:e2,:,:), X_ref(e1:e2+1), &
                                       ddc_state_new%tau_i(e1:e2,:,:) )

  end do
  !$OMP END PARALLEL DO
  
  ! loop over integration points
  pk2_i = 0.d0

  !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(iregion,Cinv,schmid0,iip,ich,islip,i,j)
  do iregion = 1,nregion

    Cinv    = all_props_dmb(iregion)%Cinv
    schmid0 = all_props_dmb(iregion)%schmid0

    do iip=region_mask(iregion,1), region_mask(iregion,2)
       
    ! assemble into internal stress
    do ich=1,Nchar
      do islip=1,Nslip
        pk2_i(iip,:,:) =  pk2_i(iip,:,:) +  &
          0.5d0*ddc_state_new%tau_i(iip,islip,ich)*(schmid0(islip,:,:) + TRANSPOSE(schmid0(islip,:,:))) ! pk2 int. conf.
      end do  !islip
    end do  !ich

    ! compute internal strain
    do i=1,3
      do j=1,3
        ddc_state_new%epsilon_i(iip,i,j) = SUM(Cinv(i,j,:,:)*pk2_i(iip,:,:))
      end do  !j
    end do  !i
  
    ! compute backstress
    do islip=1,Nslip
      ddc_state_new%tau_b(iip, islip) = SUM(pk2_i(iip,:,:)*schmid0(islip,:,:))
    end do  !islip
  end do  !iip
  end do ! iregion
  !$OMP END PARALLEL DO

  ddc_state_new%stress_i = pk2_i

  RETURN
END SUBROUTINE DDC_Update

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DDC_Compute_Backstress(n, varrho, ddc_props, dmb_props, vel, X_ref,  &
           tau_back)
  USE GlobalParams
  USE DDC
  USE DMB
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  integer, intent(in)                         :: n
  REAL(KIND=8), INTENT(IN), DIMENSION(N,Nslip,Nchar,2)  :: varrho
  REAL(KIND=8), INTENT(IN), DIMENSION(N,Nslip,Nchar)  :: vel
  TYPE(Type_DDC_PROPS), intent(in)            :: ddc_props
  TYPE(Type_DMB_PROPS), intent(inout)         :: dmb_props  
  REAL(KIND=8), INTENT(IN), DIMENSION(N+1)        :: X_ref
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), INTENT(OUT), DIMENSION(N,Nslip,Nchar)   :: tau_back
!-----------------------------------------------------------------------
!  Locals:
!-----------------------------------------------------------------------
  if(ddc_props%backstress_model=='gradient') then
    call DDC_Gradient_Backstress(n, varrho, ddc_props, X_ref, &
         tau_back)
  elseif((ddc_props%backstress_model=='zaiser').or.(ddc_props%backstress_model=='zaiser_zero')) then
    call DDC_Zaiser_Backstress(n, varrho, ddc_props, dmb_props, vel, X_ref, &
         tau_back)
!~   elseif(ddc_props%backstress_model.eq.'field') then
!~     call DDC_Field_Backstress(n, varrho, ddc_props, X_ref, &
!~          tau_back)
  elseif((ddc_props%backstress_model=='bypass').or.(ddc_props%backstress_model=='none')) then
    tau_back = 0.d0
  else  ! don't silently fall back to a default method
    call Fatal('unknown backstress_model flag: '//ddc_props%backstress_model)
  end if

  RETURN
END SUBROUTINE DDC_Compute_Backstress

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DDC_Zaiser_Backstress_noKappa(vel, ddc_props, dmb_props, tau_back_nokappa)
  USE GlobalParams
  USE DDC
  USE DMB
  use DMB_routines, only: DMB_Compute_Rotation, stroh_dislocation
  use utilities, only: LinSpace, operator(.otimes.), Trapz
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), INTENT(IN), DIMENSION(Nslip,Nchar)   :: vel
  TYPE(Type_DDC_PROPS), intent(in)            :: ddc_props
  TYPE(Type_DMB_PROPS), intent(in)            :: dmb_props
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), INTENT(OUT), DIMENSION(3,Nslip,Nchar) :: tau_back_nokappa
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: burger, tmpintegral
  INTEGER, PARAMETER :: resolution=361 !< resolve polar angle phi to 1 degree (tradeoff accuracy for speed)
  REAL(KIND=8), DIMENSION(:,:,:,:,:), allocatable :: uij
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: tau_internal
  REAL(KIND=8), DIMENSION(resolution) :: phi
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar) :: Sb, Bb
  REAL(KIND=8), DIMENSION(3,3,Nslip,Nchar) :: rotM
  REAL(KIND=8), DIMENSION(3,3,3,3) :: elC
  INTEGER :: i,j,k,l,isl,ich
  REAL(KIND=8), DIMENSION(3,Nslip)     :: slip_S !, slip_Sg !, slip_M, slip_Mg
  REAL(KIND=8), DIMENSION(3)           :: euler
  REAL(KIND=8), DIMENSION(3,3)         :: Q
  REAL(KIND=8), DIMENSION(Nslip,3,3) :: schmid0
  
  allocate(uij(3,3,resolution,Nslip,Nchar), tau_internal(resolution,Nslip,Nchar))
  schmid0 = 0.d0
  
  burger = ddc_props%burger
  slip_S = ddc_props%slip_S
  euler  = ddc_props%euler_angle
  
  rotM   = dmb_props%rotM ! rotation matrix that rotates from crystal coordinates to coordinates aligned with a disloc
  call stroh_dislocation(dmb_props,vel(:,:),resolution,uij,Sb, Bb,elC)
  
  do i=1,Nslip
    schmid0(i,:,:) = dmb_props%slip_S(:,i) .otimes. dmb_props%slip_M(:,i)
  end do
  
  tau_internal= 0.d0
  do ich=1,Nchar
    do isl=1,Nslip
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              tau_internal(:,isl,ich) = tau_internal(:,isl,ich) + schmid0(isl,i,j)*elC(i,j,k,l)*uij(k,l,:,isl,ich)
            end do
          end do
        end do
      end do
    end do
  end do
  
  
  ! get rotation matrix to rotate into global coordinates
  call DMB_Compute_Rotation(euler,Q)

  ! derived back stress from tau_internal and kappa
  call LinSpace(0.d0,2.d0*pi,resolution,phi)
  tau_back_nokappa = 0.d0
  do i=1,Nchar
    do j=1,Nslip
      ! Uij = uij/r, is then multiplied by x=r(cosphi,sinphi), so r cancels, integral is over r dr dphi, 
      ! and int dr = average dislocation separation Lbar
      ! in general, need gradient del kappa / del vec(r), but after rotating to global coordinates
      ! we only keep the x-component in our current 1D setup
      !
      ! TODO: include disloc-disloc correlations according to Zaiser et al.
      call Trapz(resolution,(cos(phi(:)))*tau_internal(:,j,i),phi,tmpintegral)
      tau_back_nokappa(1,j,i) = -burger*tmpintegral ! x-component (parallel to slip_Sg for edge)
      call Trapz(resolution,(sin(phi(:)))*tau_internal(:,j,i),phi,tmpintegral)
      tau_back_nokappa(2,j,i) = -burger*tmpintegral ! y-component (=0 for isotropic edge)
      ! rotate integral into global coordinates (will only need the x-component of that below for our 1D problem):
      tau_back_nokappa(:,j,i) = MATMUL(transpose(rotM(:,:,j,i)),tau_back_nokappa(:,j,i)) ! rotate first into crystal coordinates
      tau_back_nokappa(:,j,i) = MATMUL(Q,tau_back_nokappa(:,j,i)) ! rotate into global coordinates
    end do !Nslip
  end do !Nchar
  
  RETURN
END SUBROUTINE DDC_Zaiser_Backstress_noKappa

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DDC_Zaiser_Backstress(n, varrho, ddc_props, dmb_props, vel_state, X_ref, &
           tau_back)
  USE GlobalParams
  USE DDC
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  integer, intent(in)                         :: n
  REAL(KIND=8), INTENT(IN), DIMENSION(N,Nslip,Nchar,2)  :: varrho
  REAL(KIND=8), INTENT(IN), DIMENSION(N,Nslip,Nchar)  :: vel_state
  TYPE(Type_DDC_PROPS), intent(in)            :: ddc_props
  TYPE(Type_DMB_PROPS), intent(inout)         :: dmb_props
  REAL(KIND=8), INTENT(IN), DIMENSION(N+1)        :: X_ref
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), INTENT(OUT), DIMENSION(N,Nslip,Nchar)   :: tau_back
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: Lbar
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar) :: integral, integral_zero
  INTEGER :: el
  REAL(KIND=8), DIMENSION(N,Nslip,Nchar)   :: rho_total, kappa, grad_kappa
  REAL(KIND=8), DIMENSION(N+1,Nslip,Nchar) :: kappa_face !! djl--> remove = 0.d0
  REAL(KIND=8), DIMENSION(N) :: dX
  REAL(KIND=8), DIMENSION(Nslip,Nchar)   :: vel_zero
  REAL(KIND=8), DIMENSION(N,Nslip,Nchar)   :: vel
  
  ! decide whether to compute back stress for every velocity or just for the static disloc. field
  if (ddc_props%backstress_model=='zaiser') then
    vel = vel_state
  else
    vel = 0.d0
  end if
  vel_zero = 0.d0
  Lbar   = dmb_props%Lbar ! length scale associated with initial total disloc. density
  ! this length scale was used in gradient back stress to approximate the dislocation spacing
  ! here we might try something more sophisticated, namely the local disloc density rho_total
  
  ! compute total dislocation density and GND density fields
  rho_total(:,:,:) = varrho(:,:,:,1) + varrho(:,:,:,2)      
  kappa(:,:,:)     = varrho(:,:,:,1) - varrho(:,:,:,2)   

  ! interpolate kappa to cell edges
  kappa_face(2:N,:,:) = 0.5d0*(kappa(1:N-1,:,:) + kappa(2:N,:,:))
  kappa_face(N+1,:,:) = ( 3.d0/8.d0)*kappa(N-2,:,:) - (5.d0/4.d0)*kappa(N-1,:,:) + (15.d0/8.d0)*kappa(N,:,:)  
  kappa_face(1,:,:)     = (15.d0/8.d0)*kappa(1    ,:,:) - (5.d0/4.d0)*kappa(2    ,:,:) + ( 3.d0/8.d0)*kappa(3  ,:,:)
  
  ! compute only once for zero velocity and store within dmb_state, then retrieve at later time steps
  if (sum(abs(dmb_props%tau_back_nokappa))<1.d-15) then
    call DDC_Zaiser_Backstress_noKappa(vel_zero(:,:), ddc_props, dmb_props, integral_zero)
    dmb_props%tau_back_nokappa(:,:,:) = integral_zero
  else
    integral_zero = dmb_props%tau_back_nokappa(:,:,:)
  end if
  
  dX(1:N) = X_ref(2:N+1) - X_ref(1:N)
  do el=1,N
    ! only recompute backstress if disloc. velocity (summed over Nchar,Nslip) at el is > vcrit/1000, take zero value otherwise
    if (sum(abs(vel(el,:,:)))<1.d-3*Nslip*sum(dmb_props%wave_vel(:))) then
      integral = integral_zero
    else
      ! TODO: pre-compute & store this for N velocities between vcrit/1000 and vcrit, then interpolate (faster)
      call DDC_Zaiser_Backstress_noKappa(vel(el,:,:), ddc_props, dmb_props, integral)
    end if
  ! estimate grad(kappa) for backstress calculation (dK / dX) note: dK/dS = dK/dX * (dX/ds)
    grad_kappa(el,:,:) = (kappa_face(el+1,:,:)-kappa_face(el,:,:)) / dX(el)
!~     tau_back(el,:,:)   = (Lbar*Lbar)*integral(1,:,:)*grad_kappa(el,:,:)
    tau_back(el,:,:)   = integral(1,:,:)*grad_kappa(el,:,:)/rho_total(el,:,:) ! only need x-comp. of integral for 1D problem
  end do  !el
  
  RETURN
END SUBROUTINE DDC_Zaiser_Backstress

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DDC_Gradient_Backstress(n, varrho, ddc_props, X_ref, &
           tau_back)
  USE GlobalParams
  USE DDC
  use DMB_routines, only: DMB_Compute_Rotation
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  integer, intent(in)                         :: n
  REAL(KIND=8), INTENT(IN), DIMENSION(N,Nslip,Nchar,2)  :: varrho
  TYPE(Type_DDC_PROPS), intent(in)            :: ddc_props
  REAL(KIND=8), INTENT(IN), DIMENSION(N+1)        :: X_ref
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), INTENT(OUT), DIMENSION(N,Nslip,Nchar)   :: tau_back
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) Lrho, mu, burger
  REAL(KIND=8), DIMENSION(3,Nslip)     :: slip_S, slip_Sg
  REAL(KIND=8), DIMENSION(3)           :: euler
  REAL(KIND=8), DIMENSION(3,3)         :: Q
  REAL(KIND=8), DIMENSION(N,Nslip,Nchar)   :: rho_total, kappa, grad_kappa
  REAL(KIND=8), DIMENSION(N+1,Nslip,Nchar) :: kappa_face !! djl--> remove = 0.d0
  REAL(KIND=8), DIMENSION(N) :: dX
  INTEGER i,j,el,ich
!-----------------------------------------------------------------------
  ! unpack DDC properties
  Lrho   = ddc_props%Lrho
  mu     = ddc_props%mu
  burger = ddc_props%burger
  slip_S = ddc_props%slip_S
  euler  = ddc_props%euler_angle

  ! compute total dislocation density and GND density fields
  rho_total(:,:,:) = varrho(:,:,:,1) + varrho(:,:,:,2)      
  kappa(:,:,:)     = varrho(:,:,:,1) - varrho(:,:,:,2)   

  ! interpolate kappa to cell edges
  kappa_face(2:N,:,:) = 0.5d0*(kappa(1:N-1,:,:) + kappa(2:N,:,:))
  kappa_face(N+1,:,:) = ( 3.d0/8.d0)*kappa(N-2,:,:) - (5.d0/4.d0)*kappa(N-1,:,:) + (15.d0/8.d0)*kappa(N,:,:)  
  kappa_face(1,:,:)     = (15.d0/8.d0)*kappa(1    ,:,:) - (5.d0/4.d0)*kappa(2    ,:,:) + ( 3.d0/8.d0)*kappa(3  ,:,:)

  ! get slip vector into global coordinates
  call DMB_Compute_Rotation(euler,Q)
  do i=1,Nslip
    do j=1,3
      slip_Sg(j,i) = DOT_PRODUCT(Q(j,:),slip_S(:,i))
    end do  !j
  end do  !i

  dX(1:N) = X_ref(2:N+1) - X_ref(1:N)
  do el=1,N
  ! estimate grad(kappa) for backstress calculation (dK / dX) note: dK/dS = dK/dX * (dX/ds)
    grad_kappa(el,:,:) = (kappa_face(el+1,:,:)-kappa_face(el,:,:)) / dX(el)
    ! Comment: should generalize for anisotropic case and all characters (i.e. prefactor mubLrho will change to prefact(ich)!
    ! but this is essentially done within the DDC_Zaiser_Backstress() implementation
    do ich=1,Nchar
      tau_back(el,:,ich)   = mu*burger*Lrho*grad_kappa(el,:,ich)*slip_Sg(1,:)
    end do !ich
  end do  !el
  
  RETURN
END SUBROUTINE DDC_Gradient_Backstress

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!~ SUBROUTINE DDC_Field_Backstress(n, varrho, ddc_props, X_ref, tau_back)
!~   USE GlobalParams
!~   USE DDC
!~   IMPLICIT NONE
!~ !-----------------------------------------------------------------------
!~ !  Inputs:
!~   integer, intent(in)                         :: n
!~   REAL(KIND=8), DIMENSION(N,Nslip,Nchar,2), INTENT(IN)  :: varrho
!~   TYPE(Type_DDC_PROPS),           INTENT(IN)  :: ddc_props
!~   REAL(KIND=8), DIMENSION(N+1),       INTENT(IN)  :: X_ref
!~ !-----------------------------------------------------------------------
!~ !  Outputs:
!~   REAL(KIND=8), DIMENSION(N,Nslip,Nchar),   INTENT(OUT) :: tau_back
!~ !-----------------------------------------------------------------------
!~ !  Locals:
!~   REAL(KIND=8) :: Lrho,mu,burger
!~   INTEGER :: ii,ich
!~   REAL(KIND=8), DIMENSION(N)         :: Xip, F
!~   REAL(KIND=8), DIMENSION(N,Nslip,Nchar)   :: rho_total
!~   REAL(KIND=8), DIMENSION(N,Nslip,Nchar)   :: kappa
!~ !-----------------------------------------------------------------------
!~   Xip(1:N) = 0.5d0 * (X_ref(1:N) + X_ref(2:N+1))

!~   ! unpack DDC properties
!~   Lrho   = ddc_props%Lrho
!~   mu     = ddc_props%mu
!~   burger = ddc_props%burger

!~   ! compute total dislocation density and GND density fields
!~   rho_total(:,:,:) = varrho(:,:,:,1) + varrho(:,:,:,2)      
!~   kappa(:,:,:)     = varrho(:,:,:,1) - varrho(:,:,:,2)   

!~   do ich=1,Nchar
!~     do ii=1,Nslip
!~       call Compute_CPV(kappa(:,ii,ich),Xip,X_ref(1),X_ref(N+1),F)
!~       tau_back(:,ii,ich) = mu*burger / (2.d0*pi) * F(:)
!~     end do  !ii
!~   end do  !ich

!~   RETURN
!~ END SUBROUTINE DDC_Field_Backstress

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!~ SUBROUTINE Compute_CPV(f, x, a, b, Fout)
!~   ! Compute the Cauchy Principal Value I(x) = Int f(t)/(t-x) dt a<x<b
!~   ! for functions values f given at locations x.  
!~   USE GlobalParams
!~   use utilities, only: Fatal
!~   IMPLICIT NONE
!~ !-----------------------------------------------------------------------
!~ !  Inputs:
!~   REAL(KIND=8), DIMENSION(Nel),         INTENT(IN)  :: f
!~   REAL(KIND=8), DIMENSION(Nel),         INTENT(IN)  :: x
!~   REAL(KIND=8),                         INTENT(IN)  :: a
!~   REAL(KIND=8),                         INTENT(IN)  :: b
!~ !-----------------------------------------------------------------------
!~ !  Outputs:
!~   REAL(KIND=8), DIMENSION(Nel),         INTENT(OUT) :: Fout
!~ !-----------------------------------------------------------------------
!~ !  Locals:
!~   REAL(KIND=8), DIMENSION(Nel)         :: phit,I1,I2
!~ !  REAL(KIND=8), DIMENSION(Nel,Nslip)   :: rho_total
!~ !  REAL(KIND=8), DIMENSION(Nel,Nslip)   :: kappa
!~   REAL(KIND=8) :: fx,x_
!~   INTEGER :: ix
!~ !-----------------------------------------------------------------------
!~   call Fatal("Not Implemented! Need a composite Simpson's Rule function.")
!~   do ix=1,Nel
!~     x_ = x(ix)
!~     fx = f(ix)
!~     phit     = (f - fx) / (x-x_)
!~     if (ix .eq. 1) then
!~       phit(ix) = (f(ix+1) - f(ix)) / (x(ix+1)-x(ix))
!~     elseif (ix .eq. Nel) then
!~       phit(ix) = (f(ix) - f(ix-1)) / (x(ix)-x(ix-1))
!~     else
!~       phit(ix) = (f(ix+1) - f(ix-1)) / (x(ix+1)-x(ix-1))
!~     end if
    
!~     I1(ix) = 0.d0  ! = si.simps(phit, x=x)
!~     I2(ix) = f(ix)*log( (b-x_)/(x_-a))
!~   end do  !ix
  
!~   Fout  = I1 + I2

!~   RETURN
!~ END SUBROUTINE Compute_CPV
end module DDC_routines
