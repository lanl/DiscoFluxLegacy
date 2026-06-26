!>-----------------------------------------------------------------------
!>  Driver for 1D uniaxial deformation as well as shear deformation
!>-----------------------------------------------------------------------
! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
PROGRAM CaseA
  USE GlobalParams
  USE DMB
  USE CDT
  USE DDC
  use utilities
  use DMB_routines
  use CDT_routines
  use DDC_routines
  use OUTPUTFILES
  use readfiles
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
  CHARACTER(10)               :: simulation_type
  CHARACTER(6)                :: crystalstruct
  REAL(KIND=8)                :: impact_force, impact_velocity, impact_acc_time
  REAL(KIND=8)                :: impact_angle
  REAL(KIND=8)                :: L0
  REAL(KIND=8)                :: rhobar0
  REAL(KIND=8)                :: temperature0, specific_heat
  REAL(KIND=8)                :: C11, CT11, CP11
  REAL(KIND=8)                :: C12, CT12, CP12
  REAL(KIND=8)                :: C44, CT44, CP44
  REAL(KIND=8)                :: bulk_c1
  REAL(KIND=8)                :: bulk_c2
  REAL(KIND=8)                :: c1
  REAL(KIND=8)                :: mu
  REAL(KIND=8)                :: rho0
  REAL(KIND=8)                :: burger
  REAL(KIND=8), DIMENSION(:), allocatable :: B0, wave_vel
  REAL(KIND=8)                :: Lbar
  REAL(KIND=8)                :: g0
  REAL(KIND=8)                :: omega0
  REAL(KIND=8)                :: boltz
  REAL(KIND=8)                :: p
  REAL(KIND=8)                :: q
  REAL(KIND=8)                :: tau0
  REAL(KIND=8)                :: cdt_K
  REAL(KIND=8)                :: cdt_A
  REAL(KIND=8)                :: cdt_Y_e
  REAL(KIND=8)                :: cdt_Aint
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: Aforest !< generalization of Aint, used if input variable Aint<0
  REAL(KIND=8)                :: cdt_rhoDot_n0
  REAL(KIND=8)                :: cdt_G_n0
  REAL(KIND=8)                :: cdt_tau_n0
  REAL(KIND=8)                :: cdt_pn
  REAL(KIND=8)                :: cdt_qn

  logical               :: eos_flag
  real(KIND=8)                :: eos_v_star
  real(KIND=8)                :: eos_B_star
  real(KIND=8)                :: eos_B1
  real(KIND=8)                :: eos_rho0
  real(KIND=8)                :: eos_v0
  real(KIND=8)                :: eos_Td0
  real(KIND=8)                :: eos_a
  real(KIND=8)                :: eos_b
  real(KIND=8)                :: eos_gamma
  real(KIND=8)                :: eos_Gamma0
  real(KIND=8)                :: eos_kappa
  real(KIND=8)                :: eos_R_M

  LOGICAL               :: no_flux
  LOGICAL               :: internal_strain_flag, drag_Trho_flag
  REAL(KIND=8), DIMENSION(:,:), allocatable :: euler_angle
  REAL(KIND=8)                :: dt0
  REAL(KIND=8)                :: t_stop
  REAL(KIND=8)                :: dt_output
  REAL(KIND=8)                :: t_vel
  REAL(KIND=8)                :: gtol, tolerance=sqrt(rzero)
  REAL(KIND=8)                :: Lrho
  CHARACTER(32)         :: CDT_update
  CHARACTER(32)         :: CDTscheme
  CHARACTER(32)         :: CDTlimiter
  CHARACTER(32)         :: CDTintegrator
  CHARACTER(32)         :: CDTintFlux
  character(32)         :: drag_flag
  character(32)         :: backstress_model
  CHARACTER(32)         :: jobname
!-----------------------------------------------------------------------
  REAL(KIND=8), DIMENSION(:,:), allocatable     :: slip_M, slip_S
  REAL(KIND=8), DIMENSION(:,:,:), allocatable     :: slip_V, line_T
  REAL(KIND=8), DIMENSION(:,:,:,:), allocatable :: rotM
!~   TYPE(Type_DMB_PROPS)           :: props_dmb
  TYPE(Type_CDT_PROPS)           :: props_cdt
!~   TYPE(Type_DDC_PROPS)           :: props_ddc
  TYPE(Type_DMB_PROPS), dimension(:), allocatable :: all_props_dmb
  TYPE(Type_CDT_PROPS), dimension(:), allocatable :: all_props_cdt
  TYPE(Type_DDC_PROPS), dimension(:), allocatable :: all_props_ddc
  REAL(KIND=8), DIMENSION(4)           :: pressure
  REAL(KIND=8), DIMENSION(4)           :: free_energy
  REAL(KIND=8), DIMENSION(3)           :: entropy
  REAL(KIND=8), DIMENSION(3)           :: thermo_coeffs
  REAL(KIND=8), DIMENSION(:), allocatable  :: X_ref, lumped_mass
  REAL(KIND=8), DIMENSION(:,:), allocatable :: x_pos, x_vel, xddot
  REAL(KIND=8), DIMENSION(:), allocatable  :: cell_sizes
  REAL(KIND=8), DIMENSION(2,2)         :: matrix1
!~   REAL(KIND=8), DIMENSION(:,:), allocatable :: consistent_mass ! never used below?
  INTEGER                        :: iel
  REAL(KIND=8)                         :: cs
  REAL(KIND=8)                         :: dX
  REAL(KIND=8)                         :: dt_max, dt
  REAL(KIND=8)                         :: t_check
  REAL(KIND=8), DIMENSION(3,3)         :: eye
  TYPE(Type_DMB_STATE)           :: dmb_state_old
  TYPE(Type_DMB_STATE)           :: dmb_state_iter, dmb_state_prev
  TYPE(Type_DMB_VECTOR_BC)       :: bc_xdot, bc_force
  TYPE(Type_CDT_STATE)           :: cdt_state_old
  TYPE(Type_CDT_STATE)           :: cdt_state_i, cdt_state_ii, cdt_state_iter, cdt_state_prev
  TYPE(Type_DDC_STATE)           :: ddc_state_old
  TYPE(Type_DDC_STATE)           :: ddc_state_iter, ddc_state_prev
  !REAL(KIND=8), DIMENSION(Nel,Nslip,2) :: debug_decouple
  REAL(KIND=8), DIMENSION(:), allocatable :: flux_bc
  INTEGER                        :: iinc, iout, ipOUT
  REAL(KIND=8)                         :: t_output
  LOGICAL                        :: stop_flag
  INTEGER                        :: i,ich,iregion
  INTEGER                        :: kiter, MaxIter
  LOGICAL                        :: konverged
  REAL(KIND=8), DIMENSION(2)           :: Eddc, Edmb, Ecdt
  REAL(KIND=8)                         :: DDCtoler, DMBtoler, CDTtoler, start_time, finish_time

  ! read input
  call readcmdline(jobname)
  
  call ReadSizes(jobname)
  allocate(B0(Nchar), wave_vel(Nchar), euler_angle(Nregion,3))

  call ReadInput(jobname, simulation_type, B0, C11, C12, C44, CDT_update, CDTintegrator, &
  CDTscheme, CDTlimiter, CDTtoler, CDTintFlux, CP11, CP12, CP44, CT11, CT12, CT44, specific_heat, &
  DDCtoler, DMBtoler, L0, Lbar, Lrho, MaxIter, Boltz, bulk_c1, bulk_c2, &
  crystalstruct, burger, c1, cdt_A, cdt_Aint, cdt_G_n0, cdt_K, cdt_Y_e, &
  cdt_pn, cdt_qn, cdt_rhoDot_n0, cdt_tau_n0, dt0, dt_output, euler_angle, g0, gtol, impact_force, &
  impact_angle, impact_velocity, impact_acc_time, internal_strain_flag, ipOUT, mu, no_flux, omega0, p, &
  eos_B1, eos_B_star, eos_Gamma0, eos_R_M, eos_Td0, eos_a, eos_b, eos_flag, eos_gamma, eos_kappa, eos_rho0, &
  eos_v0, eos_v_star, q, rho0, rhobar0, t_stop, t_vel, tau0, temperature0, wave_vel, drag_flag, drag_Trho_flag, &
  backstress_model)
  
  call cpu_time(start_time)
 
  !! allocate memory for arrays that are too large to fit in the stack:
!~   allocate(consistent_mass(Nnode,Nnode))
  allocate(flux_bc(Nregion+1))
  allocate(Aforest(Nslip,Nslip,Nchar))
  allocate(slip_M(3,Nslip), slip_S(3,Nslip), rotM(3,3,Nslip,Nchar))
  allocate(slip_V(3,Nslip,Nchar), line_T(3,Nslip,Nchar))
  allocate(X_ref(Nnode), lumped_mass(Nnode), cell_sizes(Nel))
  allocate(x_pos(Nnode,3), x_vel(Nnode,3), xddot(Nnode,3))
  
  call props_cdt%allocatememory()
  allocate(all_props_dmb(Nregion), all_props_cdt(Nregion), all_props_ddc(Nregion))
  do i=1,Nregion
    call all_props_dmb(i)%allocatememory()
    call all_props_cdt(i)%allocatememory()
    call all_props_ddc(i)%allocatememory()
  end do
  call dmb_state_old%allocatememory()
!~   call dmb_state_new%allocatememory()
  call dmb_state_iter%allocatememory()
  call dmb_state_prev%allocatememory()
  call cdt_state_old%allocatememory()
!~   call cdt_state_new%allocatememory()
  call cdt_state_i%allocatememory()
  call cdt_state_ii%allocatememory()
  call cdt_state_iter%allocatememory()
  call cdt_state_prev%allocatememory()
  call ddc_state_old%allocatememory()
!~   call ddc_state_new%allocatememory()
  call ddc_state_iter%allocatememory()
  call ddc_state_prev%allocatememory()

  call generate_regions() !! initialize region_mask


!--------------------
!  0.A.1 ----> DMB Setup - set up DMB properties object
!--------------------
  !$OMP PARALLEL DO
  do iregion = 1,nregion

  call DMB_get_slip_system(crystalstruct, slip_M, slip_S, slip_V, line_T, rotM, Aforest)

!  Nslip = 12  !Defined in GlobalParams module

  all_props_dmb(iregion)%rhobar0     = rhobar0
  all_props_dmb(iregion)%T_ref       = temperature0

  all_props_dmb(iregion)%C11         = C11
  all_props_dmb(iregion)%C12         = C12
  all_props_dmb(iregion)%C44         = C44

  all_props_dmb(iregion)%bulk_c1     = bulk_c1
  all_props_dmb(iregion)%bulk_c2     = bulk_c2

  all_props_dmb(iregion)%c1          = c1
  all_props_dmb(iregion)%mu          = mu
  all_props_dmb(iregion)%burger      = burger
  all_props_dmb(iregion)%B0          = B0
  all_props_dmb(iregion)%Lbar        = Lbar
  all_props_dmb(iregion)%g0          = g0
  all_props_dmb(iregion)%omega0      = omega0
  all_props_dmb(iregion)%wave_vel    = wave_vel
  all_props_dmb(iregion)%drag_flag   = drag_flag
  all_props_dmb(iregion)%drag_Trho_flag = drag_Trho_flag
  all_props_dmb(iregion)%boltz       = boltz
  all_props_dmb(iregion)%p           = p
  all_props_dmb(iregion)%q           = q
  all_props_dmb(iregion)%euler_angle = euler_angle(iregion,:)
  all_props_dmb(iregion)%tau0        = tau0

  all_props_dmb(iregion)%slip_S      = slip_S
  all_props_dmb(iregion)%slip_M      = slip_M
  all_props_dmb(iregion)%slip_V      = slip_V
  all_props_dmb(iregion)%rotM        = rotM
  all_props_dmb(iregion)%gtol        = gtol

  all_props_dmb(iregion)%internal_strain_flag = internal_strain_flag

  all_props_dmb(iregion)%eos%eos_flag = eos_flag
  all_props_dmb(iregion)%eos%v_star   = eos_v_star
  all_props_dmb(iregion)%eos%B_star   = eos_B_star
  all_props_dmb(iregion)%eos%B1       = eos_B1
  all_props_dmb(iregion)%eos%rho0     = eos_rho0
  all_props_dmb(iregion)%eos%v0       = eos_v0
  all_props_dmb(iregion)%eos%Td0      = eos_Td0
  all_props_dmb(iregion)%eos%a        = eos_a
  all_props_dmb(iregion)%eos%b        = eos_b
  all_props_dmb(iregion)%eos%gamma    = eos_gamma
  all_props_dmb(iregion)%eos%Gamma0   = eos_Gamma0
  all_props_dmb(iregion)%eos%kappa    = eos_kappa
  all_props_dmb(iregion)%eos%R_M      = eos_R_M
  all_props_dmb(iregion)%tau_back_nokappa = 0.d0 !!initialize zaiser backstress cache to zero

  ! precompute stiffness tensor, rotation, and slip systems
  call DMB_Compute_Rotation(euler_angle(iregion,:),all_props_dmb(iregion)%Qrot)
  !! call DMB_Form_Elasticity((/C11,C12,C44/), all_props_dmb(iregion)%Qrot,all_props_dmb(iregion)%Cg)

  all_props_dmb(iregion)%elastic_constants = 0.d0
  all_props_dmb(iregion)%elastic_constants(1:3,1) = (/ C11, C12, C44/)
  all_props_dmb(iregion)%elastic_constants(1:3,2) = (/CT11,CT12,CT44/)
  all_props_dmb(iregion)%elastic_constants(1:3,3) = (/CP11,CP12,CP44/)

  call DMB_form_elasticity(all_props_dmb(iregion)%Cg, all_props_dmb(iregion)%Qrot,  &
                           all_props_dmb(iregion)%elastic_constants, 0.d0, 0.d0)

  call DMB_Invert_Tensor43(all_props_dmb(iregion)%Cg,all_props_dmb(iregion)%Cinv)

  ! get slip systems into global coordinates
  do ich=1,Nchar
    do i=1,Nslip
      all_props_dmb(iregion)%slip_Vg(:,i,ich)   = MATMUL(all_props_dmb(iregion)%Qrot, slip_V(:,i,ich))
      all_props_dmb(iregion)%slip_Sg(:,i)   = MATMUL(all_props_dmb(iregion)%Qrot, slip_S(:,i))
      all_props_dmb(iregion)%slip_Mg(:,i)   = MATMUL(all_props_dmb(iregion)%Qrot, slip_M(:,i))
      all_props_dmb(iregion)%schmid0(i,:,:) = all_props_dmb(iregion)%slip_Sg(:,i) .otimes. all_props_dmb(iregion)%slip_Mg(:,i)
    end do  !i
  end do  !ich

  call DMB_Init_EOS(all_props_dmb(iregion)%eos%Dhat)

  call DMB_Compute_EOS(all_props_dmb(iregion)%eos, 1d0, temperature0, &
                       PRESSURE, FREE_ENERGY, ENTROPY, THERMO_COEFFS)

  all_props_dmb(iregion)%eos%pressure0 = pressure(1)
  if ( abs(thermo_coeffs(1)) < 1.d-15 ) then
     thermo_coeffs(1) = specific_heat ! use input value if incomplete EOS params were given
  end if
  all_props_dmb(iregion)%specific_heat = thermo_coeffs(1)
  if (abs(debug_cycle)>0) then
    print*,"specific heat from input file:",specific_heat,"... from EOS:",thermo_coeffs(1)
  end if



  ! 0.A.2 ---->  CDT Setup - set up CDT properties object

  all_props_cdt(iregion)%update_method = CDT_update
  all_props_cdt(iregion)%scheme        = CDTscheme
  all_props_cdt(iregion)%limiter       = CDTlimiter

  all_props_cdt(iregion)%burger    = burger
  all_props_cdt(iregion)%K         = cdt_K
  all_props_cdt(iregion)%A         = cdt_A
  all_props_cdt(iregion)%Y_e       = cdt_Y_e
  all_props_cdt(iregion)%Aint      = cdt_Aint
  all_props_cdt(iregion)%Aforest   = Aforest
  all_props_cdt(iregion)%rhoDot_n0 = cdt_rhoDot_n0
  all_props_cdt(iregion)%G_n0      = cdt_G_n0
  all_props_cdt(iregion)%tau_n0    = cdt_tau_n0
  all_props_cdt(iregion)%pn        = cdt_pn
  all_props_cdt(iregion)%qn        = cdt_qn
  all_props_cdt(iregion)%boltz     = boltz

  all_props_cdt(iregion)%no_flux   = no_flux

  ! 0.A.3 ---->  DDC props setup
  all_props_ddc(iregion)%Lrho             = Lrho
  !! FIXME: why is this C44? Why not mu for edge? TODO: make this character dependent?
  all_props_ddc(iregion)%mu               = C44 ! only used in gradient backstress (not in zaiser)
  !!
  all_props_ddc(iregion)%burger           = burger
  all_props_ddc(iregion)%slip_S           = slip_S
  all_props_ddc(iregion)%euler_angle      = euler_angle(iregion,:)
  all_props_ddc(iregion)%backstress_model = backstress_model

  end do !! iregion
  !$OMP END PARALLEL DO

  props_cdt = all_props_cdt(1)

  ! 0.A.4 ---->  Global flags for switching, debuging, etc.
!  globals = Properties()
!  globals.stop = False

  ! --------------------- #
  ! 0.B DMB Solution Data # <------------------------------------------------------------------------
  ! --------------------- #

  ! 0.B.1 ---->  DMB Setup - node coordinates, initial velocities
  call LinSpace(0.d0, L0, Nnode, X_ref)
  x_pos(:,:) = 0.d0
  x_pos(:,1) = X_ref(:)
  x_vel(:,:) = 0.d0
  if (simulation_type=='shear') then
    call LinSpace(-0.5d0*impact_velocity, 0.5d0*impact_velocity, Nnode, x_vel(:,2))
  end if

  ! 0.B.2 ---->  DMB Setup - lumped system mass matrices
  lumped_mass     = 0.d0

  cell_sizes(:) = X_ref(2:)-X_ref(:Nnode-1)
  matrix1(1,:) = (/ 1.d0/3.d0, 1.d0/6.d0 /)
  matrix1(2,:) = (/ 1.d0/6.d0, 1.d0/3.d0 /)

  do iregion=1,nregion
    do iel=region_mask(iregion,1), region_mask(iregion,2)

      lumped_mass(iel:iel+1) = lumped_mass(iel:iel+1) +                          &
          0.5d0*cell_sizes(iel)*all_props_dmb(iregion)%rhobar0

!!!      uninitialized and never used, but eating huge amount of memory if Nnode is large:
!~       consistent_mass(iel:iel+1,iel:iel+1) = consistent_mass(iel:iel+1,iel:iel+1) + &
!~           all_props_dmb(iregion)%rhobar0*cell_sizes(iel)*matrix1

    end do  !iel
  end do
  ! 0.B.3 ---->  DMB Setup - check wave speed, courant
  cs = sqrt(C11 / rhobar0)
  dX = minval(cell_sizes)
  dt_max = 0.5d0*dX / cs
  t_check = L0 / cs
  if (dt0 > dt_max) then
    print*,"WARNING: time steps are too large; dt0=",dt0,", recommend dt0<", dt_max
  end if
  
  eye = Identity(3)

  call DMB_State_Init((/0.d0,0.d0,0.d0 /), X_ref, x_pos, x_vel, xddot, eye, eye, temperature0, &
       dmb_state_old)

  ! 0.B.3 ---->  DMB Setup - set up boundary conditions
  
  !initialize the values not set below
  bc_force%bc_active(1:2)   = .False.
  bc_force%dof_active       = .False.
  bc_xdot%bc_active(1:2)    = .False.
  bc_xdot%dof_active(:,:)   = .False.
  bc_force%value            = 0.d0
  bc_xdot%value             = 0.d0
  
  if (simulation_type=='shear') then
    if (abs(impact_force)<1.d-99 .or. abs(impact_velocity)>0.d0) then
      bc_xdot%bc_active(1:2)    = (/ .True., .True. /)
      bc_xdot%dof_active(1:2,1:3) = .True.
    else
      bc_force%bc_active(1:2)    = (/ .True., .True. /)
      bc_force%dof_active(1:2,1:3) = .True.
    end if
    bc_xdot%value(1,2)        = -0.5d0 * impact_velocity
    bc_xdot%value(2,2)        =  0.5d0 * impact_velocity
    bc_force%value(1,2)        = -impact_force
    bc_force%value(2,2)        =  impact_force
  elseif (simulation_type=='impact') then
    if (abs(impact_force)<1.d-99 .or. abs(impact_velocity)>0.d0) then
      bc_xdot%bc_active(1:2)    = (/ .True., .False. /)
      bc_xdot%dof_active(1,1:2) = .True.
    else
      bc_force%bc_active(1:2)    = (/ .True., .False. /)
      bc_force%dof_active(1,1:2) = .True.
    end if
    bc_xdot%value(1,1:2)      = 0.5d0 * impact_velocity *     &
                            (/ cos(impact_angle*pi/180.d0), &
                               sin(impact_angle*pi/180.d0)   /) ! prescribed velocity at impact surface
    bc_force%value(1,1:2)     = impact_force *     &
                            (/ cos(impact_angle*pi/180.d0), &
                               sin(impact_angle*pi/180.d0)   /)
  else
    call Fatal('Unknown simulation_type: '//simulation_type)
  end if
  
  ! --------------------- #
  ! 0.C CDT Solution Data # <------------------------------------------------------------------------
  ! --------------------- #

  ! 0.C.1 ---->  initialize CDT state
  call CDT_State_Init(cell_sizes, cdt_state_old)
  ! distribute total dislacation density rho0 (input parameter) equally among all characters
  ! Nchar denotes number of characters for "pos", so divide by 2 to distribute further between pos and neg
  cdt_state_old%varrho_field(:,:,:,:) = 0.5d0*rho0/(Nslip*Nchar)
  cdt_state_old%varrho_node(:,:,:,:)  = 0.5d0*rho0/(Nslip*Nchar)

  !debug_decouple(:,:,:) = 0.5d0*rho0/Nchar

  ! 0.C.3 ---->  Dislocation flux BCs

  flux_bc = rnul ! set default to "do nothing" conditions for all interfaces

  if (CDTintFlux == 'zero') then
    flux_bc = 0.d0
  elseif (CDTintFlux == 'end') THEN
    flux_bc(1)         = 0.d0
    flux_bc(Nregion+1) = 0.d0
  end if


  ! --------------------- #
  ! 0.D DDC setup Data    # <------------------------------------------------------------------------
  ! --------------------- #

  ! 0.C.2 ---->  initialize DDC backstress
  call DDC_State_Init(cell_sizes,ddc_state_old)
  ddc_state_old%tau_b = 0.d0

  !small = 1.e-6  !Defined in GlobalParams module

  ! initialize counters
  iinc     = 0
  iout     = 0
  t_output = dt_output

!--------------------------------------------------------------------------------------------------------------------#
!                                           Initialize Output Files                                                  #
!--------------------------------------------------------------------------------------------------------------------#

  call Init_All_Field_Output_Files(jobname)
  call Write_State(iinc, dmb_state_old, cdt_state_old, ddc_state_old)
  
!~   call Init_Field_File(jobname, 99    , "start", "0"    )
!~   write(99  ,*) "time"
!~   close(unit=99)

!--------------------------------------------------------------------------------------------------------------------#
!                                           Main Time Incrementation Loop                                            #
!--------------------------------------------------------------------------------------------------------------------#
  stop_flag = .False.
  do while (.not.stop_flag)

    !! write(*,'(a)') "Start of cycle"
    ! update increment counter
    iinc = iinc + 1
    ncycle = iinc
    ! compute time step size and adjust as necessary
    if (dmb_state_old%time(1)+dt0 < 1.1d0*impact_acc_time) then
      dt = 0.2d0*dt0 ! choose smaller time steps in the beginning
    else
      dt = dt0
    end if
    if ( (t_stop - dmb_state_old%time(1)) <  (1.001d0*dt) ) then
      dt = t_stop - dmb_state_old%time(1)
      stop_flag = .True.
    end if

    ! initialize iterations over sub-problems
    kiter          = 1
    konverged      = .false.
    ddc_state_iter = ddc_state_old
    dmb_state_iter = dmb_state_old
    cdt_state_iter = cdt_state_old

    ! iterate over sub-problems
    do while ( (kiter<=MaxIter).and.(.not.konverged) )

      if (CDTintegrator=="FEuler") then

        ! perform DDC update
        call DDC_Update(ddc_state_old, dmb_state_iter, cdt_state_iter, all_props_ddc, all_props_dmb, &
             ddc_state_iter)

        ! perform DMB update
        if (dmb_state_old%time(1)+dt >= t_vel) then
          bc_xdot%bc_active(:) = .False. ! remove velocity BC at t_vel to simulate release wave reflected off flyer
        end if

        call DMB_Update_State(dmb_state_old, dt, lumped_mass, bc_xdot, bc_force, all_props_dmb, cdt_state_iter, ddc_state_iter, &
                             dmb_state_iter, impact_acc_time)

        ! perform CDT update
        call CDT_Update_State(cdt_state_old, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_iter)

      elseif (CDTintegrator=="RK2") then

        ! perform DDC update
        call DDC_Update(ddc_state_old, dmb_state_iter, cdt_state_iter, all_props_ddc, all_props_dmb, &
             ddc_state_iter)

        ! perform DMB update
        if (dmb_state_old%time(1)+dt >= t_vel) then
          bc_xdot%bc_active(:) = .False. ! remove velocity BC at t_vel to simulate release wave reflected off flyer
        end if

        call DMB_Update_State(dmb_state_old, dt, lumped_mass, bc_xdot, bc_force, all_props_dmb, cdt_state_iter, ddc_state_iter, &
                             dmb_state_iter, impact_acc_time)

        ! perform CDT update
        call CDT_Update_State(cdt_state_old, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_i)
        call CDT_Update_State(cdt_state_i, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_ii)
        call CDT_Combine_State(cdt_state_old, cdt_state_ii, 0.5d0, 0.5d0, &
                              cdt_state_iter)

      else

        call Fatal('Unknown time integrator: '//CDTintegrator)

      end if

      ! check for convergence
      if (kiter==1) then
        Eddc(1) = sqrt( sum( (ddc_state_iter%epsilon_i-ddc_state_old%epsilon_i)**2 ) )
        Edmb(1) = sqrt( sum( (dmb_state_iter%cauchy-dmb_state_old%cauchy)**2 ) )
        Ecdt(1) = sqrt( sum( (cdt_state_iter%varrho_field-cdt_state_old%varrho_field)**2 ) )
        if (( Eddc(1)<tolerance ).and.( Edmb(1)<tolerance ).and.( Ecdt(1)<tolerance )) then
          konverged = .true.
        else
          dmb_state_prev = dmb_state_iter
          cdt_state_prev = cdt_state_iter
          ddc_state_prev = ddc_state_iter
          kiter = kiter + 1
        end if
      else
        Eddc(1) = sqrt( sum( (ddc_state_iter%epsilon_i-ddc_state_old%epsilon_i)**2 ) )
        Edmb(1) = sqrt( sum( (dmb_state_iter%cauchy-dmb_state_old%cauchy)**2 ) )
        Ecdt(1) = sqrt( sum( (cdt_state_iter%varrho_field-cdt_state_old%varrho_field)**2 ) )
        Eddc(2) = sqrt( sum( (ddc_state_iter%epsilon_i-ddc_state_prev%epsilon_i)**2 ) )
        Edmb(2) = sqrt( sum( (dmb_state_iter%cauchy-dmb_state_prev%cauchy)**2 ) )
        Ecdt(2) = sqrt( sum( (cdt_state_iter%varrho_field-cdt_state_prev%varrho_field)**2 ) )

        if (Eddc(1)<tolerance) then
          konverged = .true.
        elseif ( Eddc(2)/Eddc(1)<DDCtoler ) then
          konverged = .true.
        else
          konverged = .false.
          dmb_state_prev = dmb_state_iter
          cdt_state_prev = cdt_state_iter
          ddc_state_prev = ddc_state_iter
          kiter = kiter + 1
          cycle
        end if

        if (Edmb(1)<tolerance) then
          konverged = .true.
        elseif ( Edmb(2)/Edmb(1)<DMBtoler ) then
          konverged = .true.
        else
          konverged = .false.
          dmb_state_prev = dmb_state_iter
          cdt_state_prev = cdt_state_iter
          ddc_state_prev = ddc_state_iter
          kiter = kiter + 1
          cycle
        end if

        if (Ecdt(1)<tolerance) then
          konverged = .true.
        elseif ( Ecdt(2)/Ecdt(1)<CDTtoler ) then
          konverged = .true.
        else
          konverged = .false.
          dmb_state_prev = dmb_state_iter
          cdt_state_prev = cdt_state_iter
          ddc_state_prev = ddc_state_iter
          kiter = kiter + 1
        end if

      end if

    end do ! iterations

!!! Debug output
!   if ( konverged ) then
!     write(*,'(/,A,I,A)') 'Converged after ',kiter,' iterations'
!     write(*,'(A,3F,/)') 'Rel. E: ',Eddc(2)/Eddc(1), Edmb(2)/Edmb(1), Ecdt(2)/Ecdt(1)
!   else
!     write(*,'(/,A,3F)') 'Increment : ',Eddc(1), Edmb(1), Ecdt(1)
!     write(*,'(A,3F)')   'Correction: ',Eddc(2), Edmb(2), Ecdt(2)
!     write(*,'(A,3F,/)') 'Relative  : ',Eddc(2)/Eddc(1), Edmb(2)/Edmb(1), Ecdt(2)/Ecdt(1)
!   end if

!~     ddc_state_new = ddc_state_iter
!~     dmb_state_new = dmb_state_iter
!~     cdt_state_new = cdt_state_iter

    if (debug_cycle==ncycle) then
      do ich=1,Nchar
        print*,"character index",ich
        call Print_DDC_State(ddc_state_iter,ich)
        call Print_DMB_State(dmb_state_iter,ich)
        call Print_CDT_State(cdt_state_iter,ich)
      end do
    end if

!--------------------------------------------------------------------------------------------------------------------#
!                     Write Output as needed (replace with NLBeam --> paraView stuff)                                #
!--------------------------------------------------------------------------------------------------------------------#
    if ((t_output - dmb_state_iter%time(1) <= 0.5d0*dt0) .or. (stop_flag) ) then
      t_output = t_output + dt_output
      iout     = iout + 1

      write(*,110) props_cdt%scheme(1:len_trim(props_cdt%scheme)),dmb_state_iter%time(1),iinc
110   format(a,"  Time = ",es10.3,"  Cycle = ",i10)

      call Write_State(iinc, dmb_state_iter, cdt_state_iter, ddc_state_iter)
    end if

    call write_history(ipOut, dmb_state_iter) !, cdt_state_iter, ddc_state_iter)

!--------------------------------------------------------------------------------------------------------------------#
!                     Advance 'old' states for next time step                                                        #
!--------------------------------------------------------------------------------------------------------------------#

    !! write(*,'(a)') "End of cycle"
    dmb_state_old = dmb_state_iter
    cdt_state_old = cdt_state_iter
    ddc_state_old = ddc_state_iter
!    call Fatal('Done.')

  end do !while .not.stop_flag

!--------------------------------------------------------------------------------------------------------------------#
!                     Finalize Output Files                                                                          #
!--------------------------------------------------------------------------------------------------------------------#
  call Close_Field_Files()
!for var_name in output_files.keys():
!  output_files[var_name].close()
!  call Fatal('Done.')
  call cpu_time(finish_time)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
  print*,"------------------------------------------------------------"
  print*,"time: ",int(finish_time-start_time), "s"

END PROGRAM CaseA

