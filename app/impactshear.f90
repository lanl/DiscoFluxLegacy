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
PROGRAM impactshear
  use stdlib_math, only: linspace
  ! modules defined in Modules.f90:
  USE discoflux_GlobalParams
  USE discoflux_DMB, only: Type_DMB_PROPS, Type_DMB_STATE, Type_DMB_VECTOR_BC
  USE discoflux_CDT, only: Type_CDT_PROPS, Type_CDT_STATE
  USE discoflux_DDC, only: Type_DDC_PROPS, Type_DDC_STATE
  ! module defined in Utilities.f90:
  use discoflux_utilities, only: Fatal, Identity, readcmdline, operator(.otimes.)
  ! module defined in DMB.f90:
  use discoflux_DMB_routines, only: DMB_State_Init, DMB_Update_State, DMB_Init_EOS, DMB_Compute_EOS, DMB_Invert_Tensor43, &
                                    DMB_get_slip_system, DMB_form_elasticity, DMB_Compute_Rotation
  ! module defined in CDT.f90:
  use discoflux_CDT_routines, only: CDT_State_Init, CDT_Combine_State, CDT_Update_State
  ! module defined in DDC.f90:
  use discoflux_DDC_routines, only: DDC_State_Init, DDC_Update
  ! module defined in Output.f90:
  use discoflux_OUTPUTFILES, only: Init_All_Field_Output_Files, Close_Field_Files, Write_State, write_history, &
                                   Print_CDT_State, Print_DDC_State, Print_DMB_State
  ! module defined in readinputfile.f90:
  use discoflux_readfiles, only: inputdeck, ReadSizes, ReadInput
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
  type(inputdeck)                :: indat
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: Aforest !! generalization of `Aint`, used if input variable `Aint`<0
  REAL(KIND=8)                :: tolerance=sqrt(rzero)
  CHARACTER(:), allocatable :: jobname
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
  INTEGER                        :: iinc, iout
  REAL(KIND=8)                         :: t_output
  LOGICAL                        :: stop_flag
  INTEGER                        :: i,ich,iregion
  INTEGER                        :: kiter
  LOGICAL                        :: konverged
  REAL(KIND=8), DIMENSION(2)           :: Eddc, Edmb, Ecdt
  integer                        :: start_time, finish_time, countrate

  ! read input
  call readcmdline(jobname)

  call ReadSizes(jobname)
  call ReadInput(jobname, indat)

  call system_clock(start_time,countrate)

  !! allocate memory for arrays that are too large to fit in the stack:
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

  call DMB_get_slip_system(indat%crystalstruct, slip_M, slip_S, slip_V, line_T, rotM, Aforest)

  all_props_dmb(iregion)%rhobar0     = indat%rhobar0
  all_props_dmb(iregion)%T_ref       = indat%temperature0

  all_props_dmb(iregion)%C11         = indat%EC%C11
  all_props_dmb(iregion)%C12         = indat%EC%C12
  all_props_dmb(iregion)%C44         = indat%EC%C44

  all_props_dmb(iregion)%bulk_c1     = indat%DMB%bulk_c1
  all_props_dmb(iregion)%bulk_c2     = indat%DMB%bulk_c2

  all_props_dmb(iregion)%c1          = indat%DMB%c1
  all_props_dmb(iregion)%mu          = indat%DP%mu
  all_props_dmb(iregion)%burger      = indat%DP%burger
  all_props_dmb(iregion)%B0          = indat%DP%B0
  all_props_dmb(iregion)%Lbar        = indat%DP%Lbar
  all_props_dmb(iregion)%g0          = indat%DMB%g0
  all_props_dmb(iregion)%omega0      = indat%DMB%omega0
  all_props_dmb(iregion)%wave_vel    = indat%DP%wave_vel
  all_props_dmb(iregion)%drag_flag   = indat%DP%drag_flag
  all_props_dmb(iregion)%drag_Trho_flag = indat%DP%drag_Trho_flag
  all_props_dmb(iregion)%p           = indat%DMB%p
  all_props_dmb(iregion)%q           = indat%DMB%q
  all_props_dmb(iregion)%euler_angle = indat%euler_angle(iregion,:)
  all_props_dmb(iregion)%tau0        = indat%DMB%tau0

  all_props_dmb(iregion)%slip_S      = slip_S
  all_props_dmb(iregion)%slip_M      = slip_M
  all_props_dmb(iregion)%slip_V      = slip_V
  all_props_dmb(iregion)%rotM        = rotM
  all_props_dmb(iregion)%gtol        = indat%DMB%gtol

  all_props_dmb(iregion)%internal_strain_flag = indat%internal_strain_flag

  all_props_dmb(iregion)%eos%eos_flag = indat%eos%eos_flag
  all_props_dmb(iregion)%eos%v_star   = indat%eos%v_star
  all_props_dmb(iregion)%eos%B_star   = indat%eos%B_star
  all_props_dmb(iregion)%eos%B1       = indat%eos%B1
  all_props_dmb(iregion)%eos%rho0     = indat%eos%rho0
  all_props_dmb(iregion)%eos%v0       = indat%eos%v0
  all_props_dmb(iregion)%eos%Td0      = indat%eos%Td0
  all_props_dmb(iregion)%eos%a        = indat%eos%a
  all_props_dmb(iregion)%eos%b        = indat%eos%b
  all_props_dmb(iregion)%eos%gamma    = indat%eos%perc_gamma
  all_props_dmb(iregion)%eos%Gamma0   = indat%eos%perc_Gamma0
  all_props_dmb(iregion)%eos%kappa    = indat%eos%kappa
  all_props_dmb(iregion)%eos%R_M      = indat%eos%R_M
  all_props_dmb(iregion)%tau_back_nokappa = 0.d0 !!initialize zaiser backstress cache to zero

  ! precompute stiffness tensor, rotation, and slip systems
  call DMB_Compute_Rotation(indat%euler_angle(iregion,:),all_props_dmb(iregion)%Qrot)
  !! call DMB_Form_Elasticity((/indat%EC%C11,indat%EC%C12,indat%EC%C44/), all_props_dmb(iregion)%Qrot,all_props_dmb(iregion)%Cg)

  all_props_dmb(iregion)%elastic_constants = 0.d0
  all_props_dmb(iregion)%elastic_constants(1:3,1) = (/indat%EC%C11,  indat%EC%C12, indat%EC%C44 /)
  all_props_dmb(iregion)%elastic_constants(1:3,2) = (/indat%EC%CT11, indat%EC%CT12,indat%EC%CT44/)
  all_props_dmb(iregion)%elastic_constants(1:3,3) = (/indat%EC%CP11, indat%EC%CP12,indat%EC%CP44/)

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

  call DMB_Compute_EOS(all_props_dmb(iregion)%eos, 1d0, indat%temperature0, &
                       PRESSURE, FREE_ENERGY, ENTROPY, THERMO_COEFFS)

  all_props_dmb(iregion)%eos%pressure0 = pressure(1)
  if ( abs(thermo_coeffs(1)) < 1.d-15 ) then
     thermo_coeffs(1) = indat%EOS%specific_heat ! use input value if incomplete EOS params were given
  end if
  all_props_dmb(iregion)%specific_heat = thermo_coeffs(1)
  if (abs(debug_cycle)>0) then
    print*,"specific heat from input file:",indat%EOS%specific_heat,"... from EOS:",thermo_coeffs(1)
  end if



  ! 0.A.2 ---->  CDT Setup - set up CDT properties object

  all_props_cdt(iregion)%update_method = indat%CDT%update
  all_props_cdt(iregion)%scheme        = indat%CDT%scheme
  all_props_cdt(iregion)%limiter       = indat%CDT%limiter

  all_props_cdt(iregion)%burger    = indat%DP%burger
  all_props_cdt(iregion)%K         = indat%cdt%K
  all_props_cdt(iregion)%A         = indat%cdt%A
  all_props_cdt(iregion)%Y_e       = indat%cdt%Y_e
  all_props_cdt(iregion)%Aint      = indat%cdt%cdt_Aint
  all_props_cdt(iregion)%Aforest   = Aforest
  all_props_cdt(iregion)%rhoDot_n0 = indat%cdt%rhoDot_n0
  all_props_cdt(iregion)%G_n0      = indat%cdt%G_n0
  all_props_cdt(iregion)%tau_n0    = indat%cdt%tau_n0
  all_props_cdt(iregion)%pn        = indat%cdt%pn
  all_props_cdt(iregion)%qn        = indat%cdt%qn

  all_props_cdt(iregion)%no_flux   = indat%DP%no_flux

  ! 0.A.3 ---->  DDC props setup
  all_props_ddc(iregion)%Lrho             = indat%DP%Lrho
  !! FIXME: why is this C44? Why not mu for edge? TODO: make this character dependent?
  all_props_ddc(iregion)%mu               = indat%EC%C44 ! only used in gradient backstress (not in zaiser)
  !!
  all_props_ddc(iregion)%burger           = indat%DP%burger
  all_props_ddc(iregion)%slip_S           = slip_S
  all_props_ddc(iregion)%euler_angle      = indat%euler_angle(iregion,:)
  all_props_ddc(iregion)%backstress_model = indat%DDC%backstress_model

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
  X_ref = linspace(0.d0, indat%L0, Nnode)
  x_pos(:,:) = 0.d0
  x_pos(:,1) = X_ref(:)
  x_vel(:,:) = 0.d0
  xddot(:,:) = 0.d0
  if (indat%simulation_type=='shear') then
    x_vel(:,2) = linspace(-0.5d0*indat%impact_velocity, 0.5d0*indat%impact_velocity, Nnode)
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
  cs = sqrt(indat%EC%C11 / indat%rhobar0)
  dX = minval(cell_sizes)
  dt_max = 0.5d0*dX / cs
  t_check = indat%L0 / cs
  if (indat%dt0 > dt_max) then
    print*,"WARNING: time steps are too large; dt0=",indat%dt0,", recommend dt0<", dt_max
  end if

  eye = Identity(3)

  call DMB_State_Init((/0.d0,0.d0,0.d0 /), X_ref, x_pos, x_vel, xddot, eye, eye, indat%temperature0, &
       dmb_state_old)

  ! 0.B.3 ---->  DMB Setup - set up boundary conditions

  !initialize the values not set below
  bc_force%bc_active(1:2)   = .False.
  bc_force%dof_active       = .False.
  bc_xdot%bc_active(1:2)    = .False.
  bc_xdot%dof_active(:,:)   = .False.
  bc_force%value            = 0.d0
  bc_xdot%value             = 0.d0

  if (indat%simulation_type=='shear') then
    if (abs(indat%impact_force)<1.d-99 .or. abs(indat%impact_velocity)>0.d0) then
      bc_xdot%bc_active(1:2)    = (/ .True., .True. /)
      bc_xdot%dof_active(1:2,1:3) = .True.
    else
      bc_force%bc_active(1:2)    = (/ .True., .True. /)
      bc_force%dof_active(1:2,1:3) = .True.
    end if
    bc_xdot%value(1,2)        = -0.5d0 * indat%impact_velocity
    bc_xdot%value(2,2)        =  0.5d0 * indat%impact_velocity
    bc_force%value(1,2)        = -indat%impact_force
    bc_force%value(2,2)        =  indat%impact_force
  else if (indat%simulation_type=='impact') then
    if (abs(indat%impact_force)<1.d-99 .or. abs(indat%impact_velocity)>0.d0) then
      bc_xdot%bc_active(1:2)    = (/ .True., .False. /)
      bc_xdot%dof_active(1,1:2) = .True.
    else
      bc_force%bc_active(1:2)    = (/ .True., .False. /)
      bc_force%dof_active(1,1:2) = .True.
    end if
    bc_xdot%value(1,1:2)      = 0.5d0 * indat%impact_velocity *     &
                            (/ cos(indat%impact_angle*pi/180.d0), &
                               sin(indat%impact_angle*pi/180.d0)   /) ! prescribed velocity at impact surface
    bc_force%value(1,1:2)     = indat%impact_force *     &
                            (/ cos(indat%impact_angle*pi/180.d0), &
                               sin(indat%impact_angle*pi/180.d0)   /)
  else
    call Fatal('Unknown simulation_type: '//indat%simulation_type)
  end if

  ! --------------------- #
  ! 0.C CDT Solution Data # <------------------------------------------------------------------------
  ! --------------------- #

  ! 0.C.1 ---->  initialize CDT state
  call CDT_State_Init(cell_sizes, cdt_state_old)
  ! distribute total dislacation density rho0 (input parameter) equally among all characters
  ! Nchar denotes number of characters for "pos", so divide by 2 to distribute further between pos and neg
  cdt_state_old%varrho_field(:,:,:,:) = 0.5d0*indat%DP%rho0/(Nslip*Nchar)
  cdt_state_old%varrho_node(:,:,:,:)  = 0.5d0*indat%DP%rho0/(Nslip*Nchar)

  !debug_decouple(:,:,:) = 0.5d0*rho0/Nchar

  ! 0.C.3 ---->  Dislocation flux BCs

  flux_bc = rnul ! set default to "do nothing" conditions for all interfaces

  if (indat%CDT%intFlux == 'zero') then
    flux_bc = 0.d0
  else if (indat%CDT%intFlux == 'end') THEN
    flux_bc(1)         = 0.d0
    flux_bc(Nregion+1) = 0.d0
  end if


  ! --------------------- #
  ! 0.D DDC setup Data    # <------------------------------------------------------------------------
  ! --------------------- #

  ! 0.C.2 ---->  initialize DDC backstress
  call DDC_State_Init(cell_sizes,ddc_state_old)
  ddc_state_old%tau_b = 0.d0
  ddc_state_old%epsilon_i = 0.d0

  !small = 1.e-6  !Defined in GlobalParams module

  ! initialize counters
  iinc     = 0
  iout     = 0
  t_output = indat%dt_output

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
    if (dmb_state_old%time(1)+indat%dt0 < 1.1d0*indat%impact_acc_time) then
      dt = 0.2d0*indat%dt0 ! choose smaller time steps in the beginning
    else
      dt = indat%dt0
    end if
    if ( (indat%t_stop - dmb_state_old%time(1)) <  (1.001d0*dt) ) then
      dt = indat%t_stop - dmb_state_old%time(1)
      stop_flag = .True.
    end if

    ! initialize iterations over sub-problems
    kiter          = 1
    konverged      = .false.
    ddc_state_iter = ddc_state_old
    dmb_state_iter = dmb_state_old
    cdt_state_iter = cdt_state_old

    ! iterate over sub-problems
    do while ( (kiter<=indat%MaxIter).and.(.not.konverged) )

      if (indat%CDT%integrator=="FEuler") then

        ! perform DDC update
        call DDC_Update(ddc_state_old, dmb_state_iter, cdt_state_iter, all_props_ddc, all_props_dmb, &
             ddc_state_iter)

        ! perform DMB update
        if (dmb_state_old%time(1)+dt >= indat%t_vel) then
          bc_xdot%bc_active(:) = .False. ! remove velocity BC at t_vel to simulate release wave reflected off flyer
        end if

        call DMB_Update_State(dmb_state_old, dt, lumped_mass, bc_xdot, bc_force, all_props_dmb, cdt_state_iter, ddc_state_iter, &
                             dmb_state_iter, indat%impact_acc_time)

        ! perform CDT update
        call CDT_Update_State(cdt_state_old, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_iter)

      else if (indat%CDT%integrator=="RK2") then

        ! perform DDC update
        call DDC_Update(ddc_state_old, dmb_state_iter, cdt_state_iter, all_props_ddc, all_props_dmb, &
             ddc_state_iter)

        ! perform DMB update
        if (dmb_state_old%time(1)+dt >= indat%t_vel) then
          bc_xdot%bc_active(:) = .False. ! remove velocity BC at t_vel to simulate release wave reflected off flyer
        end if

        call DMB_Update_State(dmb_state_old, dt, lumped_mass, bc_xdot, bc_force, all_props_dmb, cdt_state_iter, ddc_state_iter, &
                             dmb_state_iter, indat%impact_acc_time)

        ! perform CDT update
        call CDT_Update_State(cdt_state_old, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_i)
        call CDT_Update_State(cdt_state_i, ddc_state_iter, dmb_state_iter,all_props_cdt,dt,flux_bc, &
                              cdt_state_ii)
        call CDT_Combine_State(cdt_state_old, cdt_state_ii, 0.5d0, 0.5d0, &
                              cdt_state_iter)

      else

        call Fatal('Unknown time integrator: '//indat%CDT%integrator)

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
        else if ( Eddc(2)/Eddc(1)<indat%DDC%toler ) then
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
        else if ( Edmb(2)/Edmb(1)<indat%DMB%toler ) then
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
        else if ( Ecdt(2)/Ecdt(1)<indat%CDT%toler ) then
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
    if ((t_output - dmb_state_iter%time(1) <= 0.5d0*indat%dt0) .or. (stop_flag) ) then
      t_output = t_output + indat%dt_output
      iout     = iout + 1

      write(*,110) props_cdt%scheme(1:len_trim(props_cdt%scheme)),dmb_state_iter%time(1),iinc
110   format(a,"  Time = ",es10.3,"  Cycle = ",i10)

      call Write_State(iinc, dmb_state_iter, cdt_state_iter, ddc_state_iter)
    end if

    call write_history(indat%ipOut, dmb_state_iter) !, cdt_state_iter, ddc_state_iter)

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
  call system_clock(finish_time)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
  print*,"------------------------------------------------------------"
  print*,"time: ",int(real(finish_time-start_time)/real(countrate)), "s"

END PROGRAM impactshear

