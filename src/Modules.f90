! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
MODULE GlobalParams
  IMPLICIT NONE
  public
  CHARACTER(32), PARAMETER :: prog_version="2026.04.29"
  INTEGER :: Nel   = 1500
  INTEGER :: Nnode = 1501
  INTEGER :: Nslip = 12
  INTEGER :: Nchar = 2 !< number of dislocation characters (2=edge+screw, 1=only edge)
!~   INTEGER, private :: ith_character !! temporary integer variable used only for the array constructor below
  INTEGER, DIMENSION(:), allocatable :: characters != [(ith_character,ith_character=1,Nchar)]
  INTEGER, PARAMETER :: Nsdv  = 300
  REAL(KIND=8), PARAMETER  :: small = 1.0d-6
  REAL(KIND=8), PARAMETER  :: rnul = huge(0.d0), rzero = 2.d0*tiny(0.d0)
  INTEGER            :: current_el
  INTEGER            :: debug_cycle   = 0
  INTEGER            :: debug_element = 0
  INTEGER            :: ncycle        = 0
  REAL(KIND=8), PARAMETER              :: pi = 4.d0*atan(1.d0)
  integer :: Nregion = 1
  integer, allocatable :: region_mask(:,:) !< call generate_regions() in the Driver to populate this global variable
  
  CONTAINS
    SUBROUTINE generate_regions()
      integer :: i, Nskip, ith_character
      Nskip = Nel/Nregion
      ! check consistency of Nel and Nregion:
      if (Nskip*Nregion/=Nel) then
        stop "Nel is not a multiple of Nregion."
      end if
      allocate(region_mask(Nregion,2))
      do i=1,Nregion
        region_mask(i,1) = 1 + Nskip*(i-1)
        region_mask(i,2) = Nskip*i
      end do
      
      allocate(characters(Nchar))
      characters = [(ith_character,ith_character=1,Nchar)]
      
      RETURN
    END SUBROUTINE generate_regions

END MODULE GlobalParams

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

MODULE CDT
  USE GlobalParams
  private
!-----------------------------------------------------------------------
  TYPE, public :: Type_CDT_PROPS
    CHARACTER(32)  :: update_method
    CHARACTER(32)  :: scheme
    CHARACTER(32)  :: limiter
    REAL(KIND=8)         :: burger
    REAL(KIND=8)         :: K
    REAL(KIND=8)         :: A
    REAL(KIND=8)         :: Y_e
    REAL(KIND=8)         :: Aint
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: Aforest !< generalization of Aint, used if input variable Aint<0
    REAL(KIND=8)         :: rhoDot_n0
    REAL(KIND=8)         :: G_n0
    REAL(KIND=8)         :: tau_n0
    REAL(KIND=8)         :: pn
    REAL(KIND=8)         :: qn
    REAL(KIND=8)         :: boltz
    LOGICAL        :: no_flux
    contains
      procedure :: allocatememory => CDT_PROPS_allocatememory
  END TYPE Type_CDT_PROPS

!-----------------------------------------------------------------------
  TYPE, public :: Type_CDT_STATE
    REAL(KIND=8), DIMENSION(:), allocatable          :: cell_sizes
    REAL(KIND=8), DIMENSION(:,:,:,:), allocatable :: varrho_field, varrho_node
    contains
      procedure :: allocatememory => CDT_STATE_allocatememory
  END TYPE Type_CDT_STATE
  
  CONTAINS
    SUBROUTINE CDT_PROPS_allocatememory(some_type)
      class(Type_CDT_PROPS), intent(inout) :: some_type
      allocate(some_type%Aforest(Nel,Nslip,Nchar))
    END SUBROUTINE CDT_PROPS_allocatememory
    
    SUBROUTINE CDT_STATE_allocatememory(some_type)
      class(Type_CDT_STATE), intent(inout) :: some_type
      allocate(some_type%cell_sizes(Nel))
      allocate(some_type%varrho_field(Nel,Nslip,Nchar,2), some_type%varrho_node(Nnode,Nslip,Nchar,2))
    END SUBROUTINE CDT_STATE_allocatememory

END MODULE CDT

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

MODULE DDC
  USE GlobalParams
  private
!-----------------------------------------------------------------------
  TYPE, public :: Type_DDC_PROPS
    REAL(KIND=8)                       :: Lrho
    REAL(KIND=8)                       :: mu
    REAL(KIND=8)                       :: burger
    REAL(KIND=8), DIMENSION(:,:), allocatable :: slip_S
    REAL(KIND=8), DIMENSION(3)         :: euler_angle
    CHARACTER(32)                :: backstress_model
    contains
      procedure :: allocatememory => DDC_PROPS_allocatememory
  END TYPE Type_DDC_PROPS

!-----------------------------------------------------------------------
  TYPE, public :: Type_DDC_STATE
    REAL(KIND=8), DIMENSION(:), allocatable     :: cell_sizes
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: tau_i
    REAL(KIND=8), DIMENSION(:,:), allocatable   :: tau_b
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: stress_i, epsilon_i
    contains
      procedure :: allocatememory => DDC_STATE_allocatememory
  END TYPE Type_DDC_STATE
  
  CONTAINS
    SUBROUTINE DDC_PROPS_allocatememory(some_type)
      class(Type_DDC_PROPS), intent(inout) :: some_type
      allocate(some_type%slip_S(3,Nslip))
    END SUBROUTINE DDC_PROPS_allocatememory
    
    SUBROUTINE DDC_STATE_allocatememory(some_DDC_STATE)
      class(Type_DDC_STATE), intent(inout) :: some_DDC_STATE
      allocate(some_DDC_STATE%cell_sizes(Nel))
      allocate(some_DDC_STATE%tau_i(Nel,Nslip,Nchar), some_DDC_STATE%tau_b(Nel,Nslip))
      allocate(some_DDC_STATE%stress_i(Nel,3,3), some_DDC_STATE%epsilon_i(Nel,3,3))
    END SUBROUTINE DDC_STATE_allocatememory

END MODULE DDC

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

MODULE DMB
  USE GlobalParams
  private
!-----------------------------------------------------------------------
  TYPE, public :: Type_DMB_EOS
    LOGICAL                   :: eos_flag
    REAL(KIND=8)                    :: v_star
    REAL(KIND=8)                    :: B_star
    REAL(KIND=8)                    :: B1
    REAL(KIND=8)                    :: v0
    REAL(KIND=8)                    :: rho0
    REAL(KIND=8)                    :: Td0
    REAL(KIND=8)                    :: a
    REAL(KIND=8)                    :: b
    REAL(KIND=8)                    :: gamma
    REAL(KIND=8)                    :: Gamma0
    REAL(KIND=8)                    :: kappa
    REAL(KIND=8)                    :: R_M
    REAL(KIND=8), DIMENSION(:,:), allocatable :: Dhat
    REAL(KIND=8)                    :: pressure0
    contains
      procedure :: allocatememory => DMB_EOS_allocatememory
  END TYPE Type_DMB_EOS

!-----------------------------------------------------------------------
  TYPE, public :: Type_DMB_PROPS
    REAL(KIND=8)                       :: rhobar0
    REAL(KIND=8)                       :: C11
    REAL(KIND=8)                       :: C12
    REAL(KIND=8)                       :: C44
    REAL(KIND=8)                       :: bulk_c1
    REAL(KIND=8)                       :: bulk_c2
    REAL(KIND=8)                       :: c1
    REAL(KIND=8)                       :: mu
    REAL(KIND=8)                       :: burger
    REAL(KIND=8), DIMENSION(:), allocatable :: B0, wave_vel
    REAL(KIND=8)                       :: Lbar
    REAL(KIND=8)                       :: g0
    REAL(KIND=8)                       :: omega0
    character(32)                      :: drag_flag
    REAL(KIND=8)                       :: boltz
    REAL(KIND=8)                       :: p
    REAL(KIND=8)                       :: q
    REAL(KIND=8), DIMENSION(3)         :: euler_angle
    REAL(KIND=8)                       :: tau0
    REAL(KIND=8), DIMENSION(:,:,:,:), allocatable :: rotM
    REAL(KIND=8), DIMENSION(:,:,:), allocatable   :: slip_V, slip_Vg, line_T
    REAL(KIND=8), DIMENSION(:,:), allocatable   :: slip_M, slip_S, slip_Mg, slip_Sg
    REAL(KIND=8)                       :: gtol
    LOGICAL                      :: internal_strain_flag, drag_Trho_flag
    REAL(KIND=8), dimension(3,3)       :: elastic_constants
    REAL(KIND=8), DIMENSION(3,3)       :: Qrot
    REAL(KIND=8), DIMENSION(3,3,3,3)   :: Cg
    REAL(KIND=8), DIMENSION(3,3,3,3)   :: Cinv
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: schmid0
    TYPE(Type_DMB_EOS)           :: eos
    REAL(KIND=8)                       :: specific_heat
    REAL(KIND=8)                       :: T_ref
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: tau_back_nokappa ! at vel=0
    contains
      procedure :: allocatememory => DMB_PROPS_allocatememory
  END TYPE Type_DMB_PROPS

!-----------------------------------------------------------------------
  TYPE, public :: Type_DMB_STATE
    REAL(KIND=8), DIMENSION(3)         :: time
    REAL(KIND=8), DIMENSION(:), allocatable     :: X_ref
    REAL(KIND=8), DIMENSION(:,:), allocatable   :: x, xdot, xddot
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: F
    REAL(KIND=8), DIMENSION(:,:), allocatable   :: T
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: vel, dis_acc, vel_x, vel_y, vel_z, gamma_acc
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: cauchy, Fp, Fe
    REAL(KIND=8), DIMENSION(:,:), allocatable   :: tau
    INTEGER, DIMENSION(:), allocatable          :: Nsub, i_cutback
    contains
      procedure :: allocatememory => DMB_STATE_allocatememory
  END TYPE Type_DMB_STATE

!-----------------------------------------------------------------------
  TYPE, public :: Type_DMB_FIELD_VARIABLES
    REAL(KIND=8), DIMENSION(:,:,:), allocatable :: varrho_ip
    REAL(KIND=8), DIMENSION(:), allocatable :: tau_b
    REAL(KIND=8), DIMENSION(3,3)     :: epsilon_i
    contains
      procedure :: allocatememory => DMB_FIELD_VARIABLES_allocatememory
  END TYPE Type_DMB_FIELD_VARIABLES

!-----------------------------------------------------------------------
!>  Assume BC applied to left and/or right ends of 1D mesh.  Three
!>    pieces of BC information per end of the mesh: If BC active, which
!>    degrees of freedom are affected and the BC vector itself.
  TYPE, public :: Type_DMB_VECTOR_BC
    LOGICAL, DIMENSION(2)    :: bc_active  !< True if BC applied to this end
    LOGICAL, DIMENSION(2,3)  :: dof_active !< True if DOF is affected
    REAL(KIND=8), DIMENSION(2,3)   :: value      !< BC vector value
  END TYPE Type_DMB_VECTOR_BC
  
  CONTAINS
    SUBROUTINE DMB_EOS_allocatememory(some_type)
      class(Type_DMB_EOS), intent(inout) :: some_type
      allocate(some_type%Dhat(7999,2))
    END SUBROUTINE DMB_EOS_allocatememory
    
    SUBROUTINE DMB_PROPS_allocatememory(some_type)
      class(Type_DMB_PROPS), intent(inout) :: some_type
      allocate(some_type%B0(Nchar), some_type%wave_vel(Nchar))
      allocate(some_type%rotM(3,3,Nslip,Nchar), some_type%slip_Vg(3,Nslip,Nchar))
      allocate(some_type%slip_V(3,Nslip,Nchar), some_type%line_T(3,Nslip,Nchar))
      allocate(some_type%slip_M(3,Nslip), some_type%slip_S(3,Nslip), some_type%slip_Mg(3,Nslip), some_type%slip_Sg(3,Nslip))
      allocate(some_type%schmid0(Nslip,3,3), some_type%tau_back_nokappa(3,Nslip,Nchar))
      allocate(some_type%eos%Dhat(7999,2))
    END SUBROUTINE DMB_PROPS_allocatememory
    
    SUBROUTINE DMB_STATE_allocatememory(some_type)
      class(Type_DMB_STATE), intent(inout) :: some_type
      allocate(some_type%X_ref(Nnode), some_type%F(Nel,3,3), some_type%T(Nel,1))
      allocate(some_type%x(Nnode,3), some_type%xdot(Nnode,3), some_type%xddot(Nnode,3))
      allocate(some_type%vel(Nel,Nslip,Nchar), some_type%dis_acc(Nel,Nslip,Nchar), some_type%gamma_acc(Nel,Nslip,Nchar))
      allocate(some_type%vel_x(Nel,Nslip,Nchar), some_type%vel_y(Nel,Nslip,Nchar), some_type%vel_z(Nel,Nslip,Nchar))
      allocate(some_type%cauchy(Nel,3,3), some_type%Fp(Nel,3,3), some_type%Fe(Nel,3,3))
      allocate(some_type%tau(Nel,Nslip))
      allocate(some_type%Nsub(Nel), some_type%i_cutback(Nel))
    END SUBROUTINE DMB_STATE_allocatememory
    
    SUBROUTINE DMB_FIELD_VARIABLES_allocatememory(some_type)
      class(Type_DMB_FIELD_VARIABLES), intent(inout) :: some_type
      allocate(some_type%varrho_ip(Nslip,Nchar,2), some_type%tau_b(Nslip))
    END SUBROUTINE DMB_FIELD_VARIABLES_allocatememory

END MODULE DMB

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

MODULE OUTPUT
  USE GlobalParams
  public
!-----------------------------------------------------------------------
    INTEGER :: unit_dmb_pos     = 100
    INTEGER :: unit_dmb_vel     = 110
    INTEGER :: unit_dmb_acc     = 120
    INTEGER :: unit_stress      = 130
    INTEGER, DIMENSION(:), allocatable :: unit_vel, unit_vel_x, unit_vel_y, unit_vel_z
    INTEGER, DIMENSION(:), allocatable :: unit_dis_acc , unit_rho_pos, unit_rho_neg
    INTEGER :: unit_const_soln  = 180
    INTEGER :: unit_Fp          = 190
    INTEGER :: unit_Fe          = 200
    INTEGER :: unit_Ei          = 210
    INTEGER :: unit_tau         = 220
    INTEGER :: unit_tau_back    = 230
    INTEGER :: unit_temperature = 240
    INTEGER :: unit_visar       = 250
    integer, DIMENSION(:), allocatable :: unit_gamma_acc
    integer :: unit_stress_th   = 270
    integer :: unit_velocity_th = 280
    
    
    CONTAINS
      SUBROUTINE assign_units()
        allocate(unit_vel(Nchar), unit_vel_x(Nchar), unit_vel_y(Nchar), unit_vel_z(Nchar))
        allocate(unit_dis_acc(Nchar), unit_rho_pos(Nchar), unit_rho_neg(Nchar), unit_gamma_acc(Nchar))
        unit_vel         = 140+characters(:Nchar)
        unit_vel_x       = 150+characters(:Nchar)
        unit_vel_y       = 350+characters(:Nchar)
        unit_vel_z       = 360+characters(:Nchar)
        unit_dis_acc     = 380+characters(:Nchar)
        unit_rho_pos     = 160+characters(:Nchar)
        unit_rho_neg     = 170+characters(:Nchar)
        unit_gamma_acc   = 260+characters(:Nchar)
        RETURN
      END SUBROUTINE assign_units
END MODULE OUTPUT

