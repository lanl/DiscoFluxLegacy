! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
module discoflux_readfiles
USE discoflux_GlobalParams
use discoflux_utilities, only: Fatal
implicit none

!> elastic constants and their change with pressure and temperature (used for linear extrapolation)
type :: elasticconstants
  real(8) :: C11, C12, C44 !! elastic constants (EC) for a cubic crystal at ambient conditions [MPa]
  real(8) :: CP11, CP12, CP44 !! change of (EC) with pressure at ambient conditions (dimensionless)
  real(8) :: CT11, CT12, CT44 !! change of (EC) with  temperature at ambient conditions [MPa / K]
  real(8) :: Zener=1.d0 !! Zener anisotropy ratio (will be calculated)
  real(8) :: cijk(6) = 0.d0 !! cubic TOEC; if provided, dislocation drag coeff. B0 will be caluclated
end type elasticconstants

!> dislocation properties (B0 and/or wave_vel can be provided or calculated on-the-fly)
type :: dislocationprops
  integer :: Nchar = 2 !! number of dislocation character angles to consider (2=just screw and edge, 1=just edge)
  real(8) :: burger = 3.d-7 !! Burgers vector length in mm
  real(8) :: rho0 !! dislocation density in mm^-2
  character(32) :: drag_flag = 'Austin' !! choose fct form for drag coeff.: 'Austin' (default), 'iso', or 'const '
  logical :: drag_Trho_flag = .False. !! set to .True. to make drag coeff. T and rho dependent
  !> choose between 'none' (don't update dislocation props), 'vlim' (limiting vel. only), and 'B0' (B0+vlim)
  character(32) :: update_disl = 'none'
  real(8) :: Lbar !! computed from rho0
  real(8) :: Lrho !! computed from Lbar**2
  real(8), allocatable :: B0(:) !! drag coefficient B0 [MPa-s] at low gliding velocity
  real(8), allocatable :: wave_vel(:) !! limiting velocity in mm/s
  real(8) :: mu = 0.d0 !! if this is set to zero, wave_vel and mu are computed for each disloc. character
  logical :: no_flux = .False.    !! True to turn off dislocation flux
end type dislocationprops

!> initialize some EOS model parameters to zero, if not set in input file, psi_ion =0 and/or psi_electron =0
type :: EOS_data
  logical :: eos_flag = .True.
  real(8) :: perc_gamma = 2.d0/3.d0
  real(8) :: B1 = 0.d0 ! this param. and the two below will be computed if any one of them is not set
  real(8) :: B_star = 0.d0
  real(8) :: v_star = 0.d0
  real(8) :: specific_heat = 0.d0 !! optional, only used if psi_ion=psi_electron=0
  real(8) :: perc_Gamma0 = 0.d0 !! sets psi_el=0
  real(8) :: kappa = 0.d0
  real(8) :: R_M = 0.d0 !! sets psi_ion=0
  real(8) :: Td0 = 100.d0 !! initialize a value for reference Debye temp. that will not crash the code if not set
  real(8) :: a = 0.d0
  real(8) :: b = 0.d0
  real(8) :: rho0
  real(8) :: v0 = 0.d0 !! computed from rho0
end type EOS_data

type :: CDT_data
  real(8) :: toler = 5.0d-2    !! tolerance on "relative error" in disloc. densities
  character(32) :: update = 'FEM'   !! FVM or FEM
  character(32) :: integrator = 'FEuler'  !! FEuler or RK2
  character(32) :: scheme = 'SUPG'   !! With FVM: Central, Godunov, LxF, or LLF; with FEM: Bubnov or SUPG
  character(32) :: limiter = '---'   !! With FVM: zero, Fromm, BW, LxW, MMod, SBee, or MC; with FEM: not used
  character(32) :: intFlux = 'nothing' !! allowed values: 'zero' (=old default), 'end', 'nothing'
  real(8) :: A
  real(8) :: cdt_Aint !! change to any negative value to bypass and use new source term instead
  real(8) :: G_n0
  real(8) :: K !! M = 1./K
  real(8) :: Y_e
  real(8) :: pn
  real(8) :: qn
  real(8) :: rhoDot_n0
  real(8) :: tau_n0
end type CDT_data

type :: DDC_data
  real(8) :: toler = 20.d-2    !! tolerance on "relative error" in epsilon_i
  character(32) :: backstress_model = 'gradient' !! choose between 'gradient', 'zaiser', and 'none' (='bypass')
end type DDC_data

type :: DMB_data
  real(8) :: toler = 0.5d-2    !! tolerance on "relative error" in Cauchy stress
  real(8) :: bulk_c1
  real(8) :: bulk_c2
  real(8) :: omega0
  real(8) :: p
  real(8) :: q
  real(8) :: c1
  real(8) :: g0
  real(8) :: gtol = 1.0d-03
  real(8) :: tau0
end type DMB_data

!>data structure to store the input deck information
type, public :: inputdeck
  logical :: echoinput = .True.
  character(:), allocatable :: inputfilename
  character(10) :: simulation_type = 'impact'  !! choose the driver (currently 'impact' or 'shear')
  character(6) :: crystalstruct = 'fcc' !! choose the crystal structure (currently implemented: 'fcc' and 'bcc')
  integer :: Nel = 1000         !! number of elements to simulate in this 1D simulation
  real(8), allocatable :: euler_angle(:,:) !! euler angles to define orientation of each element in the simulation
  real(8) :: L0 = 1.d0          !! length of 1D simulation box in mm
  integer :: MaxIter = 1        !! Maximum number of iterations (over sub-problems) in a time-step
  integer :: ipOut = -1         !! integration point where history variables written (default: last element)
  real(8) :: dt0 = 1.0d-11      !! default time step in s
  real(8) :: dt_output = 2.0d-9 !! write output every dt_output [s]
  real(8) :: t_stop = 5.0d-7    !! stop the simulation after this time has passed [s]
  logical :: internal_strain_flag = .False.
  real(8) :: impact_force = 0.d0    !! impact_velocity overrides impact_force
  real(8) :: impact_velocity = 0.d0
  real(8) :: impact_angle = 0.d0
  real(8) :: impact_acc_time = rzero
  real(8) :: t_vel
  real(8) :: rhobar0 !! material density in Mg / mm^3
  real(8) :: temperature0    = 300.0d0   !! initial temperature [K]
  type(elasticconstants) :: EC !! elastic constants and their change with pressure and temperature (used for linear extrapolation)
  type(dislocationprops) :: DP !! dislocation properties (B0 and/or wave_vel can be provided or calculated on-the-fly)
  type(EOS_data) :: EOS !! initialize some EOS model parameters to zero, if not set in input file, psi_ion=0 and/or psi_electron=0
  type(CDT_data) :: CDT
  type(DDC_data) :: DDC
  type(DMB_data) :: DMB
end type inputdeck

private
public :: ReadSizes, ReadInput
contains

SUBROUTINE ReadSizes(jobname)
implicit none

CHARACTER(*), INTENT(IN) :: jobname
!------------ local variables:
integer :: ios
character(6) :: crystalstruct
character(32) :: key
character(:), allocatable :: inputfilename
character(256) :: line, values
crystalstruct = 'fcc' ! default value, only used in this function to determine Nslip if not set explicitly

print*,"program version: ",prog_version
inputfilename = "input_parameters."//trim(jobname)//".dat"

open(unit=42, file=inputfilename, action="read", iostat=ios, status='old')
if (ios/=0) then
  close(unit=42)
  open(unit=42, file=jobname, action="read", iostat=ios, status='old')
  if (ios/=0) then
    close(unit=42)
    call Fatal("Error: cannot open input file, tried "//trim(inputfilename)//" and "//trim(jobname))
  end if
  print*, "reading ", jobname
else
  print*, "reading ", inputfilename
end if
do
  read(42,'(a)',iostat=ios) line
  if (ios/=0) exit
  if (line /= " ") then
    read(line,*) key,values
  else
    ! skip empty lines
    key = trim(line)
  end if
  if (key=='Nregion') read(values,*)Nregion
  if (key=='Nel') read(values,*)Nel
  if (key=='Nchar') read(values,*)Nchar
  if (key=='crystalstruct') read(values,*)crystalstruct
end do ! read file
close(unit=42)
Nnode = Nel+1
if(crystalstruct == "fcc") then
  Nslip = 12
elseif(crystalstruct == "bcc") then
  Nslip = 48
else
  call Fatal('Bad structure type: crystal symmetry '//crystalstruct//' unkown/not implemented')
end if
print*,"Nregion ",Nregion," Nel ",Nel," Nchar ",Nchar

END SUBROUTINE ReadSizes


subroutine ReadInput(jobname,indat)

use dislocdyn_dislocations, only: disloc, phonondrag
!$   Use omp_lib
implicit none

CHARACTER(*), INTENT(IN) :: jobname
!------------ local variables:
type(disloc) :: disl
integer :: ios, j, j_eul
character(32) :: key
REAL(KIND=8) :: value1, lat_a, Millerb(3), Millern0(3)
real(8), allocatable :: drag(:,:), vlim(:,:)
character(256) :: line, values
!--------------------------- outputs:
type(inputdeck), intent(out) :: indat

allocate(indat%DP%B0(Nchar), indat%DP%wave_vel(Nchar), indat%euler_angle(Nregion,3))
do j=1,Nregion
  indat%euler_angle(j,:) = (/0d0, 0.d0, 0.d0/)
end do
j_eul = 1
indat%DP%wave_vel(:) = 0.d0
indat%DP%B0 = 3.d-11 ! default value for low velocity drag coefficient is a rough order-of-magnitude estimate

indat%inputfilename = "input_parameters."//trim(jobname)//".dat"
!~ print*, "reading ", indat%inputfilename

open(unit=42, file=indat%inputfilename, action="read", iostat=ios, status='old')
if (ios/=0) then
  close(unit=42)
  open(unit=42, file=jobname, action="read", iostat=ios, status='old')
  if (ios/=0) then
    close(unit=42)
    call Fatal("Error: cannot open input file, tried "//trim(indat%inputfilename)//" and "//trim(jobname))
  end if
end if
if ( (Nchar>999) .or. (Nchar<1) ) then
  call Fatal("Error, Nchar must be between 1 and 999.")
end if
do
  read(42,'(a)',iostat=ios) line
  if (ios/=0) exit
  if (line /= " ") then
    read(line,*) key,values
  else
    ! skip empty lines
    key = trim(line)
  end if
  if (key=='echoinput') read(values,*)indat%echoinput
  if (key=='simulation_type') read(values,*)indat%simulation_type
  if (key=='Nslip') read(values,*)Nslip
  if (key=='crystalstruct') read(values,*)indat%crystalstruct
  if (key=='B0') read(line,*) key,indat%DP%B0(1:Nchar)
  if (key=='cijk') read(line,*) key,indat%EC%cijk(1:6)
  if (key=='C11') read(values,*)indat%EC%C11
  if (key=='C12') read(values,*)indat%EC%C12
  if (key=='C44') read(values,*)indat%EC%C44
  if (key=='CDT_update') read(values,*)indat%CDT%update
  if (key=='CDTintegrator') read(values,*)indat%CDT%integrator
  if (key=='CDTlimiter') read(values,*)indat%CDT%limiter
  if (key=='CDTscheme') read(values,*)indat%CDT%scheme
  if (key=='CDTtoler') read(values,*)indat%CDT%toler
  if (key=='CDTintFlux') read(values,*)indat%CDT%intFlux
  if (key=='backstress_model') read(values,*)indat%DDC%backstress_model
  if (key=='CP11') read(values,*)indat%EC%CP11
  if (key=='CP12') read(values,*)indat%EC%CP12
  if (key=='CP44') read(values,*)indat%EC%CP44
  if (key=='CT11') read(values,*)indat%EC%CT11
  if (key=='CT12') read(values,*)indat%EC%CT12
  if (key=='CT44') read(values,*)indat%EC%CT44
  if (key=='specific_heat') read(values,*)indat%EOS%specific_heat
  
  if (key=='DDCtoler') read(values,*)indat%DDC%toler
  if (key=='DMBtoler') read(values,*)indat%DMB%toler
  if (key=='L0') read(values,*)indat%L0
  if (key=='boltz') read(values,*)Boltz
  if (key=='bulk_c1') read(values,*)indat%DMB%bulk_c1
  if (key=='bulk_c2') read(values,*)indat%DMB%bulk_c2
  if (key=='burger') read(values,*)indat%DP%burger
  if (key=='c1') read(values,*)indat%DMB%c1
  if (key=='cdt_A') read(values,*)indat%cdt%A
  if (key=='cdt_Aint') read(values,*)indat%cdt%cdt_Aint
  if (key=='cdt_G_n0') read(values,*)indat%cdt%G_n0
  if (key=='cdt_K') read(values,*)indat%cdt%K
  if (key=='cdt_Y_e') read(values,*)indat%cdt%Y_e
  if (key=='cdt_pn') read(values,*)indat%cdt%pn
  if (key=='cdt_qn') read(values,*)indat%cdt%qn
  if (key=='cdt_rhoDot_n0') read(values,*)indat%cdt%rhoDot_n0
  if (key=='cdt_tau_n0') read(values,*)indat%cdt%tau_n0
  if (key=='dt0') read(values,*)indat%dt0
  if (key=='dt_output') read(values,*)indat%dt_output
  
  if (key=='g0') read(values,*)indat%DMB%g0
  if (key=='gtol') read(values,*)indat%DMB%gtol
  if (key=='impact_angle') read(values,*)indat%impact_angle
  if (key=='impact_velocity') read(values,*)indat%impact_velocity
  if (key=='impact_force') read(values,*)indat%impact_force
  if (key=='impact_acc_time') read(values,*)indat%impact_acc_time
  if (key=='internal_strain_flag') read(values,*)indat%internal_strain_flag
  if (key=='mu') read(values,*)indat%DP%mu
  if (key=='no_flux') read(values,*)indat%DP%no_flux
  if (key=='omega0') read(values,*)indat%DMB%omega0
  if (key=='p') read(values,*)indat%DMB%p
  if (key=='props_dmb_perc_eos_perc_B1'.or.key=='eos_B1') read(values,*)indat%eos%B1
  if (key=='props_dmb_perc_eos_perc_B_star'.or.key=='eos_B_star') read(values,*)indat%eos%B_star
  if (key=='props_dmb_perc_eos_perc_Gamma0'.or.key=='eos_Gamma0') read(values,*)indat%eos%perc_Gamma0
  if (key=='props_dmb_perc_eos_perc_R_M'.or.key=='eos_R_M') read(values,*)indat%eos%R_M
  if (key=='props_dmb_perc_eos_perc_Td0'.or.key=='eos_Td0') read(values,*)indat%eos%Td0
  if (key=='props_dmb_perc_eos_perc_a'.or.key=='eos_a') read(values,*)indat%eos%a
  if (key=='props_dmb_perc_eos_perc_b'.or.key=='eos_b') read(values,*)indat%eos%b
  if (key=='props_dmb_perc_eos_perc_eos_flag'.or.key=='eos_flag') read(values,*)indat%eos%eos_flag
  if (key=='props_dmb_perc_eos_perc_gamma'.or.key=='eos_gamma') read(values,*)indat%eos%perc_gamma
  if (key=='props_dmb_perc_eos_perc_kappa'.or.key=='eos_kappa') read(values,*)indat%eos%kappa
  if (key=='props_dmb_perc_eos_perc_rho0'.or.key=='eos_rho0') read(values,*)indat%eos%rho0
  if (key=='props_dmb_perc_eos_perc_v0'.or.key=='eos_v0') read(values,*)indat%eos%v0 ! computed from eos_perc_rho0 below
  if (key=='props_dmb_perc_eos_perc_v_star'.or.key=='eos_v_star') read(values,*)indat%eos%v_star
  
  if (key=='q') read(values,*)indat%DMB%q
  if (key=='rho0') read(values,*)indat%DP%rho0
  if (key=='rhobar0') read(values,*)indat%rhobar0
  if (key=='t_stop') read(values,*)indat%t_stop
  if (key=='t_vel') read(values,*)indat%t_vel
  if (key=='tau0') read(values,*)indat%DMB%tau0
  if (key=='temperature0') read(values,*)indat%temperature0
  if (key=='wave_vel') read(line,*) key,indat%DP%wave_vel(1:Nchar) ! computed from mu below, if mu is given explicitly!
  if (key=='all_euler_angles') then
    !! use this keyword to read all euler angles from one input line
    read(line,*)key,(indat%euler_angle(j,1),indat%euler_angle(j,2),indat%euler_angle(j,3), j=1,Nregion)
  end if
  if ((key=='euler_angle') .and. (j_eul<=Nregion)) then
    !! use this keyword to read one set of euler angles, use multiple entries in order to read in for all Nregion
    read(line,*)key,indat%euler_angle(j_eul,1),indat%euler_angle(j_eul,2),indat%euler_angle(j_eul,3)
    j_eul = j_eul + 1
  end if
  if (key=='drag_flag') read(values,*)indat%DP%drag_flag
  if (key=='drag_Trho_flag') read(values,*)indat%DP%drag_Trho_flag
  if (key=='debug_cycle') then
    read(values,*)debug_cycle
    print*,"setting debug_cycle=",debug_cycle
  end if
  
  if (key=='MaxIter') then
    read(values,*)value1
    indat%MaxIter = int(value1)
  elseif (key=='ipOut') then
    read(values,*)value1
    indat%ipOut = int(value1)
  end if
end do ! read file
close(unit=42)

!! if only a subset of euler angles was set explicitly, fill the rest with the last value that was read:
!! (useful for single crystal calcs. where Nregion>1 was set to take advantage of openmp parallelization of loops over Nregion)
if ((j_eul>1) .and. (j_eul<=Nregion)) then
  do j=j_eul,Nregion
    indat%euler_angle(j,:) = indat%euler_angle(j_eul-1,:)
  end do
end if

! compute some derived values:
indat%DP%Lbar = 1.d0 / sqrt(indat%DP%rho0)
indat%DP%Lrho = indat%DP%Lbar**2
if (indat%DP%mu > 0.d0) then
  indat%DP%wave_vel(:) = sqrt(indat%DP%mu / indat%rhobar0) ! mm / s
end if
indat%EC%Zener = 2*indat%EC%C44/(indat%EC%C11-indat%EC%C12)

! calculate character dependent wave_vel, if not provided by the user:
if (abs(indat%DP%wave_vel(1)) < 1.d-15) then
    ! dislocdynlib needs SI units:
    if (indat%crystalstruct=='fcc') then
      lat_a = indat%DP%burger*sqrt(2.d0)*1.d-3
      Millerb=[0.5d0,0.5d0,0.d0]
      Millern0=[-1.d0,1.d0,-1.d0]
    elseif (indat%crystalstruct=='bcc') then
      lat_a = indat%DP%burger*sqrt(3.d0)*0.5d-3
      Millerb=[0.5d0,-0.5d0,0.5d0]
      Millern0=[1.d0,1.d0,0.d0]
      ! TODO: enlarge B0 and wave_vel so that they hold values for each slip system
      ! for bcc we need B0 and wave_vel for each of 110, 112, and 123 (only 110 implemented here for now)
    end if
    disl = disloc(sym=indat%crystalstruct,metal="mat",rho=indat%rhobar0*1.d12,lat_a = [lat_a,0.d0,0.d0])
    disl%cij = [indat%EC%C11*1.d6, indat%EC%C12*1.d6, indat%EC%C44*1.d6]
    disl%cijk = indat%EC%cijk*1.d6
    disl%ntheta = Nchar
    call disl%init(Millerb=Millerb,Millern0=Millern0) ! infers mat%burgers from Millerb
    allocate(vlim(Nchar,3))
    call disl%computevcrit(vlim)
    ! convert wave_vel back from SI units:
    indat%DP%wave_vel = vlim(:,1)*1.d3
    if (abs(indat%EC%cijk(1)) > 1.d0) then
      call phonondrag(drag,disl,[0.d0])
      indat%DP%B0 = drag(:,1)*1.d-9
    end if
end if

! fall-back for bulk modulus Bstar, its pressure derivative B1 and reference volume vstar:
if ( abs(indat%eos%B1*indat%eos%B_star*indat%eos%v_star) < 1.d-15 ) then
   ! print*,"Warning: missing some eos values, computing B_star, B1, and v_star from SOEC and rhobar0"
   indat%eos%v_star = 1.d0/indat%rhobar0 ! reference volume at ambient conditions
   indat%eos%B_star = (indat%EC%C11 + 2.d0*indat%EC%C12)/3.d0 ! bulk modulus from cubic SOEC at v_star
   indat%eos%B1 = (indat%EC%CP11 + 2.d0*indat%EC%CP12)/3.d0 ! pressure derivative of bulk modulus from CPij at v_star
end if

if ( abs(indat%eos%rho0) < 1.d-15) then
   indat%eos%rho0 = indat%rhobar0
end if
if ( abs(indat%eos%v0) < 1.d-15) then
   indat%eos%v0       = 1.d0 / indat%eos%rho0 !((Mg/mm^3)
end if

if ((indat%crystalstruct == "fcc") .and. (Nslip>12)) then
  call Fatal('Nslip must be <= 12 for fcc')
elseif ((indat%crystalstruct == "bcc") .and. (Nslip>48)) then
  call Fatal('Nslip must be <= 48 for fcc')
end if

if (indat%impact_acc_time < rzero) indat%impact_acc_time=rzero

if (indat%ipOut<0) then
  indat%ipOut = Nel + indat%ipOut + 1 ! allow python-like syntax where -1 means the last element of the array
end if

if (indat%echoinput) then
  print*,"running '",trim(indat%simulation_type),"' simulation using:"
  print*,"crystalstruct ",indat%crystalstruct, " Nslip ",Nslip
  print*,"B0",indat%DP%B0
  print*,"C11",indat%EC%C11
  print*,"C12",indat%EC%C12
  print*,"C44",indat%EC%C44
  print*,"Zener ratio",indat%EC%Zener
  print*,"cijk",indat%EC%cijk
  print*,"CDT_update ",indat%CDT%update
  print*,"CDTintegrator ",indat%CDT%integrator
  print*,"CDTlimiter ",indat%CDT%limiter
  print*,"CDTscheme ",indat%CDT%scheme
  print*,"CDTtoler",indat%CDT%toler
  print*,"CDTintFlux ",indat%CDT%intFlux
  print*,"CP11",indat%EC%CP11
  print*,"CP12",indat%EC%CP12
  print*,"CP44",indat%EC%CP44
  print*,"CT11",indat%EC%CT11
  print*,"CT12",indat%EC%CT12
  print*,"CT44",indat%EC%CT44
  print*,"specific_heat",indat%EOS%specific_heat

  print*,"DDCtoler",indat%DDC%toler
  print*,"DMBtoler",indat%DMB%toler
  print*,"L0",indat%L0
  print*,"Lbar",indat%DP%Lbar
  print*,"Lrho=Lbar**2",indat%DP%Lrho
  print*,"MaxIter",indat%MaxIter
  print*,"boltz",Boltz
  print*,"bulk_c1",indat%DMB%bulk_c1
  print*,"bulk_c2",indat%DMB%bulk_c2
  print*,"burger",indat%DP%burger
  print*,"c1",indat%DMB%c1
  print*,"cdt_A",indat%cdt%A
  print*,"cdt_Aint",indat%cdt%cdt_Aint
  print*,"cdt_G_n0",indat%cdt%G_n0
  print*,"cdt_K",indat%cdt%K
  print*,"cdt_Y_e",indat%cdt%Y_e
  print*,"cdt_pn",indat%cdt%pn
  print*,"cdt_qn",indat%cdt%qn
  print*,"cdt_rhoDot_n0",indat%cdt%rhoDot_n0
  print*,"cdt_tau_n0",indat%cdt%tau_n0
  print*,"dt0",indat%dt0
  print*,"dt_output",indat%dt_output

  print*,"g0",indat%DMB%g0
  print*,"gtol",indat%DMB%gtol
  print*,"impact_angle",indat%impact_angle
  print*,"impact_force",indat%impact_force
  print*,"impact_velocity",indat%impact_velocity
  print*,"impact_acc_time",indat%impact_acc_time
  print*,"internal_strain_flag",indat%internal_strain_flag
  print*,"ipOut",indat%ipOut
  print*,"mu",indat%DP%mu
  print*,"no_flux",indat%DP%no_flux
  print*,"omega0",indat%DMB%omega0
  print*,"p",indat%DMB%p
  print*,"eos_B1",indat%eos%B1
  print*,"eos_B_star",indat%eos%B_star
  print*,"eos_Gamma0",indat%eos%perc_Gamma0
  print*,"eos_R_M",indat%eos%R_M
  print*,"eos_Td0",indat%eos%Td0
  print*,"eos_a",indat%eos%a
  print*,"eos_b",indat%eos%b
  print*,"eos_flag",indat%eos%eos_flag
  print*,"eos_gamma",indat%eos%perc_gamma
  print*,"eos_kappa",indat%eos%kappa
  print*,"eos_rho0",indat%eos%rho0
  print*,"eos_v0",indat%eos%v0
  print*,"eos_v_star",indat%eos%v_star
  print*,"q",indat%DMB%q
  print*,"rho0",indat%DP%rho0
  print*,"rhobar0",indat%rhobar0
  print*,"t_stop",indat%t_stop
  print*,"t_vel",indat%t_vel
  print*,"tau0",indat%DMB%tau0
  print*,"temperature0",indat%temperature0
  print*,"wave_vel",indat%DP%wave_vel
  print*,"drag_flag ",indat%DP%drag_flag
  print*,"drag_Trho_flag ",indat%DP%drag_Trho_flag
  print*,"backstress_model ",indat%DDC%backstress_model
  print*,"list of euler_angles:"
  do j_eul=1,Nregion
    print*,j_eul,indat%euler_angle(j_eul,:)
  end do
  print*,"done reading input"
end if

!$    print*, 'OpenMP parallelization enabled:'
!$    print*, 'using ',omp_get_max_threads(),' of ',omp_get_num_procs(),' processors'

END SUBROUTINE ReadInput
end module discoflux_readfiles
