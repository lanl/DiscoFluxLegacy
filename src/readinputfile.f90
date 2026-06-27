! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
module readfiles
public
contains
SUBROUTINE ReadSizes(jobname)

USE GlobalParams
use utilities, only: Fatal
implicit none

CHARACTER(32), INTENT(IN) :: jobname
!------------ local variables:
integer :: ios
character(6) :: crystalstruct
character(32) :: key
character(54) :: inputfilename
character(256) :: line, values
crystalstruct = 'fcc' ! default value, only used in this function to determine Nslip if not set explicitly

print*,"program version: ",prog_version
inputfilename = "input_parameters."//trim(jobname)//".dat"
print*, "reading ", inputfilename

open(unit=42, file=inputfilename, action="read", iostat=ios, status='old')
if (ios/=0) then
  close(unit=42)
  print*, "failed, now reading ", jobname
  open(unit=42, file=jobname, action="read", iostat=ios, status='old')
  if (ios/=0) then
    close(unit=42)
    call Fatal("Error: cannot open input file, tried "//trim(inputfilename)//" and "//jobname)
  end if
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

SUBROUTINE ReadInput(jobname, simulation_type, B0, C11, C12, C44, CDT_update, CDTintegrator, &
  CDTscheme, CDTlimiter, CDTtoler, CDTintFlux, CP11, CP12, CP44, CT11, CT12, CT44, specific_heat, &
  DDCtoler, DMBtoler, L0, Lbar, Lrho, MaxIter, Boltz, bulk_c1, bulk_c2, &
  crystalstruct, burger, c1, cdt_A, cdt_Aint, cdt_G_n0, cdt_K, cdt_Y_e, &
  cdt_pn, cdt_qn, cdt_rhoDot_n0, cdt_tau_n0, dt0, dt_output, euler_angle, g0, gtol, impact_force, &
  impact_angle, impact_velocity, impact_acc_time, internal_strain_flag, ipOut, mu, no_flux, omega0, p, &
  eos_perc_B1, eos_perc_B_star, eos_perc_Gamma0, eos_perc_R_M, eos_perc_Td0, eos_perc_a, &
  eos_perc_b, eos_perc_eos_flag, eos_perc_gamma, eos_perc_kappa, eos_perc_rho0, eos_perc_v0, &
  eos_perc_v_star, q, rho0, rhobar0, t_stop, t_vel, tau0, temperature0, wave_vel, drag_flag, drag_Trho_flag, &
  backstress_model)

USE GlobalParams
use utilities, only: Fatal
!$   Use omp_lib
implicit none

CHARACTER(32), INTENT(IN) :: jobname
!------------ local variables:
integer :: ios, j, ich, j_eul
character(32) :: key
character(54) :: inputfilename
REAL(KIND=8) :: value1, Zener
character(256) :: line, values
LOGICAL :: echoinput
!--------------------------- outputs:
REAL(KIND=8), INTENT(OUT) :: B0(Nchar), C11, C12, C44, CDTtoler, CP11, CP12, CP44, CT11, CT12, CT44, specific_heat
character(32), INTENT(OUT) :: CDT_update, CDTintegrator, CDTscheme, CDTlimiter, CDTintFlux, drag_flag, backstress_model
CHARACTER(6), INTENT(OUT) :: crystalstruct
CHARACTER(10), INTENT(OUT) :: simulation_type
REAL(KIND=8), INTENT(OUT) :: DDCtoler, DMBtoler, L0, Lbar, Lrho, Boltz, bulk_c1, bulk_c2
INTEGER, INTENT(OUT) :: MaxIter, ipOut
REAL(KIND=8), INTENT(OUT) :: burger, c1, cdt_A, cdt_Aint, cdt_G_n0, cdt_K, cdt_Y_e
REAL(KIND=8), INTENT(OUT) :: cdt_pn, cdt_qn, cdt_rhoDot_n0, cdt_tau_n0, dt0, dt_output, g0, gtol
REAL(KIND=8), INTENT(OUT) :: impact_angle, impact_force, impact_velocity, impact_acc_time, mu, omega0, p
REAL(KIND=8), INTENT(OUT) :: eos_perc_B1, eos_perc_B_star, eos_perc_Gamma0, eos_perc_R_M
REAL(KIND=8), INTENT(OUT) :: eos_perc_Td0, eos_perc_a, eos_perc_b, eos_perc_gamma, eos_perc_kappa
REAL(KIND=8), INTENT(OUT) :: eos_perc_rho0, eos_perc_v0, eos_perc_v_star
REAL(KIND=8), INTENT(OUT) :: q, rho0, rhobar0, t_stop, t_vel, tau0, temperature0, wave_vel(Nchar)
LOGICAL, INTENT(OUT) :: internal_strain_flag, no_flux, eos_perc_eos_flag, drag_Trho_flag
REAL(KIND=8), DIMENSION(Nregion,3), INTENT(OUT) :: euler_angle

! SET SOME DEFAULT VALUES
echoinput = .True.
simulation_type = 'impact'
crystalstruct = 'fcc'
Zener=1.d0
impact_force=0.d0    ! impact_velocity overrides impact_force
impact_velocity=0.d0
impact_acc_time=rzero
CDTtoler = 5.0d-2    ! tolerance on "relative error" in disloc. densities
CDT_update = 'FEM'   ! FVM or FEM
CDTintegrator = 'FEuler'  ! FEuler or RK2
CDTscheme = 'SUPG'   ! With FVM: Central, Godunov, LxF, or LLF; with FEM: Bubnov or SUPG
CDTlimiter = '---'   ! With FVM: zero, Fromm, BW, LxW, MMod, SBee, or MC; with FEM: not used
CDTintFlux = 'nothing' ! allowed values: 'zero' (=old default), 'end', 'nothing'
DDCtoler = 20.d-2    ! tolerance on "relative error" in epsilon_i
DMBtoler = 0.5d-2    ! tolerance on "relative error" in Cauchy stress
L0 = 1.d0 ! mm
MaxIter = 1          ! Maximum number of iterations (over sub-problems) in a time-step
ipOut = -1          ! integration point where history variables written (default: last element)
dt0 = 1.0d-11 ! s
dt_output = 2.0d-9  ! s
gtol = 1.0d-03
t_stop = 5.0d-7  ! s
internal_strain_flag = .False.
no_flux = .False.    !True to turn off dislocation flux
eos_perc_eos_flag = .True.
eos_perc_gamma = 2.d0/3.d0
eos_perc_B1 = 0.d0 ! this param. and the two below will be computed if any one of them is not set
eos_perc_B_star = 0.d0
eos_perc_v_star = 0.d0
! initialize some EOS model parameters to zero, if not set in input file, psi_ion =0 and/or psi_electron =0
specific_heat = 0.d0 ! optional, only used if psi_ion=psi_electron=0
eos_perc_Gamma0 = 0.d0 ! sets psi_el=0
eos_perc_kappa = 0.d0
eos_perc_R_M = 0.d0 ! sets psi_ion=0
eos_perc_Td0 = 100.d0 ! initialize a value for reference Debye temp. that will not crash the code if not set
eos_perc_a = 0.d0
eos_perc_b = 0.d0
do j=1,Nregion
  euler_angle(j,:) = (/0d0, 0.d0, 0.d0/)
end do
j_eul = 1
Boltz = 1.380649d-20 ! mJ/K
mu = 0.d0 ! if this is set to zero, wave_vel and mu are computed for each disloc. character
wave_vel(:) = 0.d0
temperature0    = 300.0d0   ! K
drag_flag = 'Austin' ! choose fct form for drag coeff.: 'Austin' (default), 'iso', or 'const '
drag_Trho_flag = .False. ! set to .True. to make drag coeff. T and rho dependent
B0 = 3.d-11 ! default value for low velocity drag coefficient is a rough order-of-magnitude estimate
backstress_model = 'gradient'

inputfilename = "input_parameters."//trim(jobname)//".dat"
!~ print*, "reading ", inputfilename

open(unit=42, file=inputfilename, action="read", iostat=ios, status='old')
if (ios/=0) then
  close(unit=42)
!~   print*, "failed, now reading ", jobname
  open(unit=42, file=jobname, action="read", iostat=ios, status='old')
  if (ios/=0) then
    close(unit=42)
    call Fatal("Error: cannot open input file, tried "//trim(inputfilename)//" and "//jobname)
  end if
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
  if (key=='echoinput') read(values,*)echoinput
  if (key=='simulation_type') read(values,*)simulation_type
  if (key=='Nslip') read(values,*)Nslip
  if (key=='crystalstruct') read(values,*)crystalstruct
  if (key=='B0') read(line,*) key,B0(1:Nchar)
  if (key=='C11') read(values,*)C11
  if (key=='C12') read(values,*)C12
  if (key=='C44') read(values,*)C44
  if (key=='CDT_update') read(values,*)CDT_update
  if (key=='CDTintegrator') read(values,*)CDTintegrator
  if (key=='CDTlimiter') read(values,*)CDTlimiter
  if (key=='CDTscheme') read(values,*)CDTscheme
  if (key=='CDTtoler') read(values,*)CDTtoler
  if (key=='CDTintFlux') read(values,*)CDTintFlux
  if (key=='backstress_model') read(values,*)backstress_model
  if (key=='CP11') read(values,*)CP11
  if (key=='CP12') read(values,*)CP12
  if (key=='CP44') read(values,*)CP44
  if (key=='CT11') read(values,*)CT11
  if (key=='CT12') read(values,*)CT12
  if (key=='CT44') read(values,*)CT44
  if (key=='specific_heat') read(values,*)specific_heat
  
  if (key=='DDCtoler') read(values,*)DDCtoler
  if (key=='DMBtoler') read(values,*)DMBtoler
  if (key=='L0') read(values,*)L0
!  if (key=='Lbar') read(values,*)Lbar ! computed from rho0 below
!  if (key=='Lrho') read(values,*)Lrho ! computed from Lbar**2 below
  if (key=='boltz') read(values,*)Boltz
  if (key=='bulk_c1') read(values,*)bulk_c1
  if (key=='bulk_c2') read(values,*)bulk_c2
  if (key=='burger') read(values,*)burger
  if (key=='c1') read(values,*)c1
  if (key=='cdt_A') read(values,*)cdt_A
  if (key=='cdt_Aint') read(values,*)cdt_Aint
  if (key=='cdt_G_n0') read(values,*)cdt_G_n0
  if (key=='cdt_K') read(values,*)cdt_K
  if (key=='cdt_Y_e') read(values,*)cdt_Y_e
  if (key=='cdt_pn') read(values,*)cdt_pn
  if (key=='cdt_qn') read(values,*)cdt_qn
  if (key=='cdt_rhoDot_n0') read(values,*)cdt_rhoDot_n0
  if (key=='cdt_tau_n0') read(values,*)cdt_tau_n0
  if (key=='dt0') read(values,*)dt0
  if (key=='dt_output') read(values,*)dt_output
  
  if (key=='g0') read(values,*)g0
  if (key=='gtol') read(values,*)gtol
  if (key=='impact_angle') read(values,*)impact_angle
  if (key=='impact_velocity') read(values,*)impact_velocity
  if (key=='impact_force') read(values,*)impact_force
  if (key=='impact_acc_time') read(values,*)impact_acc_time
  if (key=='internal_strain_flag') read(values,*)internal_strain_flag
  if (key=='mu') read(values,*)mu
  if (key=='no_flux') read(values,*)no_flux
  if (key=='omega0') read(values,*)omega0
  if (key=='p') read(values,*)p
  if (key=='props_dmb_perc_eos_perc_B1'.or.key=='eos_B1') read(values,*)eos_perc_B1
  if (key=='props_dmb_perc_eos_perc_B_star'.or.key=='eos_B_star') read(values,*)eos_perc_B_star
  if (key=='props_dmb_perc_eos_perc_Gamma0'.or.key=='eos_Gamma0') read(values,*)eos_perc_Gamma0
  if (key=='props_dmb_perc_eos_perc_R_M'.or.key=='eos_R_M') read(values,*)eos_perc_R_M
  if (key=='props_dmb_perc_eos_perc_Td0'.or.key=='eos_Td0') read(values,*)eos_perc_Td0
  if (key=='props_dmb_perc_eos_perc_a'.or.key=='eos_a') read(values,*)eos_perc_a
  if (key=='props_dmb_perc_eos_perc_b'.or.key=='eos_b') read(values,*)eos_perc_b
  if (key=='props_dmb_perc_eos_perc_eos_flag'.or.key=='eos_flag') read(values,*)eos_perc_eos_flag
  if (key=='props_dmb_perc_eos_perc_gamma'.or.key=='eos_gamma') read(values,*)eos_perc_gamma
  if (key=='props_dmb_perc_eos_perc_kappa'.or.key=='eos_kappa') read(values,*)eos_perc_kappa
  if (key=='props_dmb_perc_eos_perc_rho0'.or.key=='eos_rho0') read(values,*)eos_perc_rho0
  if (key=='props_dmb_perc_eos_perc_v0'.or.key=='eos_v0') read(values,*)eos_perc_v0
!  if (key=='props_dmb_perc_eos_perc_v0'.or.key=='eos_v0') read(values,*)eos_perc_v0 ! computed from eos_perc_rho0 below
  if (key=='props_dmb_perc_eos_perc_v_star'.or.key=='eos_v_star') read(values,*)eos_perc_v_star
  
  if (key=='q') read(values,*)q
  if (key=='rho0') read(values,*)rho0
  if (key=='rhobar0') read(values,*)rhobar0
  if (key=='t_stop') read(values,*)t_stop
  if (key=='t_vel') read(values,*)t_vel
  if (key=='tau0') read(values,*)tau0
  if (key=='temperature0') read(values,*)temperature0
  if (key=='wave_vel') read(line,*) key,wave_vel(1:Nchar) ! computed from mu below, if mu is given explicitly!
  if (key=='all_euler_angles') then
    !! use this keyword to read all euler angles from one input line
    read(line,*)key,(euler_angle(j,1),euler_angle(j,2),euler_angle(j,3), j=1,Nregion)
  end if
  if ((key=='euler_angle') .and. (j_eul<=Nregion)) then
    !! use this keyword to read one set of euler angles, use multiple entries in order to read in for all Nregion
    read(line,*)key,euler_angle(j_eul,1),euler_angle(j_eul,2),euler_angle(j_eul,3)
    j_eul = j_eul + 1
  end if
  if (key=='drag_flag') read(line,*)key,drag_flag
  if (key=='drag_Trho_flag') read(line,*)key,drag_Trho_flag
  if (key=='debug_cycle') then
    read(values,*)debug_cycle
    print*,"setting debug_cycle=",debug_cycle
  end if
  
  if (key=='MaxIter') then
    read(values,*)value1
    MaxIter = int(value1)
  elseif (key=='ipOut') then
    read(values,*)value1
    ipOut = int(value1)
  end if
end do ! read file
close(unit=42)

!! if only a subset of euler angles was set explicitly, fill the rest with the last value that was read:
!! (useful for single crystal calcs. where Nregion>1 was set to take advantage of openmp parallelization of loops over Nregion)
if ((j_eul>1) .and. (j_eul<=Nregion)) then
  do j=j_eul,Nregion
    euler_angle(j,:) = euler_angle(j_eul-1,:)
  end do
end if

! compute some derived values:
Lbar = 1.d0 / sqrt(rho0)   !((Mg/mm^3)  
Lrho = Lbar**2
if (mu > 0.d0) then
  wave_vel(:) = sqrt(mu / rhobar0) ! mm / s
end if
Zener = 2*C44/(C11-C12)

! calculate an estimate for character dependent wave_vel, if not provided by the user (WARNING: implemented only for fcc):
if ( abs(wave_vel(1)) < 1.d-15) then
  wave_vel(1) = sqrt((3.d0*C44*(C11-C12))/(rhobar0*2.d0*(C44+C11-C12))) ! analytic solution for vcrit of fcc screw
  wave_vel(Nchar) = min(sqrt((C11-C12)/(2.d0*rhobar0)),sqrt(C44/rhobar0)) ! analytic solution for vcrit of fcc edge
  do ich=2, Nchar-1
    !TODO: generalize! (for now we just interpolate between vcrit for edge and screw)
    wave_vel(ich) = max(wave_vel(1),wave_vel(Nchar)) - abs(wave_vel(1) - wave_vel(Nchar))*(ich-1.d0)/(Nchar-1.d0)
  end do !ich
end if

! fall-back for bulk modulus Bstar, its pressure derivative B1 and reference volume vstar:
if ( abs(eos_perc_B1*eos_perc_B_star*eos_perc_v_star) < 1.d-15 ) then
   ! print*,"Warning: missing some eos values, computing B_star, B1, and v_star from SOEC and rhobar0"
   eos_perc_v_star = 1.d0/rhobar0 ! reference volume at ambient conditions
   eos_perc_B_star = (C11 + 2.d0*C12)/3.d0 ! bulk modulus from cubic SOEC at v_star
   eos_perc_B1 = (CP11 + 2.d0*CP12)/3.d0 ! pressure derivative of bulk modulus from CPij at v_star
end if

if ( abs(eos_perc_rho0) < 1.d-15) then
   eos_perc_rho0 = rhobar0
end if
if ( abs(eos_perc_v0) < 1.d-15) then
   eos_perc_v0       = 1.d0 / eos_perc_rho0 !((Mg/mm^3)
end if

if ((crystalstruct == "fcc") .and. (Nslip>12)) then
  call Fatal('Nslip must be <= 12 for fcc')
elseif ((crystalstruct == "bcc") .and. (Nslip>48)) then
  call Fatal('Nslip must be <= 48 for fcc')
end if

if (impact_acc_time < rzero) impact_acc_time=rzero

if (ipOut<0) then
  ipOut = Nel + ipOut + 1 ! allow python-like syntax where -1 means the last element of the array
end if

if (echoinput) then
  print*,"running '",trim(simulation_type),"' simulation using:"
  print*,"crystalstruct ",crystalstruct, " Nslip ",Nslip
  print*,"B0",B0
  print*,"C11",C11
  print*,"C12",C12
  print*,"C44",C44
  print*,"Zener ratio",Zener
  print*,"CDT_update ",CDT_update
  print*,"CDTintegrator ",CDTintegrator
  print*,"CDTlimiter ",CDTlimiter
  print*,"CDTscheme ",CDTscheme
  print*,"CDTtoler",CDTtoler
  print*,"CDTintFlux ",CDTintFlux
  print*,"CP11",CP11
  print*,"CP12",CP12
  print*,"CP44",CP44
  print*,"CT11",CT11
  print*,"CT12",CT12
  print*,"CT44",CT44
  print*,"specific_heat",specific_heat

  print*,"DDCtoler",DDCtoler
  print*,"DMBtoler",DMBtoler
  print*,"L0",L0
  print*,"Lbar",Lbar
  print*,"Lrho=Lbar**2",Lrho
  print*,"MaxIter",MaxIter
  print*,"boltz",Boltz
  print*,"bulk_c1",bulk_c1
  print*,"bulk_c2",bulk_c2
  print*,"burger",burger
  print*,"c1",c1
  print*,"cdt_A",cdt_A
  print*,"cdt_Aint",cdt_Aint
  print*,"cdt_G_n0",cdt_G_n0
  print*,"cdt_K",cdt_K
  print*,"cdt_Y_e",cdt_Y_e
  print*,"cdt_pn",cdt_pn
  print*,"cdt_qn",cdt_qn
  print*,"cdt_rhoDot_n0",cdt_rhoDot_n0
  print*,"cdt_tau_n0",cdt_tau_n0
  print*,"dt0",dt0
  print*,"dt_output",dt_output

  print*,"g0",g0
  print*,"gtol",gtol
  print*,"impact_angle",impact_angle
  print*,"impact_force",impact_force
  print*,"impact_velocity",impact_velocity
  print*,"impact_acc_time",impact_acc_time
  print*,"internal_strain_flag",internal_strain_flag
  print*,"ipOut",ipOut
  print*,"mu",mu
  print*,"no_flux",no_flux
  print*,"omega0",omega0
  print*,"p",p
  print*,"eos_B1",eos_perc_B1
  print*,"eos_B_star",eos_perc_B_star
  print*,"eos_Gamma0",eos_perc_Gamma0
  print*,"eos_R_M",eos_perc_R_M
  print*,"eos_Td0",eos_perc_Td0
  print*,"eos_a",eos_perc_a
  print*,"eos_b",eos_perc_b
  print*,"eos_flag",eos_perc_eos_flag
  print*,"eos_gamma",eos_perc_gamma
  print*,"eos_kappa",eos_perc_kappa
  print*,"eos_rho0",eos_perc_rho0
  print*,"eos_v0",eos_perc_v0
  print*,"eos_v_star",eos_perc_v_star
  print*,"q",q
  print*,"rho0",rho0
  print*,"rhobar0",rhobar0
  print*,"t_stop",t_stop
  print*,"t_vel",t_vel
  print*,"tau0",tau0
  print*,"temperature0",temperature0
  print*,"wave_vel",wave_vel
  print*,"drag_flag ",drag_flag
  print*,"drag_Trho_flag ",drag_Trho_flag
  print*,"backstress_model ",backstress_model
  print*,"list of euler_angles:"
  do j_eul=1,Nregion
    print*,j_eul,euler_angle(j_eul,:)
  end do
  print*,"done reading input"
end if

!$    print*, 'OpenMP parallelization enabled:'
!$    print*, 'using ',omp_get_max_threads(),' of ',omp_get_num_procs(),' processors'

END SUBROUTINE ReadInput
end module readfiles
