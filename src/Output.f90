! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
module OUTPUTFILES
public
contains
!-----------------------------------------------------------------------
!>  Opens all output files.
SUBROUTINE Init_All_Field_Output_Files(jobname)
!-----------------------------------------------------------------------
  USE OUTPUT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  CHARACTER*(*) jobname    !<job name
  CHARACTER(1) str_ich
  INTEGER :: ich
!-----------------------------------------------------------------------
!  Locals:
!-----------------------------------------------------------------------
  call assign_units()
  
  call Init_Field_File(jobname, unit_dmb_pos    , "node", "dmb_pos"    )
  call Init_Field_File(jobname, unit_dmb_vel    , "node", "dmb_vel"    )
  call Init_Field_File(jobname, unit_dmb_acc    , "node", "dmb_acc"    )
  call Init_Field_File(jobname, unit_stress     , "ip",   "stress"     )
  do ich=1,Nchar
    write(str_ich, '(i1)' )ich
    call Init_Field_File(jobname, unit_vel(ich)        , "ip",   "dis_vel." // str_ich  )
    call Init_Field_File(jobname, unit_vel_x(ich)      , "ip",   "dis_vel_x." // str_ich  )
    call Init_Field_File(jobname, unit_vel_y(ich)      , "ip",   "dis_vel_y." // str_ich  )
    call Init_Field_File(jobname, unit_vel_z(ich)      , "ip",   "dis_vel_z." // str_ich  )
    call Init_Field_File(jobname, unit_dis_acc(ich)    , "ip",   "dis_acc." // str_ich  )
    call Init_Field_File(jobname, unit_rho_pos(ich)    , "ip",   "rho_pos." // str_ich    )
    call Init_Field_File(jobname, unit_rho_neg(ich)    , "ip",   "rho_neg." // str_ich    )
    call Init_Field_File(jobname, unit_gamma_acc(ich)  , "node", "gamma_acc." // str_ich)
  end do
  call Init_Field_File(jobname, unit_const_soln , "ip",   "const_soln" )
  call Init_Field_File(jobname, unit_Fp         , "ip",   "Fp"         )
  call Init_Field_File(jobname, unit_Fe         , "ip",   "Fe"         )
  call Init_Field_File(jobname, unit_Ei         , "ip",   "Ei"         )
  call Init_Field_File(jobname, unit_tau        , "ip",   "tau"       )
  call Init_Field_File(jobname, unit_tau_back   , "ip",   "tau_back")
  call Init_Field_File(jobname, unit_temperature, "ip",   "T")
  call Init_Field_File(jobname, unit_visar      , "th",   "visar")
  call Init_Field_File(jobname, unit_stress_th  , "th",   "stress")  
  call Init_Field_File(jobname, unit_velocity_th, "th",   "vel")


  write(unit_dmb_pos    ,100) "position", "inc", "time", "node id", "x1 (pos)", "x2 (pos)", "x3 (pos)"
  write(unit_dmb_vel    ,100) "velocity", "inc", "time", "node id", "v1 (vel)", "v2 (vel)", "v3 (vel)"
  write(unit_dmb_acc    ,100) "accel", "inc", "time", "node id", "a1 (acc)", "a2 (acc)", "a3 (acc)"
  write(unit_stress     ,101) "inc", "time", "ip id",   "sigma-11", "sigma-22",  "sigma-33", "sigma-23", "sigma-31", "sigma-12"
  do ich=1,Nchar
    write(unit_vel(ich)        ,102) "inc", "time", "ip id",   "vel-1",       "vel-2",     "vel-3",  "..."
    write(unit_vel_x(ich)      ,103) "inc", "time", "ip id",   "vel-1",       "vel-2",     "vel-3",  "..."
    write(unit_vel_y(ich)      ,103) "inc", "time", "ip id",   "vel-1",       "vel-2",     "vel-3",  "..."
    write(unit_vel_z(ich)      ,103) "inc", "time", "ip id",   "vel-1",       "vel-2",     "vel-3",  "..."
    write(unit_dis_acc(ich)    ,103) "inc", "time", "ip id",   "acc-1",       "acc-2",     "acc-3",  "..."
    write(unit_rho_pos(ich)    ,104) "inc", "time", "ip id",   "dd(+)-1",     "dd(+)-2",   "dd(+)-3",  "..."
    write(unit_rho_neg(ich)    ,105) "inc", "time", "ip id",   "dd(-)-1",     "dd(-)-2",   "dd(-)-3",  "..."
    write(unit_gamma_acc(ich)  ,114) "inc", "time", "ip id",   "gamma-1", "gamma-2", "gamma-3", "..."
  end do
  write(unit_Fp         ,106) "inc", "time", "ip id",   "Fp-11", "Fp-12", "Fp-13", "Fp-21", "Fp-22", "Fp-23", "Fp-31", &
    "Fp-32", "Fp-33"
  write(unit_Fe         ,107) "inc", "time", "ip id",   "Fe-11", "Fe-12", "Fe-13", "Fe-21", "Fe-22", "Fe-23", "Fe-31", &
    "Fe-32", "Fe-33"
  write(unit_tau        ,108) "inc", "time", "ip id",   "tau-1",        "tau2",     "tau-3",  "..."
  write(unit_tau_back   ,111) "inc", "time", "ip id",   "tau_b-1",      "tau_b2",   "tau_b-3",  "..."
  write(unit_const_soln ,109) "inc", "time", "ip id",   "subinc",      "cutback"
  write(unit_temperature,110) "inc", "time", "ip id",   "temperature"
  write(unit_Ei         ,112) "inc", "time", "ip id",   "Ei-11", "Ei-22", "Ei-33", "Ei-23", "Ei-31", "Ei-12"
  write(unit_visar      ,113) "time", "node id", "xdot1", "xdot2", "xdot3"
  write(unit_velocity_th,113) "time", "node id", "xdot1", "xdot2", "xdot3"  
  write(unit_stress_th  ,115) "time", "ip id",   "sigma-11", "sigma-22",  "sigma-33", "sigma-23", "sigma-31", "sigma-12"



100 format("## dmb_dof field file contains ",a," of nodes (at nodes)",/,  &
   "##",a7,a13,a8,3a17,/,  &
   "## ------   ---------- -------",3(3x, 14('-')))

101 format("## stress field file contains components of cauchy stress at integration points",/,  &
   "##",a7,a13,a8,6a17,/,  &
   "## ------   ---------- -------",6(3x, 14('-')))

102 format("## vel field file contains dislocation velocities on each slip system (vel-sys#) at integration points",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

103 format("## vel_x field file contains x-component of dislocation velocities on each slip system (vel-sys#) at IP",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

104 format("## rho_pos field file contains dislocation densities of + character on each slip system (dd-sys#) at IP",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

105 format("## rho_neg field file contains dislocation densities of - character on each slip system (dd-sys#) at IP",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

106 format("## Fp field file contains coefficients of plastic def. gradient at integration points",/,  &
   "##",a7,a13,a8,9a17,/,  &
   "## ------   ---------- -------",9(3x, 14('-')))

107 format("## Fe field file contains coefficients of elastic def. gradient at integration points",/,  &
   "##",a7,a13,a8,9a17,/,  &
   "## ------   ---------- -------",9(3x, 14('-')))

108 format("## tau field file contains projection of cauchy stress on slip systems (no backstress) at IP",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

109 format("## const_soln field file contains number of sub-increments, and total number of constitutive cutbacks at IP",/,  &
   "##",a7,a13,a8,2a17,/,  &
   "## ------   ---------- -------",2(3x, 14('-')))

110 format("## temp field files contains temperature at integration points",/,  &
   "##",a7,a13,a8,a17,/,  &
   "## ------   ---------- -------",1(3x, 14('-')))

111 format("## tau_back field file contains backstress on slip systems at integration points",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-'))) 

112 format("## Ei field file contains internal strain due to dislocation configuration",/,  &
   "##",a7,a13,a8,6a17,/,  &
   "## ------   ---------- -------",6(3x, 14('-'))) 

113 format("## free surface velocity, i.e. components of velocity vector at the right boundary",/,  &
   "##",a9,a8,3a22,/,  &
   "## -------- -------",3(3x, 19('-'))) 

114 format("## accumulated magnitudes of slip on each slip system at IPs",/,  &
   "##",a7,a13,a8,4a17,/,  &
   "## ------   ---------- -------",4(3x, 14('-')))

115 format("## stress history file contains components of cauchy stress at each time",/,  &
   "##",a9,a8,6a22,/,  &
   "## -------- -------",6(3x, 19('-')))

  RETURN
END SUBROUTINE Init_All_Field_Output_Files

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Opens one output file.
SUBROUTINE Init_Field_File(jobname, unitnumber, loc, vrb)
!-----------------------------------------------------------------------
  USE OUTPUT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  CHARACTER(32) :: jobname      !< job name
  INTEGER       :: unitnumber   !< unit number
  CHARACTER*(*) :: loc          !< location of value in mesh
  CHARACTER*(*) :: vrb          !< Variable name
!-----------------------------------------------------------------------
!  Locals:
  CHARACTER(128) :: filename
!-----------------------------------------------------------------------
  filename = jobname(1:len_trim(jobname)) // "." // loc // "." // vrb // ".F90txt"
  open(unit=unitnumber, file=filename)

  RETURN
END SUBROUTINE Init_Field_File

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Closes all output files.
SUBROUTINE Close_Field_Files()
!-----------------------------------------------------------------------
  USE OUTPUT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i
!-----------------------------------------------------------------------
  do i=10,27
    close(unit=i)
  end do  !i

  RETURN
END SUBROUTINE Close_Field_Files

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE Write_State(iinc, dmb_state_new, cdt_state_new, ddc_state_new)
!-----------------------------------------------------------------------
  USE GlobalParams
  USE OUTPUT
  USE DMB
  USE CDT
  USE DDC
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  INTEGER, intent(in)               :: iinc
  TYPE(Type_DMB_STATE), intent(in)  :: dmb_state_new
  TYPE(Type_CDT_STATE), intent(in)  :: cdt_state_new
  TYPE(Type_DDC_STATE), intent(in)  :: ddc_state_new
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: t_out
  REAL(KIND=8), allocatable, DIMENSION(:,:,:) :: cauchy, Ei
  INTEGER :: ip,i,slip,j,ich
!-----------------------------------------------------------------------
  allocate(cauchy(Nel,3,3), Ei(Nel,3,3))
  
  t_out = dmb_state_new%time(1)

  cauchy = dmb_state_new%cauchy
  do ip=1,Nel
    write(unit_stress,100) iinc, t_out, ip, (cauchy(ip,i,i),i=1,3),cauchy(ip,2,3),cauchy(ip,3,1),cauchy(ip,1,2)
  end do  !ip

  do ip=1,Nnode
    write(unit_dmb_pos,  100) iinc, t_out, ip, dmb_state_new%x(ip,:)
    write(unit_dmb_vel,  100) iinc, t_out, ip, dmb_state_new%xdot(ip,:)
    write(unit_dmb_acc,  100) iinc, t_out, ip, dmb_state_new%xddot(ip,:)
  end do  !ip

!--- begin loop over characters
  do ich=1,Nchar
  do ip=1,Nel
    write(unit_vel(ich),100) iinc, t_out, ip, (dmb_state_new%vel(ip,slip,ich),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_vel_x(ich),100) iinc, t_out, ip, (dmb_state_new%vel_x(ip,slip,ich),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_vel_y(ich),100) iinc, t_out, ip, (dmb_state_new%vel_y(ip,slip,ich),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_vel_z(ich),100) iinc, t_out, ip, (dmb_state_new%vel_z(ip,slip,ich),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_dis_acc(ich),100) iinc, t_out, ip, (dmb_state_new%dis_acc(ip,slip,ich),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_rho_pos(ich),100) iinc, t_out, ip, (cdt_state_new%varrho_field(ip,slip,ich,1),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_rho_neg(ich),100) iinc, t_out, ip, (cdt_state_new%varrho_field(ip,slip,ich,2),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_gamma_acc(ich),100) iinc, t_out, ip, (dmb_state_new%gamma_acc(ip,slip,ich),slip=1,Nslip)
  end do  !ip
  
  write(unit_vel_x(ich),110)
  write(unit_vel_y(ich),110)
  write(unit_vel_z(ich),110)
  write(unit_vel(ich),110)
  write(unit_dis_acc(ich),110)
  write(unit_rho_pos(ich),110)
  write(unit_rho_neg(ich),110)

  end do !ich
!--- end loop over characters

  do ip=1,Nel
    write(unit_Fp,100) iinc, t_out, ip, ((dmb_state_new%Fp(ip,i,j),j=1,3),i=1,3)
  end do  !ip

  do ip=1,Nel
    write(unit_Fe,100) iinc, t_out, ip, ((dmb_state_new%Fe(ip,i,j),j=1,3),i=1,3)
  end do  !ip

  Ei = ddc_state_new%epsilon_i
  do ip=1,Nel
    write(unit_Ei,100) iinc, t_out, ip, (Ei(ip,j,j),j=1,3),Ei(ip,2,3),Ei(ip,3,1),Ei(ip,1,2)
  end do  !ip

  do ip=1,Nel
    write(unit_temperature,100) iinc, t_out, ip, dmb_state_new%T(ip,1)
  end do  !ip

  do ip=1,Nel
    write(unit_const_soln,100) iinc, t_out, ip, dble(dmb_state_new%Nsub(ip)), dble(dmb_state_new%i_cutback(ip))
  end do  !ip
  
  do ip=1,Nel
    write(unit_tau,100) iinc, t_out, ip, (dmb_state_new%tau(ip,slip),slip=1,Nslip)
  end do  !ip

  do ip=1,Nel
    write(unit_tau_back,100) iinc, t_out, ip, (ddc_state_new%tau_b(ip,slip),slip=1,Nslip)
  end do  !ip

  write(unit_stress,110)
  write(unit_dmb_pos,110)
  write(unit_dmb_vel,110)
  write(unit_dmb_acc,110)
  write(unit_Fp,110)
  write(unit_Fe,110)
  write(unit_Ei,110)
  write(unit_temperature,110)
  write(unit_const_soln,110)
  write(unit_tau,110)
  write(unit_tau_back,110)


!100 format(i9,2x,es11.4,i8,20(2x,es25.18))
!! dnb: do we really need so many digits in all our output files?
100 format(i9,2x,es11.4,i8,48(2x,es16.9)) !! 48=max Nslip on same line
110 format(//)

  RETURN
END SUBROUTINE Write_State

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------


subroutine write_history(ipOut, dmb_state_new) !, cdt_state_new, ddc_state_new)
!-----------------------------------------------------------------------
  USE GlobalParams
  USE OUTPUT
  USE DMB
  USE CDT
  USE DDC
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_STATE), intent(in)  :: dmb_state_new
!~   TYPE(Type_CDT_STATE)           :: cdt_state_new
!~   TYPE(Type_DDC_STATE)           :: ddc_state_new
  integer, intent(in   )         :: ipOut
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8)  :: t_out, cauchy(3,3)
  integer :: i 
!-----------------------------------------------------------------------
  t_out = dmb_state_new%time(1)
  
  cauchy = dmb_state_new%cauchy(ipOut,:,:)
  
  write(unit_visar,     100)  t_out, Nnode, dmb_state_new%xdot(Nnode,:)  
  write(unit_velocity_th, 100) t_out, (ipOut+1), dmb_state_new%xdot(ipOut+1, :)
  write(unit_stress_th, 100)  t_out, ipOut, (cauchy(i,i),i=1,3),cauchy(2,3),cauchy(3,1),cauchy(1,2)
  
  
100 format(es11.4,i8,20(2x,es20.13))

  return 
  
end subroutine write_history


!>  Print DDC state to terminal.
SUBROUTINE Print_DDC_State(ddc_state,ich)
!-----------------------------------------------------------------------
  USE GlobalParams
  USE OUTPUT
  USE DDC
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DDC_STATE), intent(in)         :: ddc_state
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i,ich
  REAL(KIND=8), DIMENSION(Nel)       :: cell_sizes
  REAL(KIND=8), allocatable, DIMENSION(:,:,:) :: tau_i
  REAL(KIND=8), allocatable, DIMENSION(:,:) :: tau_b
  REAL(KIND=8), allocatable, DIMENSION(:,:,:)   :: stress_i
  REAL(KIND=8), allocatable, DIMENSION(:,:,:)   :: epsilon_i
!-----------------------------------------------------------------------
  allocate(tau_i(Nel,Nslip,Nchar), tau_b(Nel,Nslip))
  allocate(stress_i(Nel,3,3), epsilon_i(Nel,3,3))
  
  !DDC state members:
  cell_sizes = ddc_state%cell_sizes  !< (Ncell)
  tau_i      = ddc_state%tau_i       !< (Ncell, Nslip, Nchar)
  tau_b      = ddc_state%tau_b       !< (Ncell, Nslip)
  stress_i   = ddc_state%stress_i    !< (Ncell, 3,3)
  epsilon_i  = ddc_state%epsilon_i   !< (Ncell, 3,3)

  write(*,'(a)') "DDC State:"
  write(*,'(a)') "  cell_sizes:"
  do i=1,Nel
    write(*,'(i2,2x,es28.16)') i,cell_sizes(i)
  end do
  write(*,'(a)') "  tau_i:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,tau_i(i,:,ich)
  end do
  write(*,'(a)') "  tau_b:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,tau_b(i,:)
  end do
  write(*,'(a)') "  stress_i:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,stress_i(i,:,1)
    write(*,'(   4x,3es28.16)')   stress_i(i,:,2)
    write(*,'(   4x,3es28.16)')   stress_i(i,:,3)
  end do
  write(*,'(a)') "  epsilon_i:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,epsilon_i(i,:,1)
    write(*,'(   4x,3es28.16)')   epsilon_i(i,:,2)
    write(*,'(   4x,3es28.16)')   epsilon_i(i,:,3)
  end do

  RETURN
END SUBROUTINE Print_DDC_State

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Print DMB state to terminal.
SUBROUTINE Print_DMB_State(dmb_state,ich)
!-----------------------------------------------------------------------
  USE GlobalParams
  USE OUTPUT
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_STATE), intent(in)   :: dmb_state
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i,ich
  REAL(KIND=8), DIMENSION(Nnode,3)   :: x
  REAL(KIND=8), DIMENSION(Nnode)     :: X_ref
  REAL(KIND=8), DIMENSION(Nnode,3)   :: xdot
  REAL(KIND=8), DIMENSION(Nnode,3)   :: xddot
  REAL(KIND=8), DIMENSION(Nel,1)     :: T
  REAL(KIND=8), allocatable, DIMENSION(:,:,:)   :: F, vel, vel_x, cauchy, Fp, Fe
  REAL(KIND=8), allocatable, DIMENSION(:,:)   :: tau
  INTEGER, DIMENSION(Nel)      :: Nsub
  INTEGER, DIMENSION(Nel)      :: i_cutback
!-----------------------------------------------------------------------
  allocate(F(Nel,3,3), cauchy(Nel,3,3), Fp(Nel,3,3), Fe(Nel,3,3))
  allocate(vel(Nel,Nslip,Nchar), vel_x(Nel,Nslip,Nchar))
  allocate(tau(Nel,Nslip))
  
  !DMB state members:
  x         = dmb_state%x
  X_ref     = dmb_state%X_ref
  xdot      = dmb_state%xdot
  xddot     = dmb_state%xddot
  F         = dmb_state%F
  T         = dmb_state%T
  vel       = dmb_state%vel
  vel_x     = dmb_state%vel_x
  cauchy    = dmb_state%cauchy
  Fp        = dmb_state%Fp
  Fe        = dmb_state%Fe
  tau       = dmb_state%tau
  Nsub      = dmb_state%Nsub
  i_cutback = dmb_state%i_cutback

  write(*,'(a)') "DMB State:"
  write(*,'(a)') "  x:"
  do i=1,Nnode
    write(*,'(i2,2x,3es28.16)') i,x(i,:)
  end do
  write(*,'(a)') "  X_ref:"
  do i=1,Nnode
    write(*,'(i2,2x,es28.16)') i,X_ref(i)
  end do
  write(*,'(a)') "  xdot:"
  do i=1,Nnode
    write(*,'(i2,2x,3es28.16)') i,xdot(i,:)
  end do
  write(*,'(a)') "  xddot:"
  do i=1,Nnode
    write(*,'(i2,2x,3es28.16)') i,xddot(i,:)
  end do
  write(*,'(a)') "  F:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,F(i,:,1)
    write(*,'(   4x,3es28.16)')   F(i,:,2)
    write(*,'(   4x,3es28.16)')   F(i,:,3)
  end do
  write(*,'(a)') "  T:"
  do i=1,Nel
    write(*,'(i2,2x,es28.16)') i,T(i,1)
  end do
  write(*,'(a)') "  vel:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,vel(i,:,ich)
  end do
  write(*,'(a)') "  vel_x:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,vel_x(i,:,ich)
  end do
  write(*,'(a)') "  cauchy:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,cauchy(i,:,1)
    write(*,'(   4x,3es28.16)')   cauchy(i,:,2)
    write(*,'(   4x,3es28.16)')   cauchy(i,:,3)
  end do
  write(*,'(a)') "  Fp:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,Fp(i,:,1)
    write(*,'(   4x,3es28.16)')   Fp(i,:,2)
    write(*,'(   4x,3es28.16)')   Fp(i,:,3)
  end do
  write(*,'(a)') "  Fe:"
  do i=1,Nel
    write(*,'(i2,2x,3es28.16)') i,Fe(i,:,1)
    write(*,'(   4x,3es28.16)')   Fe(i,:,2)
    write(*,'(   4x,3es28.16)')   Fe(i,:,3)
  end do
  write(*,'(a)') "  tau:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,tau(i,:)
  end do
  write(*,'(a)') "  Nsub:"
  do i=1,Nel
    write(*,'(i2,2x,i5)') i,Nsub(i)
  end do
  write(*,'(a)') "  i_cutback:"
  do i=1,Nel
    write(*,'(i2,2x,i5)') i,i_cutback(i)
  end do

  RETURN
END SUBROUTINE Print_DMB_State

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Print CDT state to terminal.
SUBROUTINE Print_CDT_State(cdt_state,ich)
!-----------------------------------------------------------------------
  USE GlobalParams
  USE OUTPUT
  USE CDT
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_CDT_STATE), intent(in)     :: cdt_state
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i,ich
  REAL(KIND=8), DIMENSION(Nel)           :: cell_sizes
  REAL(KIND=8), allocatable, DIMENSION(:,:,:,:)   :: varrho_field, varrho_node
!-----------------------------------------------------------------------
  allocate(varrho_field(Nel,Nslip,Nchar,2), varrho_node(Nnode,Nslip,Nchar,2))
  
  !CDT state members:
  cell_sizes   = cdt_state%cell_sizes   !(Nel)
  varrho_field = cdt_state%varrho_field !(Nel,   Nslip, Nchar, 2)
  varrho_node  = cdt_state%varrho_node  !(Nnode, Nslip, Nchar, 2)

  write(*,'(a)') "CDT State:"
  write(*,'(a)') "  cell_sizes:"
  do i=1,Nel
    write(*,'(i2,2x,es28.16)') i,cell_sizes(i)
  end do
  write(*,'(a)') "  varrho_field 1:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,varrho_field(i,:,ich,1)
  end do
  write(*,'(a)') "  varrho_field 2:"
  do i=1,Nel
    write(*,'(i2,2x,12es28.16)') i,varrho_field(i,:,ich,2)
  end do
  write(*,'(a)') "  varrho_node 1:"
  do i=1,Nnode
    write(*,'(i2,2x,12es28.16)') i,varrho_node(i,:,ich,1)
  end do
  write(*,'(a)') "  varrho_node 2:"
  do i=1,Nnode
    write(*,'(i2,2x,12es28.16)') i,varrho_node(i,:,ich,2)
  end do

  RETURN
END SUBROUTINE Print_CDT_State
end module OUTPUTFILES
