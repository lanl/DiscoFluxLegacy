! © 2026. Triad National Security, LLC. All rights reserved.
! This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
! which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration.
! All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
! Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up,
! irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform
! publicly and display publicly, and to permit others to do so.
module checks
contains
SUBROUTINE testequal(A,B,N,M,string,tolerance,count_pass,count_fail)
!  Check if two NxM matrices are equal (with some tolerance), and print 'test <string> PASSED/FAILED'
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: N, M
  REAL(KIND=8), INTENT(IN)  :: tolerance
  REAL(KIND=8), DIMENSION(N,M), INTENT(IN)  :: A, B
  character(*), INTENT(IN) :: string
  INTEGER :: count_pass,count_fail
  logical :: equal
  
  equal = all(abs(A-B) <= tolerance)
  
  if (equal) then
    count_pass = count_pass+1
    print*,"test " // string//": "//char(9)//"PASSED"
  else
    count_fail = count_fail+1
    print*,"test " // string//": "//char(9)//"FAILED"
  end if
  
  RETURN
END SUBROUTINE testequal

SUBROUTINE testequalarray(A,B,N,string,tolerance,count_pass,count_fail)
!  Check if two arrays are equal (with some tolerance), and print 'test <string> PASSED/FAILED'
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: N
  REAL(KIND=8), INTENT(IN)  :: tolerance
  REAL(KIND=8), DIMENSION(N), INTENT(IN)  :: A, B
  character(*), INTENT(IN) :: string
  INTEGER :: count_pass,count_fail
  logical :: equal
  
  equal = all(abs(A-B) <= tolerance)
  
  if (equal) then
    count_pass = count_pass+1
    print*,"test " // string//": "//char(9)//"PASSED"
  else
    count_fail = count_fail+1
    print*,"test " // string//": "//char(9)//"FAILED"
  end if
  
  RETURN
END SUBROUTINE testequalarray

SUBROUTINE testzero(A,string,tolerance,count_pass,count_fail)
! check if A=0 (with some tolerance), and print 'test <string> PASSED/FAILED'
  IMPLICIT NONE
  
  REAL(KIND=8), INTENT(IN)  :: tolerance, A
  character(*), INTENT(IN) :: string
  INTEGER :: count_pass,count_fail
  call testequalarray([A],[0.d0],1,string,tolerance,count_pass,count_fail)
  RETURN
END SUBROUTINE testzero
end module checks

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

PROGRAM RunTests
  use stdlib_linalg, only: expm
  USE GlobalParams
  USE DMB
  USE CDT
  USE DDC
  use DMB_routines, only: DMB_Compute_Rotation, DMB_form_elasticity, DMB_get_slip_system, DMB_Invert_Tensor26, DMB_Invert_Tensor43,&
                        stroh_dislocation, operator(.voigt.), operator(.unvoigt.)
  use CDT_routines
  use DDC_routines
  use utilities
  use checks

  IMPLICIT NONE
  
  CHARACTER(32) :: PASSED, FAILED
  INTEGER, PARAMETER :: resolution=361 ! resolve polar angle phi to 1 degree (tradeoff accuracy for speed)
  REAL(KIND=8) :: num, nu_poisson, rho0, Fwant, tmpintegral
  Real(KIND=8), DIMENSION(5) :: array1,array2
  REAL(KIND=8), DIMENSION(3)   :: euler, unit_y, unit_z
  REAL(KIND=8), DIMENSION(3,3) :: eye, A, B, rot, Q
  REAL(KIND=8), DIMENSION(6,6) :: eye6, mat2, A6, B6
  REAL(KIND=8), DIMENSION(3,3,3,3) :: A3, B3, elC, Cg
  REAL(KIND=8), allocatable, DIMENSION(:) :: x, func
  REAL(KIND=8), allocatable, DIMENSION(:) :: Integral
  REAL(KIND=8), allocatable, DIMENSION(:,:,:,:,:) :: uij
  REAL(KIND=8), DIMENSION(:), allocatable :: zeros_nslip
  REAL(KIND=8), DIMENSION(resolution) :: zeros, strainenergydensity, phi
  REAL(KIND=8), DIMENSION(:,:,:), allocatable :: Sb, Bb, tau_back_nokappa
  REAL(KIND=8), DIMENSION(:), allocatable       :: X_ref
  REAL(KIND=8), allocatable, DIMENSION(:,:,:,:)  :: varrho
  REAL(KIND=8), allocatable, DIMENSION(:,:,:)   :: tau_back1, tau_back2, rho_total
  REAL(KIND=8), allocatable, DIMENSION(:,:) :: zeros_threenslip, vel, tline, zeros2
  INTEGER :: i, j, k, l, count_fail, count_pass, start_time, finish_time, countrate
  TYPE(Type_DMB_PROPS)           :: props_dmb
  TYPE(Type_DMB_STATE)           :: dmb_state
  TYPE(Type_CDT_PROPS)           :: props_cdt
  TYPE(Type_DDC_PROPS)           :: ddc_props
  
  PASSED = char(9)//"PASSED"; FAILED = char(9)//"FAILED"
  count_fail=0; count_pass=0
  call system_clock(start_time,countrate)
  ! allocate memory
  allocate(x(10000), func(10000))
  allocate(integral(9999), zeros_nslip(Nslip), X_ref(Nnode))
  allocate(uij(3,3,resolution,Nslip,Nchar))
  allocate(varrho(Nel,Nslip,Nchar,2), zeros_threenslip(Nslip,3))
  allocate(tau_back1(Nel,Nslip,Nchar), tau_back2(Nel,Nslip,Nchar), rho_total(Nel,Nslip,Nchar))
  allocate(Sb(3,Nslip,Nchar), Bb(3,Nslip,Nchar), tau_back_nokappa(3,Nslip,Nchar))
  allocate(vel(Nslip,Nchar), tline(Nslip,3), zeros2(resolution,Nslip))
  
  call props_dmb%allocatememory()
  call dmb_state%allocatememory()
  call props_cdt%allocatememory()
  call ddc_props%allocatememory()
  
  zeros = 0.d0
  zeros2 = 0.d0
  zeros_nslip = 0.d0
  zeros_threenslip = 0.d0
  euler = [0d0, 0.d0, 0.d0]
  unit_y=[0.d0,1.d0,0.d0]
  unit_z=[0.d0,0.d0,1.d0]
  rho0 = 10000000.0d0

  print*,"running tests for program version: ",prog_version
  print*,"compiled with these defaults:"
  print*,"# of elements, Nel = ",Nel
  print*,"# of regions, Nregion = ",Nregion
  print*,"# of slip systems, Nslip = ",Nslip
  print*,"# of dislocation character angles, Nchar = ",Nchar
  print*,"------------------------------------------------------------"
  
  call generate_regions() ! initialize region_mask
  call testzero(dble(sum(region_mask(1,:)+region_mask(Nregion,:))-(2*(1+Nel))),"region_mask",0.d0,count_pass,count_fail)
  
  ! TestVoigt, initialize A6 with random numbers, then symmetrize before running the test
  CALL RANDOM_NUMBER(A6)
  do i=1,6
    do j=1,i
      A6(j,i) = A6(i,j)
    end do
  end do
  A3 = .unvoigt.A6
  B6 = .voigt.A3
  call testequal(A6,B6,6,6,"Voigt",1.d-15,count_pass,count_fail)
  
  ! test Det3x3
  eye = Identity(3)
  num = .det. eye
  call testzero(num-1.d0,"Det3x3",0.d0,count_pass,count_fail)
  
  ! test DetInv3x3
  CALL RANDOM_NUMBER(A)
  call DetInv3x3(A,num,B)
  B = MATMUL(A,B)
  if (abs(num - .det. A)<1.d-99) then
    call testequal(eye,B,3,3,"DetInv3x3",1.d-12,count_pass,count_fail)
  else
    count_fail = count_fail+1
    print*,"test DetInv3x3: "//FAILED
  end if
  
  ! test DMB_Invert_Tensor26
  call DMB_Invert_Tensor26(A6,B6)
  eye6 = Identity(6)
  B6 = MATMUL(A6,B6)
  call testequal(eye6,B6,6,6,"DMB_Invert_Tensor26",1.d-12,count_pass,count_fail)
  
  
  ! test DMB_Invert_Tensor43
  A3 = .unvoigt.A6
  call DMB_Invert_Tensor43(A3,B3)
  B6 = .voigt.B3
  mat2 = Identity(6)
  do i=4,6
    mat2(i,i) = 2.d0
  end do
  B6 = MATMUL(A6,MATMUL(mat2,MATMUL(B6,mat2)))
  call testequal(eye6,B6,6,6,"DMB_Invert_Tensor43",1.d-12,count_pass,count_fail)
  
  ! test linspace
  call LinSpace(0.2d0,1.d0,5,array1)
  array2 = [0.2d0,0.4d0,0.6d0,0.8d0,1.d0]
  call testequalarray(array1,array2,5,"LinSpace",1.d-15,count_pass,count_fail)
  
  ! test CumTrapz / Trapz
  call LinSpace(1.d-15, 100.d0, 10000, x)
  func(:) = sin(pi*x(:)/25.d0) + x(:)**3 / (exp(x(:))-1.d0)
  call CumTrapz(10000,func,x,Integral)
  call testzero(pi**4/15.d0 - integral(9999),"CumTrapz",1.d-6,count_pass,count_fail)
  call Trapz(10000,func,x,tmpintegral)
  call testzero(tmpintegral - integral(9999),"Trapz",1.d-13,count_pass,count_fail)
  
  ! test Invert_Tensor
  CALL RANDOM_NUMBER(A)
  call Invert_Tensor(A,3,B)
  call testequal(eye,matmul(A,B),3,3,"Invert_Tensor",1.d-12,count_pass,count_fail)
  
  ! test ExpM
  A=0.d0
  A = ExpM(A,order=6)
  B = ExpM(1.123d0*eye,order=6)
  call testequal(exp(1.123d0)*eye,matmul(A,B),3,3,"ExpM, trivial cases",1.d-13,count_pass,count_fail)
  call RANDOM_NUMBER(A)
  B = ExpM(1.d-6*A)
  A = ExpM(-1.d-6*A)
  call testequal(A+B,2*eye,3,3,"ExpM, small argument",1.d-10,count_pass,count_fail)
  
  ! test rotation matrix with euler angles
  call DMB_Compute_Rotation(euler,Q)
  call testequal(Q,eye,3,3,"DMB_Compute_Rotation, eulers=0",1.d-15,count_pass,count_fail)
  call DMB_Compute_Rotation([pi/4.d0,0.d0,0.d0],Q)
  
  ! test interpolation and find interval:
!   print*,array1
!   print*,sin(array1)
  call FindInterval(array1, 5, 0.53d0, i)
  call Interp1D(0.3d0, 5, array1, sin(array1), Fwant)
  call testzero(0.5d0*(sin(array1(1))+sin(array1(2)))-Fwant + 2 - i,"Interp1D, FindInterval",1.d-15,count_pass,count_fail)
  
  !*********************************************************************
  ! now run some more advanced tests
  !*********************************************************************
  call DMB_get_slip_system('fcc', props_dmb%slip_M, props_dmb%slip_S, props_dmb%slip_V, props_dmb%line_T, props_dmb%rotM, &
                            props_cdt%Aforest)
!   do i=1,Nslip
!     print*,props_cdt%Aforest(9,i,:)
!   end do
!   print*,sum(props_cdt%Aforest(9,:,1))/12.d0
!   print*,sum(props_cdt%Aforest(9,:,2))/12.d0

  ! test RotAlign
  do i=1,min(Nslip,12)
  rot = props_dmb%rotM(:,:,i,1)
  tline(i,:) = matmul(rot,props_dmb%slip_M(:,i))+matmul(rot,props_dmb%line_T(:,i,1))
!   call testequalarray(tline(i,:),unit_y+unit_z,3,"RotAlign, screw", 1.d-15,count_pass,count_fail)
  tline(i,:) = tline(i,:) - (unit_y+unit_z)
  end do
  if (Nchar>1) then
    call testequal(tline,zeros_threenslip,Nslip,3,"RotAlign, screw", 1.d-15,count_pass,count_fail)
  end if
  do i=1,min(Nslip,12)
  rot = props_dmb%rotM(:,:,i,Nchar)
  tline(i,:) = matmul(rot,props_dmb%slip_M(:,i))+matmul(rot,props_dmb%line_T(:,i,Nchar))
!   call testequalarray(tline(i,:),unit_y+unit_z,3,"RotAlign, edge", 1.d-15,count_pass,count_fail)
  tline(i,:) = tline(i,:) - (unit_y+unit_z)
  end do
  call testequal(tline,zeros_threenslip,Nslip,3,"RotAlign, edge", 1.d-15,count_pass,count_fail)
  
  ! test Cu
  props_dmb%C11     = 168300.0d0    
  props_dmb%C12     = 121200.0d0
  props_dmb%C44     = 75700.0d0
  props_dmb%elastic_constants = 0.d0
  props_dmb%elastic_constants(1,1) = props_dmb%C11
  props_dmb%elastic_constants(2,1) = props_dmb%C12
  props_dmb%elastic_constants(3,1) = props_dmb%C44
  props_dmb%burger  = 3d-07
  props_dmb%rhobar0 = 8.96d-09
  props_dmb%wave_vel(1) = 2039915.1008740345d0
  props_dmb%wave_vel(Nchar) = 1605654.0723331412d0
  vel = 1500.d3
    
  call stroh_dislocation(props_dmb,vel,resolution,uij,Sb, Bb,elC)
  if (abs(sum(Sb(:,9,1))) < 1.d-15 .or. Nchar==1) then
    call testzero(Sb(1,9,Nchar)+Sb(2,9,Nchar),"Sb",1.d-15,count_pass,count_fail)
  else
    count_fail = count_fail+1
    print*,"test Sb: "//FAILED
  end if
  call testzero((Bb(1,9,Nchar)-Bb(2,9,Nchar)+Bb(3,9,Nchar))/props_dmb%C44,"Bb",1.d-15,count_pass,count_fail)
!   do i=1,Nslip
!   print*,"testing trace uij_screw for slip system",i
  if (Nchar > 1) then
    call testequalarray(uij(1,1,:,:,1)+uij(2,2,:,:,1)+uij(3,3,:,:,1),zeros2,resolution,"trace uij_screw",&
                        1.d-15,count_pass,count_fail)
  end if
  
  vel = 0.d0
  call stroh_dislocation(props_dmb,vel,resolution,uij,Sb, Bb,elC)
!   print*,Sb(:,9,1)
!   print*,Sb(:,9,Nchar)
!   print*,Bb(:,9,1)/props_dmb%C44
!   print*,Bb(:,9,Nchar)/props_dmb%C44
  call testequalarray(Sb(:,9,Nchar),[-0.01656451d0,  0.01656451d0, -0.02989161d0],3, &
                      "Sb (fcc Cu): validate against python implementation",1.d-8,count_pass,count_fail)
  if (Nchar > 1) then
    call testequalarray(Bb(:,9,1)/props_dmb%C44,[0.06277008d0,  0.06277008d0, 0.d0],3, &
                      "Bb (fcc Cu) screw: validate against python implementation",1.d-8,count_pass,count_fail)
  end if
  call testequalarray(Bb(:,9,Nchar)/props_dmb%C44,[0.110293506d0,  0.110293506d0, 0.d0],3, &
                      "Bb (fcc Cu) edge: validate against python implementation",1.d-8,count_pass,count_fail)

  ! test form_elasticity
  call DMB_form_elasticity(Cg, props_dmb%rotM(:,:,3,Nchar), props_dmb%elastic_constants, 0.d0, 0.d0)
  A6 = .voigt.Cg
!   print*,A6
  call testzero(sum(A6(:3,:3))-3*props_dmb%C11-6*props_dmb%C12,"form_elasticity",1.d-6,count_pass,count_fail)
  
  !test isotropic limit
  nu_poisson = 1.d0/3.d0
!   nu_poisson = 0.35d0
  props_dmb%C44     = 1.0d0
  props_dmb%C12     = props_dmb%C44*2.d0*nu_poisson/(1.d0-2.d0*nu_poisson)
  props_dmb%C11     = props_dmb%C12+2.d0*props_dmb%C44
  props_dmb%rhobar0 = 1.d0
  props_dmb%wave_vel = 1.d0
  vel = 0.d0
  
  call stroh_dislocation(props_dmb,vel,resolution,uij,Sb, Bb,elC)
  
!   print*,Sb(:,9,1)
!   print*,Sb(:,9,Nchar)
!   print*,Bb(:,9,1)/props_dmb%C44
!   print*,Bb(:,9,Nchar)/props_dmb%C44
  if ((abs(Sb(1,9,Nchar) + props_dmb%C44/(props_dmb%C12+2.d0*props_dmb%C44)/sqrt(3.d0)/2.d0/pi) < 1.d-15)  .and. &
     (abs(Bb(1,9,Nchar)/props_dmb%C44 - (props_dmb%C12+props_dmb%C44)/(props_dmb%C12+2.d0*props_dmb%C44)/pi/sqrt(2.d0)) &
                    < 1.d-15)) then
    if ((Nchar<2) .OR. (abs(Bb(1,9,1)/props_dmb%C44 - 0.5d0/(pi*sqrt(2.d0))) < 1.d-15)) then
      count_pass = count_pass+1
      print*,"test Sb, Bb isotropic: "//PASSED
    else
      count_fail = count_fail+1
      print*,"test Bb isotropic: "//FAILED
    end if
  else
    count_fail = count_fail+1
    print*,"test Sb, Bb isotropic: "//FAILED
  end if
  strainenergydensity = 0.d0
  do i=1,3
    do j=1,3
      do k=1,3
        do l=1,3
          strainenergydensity(:) = strainenergydensity(:) + 0.5d0*elC(i,j,k,l)*uij(i,j,:,9,Nchar)*uij(k,l,:,9,Nchar)
        end do
      end do
    end do
  end do
  call LinSpace(0.d0,2.d0*pi,resolution,phi)
  call Trapz(resolution, strainenergydensity, phi, tmpintegral)
  call testzero(tmpintegral-props_dmb%C44/(4.d0*pi*(1-nu_poisson)),&
                      "isotropic strainenergy",1.d-15,count_pass,count_fail)
  
                      
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
  ! test new backstress
  ddc_props%slip_S = props_dmb%slip_S
  ddc_props%Lrho = 1.d0 / rho0
  ddc_props%mu = 1.d0
  ddc_props%burger = props_dmb%burger
  ddc_props%euler_angle = euler
  props_dmb%Lbar = 1.d0 / sqrt(rho0)

  call DDC_Zaiser_Backstress_noKappa(vel, ddc_props, props_dmb, tau_back_nokappa)
  do i=1,Nchar
    do j=1,Nslip
      ! undo rotation made in DDC_Zaiser_Backstress_noKappa(), i.e. align with dislocation
      tau_back_nokappa(:,j,i) = -MATMUL(props_dmb%rotM(:,:,j,i),tau_back_nokappa(:,j,i))/(props_dmb%burger*(props_dmb%Lbar**2))
    end do
  end do
!   print*,tau_back_nokappa(:,9,1)/(rho0*props_dmb%C44),tau_back_nokappa(:,9,Nchar)/(rho0*props_dmb%C44)
!   print*,MATMUL(transpose(props_dmb%rotM(:,:,9,Nchar)),[1.d0,0.d0,0.d0])
!   print*,props_dmb%slip_S(:,9)
!   print*,props_dmb%slip_V(:,9,2)
!   print*,MATMUL(transpose(props_dmb%rotM(:,:,9,1)),[1.d0,0.d0,0.d0])
!   print*,props_dmb%slip_V(:,9,1)
!   print*,MATMUL(transpose(props_dmb%rotM(:,:,9,Nchar)),[0.d0,1.d0,0.d0])
!   print*,MATMUL(transpose(props_dmb%rotM(:,:,9,1)),[0.d0,1.d0,0.d0])
!   print*,props_dmb%slip_M(:,9)
  if (Nchar>1) then
    call testequal(0.75d0*tau_back_nokappa(:,:,1)/(rho0*props_dmb%C44),tau_back_nokappa(:,:,Nchar)/(rho0*props_dmb%C44), &
                      3, Nslip, "isotropic backstress, screw vs edge",1.d-12,count_pass,count_fail)
  end if

  varrho(:,:,:,:) = 0.5d0*rho0/(Nslip*Nchar)
  varrho(123,:,:,2) = 0.4d0*rho0/(Nslip*Nchar)
  varrho(125,:,:,2) = 0.6d0*rho0/(Nslip*Nchar)
  rho_total(:,:,:) = varrho(:,:,:,1) + varrho(:,:,:,2)
  dmb_state%vel = 1.d0
  call LinSpace(0.d0, 1.d0, Nnode, X_ref)
  call DDC_Gradient_Backstress(Nel,varrho, ddc_props, X_ref, tau_back1)
!   print*,tau_back1(124,9,1)
!   print*,tau_back1(124,9,2)
  ddc_props%backstress_model = 'zaiser_zero' ! ignore dmb_state%vel and compute for vel=0 only:
  call DDC_Zaiser_Backstress(Nel,varrho, ddc_props, props_dmb, dmb_state%vel, X_ref, tau_back2)
!   print*,tau_back2(124,9,1)
!   print*,tau_back2(124,9,2)
  if (Nchar>1) then
    call testzero(tau_back1(124,9,1)*0.5d0*props_dmb%slip_V(1,9,1)/props_dmb%slip_S(1,9) - &
      tau_back2(124,9,1)*rho_total(124,9,1)/rho0 + &
      tau_back1(124,9,Nchar)*0.5d0*0.75d0 - tau_back2(124,9,Nchar)*rho_total(124,9,Nchar)/rho0, &
      "compare gradient and zaiser back stress",1.d-15,count_pass,count_fail)
  else
    call testzero(tau_back1(124,9,Nchar)*0.5d0*0.75d0 - tau_back2(124,9,Nchar)*rho_total(124,9,Nchar)/rho0, &
      "compare gradient and zaiser back stress",1.d-15,count_pass,count_fail)
  end if
  
  ! test bcc:
  call DMB_get_slip_system('bcc', props_dmb%slip_M, props_dmb%slip_S, props_dmb%slip_V, props_dmb%line_T, props_dmb%rotM, &
                            props_cdt%Aforest)
  ! test Fe
  props_dmb%C11     = 226000.0d0    
  props_dmb%C12     = 140000.0d0
  props_dmb%C44     = 116000.0d0
  props_dmb%elastic_constants = 0.d0
  props_dmb%elastic_constants(1,1) = props_dmb%C11
  props_dmb%elastic_constants(2,1) = props_dmb%C12
  props_dmb%elastic_constants(3,1) = props_dmb%C44
  props_dmb%burger  = 0.d0
  props_dmb%rhobar0 = 0.d0
  props_dmb%wave_vel(1) = 0.d0
  props_dmb%wave_vel(Nchar) = 0.d0
  vel = 0.d0
  call stroh_dislocation(props_dmb,vel,resolution,uij,Sb, Bb,elC)
  
  call testequalarray(Sb(:,9,Nchar),[-3.24487029d-02,-3.24487029d-02, 0.d0],3, &
                      "Sb (bcc Fe): validate against python implementation",1.d-8,count_pass,count_fail)
  if (Nchar > 1) then
    call testequalarray(Bb(:,9,1)/props_dmb%C44,[0.04651684d0, -0.04651684d0,  0.04651684d0],3, &
                      "Bb (bcc Fe) screw: validate against python implementation",1.d-8,count_pass,count_fail)
  end if
  call testequalarray(Bb(:,9,Nchar)/props_dmb%C44,[0.08475254d0, -0.08475254d0,  0.10570937d0],3, &
                      "Bb (bcc Fe) edge: validate against python implementation",1.d-8,count_pass,count_fail)
  
  
  call system_clock(finish_time)
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
  print*,"------------------------------------------------------------"
  print*,"SUMMARY:", count_pass," passed and ",count_fail," failed"
  print*,"time: ",int(1000.d0*real(finish_time-start_time)/real(countrate)), "ms"
  
  ! check this last line is working:
  count_fail = count_fail+1
  if (count_fail>0) error stop "some tests failed"
END PROGRAM RunTests
