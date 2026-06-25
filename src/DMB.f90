module DMB_routines
private
interface operator(.voigt.)
  module procedure Voigt43to26
end interface
interface operator(.unvoigt.)
  module procedure Voigt26to43
end interface
public :: DMB_Init_EOS, DMB_Compute_EOS, DMB_Compute_Rotation, DMB_Invert_Tensor26, DMB_Invert_Tensor43, DMB_State_Init, &
          DMB_Update_State, DMB_form_elasticity, DMB_get_slip_system, Voigt26to43, Voigt43to26, stroh_dislocation, &
          operator(.voigt.), operator(.unvoigt.)
contains
SUBROUTINE DMB_Init_EOS(Dhat)
  use utilities, only: LinSpace, CumTrapz
  IMPLICIT NONE

!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(7999,2), INTENT(OUT) :: Dhat
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(8000) :: x, func!, bigI
  REAL(KIND=8), DIMENSION(7999) :: Integral
!-----------------------------------------------------------------------
  call LinSpace(0.01d0, 20.d0, 8000, X)

  func(:) = x(:)**3 / (exp(x(:))-1.d0)

  call CumTrapz(8000,func,x,Integral)

  Dhat(1:7999,1) = x(2:8000)
  Dhat(1:7999,2) = Integral(1:7999)*3.d0/x(2:8000)**3

END SUBROUTINE DMB_Init_EOS

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DMB_Compute_EOS(eos_props, detF, temperature, &
           PRESSURE, FREE_ENERGY, ENTROPY, THERMO_COEFFS)
!-----------------------------------------------------------------------
  use GlobalParams
  USE DMB
  use utilities, only: Interp1D
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_EOS)           :: eos_props
  REAL(KIND=8)                       :: detF
  REAL(KIND=8)                       :: temperature
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(4) :: pressure
  REAL(KIND=8), DIMENSION(4) :: free_energy
  REAL(KIND=8), DIMENSION(3) :: entropy
  REAL(KIND=8), DIMENSION(3) :: thermo_coeffs

  INTENT(IN) eos_props, detF, temperature
  INTENT(OUT) pressure, free_energy, entropy, thermo_coeffs
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: third, rho0, v_star, Bstar, B1, Td0, a, b, gamma, Gamma0
  REAL(KIND=8) :: kappa, R_M, v0
  REAL(KIND=8), DIMENSION(:,:), allocatable :: Dhat_array
  REAL(KIND=8) :: v, X, expX, grun, p_static, p_elect, temp_debye, temp_ratio
  REAL(KIND=8) :: exp_n, exp_p, exp_pm1, exp_nm1, exp_n1m
  REAL(KIND=8) :: Dhat, p_ion, s_elect, s_ion, psi_static
  REAL(KIND=8) :: psi_ion, psi_elect, term1, term2, term3, specific_heat
  REAL(KIND=8) :: te_coupling, bulk_mod, p_total, psi_total, s_total
!-----------------------------------------------------------------------
  allocate(Dhat_array(7999,2))
  third = 1.d0 / 3.d0
  
  ! assign properties       
  rho0   = eos_props%rho0
  v_star = eos_props%v_star   
  Bstar  = eos_props%B_star
  B1     = eos_props%B1    
  Td0    = eos_props%Td0  
  a      = eos_props%a        
  b      = eos_props%b     
  gamma  = eos_props%gamma 
  Gamma0 = eos_props%Gamma0
  kappa  = eos_props%kappa 
  R_M    = eos_props%R_M   
  v0     = eos_props%v0

  Dhat_array = eos_props%Dhat
  
  v      = detF / rho0
  X      = 1.5d0*(B1-1.d0)*( (v/v_star)**(third) - 1.d0)
  expX   = exp(-X)
  grun   = gamma/v + a/v0 + b*v/(v0**2)
        
  ! compute pressures
  p_static = -3.d0*Bstar*((v_star/v)**third - (v_star/v)**(2.d0*third))*expX
  p_elect  =  0.5d0*(kappa/v)*Gamma0*(v/v0)**kappa*temperature**2
        
  temp_debye = Td0*(v/v0)**(-gamma)*exp(-a*((v/v0) - 1.d0) - 0.5d0*b*( (v/v0)**2 - 1.d0) )
  temp_ratio = temp_debye / temperature
        
  exp_n   = exp(-temp_ratio)
  exp_p   = exp( temp_ratio)
  exp_pm1 = exp_p - 1.d0
  exp_nm1 = exp_n - 1.d0
  exp_n1m = 1.d0 - exp_n
  
  call Interp1D(temp_ratio, 7999, Dhat_array(:,1), Dhat_array(:,2), Dhat)
  
  p_ion = 3.d0*R_M*grun*temp_debye*(3.d0/8.d0 + exp_n / exp_n1m - 1.d0 / exp_pm1 + Dhat/temp_ratio)
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,*) "rho0=",rho0       
    write(*,*) "v_star=",v_star       
    write(*,*) "Bstar=",Bstar       
    write(*,*) "B1=",B1       
    write(*,*) "Td0=",Td0       
    write(*,*) "a=",a       
    write(*,*) "b=",b       
    write(*,*) "gamma=",gamma       
    write(*,*) "Gamma0=",Gamma0       
    write(*,*) "kappa=",kappa       
    write(*,*) "R_M=",R_M       
    write(*,*) "v0=",v0       
    write(*,*) "v=",v       
    write(*,*) "X=",X       
    write(*,*) "expX=",expX       
    write(*,*) "p_static=",p_static       
    write(*,*) "p_elect=",p_elect       
    write(*,*) "grun=",grun       
    write(*,*) "temp_debye=",temp_debye       
    write(*,*) "exp_n=",exp_n       
    write(*,*) "exp_p=",exp_p       
    write(*,*) "Dhat=",Dhat       
    write(*,*) "temp_ratio=",temp_ratio       
  end if

  ! compute entropies
  s_elect = Gamma0*temperature*(v/v0)**kappa
  s_ion   = -R_M*(3.d0*log(exp_n1m) + 3.d0*temp_ratio*( 1.d0/exp_pm1 - exp_n/exp_n1m ) - 4.d0*Dhat)
    
  !compute free energies
  psi_static = 4.d0*v_star*Bstar/((B1-1.d0)**2) * (1.d0-(1.d0+X)*expX )
  psi_ion    = R_M*( (9.d0/8.d0)*temp_debye + 3.d0*temperature*log(exp_n1m)-temperature*Dhat)
  psi_elect  = -0.5d0*Gamma0*((v/v0)**kappa)*temperature**2
      
  !compute specific heat
  term1 = Gamma0*temperature*(v/v0)**kappa
  term2 = temp_ratio**2 * exp_n/(exp_nm1**2) 
  term3 = temp_ratio*(1.d0/exp_pm1)*(3.d0 + temp_ratio*(exp_p/exp_pm1 ) )
  specific_heat = term1 + 3.d0*R_M*(term2 - term3 + 4.d0*Dhat)
        
  ! compute thermoelastic coupling coefficient
  term1 = 3.d0*R_M*grun*temp_ratio**2 *( exp_n/(exp_nm1**2) &
           - exp_p/(exp_pm1**2) - 3.d0 / exp_pm1 / temp_ratio + 4.d0*Dhat/( (temp_ratio)**2) )
  term2 = (kappa/v)*Gamma0*temperature*(v/v0)**kappa
  te_coupling = -(term1+term2)
        
  ! compute isothermal bulk modulus
  term1 = -Bstar*expX/v*(                                                &
        (v_star/v)**third *( 0.5d0*B1 - 1.5d0 + 2.d0*(v_star/v)**third)  &
          - 0.5d0*(B1 - 1.d0)      )
        
  term2 = -3.d0*R_M*temp_debye*grun**2 * (                                        &
        (3.d0/8.d0 - exp_n/(exp_nm1)*(1.d0 + temp_ratio/(exp_nm1) ) +   &
        (1.d0 / (exp_pm1))*(2.d0 + temp_ratio*exp_p/(exp_pm1))           &
        -3.d0*Dhat/temp_ratio) )
   
  term3 = kappa*(kappa - 1.d0)/(2.d0*v**2)*Gamma0*temperature**2*(v/v0)**kappa
        
  bulk_mod = -v*(term1+term2+term3)
        
  ! populate return arrays
  p_total          = p_static + p_ion + p_elect 
  pressure(:)      = (/ p_total, p_static, p_ion, p_elect /)
        
  psi_total        = psi_static + psi_ion + psi_elect
  free_energy(:)   = (/ psi_total, psi_static, psi_ion, psi_elect /)
        
  s_total          = s_ion + s_elect
  entropy(:)       = (/ s_total, s_ion, s_elect /)
  
  thermo_coeffs(:) = (/ specific_heat, te_coupling, bulk_mod /)

!  write(*,*) "--------------------"
!  write(*,10) "pressure      = ",pressure(1:4)
!  write(*,10) "free_energy   = ",free_energy(1:4)
!  write(*,10) "entropy       = ",entropy(1:3)
!  write(*,10) "thermo_coeffs = ",thermo_coeffs(1:3)
!10 format(a,5es21.13)

END SUBROUTINE DMB_Compute_EOS

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DMB_State(times, state)
  USE GlobalParams
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3),     INTENT(IN)    :: times
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_DMB_STATE),     INTENT(INOUT)   :: state
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i
!-----------------------------------------------------------------------
!~   call state%allocatememory()

  state%time(:)        = times
  state%X_ref(:)       = 0.d0
  state%x(:,:)         = 0.d0
  state%xdot(:,:)      = 0.d0
  state%xddot(:,:)     = 0.d0
  state%F(:,:,:)       = 0.d0
  state%T(:,1)         = 0.d0
  state%vel(:,:,:)       = 0.d0
  state%dis_acc(:,:,:)  = 0.d0
  state%vel_x(:,:,:)     = 0.d0
  state%vel_y(:,:,:)     = 0.d0
  state%vel_z(:,:,:)     = 0.d0
  state%cauchy(:,:,:)  = 0.d0
  state%Fp(:,:,:)      = 0.d0
  state%Fe(:,:,:)      = 0.d0
  state%tau(:,:)       = 0.d0
  state%gamma_acc(:,:,:) = 0.d0
  state%Nsub(:)        = 0
  state%i_cutback(:)   = 0

  do i=1,3
    state%F(:,i,i) = 1.d0
  end do  !i

END SUBROUTINE DMB_State


!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DMB_State_Init(times, X_ref, x, xdot, xddot, Fp, Fe, T, &
           DMB_State)
  USE GlobalParams
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3),       INTENT(IN)    :: times
  REAL(KIND=8), DIMENSION(Nnode),   INTENT(IN)    :: X_ref
  REAL(KIND=8), DIMENSION(Nnode,3), INTENT(IN)    :: x
  REAL(KIND=8), DIMENSION(Nnode,3), INTENT(IN)    :: xdot
  REAL(KIND=8), DIMENSION(Nnode,3), INTENT(IN)    :: xddot
  REAL(KIND=8), DIMENSION(3,3),     INTENT(IN)    :: Fp
  REAL(KIND=8), DIMENSION(3,3),     INTENT(IN)    :: Fe
  REAL(KIND=8),                     INTENT(IN)    :: T
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_DMB_STATE),       INTENT(OUT)   :: DMB_State
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i
!-----------------------------------------------------------------------
  call DMB_State%allocatememory()
  DMB_State%time(:)       = times
  DMB_State%X_ref(:)      = X_ref
  DMB_State%x(:,:)        = x
  DMB_State%xdot(:,:)     = xdot
  DMB_State%xddot(:,:)    = xddot
  DMB_State%F(:,:,:)      = 0.d0
  DMB_State%T(:,1)        = T
  DMB_State%gamma_acc(:,:,:)= 0.d0
  DMB_State%vel(:,:,:)      = 0.d0
  DMB_State%dis_acc(:,:,:)  = 0.d0
  DMB_State%vel_x(:,:,:)    = 0.d0
  DMB_State%vel_y(:,:,:)    = 0.d0
  DMB_State%vel_z(:,:,:)    = 0.d0
  DMB_State%cauchy(:,:,:) = 0.d0
  DMB_State%Fp(:,:,:)     = 0.d0
  DMB_State%Fe(:,:,:)     = 0.d0
  DMB_State%Nsub(:) = 0
  DMB_State%i_cutback = 0

  do i=1,3
    DMB_State%F(:,i,i) = 1.d0
  end do  !i

  do i=1,Nel
    DMB_State%Fp(i,:,:) = Fp(:,:)
    DMB_State%Fe(i,:,:) = Fe(:,:)
  end do  !i

END SUBROUTINE DMB_State_Init

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!> compute "active" rotation matrix using Bunge euler angles
!> "active" rotation matrix is defined such that E_i = R_ij * e_j where 
!> e_i are the crystal base vectors and E_j are the global or "sample" base vectors
SUBROUTINE DMB_Compute_Rotation(euler_angle, rot)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3),           INTENT(IN)    :: euler_angle
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3),         INTENT(OUT)   :: rot
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) phi1, Phi, phi2, cPhi, sPhi, c1, s1, c2, s2
!-----------------------------------------------------------------------

  phi1 = euler_angle(1)
  Phi  = euler_angle(2)
  phi2 = euler_angle(3)
  
  cPhi = cos(Phi)
  sPhi = sin(Phi)

  c1   = cos(phi1)
  s1   = sin(phi1)

  c2   = cos(phi2)
  s2   = sin(phi2)

  ! djl: transposed the rotation matrix...
  ! makes no difference for euler_angle = (/0, 0, 0/)
  rot(1, 1:3) = (/ c1*c2 - s1*s2*cPhi,  -c1*s2 - s1*c2*cPhi,   s1*sPhi /)
  rot(2, 1:3) = (/ s1*c2 + c1*s2*cPhi,  -s1*s2 + c1*c2*cPhi,  -c1*sPhi /)
  rot(3, 1:3) = (/            s2*sPhi,              c2*sPhi,      cPhi /)  
  
END SUBROUTINE DMB_Compute_Rotation

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!> form a fourth order elasticity tensor from elastic constants (in crystal basis) 
!> and rotate into global basis
subroutine DMB_form_elasticity(Cg, rotMat, elastic_constants, temperature, pressure)

implicit none

REAL(KIND=8), intent(in   ) :: rotMat(3,3)
REAL(KIND=8), intent(in   ) :: elastic_constants(3,3)
REAL(KIND=8), intent(in   ) :: temperature
REAL(KIND=8), intent(in   ) :: pressure

REAL(KIND=8), intent(  out) :: Cg(3,3,3,3)

REAL(KIND=8), parameter :: zero=0.d0
REAL(KIND=8)            :: Cc(3,3,3,3), C0(3), CT(3), CP(3)
REAL(KIND=8)            :: T_, P_

integer           :: q,r,s,t !, i,j,k,l

C0 = elastic_constants(:,1)
CT = elastic_constants(:,2)
CP = elastic_constants(:,3)

T_ = temperature
p_ = pressure

Cc=zero; Cg = zero
Cc(1,1,1,1) = C0(1) + CT(1)*T_ + CP(1)*P_ !! C11
Cc(2,2,2,2) = C0(1) + CT(1)*T_ + CP(1)*P_ !! C22
Cc(3,3,3,3) = C0(1) + CT(1)*T_ + CP(1)*P_ !! C33

Cc(1,1,2,2) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C12
Cc(1,1,3,3) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C13
Cc(2,2,3,3) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C23

Cc(2,2,1,1) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C12
Cc(3,3,1,1) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C13
Cc(3,3,2,2) = C0(2) + CT(2)*T_ + CP(2)*P_ !! C23

Cc(2,3,2,3) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C44
Cc(2,3,3,2) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C44
Cc(3,2,2,3) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C44
Cc(3,2,3,2) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C44

Cc(3,1,3,1) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C55
Cc(3,1,1,3) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C55
Cc(1,3,3,1) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C55
Cc(1,3,1,3) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C55

Cc(2,1,2,1) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C66
Cc(2,1,1,2) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C66
Cc(1,2,2,1) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C66
Cc(1,2,1,2) = C0(3) + CT(3)*T_ + CP(3)*P_ !! C66    

!~ do q=1,3; do r=1,3; do s=1,3; do t=1,3
!~     do i=1,3; do j=1,3; do k=1,3; do l=1,3
!~         Cg(q,r,s,t) = Cg(q,r,s,t) + rotMat(q,i)*rotMat(r,j)*rotMat(s,k)*rotMat(t,l)*Cc(i,j,k,l)
!~     end do; end do; end do; end do
!~ end do; end do; end do; end do      

! faster than the nested loops above:
do s=1,3; do r=1,3; do q=1,3
  Cg(q,r,s,:) = MATMUL(rotMat(:,:),Cc(q,r,s,:))
end do; end do; end do
do t=1,3; do r=1,3; do q=1,3
  Cg(q,r,:,t) = MATMUL(rotMat(:,:),Cg(q,r,:,t))
end do; end do; end do
do t=1,3; do s=1,3; do q=1,3
  Cg(q,:,s,t) = MATMUL(rotMat(:,:),Cg(q,:,s,t))
end do; end do; end do
do t=1,3; do s=1,3; do r=1,3
  Cg(:,r,s,t) = MATMUL(rotMat(:,:),Cg(:,r,s,t))
end do; end do; end do

end subroutine DMB_form_elasticity



!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Convert a 3x3x3x3 symmetric 4th order tensor to a 6x6 2nd order tensor
pure function Voigt43to26(A4) result(A2)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3,3,3), INTENT(IN)  :: A4
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(6,6) :: A2
!-----------------------------------------------------------------------
!  Locals:
  INTEGER, DIMENSION(6) :: iind,jind
  INTEGER :: iv,jv
!-----------------------------------------------------------------------
  iind = (/1, 2, 3, 2, 3, 1/)
  jind = (/1, 2, 3, 3, 1, 2/)
  do iv=1,6
    do jv=1,6
      a2(iv,jv) = a4(iind(iv),jind(iv),iind(jv),jind(jv)) 
    end do  !jv
  end do  !iv

end function Voigt43to26

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Convert a 6x6 2nd order tensor to a 3x3x3x3 symmetric 4th order tensor
pure function Voigt26to43(A2) result(A4)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(6,6),     INTENT(IN) :: A2
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3,3,3)  :: A4
!-----------------------------------------------------------------------
!  Locals:
!  INTEGER, DIMENSION(6) :: iind,jind
  INTEGER :: iv,jv, i, j, k, l
!-----------------------------------------------------------------------
  ! djl: the following does not populate entire 3x3x3x3
  ! iind = (/1, 2, 3, 2, 3, 1/)
  ! jind = (/1, 2, 3, 3, 1, 2/)
  ! do iv=1,6
  !   do jv=1,6
  !     a4(iind(iv),jind(iv),iind(jv),jind(jv)) = a2(iv,jv)
  !   end do  !jv
  ! end do  !iv
  
  ! djl: add block below to fully populate 4th order tensor
  do i = 1,3
    do j = 1,3
      if (i==j) then
        iv = i
      else
        iv = 9-(i+j)
      end if 
      
      do k = 1,3
        do l = 1,3
          if (k==l) then
            jv = k
          else
            jv = 9-(k+l)
          end if
	  
          a4(i,j,k,l) = a2(iv,jv)
	  
        end do
      end do
    end do
  end do 

end function Voigt26to43

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Invert a 3x3x3x3 symmetric tensor.
SUBROUTINE DMB_Invert_Tensor43(A,Ainv)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3,3,3), INTENT(IN)  :: A
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3,3,3), INTENT(OUT) :: Ainv
!-----------------------------------------------------------------------
!  Locals:
!  REAL(KIND=8), DIMENSION(6)   :: Winv = (/1.d0,1.d0,1.d0,0.25d0,0.25d0,0.25d0/)
  REAL(KIND=8), DIMENSION(6,6) :: Amat,Aminv
!  INTEGER :: i,j
!-----------------------------------------------------------------------
! Convert to Voigt form
  Amat = .voigt. A

!  Aminv  = nmp.dot( Winv, nmp.dot(nmp.linalg.inv(Amat), Winv))
!
!    but Winv is a diagonal matrix, so I square the original Winv
!    values and scale the inverse columns by it.

! Invert the 6x6
  call DMB_Invert_Tensor26(Amat,Aminv)

  ! Scale the columns by Winv
  !  do j=1,6
  !   do i=1,6
  !     Aminv(i,j) = Aminv(i,j)*Winv(j)
  !   end do  !i
  !  end do  !j
  
  ! djl --> above doesn't scale terms correctly for off-diagonal blocks
  ! shouldn't make a difference for cubic symmetry...
  Aminv(1:3, 4:6) = 0.50d0*Aminv(1:3, 4:6)  ! top right block
  Aminv(4:6, 1:3) = 0.50d0*Aminv(4:6, 1:3)  ! bot left block
  Aminv(4:6, 4:6) = 0.25d0*Aminv(4:6, 4:6)  ! bot right block

! Convert from Voigt form
  Ainv = .unvoigt. Aminv

END SUBROUTINE DMB_Invert_Tensor43

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Invert a 6x6 tensor.
SUBROUTINE DMB_Invert_Tensor26(A,Ainv)
  use utilities, only: Fatal
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(6,6), INTENT(IN)  :: A
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(6,6), INTENT(OUT) :: Ainv
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(36) :: work
  INTEGER, DIMENSION(6) :: ipiv
  INTEGER :: n, nn
  INTEGER :: info
!-----------------------------------------------------------------------
  n = 6; nn=36
  Ainv = A

  call DGETRF(n,n,Ainv,n,ipiv,info)

  if(info/=0) then
    write(*,*) 'Invert_Tensor26(): Lapack DGETRF() failed.'
    call Fatal('Bad matrix.')
  end if

  call DGETRI(n,Ainv,n,ipiv,work,nn,info)

  if(info/=0) then
    write(*,*) 'Invert_Tensor26(): Lapack DGETRI() failed.'
    call Fatal('Bad matrix.')
  end if

END SUBROUTINE DMB_Invert_Tensor26

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Populate the slip system normals and directions.
SUBROUTINE DMB_get_slip_system(name, slip_M_out, slip_S_out, slip_V, line_T, rotM, Aforest)
  USE GlobalParams
!-----------------------------------------------------------------------
  use utilities, only: Fatal, operator(.cross.), RotAlign
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  CHARACTER*(*), INTENT(IN) :: name    !< Structure
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,Nslip), INTENT(OUT) :: slip_M_out, slip_S_out  !< Slip normals M, slip directions S
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar), INTENT(OUT) :: slip_V, line_T  !< velocity directions slip_V and line directions line_T
  REAL(KIND=8), DIMENSION(3,3,Nslip,Nchar), INTENT(OUT) :: rotM
  REAL(KIND=8), DIMENSION(Nslip,Nslip,Nchar), INTENT(OUT) :: Aforest !< slip_S is char independent, hence alpha=Nslip, beta=Nslip+Nchar
  !! rotM = rotation matrices which will align slip_M and line sense with Cartesian directions y and z

!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(3,max(48,Nslip)) :: slip_M, slip_S  !Slip normals M, slip directions S (local copy)
  REAL(KIND=8) s, m, O, x, y
  REAL(KIND=8), DIMENSION(3,Nslip) :: crossMB
  ! arrays of cos/sin of character angles:
  REAL(KIND=8), DIMENSION(Nchar) :: coschar,sinchar 
  INTEGER :: islip, ich, jslip
!-----------------------------------------------------------------------
  ! screw
  coschar(1) = 1.d0
  sinchar(1) = 0.d0
  ! edge
  coschar(Nchar) = 0.d0
  sinchar(Nchar) = 1.d0
  ! other characters
  do ich=2,Nchar-1
    ! ich-1 because fortran indices are 1-based; 1.d0 to force int to real conversion before dividing
    coschar(ich) = cos(0.5d0*pi*(ich-1.d0)/(Nchar-1.d0))
    sinchar(ich) = sin(0.5d0*pi*(ich-1.d0)/(Nchar-1.d0))
  end do !ich
  
  slip_M = 0.d0
  slip_S = 0.d0
  
  if(name == "fcc") then
    s = 1.d0 / sqrt(2.d0)
    m = 1.d0 / sqrt(3.d0)
    O = 0.d0
!--------------------
!  slip normals
!--------------------
    slip_M(:,1)  = (/m, m, m/)
    slip_M(:,2)  = (/m, m, m/)
    slip_M(:,3)  = (/m, m, m/)

    slip_M(:,4)  = (/-m, m, m/)
    slip_M(:,5)  = (/-m, m, m/)
    slip_M(:,6)  = (/-m, m, m/)
  
    slip_M(:,7)  = (/m, -m, m/)
    slip_M(:,8)  = (/m, -m, m/)
    slip_M(:,9)  = (/m, -m, m/)
  
    slip_M(:,10) = (/-m, -m, m/)
    slip_M(:,11) = (/-m, -m, m/)
    slip_M(:,12) = (/-m, -m, m/)
!--------------------
!  slip directions
!--------------------
    slip_S(:,1)  = (/  s, -s,  O/)
    slip_S(:,2)  = (/ -s,  O,  s/)
    slip_S(:,3)  = (/  O,  s, -s/)
  
    slip_S(:,4)  = (/  s,  O,  s/)
    slip_S(:,5)  = (/ -s, -s,  O/)
    slip_S(:,6)  = (/  O,  s, -s/)
  
    slip_S(:,7)  = (/ -s,  O,  s/)
    slip_S(:,8)  = (/  O, -s, -s/)
    slip_S(:,9)  = (/  s,  s,  O/)
  
    slip_S(:,10) = (/ -s,  s,  O/)
    slip_S(:,11) = (/  s,  O,  s/)
    slip_S(:,12) = (/  O, -s, -s/)
  
  elseif(name == "bcc") then
    s = 1.d0 / sqrt(2.d0)
    m = 1.d0 / sqrt(3.d0)
    O = 0.d0
    x = 1.d0 / sqrt(6.d0)
    y = 1.d0 / sqrt(14.d0)
!--------------------
!  slip normals
!--------------------
    ! {110} types:
    slip_M(:,1)  = (/  s, -s,  O/)
    slip_M(:,2)  = (/ -s,  O,  s/)
    slip_M(:,3)  = (/  O,  s, -s/)
  
    slip_M(:,4)  = (/  s,  O,  s/)
    slip_M(:,5)  = (/ -s, -s,  O/)
    slip_M(:,6)  = (/  O,  s, -s/)
  
    slip_M(:,7)  = (/ -s,  O,  s/)
    slip_M(:,8)  = (/  O, -s, -s/)
    slip_M(:,9)  = (/  s,  s,  O/)
  
    slip_M(:,10) = (/ -s,  s,  O/)
    slip_M(:,11) = (/  s,  O,  s/)
    slip_M(:,12) = (/  O, -s, -s/)
    
    ! {112} types:
    slip_M(:,13) = (/ x, x, 2.d0*x /)
    slip_M(:,14) = (/ x, x, -2.d0*x /)
    slip_M(:,15) = (/ x, -x, 2.d0*x /)
    slip_M(:,16) = (/ -x, x, 2.d0*x /)
    
    slip_M(:,17) = (/ x, 2.d0*x, x /)
    slip_M(:,18) = (/ x, 2.d0*x, -x /)
    slip_M(:,19) = (/ x, -2.d0*x, x /)
    slip_M(:,20) = (/ -x, 2.d0*x, x /)
    
    slip_M(:,21) = (/ 2.d0*x, x, x /)
    slip_M(:,22) = (/ 2.d0*x, x, -x /)
    slip_M(:,23) = (/ 2.d0*x, -x, x /)
    slip_M(:,24) = (/ -2.d0*x, x, x /)
    
    ! {123} types:
    slip_M(:,25) = (/ y, 2.d0*y, 3.d0*y/)
    slip_M(:,26) = (/ y, 2.d0*y, -3.d0*y/)
    slip_M(:,27) = (/ y, -2.d0*y, 3.d0*y/)
    slip_M(:,28) = (/ -y, 2.d0*y, 3.d0*y/)
    
    slip_M(:,29) = (/ 2.d0*y, y, 3.d0*y/)
    slip_M(:,30) = (/ 2.d0*y, y, -3.d0*y/)
    slip_M(:,31) = (/ 2.d0*y, -y, 3.d0*y/)
    slip_M(:,32) = (/ -2.d0*y, y, 3.d0*y/)
    
    slip_M(:,33) = (/ 2.d0*y, 3.d0*y, y/)
    slip_M(:,34) = (/ 2.d0*y, 3.d0*y, -y/)
    slip_M(:,35) = (/ 2.d0*y, -3.d0*y, y/)
    slip_M(:,36) = (/ -2.d0*y, 3.d0*y, y/)
    
    slip_M(:,37) = (/ y, 3.d0*y, 2.d0*y/)
    slip_M(:,38) = (/ y, 3.d0*y, -2.d0*y/)
    slip_M(:,39) = (/ y, -3.d0*y, 2.d0*y/)
    slip_M(:,40) = (/ -y, 3.d0*y, 2.d0*y/)
    
    slip_M(:,41) = (/ 3.d0*y, y, 2.d0*y/)
    slip_M(:,42) = (/ 3.d0*y, y, -2.d0*y/)
    slip_M(:,43) = (/ 3.d0*y, -y, 2.d0*y/)
    slip_M(:,44) = (/ -3.d0*y, y, 2.d0*y/)
    
    slip_M(:,45) = (/ 3.d0*y, 2.d0*y, y/)
    slip_M(:,46) = (/ 3.d0*y, 2.d0*y, y/)
    slip_M(:,47) = (/ 3.d0*y, 2.d0*y, y/)
    slip_M(:,48) = (/ 3.d0*y, 2.d0*y, y/)
        
!--------------------
!  slip directions
!--------------------
    ! {110} types:
    slip_S(:,1)  = (/m, m, m/)
    slip_S(:,2)  = (/m, m, m/)
    slip_S(:,3)  = (/m, m, m/)

    slip_S(:,4)  = (/-m, m, m/)
    slip_S(:,5)  = (/-m, m, m/)
    slip_S(:,6)  = (/-m, m, m/)
  
    slip_S(:,7)  = (/m, -m, m/)
    slip_S(:,8)  = (/m, -m, m/)
    slip_S(:,9)  = (/m, -m, m/)
  
    slip_S(:,10) = (/-m, -m, m/)
    slip_S(:,11) = (/-m, -m, m/)
    slip_S(:,12) = (/-m, -m, m/)
    
    ! {112} types:
    slip_S(:,13) = (/ m, m, -m/)
    slip_S(:,14) = (/ m, m, m/)
    slip_S(:,15) = (/ m, -m, -m/)
    slip_S(:,16) = (/ -m, m, -m/)
    
    slip_S(:,17) = (/ m, -m, m/)
    slip_S(:,18) = (/ m, -m, -m/)
    slip_S(:,19) = (/ m, m, m/)
    slip_S(:,20) = (/ -m, -m, m/)
    
    slip_S(:,21) = (/ -m, m, m/)
    slip_S(:,22) = (/ -m, m, -m/)
    slip_S(:,23) = (/ -m, -m, m/)
    slip_S(:,24) = (/ m, m, m/)
    
    ! {123} types:
    slip_S(:,25) = (/ m, m, -m/)
    slip_S(:,26) = (/ m, m, m/)
    slip_S(:,27) = (/ m, -m, -m/)
    slip_S(:,28) = (/ -m, m, -m/)
    
    slip_S(:,29) = (/ m, m, -m/)
    slip_S(:,30) = (/ m, m, m/)
    slip_S(:,31) = (/ m, -m, -m/)
    slip_S(:,32) = (/ -m, m, -m/)
    
    slip_S(:,33) = (/ m, -m, m/)
    slip_S(:,34) = (/ m, -m, -m/)
    slip_S(:,35) = (/ m, m, m/)
    slip_S(:,36) = (/ -m, -m, m/)
    
    slip_S(:,37) = (/ m, -m, m/)
    slip_S(:,38) = (/ m, -m, -m/)
    slip_S(:,39) = (/ m, m, m/)
    slip_S(:,40) = (/ -m, -m, m/)
    
    slip_S(:,41) = (/ -m, m, m/)
    slip_S(:,42) = (/ -m, m, -m/)
    slip_S(:,43) = (/ -m, -m, m/)
    slip_S(:,44) = (/ m, m, m/)
    
    slip_S(:,45) = (/ -m, m, m/)
    slip_S(:,46) = (/ -m, m, -m/)
    slip_S(:,47) = (/ -m, -m, m/)
    slip_S(:,48) = (/ m, m, m/)
  
  else
    call Fatal('Bad structure type: crystal symmetry '//name//' unkown/not implemented')
  end if
  
  slip_M_out(:,:) = slip_M(:,:Nslip)
  slip_S_out(:,:) = slip_S(:,:Nslip)
!--------------------    
!  disloc. vel directions (depend on dislocation character, edge=parallel to burgers)
!  slip_V = line(char) x slip_M (line is normal/parallel to burgers for edge/screw)
    do islip=1,Nslip
      crossMB(:,islip) = slip_M(:,islip) .cross. slip_S(:,islip)
    end do !islip
    do ich=1,Nchar
      do islip=1,Nslip
        ! present char angle theta is minus theta from daniels linetension paper
        line_T(:,islip,ich) = coschar(ich)*slip_S(:,islip) + sinchar(ich)*crossMB(:,islip)
        slip_V(:,islip,ich) = - coschar(ich)*crossMB(:,islip) + sinchar(ich)*slip_S(:,islip)
        call RotAlign(slip_M(:,islip),line_T(:,islip,ich),(/0.d0,1.d0,0.d0/),(/0.d0,0.d0,1.d0/),rotM(:,:,islip,ich))
      end do !islip
    end do !ich
    
    ! compute matrix Aforest that will be used to determine the dislocation source term of CDT subsystem
    Aforest = 0.d0
    do islip=1,Nslip
      do jslip=1,Nslip
        do ich=1,Nchar
          ! should always be positive
          Aforest(islip,jslip,ich) = abs(DOT_PRODUCT(slip_M(:,islip),line_T(:,jslip,ich)))
        end do !ich
      end do !jslip
    end do !islip

END SUBROUTINE DMB_get_slip_system

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE xtal_Update_Stress(F_new, Fp_old, F_i, Cg, slip_Sg, slip_Mg, dgamma, &
           cauchy_new, Fe_new, Fp_new)
  use stdlib_linalg, only: expm
  USE GlobalParams
  USE DMB
  use utilities, only: Identity, DetInv3x3, operator(.det.)!, ExpM
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fp_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_i
  REAL(KIND=8), DIMENSION(3,3,3,3),       INTENT(IN)    :: Cg
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Sg
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Mg
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(IN)    :: dgamma
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: cauchy_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fe_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fp_new
!-----------------------------------------------------------------------
!  Locals:
  INTEGER :: i,j,k,ich
  REAL(KIND=8), DIMENSION(3,3) :: dE, exp_dE, Fbar_inv, Fbar, eye, Ee, pk2
  REAL(KIND=8) :: D

!-----------------------------------------------------------------------
  dE = 0.d0
  eye = Identity(3)
  do ich=1,Nchar
    do i=1,Nslip
      do j=1,3
        do k=1,3
          dE(k,j) = dE(k,j) + dgamma(i,ich)*slip_Sg(k,i)*slip_Mg(j,i)
        end do  !k
      end do  !j
    end do !i
  end do !ich


  exp_dE = ExpM(dE,order=6)
  Fp_new = MATMUL(exp_dE, Fp_old)
  
  ! djl --> just as a check 
  !Fp_new  = matmul( (eye + dE), Fp_old )  ! 1st-Order Pade approx
  
  Fbar = MATMUL(F_i, Fp_new)
  call DetInv3x3(Fbar,D,Fbar_inv)
  Fe_new = MATMUL(F_new, Fbar_inv)
 
  ! different strain measure could go here
  if (.True.) then
    Ee = 0.5d0*(MATMUL(transpose(Fe_new),Fe_new) - eye)
    do j=1,3
      do i=1,3
        pk2(i,j) = sum(Cg(i,j,:,:)*Ee(:,:))
      end do  !i
    end do  !j
  
    cauchy_new = MATMUL(MATMUL(Fe_new, pk2), TRANSPOSE(Fe_new)) / (.det.Fe_new)

  else
    Ee = 0.5d0*(transpose(Fe_new) + Fe_new) - eye
    do j=1,3
      do i=1,3
        cauchy_new(i,j) = sum(Cg(i,j,:,:)*Ee(:,:))
      end do  !i
    end do  !j
   
  end if

  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "xtal_update_stress()"
    write(*,'("F_new:",3(/,3es28.16))') ((F_new(i,j),j=1,3),i=1,3)
    write(*,'("Fp_old:",3(/,3es28.16))') ((Fp_old(i,j),j=1,3),i=1,3)
    write(*,'("F_i:",3(/,3es28.16))') ((F_i(i,j),j=1,3),i=1,3)
    write(*,'("dE:",3(/,3es28.16))') ((dE(i,j),j=1,3),i=1,3)
    write(*,'("Fbar:",3(/,3es28.16))') ((Fbar(i,j),j=1,3),i=1,3)
    write(*,*)
  end if

END SUBROUTINE xtal_Update_Stress

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE Compute_Velocities(sigma, tau_b, temp, Fe, slip_Sg, slip_Mg, rho,  props, &
           detF, pressure, vel, tau)
  
  USE GlobalParams
  USE DMB
  use utilities, only: DetInv3x3
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: sigma
  REAL(KIND=8), DIMENSION(Nslip),         INTENT(IN)    :: tau_b
  REAL(KIND=8),                           INTENT(IN)    :: detF, pressure
  REAL(KIND=8), DIMENSION(1),             INTENT(IN)    :: temp
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fe
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Sg
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Mg
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(IN)    :: rho
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: vel
  REAL(KIND=8), DIMENSION(Nslip),         INTENT(OUT)   :: tau
!-----------------------------------------------------------------------
!  Locals:
!  REAL(KIND=8) :: small2 = 1.d-10!, exp_lim = 2.d2
  REAL(KIND=8), DIMENSION(3,3)  :: FeInv
  REAL(KIND=8) :: detFe !, tau_a
  INTEGER :: i
  REAL(KIND=8), DIMENSION(3) :: S, M
  REAL(KIND=8), DIMENSION(Nslip) :: tau_eff
!-----------------------------------------------------------------------

  ! Compute resolved shear stress
  call DetInv3x3(Fe,detFe,FeInv)

  do i=1,Nslip
    S      = MATMUL(Fe              , slip_Sg(:,i))
    M      = MATMUL(transpose(FeInv), slip_Mg(:,i))
    tau(i) = detFe * dot_product( S, MATMUL(sigma, M))
  end do  !i    

  ! tau_a = 1.1d0 ! test with athermal barrier to dislocation motion
  ! tau_eff(:) = sign(abs(abs(tau(:) - tau_b(:)) - tau_a),tau(:) - tau_b(:))
  tau_eff(:) = (tau(:) - tau_b(:)) ! * (sign(0.5d0,abs(tau(:) - tau_b(:)) - tau_a)+0.5d0)
      
  ! mean velocity is combination
  call Dislocation_Velocity(tau_eff, temp, rho,  props, detF, pressure, &
                            vel)
  !write(*,'("tau:",/,12es28.16)') tau(:)
  !write(*,'("vel:",/,12es28.16)') vel(:)
  ! debug velocity
  ! vel(:) = 0.d0

END SUBROUTINE Compute_Velocities

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE Dislocation_Velocity(tau_eff, temp, rho,  props, detF, pressure, &
           vel)

  USE GlobalParams
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(Nslip),         INTENT(IN)    :: tau_eff
  REAL(KIND=8),                           INTENT(IN)    :: detF, pressure
  REAL(KIND=8), DIMENSION(1),             INTENT(IN)    :: temp
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(IN)    :: rho
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: vel
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: small2, exp_lim, muP, muT, deltamu
  REAL(KIND=8), DIMENSION(3) :: CP, CT
  REAL(KIND=8) :: c1, brgr, B0(Nchar), Lbar, g0, omega0, boltz, p, q, tau0
  REAL(KIND=8) :: inner, deltaG, exp_arg, rhobar0
  REAL(KIND=8), DIMENSION(Nslip) :: tau_eff2, sign_tau
  REAL(KIND=8), DIMENSION(Nslip,Nchar) :: rho_eff, slip_r
  REAL(KIND=8), DIMENSION(Nslip,Nchar) :: t_wait, t_run, vel_run
  REAL(KIND=8), DIMENSION(Nchar) :: deltaG0, xi0, vcrit, mu
  character(32) :: drag_flag
  LOGICAL :: drag_Trho_flag
  INTEGER :: i, ich
!-----------------------------------------------------------------------
  small2 = 1.d-10; exp_lim = 2.d2
  t_wait = 0.d0
  t_run = 0.d0
  vel_run = 0.d0

  ! Assign properties
  drag_Trho_flag = props%drag_Trho_flag ! decide whether to include T,rho dependence in drag coeff.
  c1          = props%c1
  Lbar        = props%Lbar    
  g0          = props%g0      
  omega0      = props%omega0
  boltz       = props%boltz   
  p           = props%p       
  q           = props%q
  tau0        = props%tau0
  drag_flag   = props%drag_flag

  if (drag_Trho_flag) then
    CP = props%elastic_constants(:,3)
    CT = props%elastic_constants(:,2)
    muP = (CP(1)-CP(2)+2.d0*CP(3))/4.d0 ! rough estimate: average change of shear moduli with P
    muT = (CT(1)-CT(2)+2.d0*CT(3))/4.d0 ! rough estimate: average change of shear moduli with T
    deltamu = muT*(temp(1)-props%T_ref) + muP*pressure
    rhobar0 = props%rhobar0/detF
    brgr    = props%burger*detF**(1./3.)
    mu(:)   = rhobar0*props%wave_vel(:)**2 + deltamu
    vcrit   = sqrt(mu/rhobar0)
    B0      = props%B0*(temp(1)/props%T_ref)*sqrt(1.d0 - deltamu/mu)/(detF**(7./6.))
  else
    rhobar0 = props%rhobar0
    brgr    = props%burger
    vcrit   = props%wave_vel
    mu(:)   = rhobar0*vcrit(:)**2
    B0      = props%B0
  end if

  !vcrit = c_S, brgr is Burgers vector magnitude, B0 is the drag coeff. at zero velocity
  ! (set in input file): can generalize to B0(T,rho)!
  xi0(:)         = 0.5d0*B0(:)*vcrit(:) / brgr
  deltaG0     = g0*mu*brgr**3
  !print*, vcrit, mu, deltaG0

  ! Compute slip resistance
  rho_eff(:,:) = rho(:,:)  ! can add interaction matrix here
  do ich=1,Nchar
    slip_r(:,ich)  = c1*mu(ich)*brgr*sqrt(rho_eff(:,ich)) + tau0
  end do
  sign_tau(:) = sign(1.d0,tau_eff(:))
  tau_eff2(:) = abs(tau_eff(:))

  ! Instantaneous velocity (drag limited + relativistic effect)     
  t_wait(:,:)  = 0.d0
  t_run(:,:)   = 0.d0
  vel_run(:,:) = 0.d0
  vel(:,:)     = 0.d0
  
  if (drag_flag=='iso') then
    do ich=1,Nchar
      do i=1,Nslip
        if (tau_eff2(i) > small2) then
          ! inspired by asymptotic form of edge disloc. in isotropic limit of B=B0/sqrt(1+beta**2):
          ! but here beta=vel_run/vcrit and vcrit may be character dependent
          vel_run(i,ich) = vcrit(ich)/(( (2*xi0(ich)/tau_eff2(i))**2 + 1.d0 )**0.5d0)
        end if
      end do  !i
    end do  !i
  elseif (drag_flag=='const') then
    do ich=1,Nchar
      do i=1,Nslip
        if (tau_eff2(i) > small2) then
          ! constant B:
          vel_run(i,ich) = brgr*tau_eff2(i)/B0(ich)
        end if
      end do  !i
    end do  !i
  else !default is 'Austin':
    do ich=1,Nchar
      do i=1,Nslip
        if (tau_eff2(i) > small2) then
          ! R. Austin's relativistic form of B=B0/(1+beta**2), with beta=vel_run/vcrit, leads to this:
          vel_run(i,ich) = vcrit(ich)*(( (xi0(ich)/tau_eff2(i))**2 + 1.d0 )**0.5d0 - (xi0(ich)/tau_eff2(i)))
        end if
      end do  !i
    end do  !i
  end if ! choose drag

  ! Wait / run times (Kocks thermal activiation)    
  do ich=1,Nchar
  do i=1,Nslip
    inner = 1.d0 - ( tau_eff2(i) / slip_r(i,ich) )**p
    if (inner > 0.d0) then
      deltaG  = deltaG0(ich)*( inner**q )
      exp_arg = deltaG / (boltz*temp(1))
      if (exp_arg < exp_lim) then
        t_wait(i,ich) = (exp(exp_arg) - 1.d0) / omega0
      else
        vel_run(i,ich) = 0.d0 ! i.e., t_wait = inf
      end if
    else
      t_wait(i,ich) = 0.d0
    end if
    if (vel_run(i,ich) > small2) then
      t_run(i,ich) = Lbar / vel_run(i,ich)
    end if
    
    if (vel_run(i,ich) > small2) then
      ! mean velocity is combination
      vel(i,ich) = sign_tau(i)*Lbar / (t_wait(i,ich) + t_run(i,ich))
    end if
  end do  !i
  end do !ich

END SUBROUTINE Dislocation_Velocity

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!> returns gamma_dot using 4th order RK scheme
SUBROUTINE RK_Rates(dt, temp, F_old, F_new, Fe_old, Fp_old, cauchy_old, J_new, pressure, &
           field_variables, slip_Sg, slip_Mg, Cg, props, &
           gamma_dot, error)
  
  USE GlobalParams
  USE DMB
  use utilities, only: Identity
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8),                           INTENT(IN)    :: dt, J_new, pressure
  REAL(KIND=8), DIMENSION(1),             INTENT(IN)    :: temp
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fe_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fp_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: cauchy_old
  TYPE(Type_DMB_FIELD_VARIABLES),   INTENT(IN)    :: field_variables
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Sg
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Mg
  REAL(KIND=8), DIMENSION(3,3,3,3),       INTENT(IN)    :: Cg
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: gamma_dot
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: error
!-----------------------------------------------------------------------
!  Locals:
  LOGICAL :: internal_strain_flag
  REAL(KIND=8) :: brgr
  REAL(KIND=8), DIMENSION(Nslip,Nchar,2)  :: varrho_ip
  REAL(KIND=8), DIMENSION(Nslip)          :: tau_b, tau
  REAL(KIND=8), DIMENSION(Nslip,Nchar)    :: varrho_mobile_sum, vel
  REAL(KIND=8), DIMENSION(3,3)      :: epsilon_i
  REAL(KIND=8), DIMENSION(Nslip,Nchar)    :: gamma_dot_1, gamma_dot_2, gamma_dot_3, gamma_dot_4, gamma_dot_5
  REAL(KIND=8), DIMENSION(Nslip,Nchar)    :: dgamma, gamma_dot_third
  REAL(KIND=8), DIMENSION(3,3)      :: F_i, F_mid, cauchy, F, Fe, Fp!, eye, pk2
!  INTEGER :: i,j
!-----------------------------------------------------------------------

  ! get properties
  brgr                 = props%burger
  internal_strain_flag = props%internal_strain_flag

  varrho_ip = field_variables%varrho_ip
  tau_b     = field_variables%tau_b
  epsilon_i = field_variables%epsilon_i

  varrho_mobile_sum(:,:) = varrho_ip(:,:,1) + varrho_ip(:,:,2)

  F_i = Identity(3)
  if(internal_strain_flag) F_i = F_i + epsilon_i

  F_mid(:,:) = 0.5d0*(F_old(:,:) + F_new(:,:)) ! approximation of above !

  !==============================================================================!
  !       RK Step 1                                                              !
  !==============================================================================!

  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "rk_rates(), step 1"
  end if
  call Compute_Velocities(cauchy_old, tau_b, temp, Fe_old, slip_Sg, slip_Mg, varrho_mobile_sum, props, &
                          J_new, pressure, vel, tau)

  gamma_dot_1(:,:)  = varrho_mobile_sum(:,:)*vel(:,:)*brgr
  !write(*,'("RK2_pre, gamma_dot_1:",/,12es28.16)') gamma_dot_1(:)

  !====================================================================================!
  !       RK Step 2: F=F_mid, dgamma = gdot1*dt/2 --> sig2, vel2, gdot2     !
  !====================================================================================!
      
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "rk_rates(), step 2"
  end if
  dgamma(:,:) = 0.5d0*dt*gamma_dot_1(:,:)
  F(:,:)    = F_mid(:,:)
  !write(*,'("dt=",es28.16)') dt
  !write(*,'("RK2_pre, dgamma:",/,12es28.16)') dgamma(:)
  !write(*,'("RK2_pre, F:",3(/,3es28.16))') ((F(i,j),j=1,3),i=1,3)
  
  call xtal_Update_Stress(F, Fp_old, F_i, Cg, slip_Sg, slip_Mg, dgamma, &
                          cauchy, Fe, Fp)

  !write(*,'("RK2_post, pk2:",3(/,3es28.16))') ((pk2(i,j),j=1,3),i=1,3)
  !write(*,'("RK2, Fe:",3(/,3es28.16))') ((Fe(i,j),j=1,3),i=1,3)
  !write(*,'("RK2, Fp:",3(/,3es28.16))') ((Fp(i,j),j=1,3),i=1,3)
  !write(*,'("RK2, cauchy:",3(/,3es28.16))') ((cauchy(i,j),j=1,3),i=1,3)

  call Compute_Velocities(cauchy, tau_b, temp, Fe, slip_Sg, slip_Mg, varrho_mobile_sum,  props, &
                          J_new, pressure, vel, tau)
  gamma_dot_2(:,:)  = varrho_mobile_sum(:,:)*vel(:,:)*brgr

  !====================================================================================!
  !       RK Step 3: F=F_mid, dgamma = gdot2*dt/2 --> sig3, vel3, rho_dot_3, gdot3     !
  !====================================================================================!
      
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "rk_rates(), step 3"
  end if
  dgamma(:,:) = 0.5d0*dt*gamma_dot_2(:,:)
  F(:,:)    = F_mid(:,:)

  call xtal_Update_Stress(F, Fp_old, F_i, Cg, slip_Sg, slip_Mg, dgamma, &
                          cauchy, Fe, Fp)
      
  call Compute_Velocities(cauchy, tau_b, temp, Fe, slip_Sg, slip_Mg, varrho_mobile_sum,  props, &
                          J_new, pressure, vel, tau)
			  
  gamma_dot_3(:,:)  = varrho_mobile_sum(:,:)*vel(:,:)*brgr
      
  !====================================================================================!
  !       RK Step 4: F=F_new, dgamma = gdot3*dt --> sig4, vel4, rho_dot_4, gdot4       !
  !====================================================================================!
      
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "rk_rates(), step 4"
  end if
  dgamma(:,:) = dt*gamma_dot_3(:,:)
  F(:,:)    = F_new(:,:)

  call xtal_Update_Stress(F, Fp_old, F_i, Cg, slip_Sg, slip_Mg, dgamma, &
                          cauchy, Fe, Fp)
  
  call Compute_Velocities(cauchy, tau_b, temp, Fe, slip_Sg, slip_Mg, varrho_mobile_sum,  props, &
                          J_new, pressure, vel, tau)
  gamma_dot_4(:,:)  = varrho_mobile_sum(:,:)*vel(:,:)*brgr


  !====================================================================================!
  !       RK Step 5: F=F_new, (this is for third order RK to estimate error            !
  !====================================================================================!
      
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'(a)') "rk_rates(), step 5"
  end if
  dgamma(:,:) = dt*( 2.d0*gamma_dot_2(:,:) - gamma_dot_1(:,:))
  F      = F_new
  
  call xtal_Update_Stress(F, Fp_old, F_i, Cg, slip_Sg, slip_Mg, dgamma, &
                          cauchy, Fe, Fp)
  
  call Compute_Velocities(cauchy, tau_b, temp, Fe, slip_Sg, slip_Mg, varrho_mobile_sum,  props, &
                          J_new, pressure, vel, tau)
			  
  gamma_dot_5(:,:)  = varrho_mobile_sum(:,:)*vel(:,:)*brgr

  !====================================================================================!
  !       Final update of evolving vars based on DEs                                   !
  !====================================================================================!        
      
  gamma_dot(:,:)    = (gamma_dot_1(:,:) + 2.d0*(gamma_dot_2(:,:) + gamma_dot_3(:,:)) + gamma_dot_4(:,:)) / (6.d0)
  
  gamma_dot_third(:,:) = (gamma_dot_1(:,:) + 4.d0*gamma_dot_2(:,:)  + gamma_dot_5(:,:)) / (6.d0)
 
  error(:,:) = abs((gamma_dot(:,:) - gamma_dot_third(:,:))) / (abs(gamma_dot(:,:)) + 1.d-3) ! error norm

END SUBROUTINE RK_Rates

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE RK_SubIncrement(dt, T_old,  F_old, F_new, J_new, pressure, Fe_old, &
               Fp_old, cauchy_old, field_variables, slip_Sg, slip_Mg, Cmat, props, &
               cauchy, Fe_new, Fp_new, dgamma, Nsub, icutback)
  
  USE GlobalParams
  USE DMB
  use utilities, only: Identity
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8),                           INTENT(IN)    :: dt, J_new, pressure
  REAL(KIND=8), DIMENSION(1),             INTENT(IN)    :: T_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fe_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fp_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: cauchy_old
  TYPE(Type_DMB_FIELD_VARIABLES),   INTENT(IN)    :: field_variables
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Sg
  REAL(KIND=8), DIMENSION(3,Nslip),       INTENT(IN)    :: slip_Mg
  REAL(KIND=8), DIMENSION(3,3,3,3),       INTENT(IN)    :: Cmat
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: cauchy
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fe_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fp_new
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: dgamma
  INTEGER,                          INTENT(OUT)   :: Nsub
  INTEGER,                          INTENT(OUT)   :: icutback
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(Nslip,Nchar,2) :: varrho_ip
!  REAL(KIND=8), DIMENSION(Nslip,Nchar)   :: varrho_mobile_sum
  REAL(KIND=8) gtol, t_prime, t_final, dt_prime, dt_fac
  REAL(KIND=8), DIMENSION(3,3) :: F_i, epsilon_i, F_dot, F0, cauchy0, Fe0, Fp0, Fp, F1, pk2
  LOGICAL :: flag_continue, flag_last_step
  REAL(KIND=8), DIMENSION(Nslip) :: tau_b
  REAL(KIND=8), DIMENSION(Nslip,Nchar) :: gamma_dot, gdot_error, dgamma_total
  INTEGER :: i,j
!-----------------------------------------------------------------------
  gtol = props%gtol

  t_prime = 0.d0
  t_final = dt

  ! compute def grad. associated with dislocation configuration
  
  !   -- note, this should probably be incremented from 
  !      "old" value throughout the subincrements. would 
  !      require more wiring than I want to do right now.
  varrho_ip = field_variables%varrho_ip
  tau_b     = field_variables%tau_b
  epsilon_i = field_variables%epsilon_i

  F_i = Identity(3)
  if(props%internal_strain_flag) F_i = F_i + epsilon_i
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'("rk_subincrement()")')
    write(*,'("epsilon_i:",3(/,3es28.16))') ((epsilon_i(i,j),j=1,3),i=1,3)
    write(*,'("F_i:",3(/,3es28.16))') ((F_i(i,j),j=1,3),i=1,3)
  end if

  ! Compute constant deformation rate
  F_dot(:,:) = (F_new(:,:)-F_old(:,:))/dt
  ! L     = sl.logm(matmul(F_new, k_inv3x3(F_old)) ) / dt
  
  ! Initialize state
  dgamma_total = 0.d0
  dt_prime     = dt
  F0(:,:)      = F_old(:,:)
  cauchy0(:,:) = cauchy_old (:,:)
  Fe0(:,:)     = Fe_old (:,:)
  Fp0(:,:)     = Fp_old(:,:)
  F1(:,:)      = F_new(:,:)
  !write(*,'("F_dot:",3(/,3es28.16))') ((F_dot(i,j),j=1,3),i=1,3)
  !write(*,'("F0:",3(/,3es28.16))') ((F0(i,j),j=1,3),i=1,3)
  !write(*,'("cauchy0:",3(/,3es28.16))') ((cauchy0(i,j),j=1,3),i=1,3)
  !write(*,'("Fe0:",3(/,3es28.16))') ((Fe0(i,j),j=1,3),i=1,3)
  !write(*,'("Fp0:",3(/,3es28.16))') ((Fp0(i,j),j=1,3),i=1,3)
  !write(*,'("F1:",3(/,3es28.16))') ((F1(i,j),j=1,3),i=1,3)
  
  ! Loop until full time step has been accomplished (or gives up)
  icutback       = 0
  flag_continue  = .True.
  flag_last_step = .True.
  Nsub           = 0
  dt_fac         = 1.2d0
  do while(flag_continue)

    if(current_el==1.and.ncycle==debug_cycle) then
      write(*,'("rk_subincrement(), iterating, Nsub=",i5,", icutback=",i5)') Nsub,icutback
    end if
    ! update F1
    if (flag_last_step) then
      F1(:,:) = F_new(:,:)
    else
      F1(:,:) = F0(:,:) + F_dot(:,:)*dt_prime
    end if

    call RK_Rates(dt_prime, T_old,  F0, F1, Fe0, Fp0, cauchy0, J_new, pressure, &
                  field_variables, slip_Sg, slip_Mg, Cmat, props, &
                  gamma_dot, gdot_error)
    if(current_el==1.and.ncycle==debug_cycle) then
      write(*,'("gamma_dot:",/,12es28.16)') (gamma_dot(i,Nchar),i=1,Nslip)! writing edge TODO
      write(*,'("Nsub,icutback,dt,t_prime,dt_prime=",2i5,3es10.3)') Nsub,icutback,dt,t_prime,dt_prime
    end if

    if (maxval(gdot_error(:,:)) > gtol) then
      ! reduce time step size and try again
      dt_prime       = 0.2d0*dt_prime
      flag_last_step = .False.
      icutback       = icutback + 1
      if(current_el==1.and.ncycle==debug_cycle) then
        write(*,'(a)') "Cutting back..."
      end if
      
    else
      ! check if this was the last step
      if (flag_last_step) then
        flag_continue = .False.
        if(current_el==1.and.ncycle==debug_cycle) then
          write(*,'(a)') "Last step..."
        end if
      end if
      
      ! update state and continue on

      dgamma(:,:)    = gamma_dot(:,:)*dt_prime
      dgamma_total = dgamma_total + dgamma
      
      F0(:,:)      = F1(:,:)
      !write(*,'("dgamma:",/,12es28.16)') (dgamma(j),j=1,Nslip)
      !write(*,'("F0:",3(/,3es28.16))') ((F0(i,j),j=1,3),i=1,3)
      
      ! update time
      if (t_prime+ 1.01d0*dt_prime >= t_final) then
        dt_prime       = t_final-t_prime
        t_prime        = t_final
        flag_last_step = .True.
        
      else
        if(t_prime+dt_fac*dt_prime >= t_final) then
          dt_prime       = t_final - t_prime
          t_prime        = t_final
          flag_last_step = .True.
        else
          dt_prime = dt_fac * dt_prime
          t_prime  = t_prime + dt_prime
          if(current_el==1.and.ncycle==debug_cycle) then
            write(*,'(a,es10.3)') "Sub-cycling..."
          end if
        end if
      end if
        
      ! update Fe, Fp, cauchy
      call xtal_Update_Stress(F0, Fp0, F_i, Cmat, slip_Sg, slip_Mg, dgamma, &
                              cauchy0, Fe0, Fp)
      Fp0 = Fp

      if(current_el==debug_element.and.ncycle==debug_cycle) then
        write(*,'("pk2:",3(/,3es28.16))') ((pk2(i,j),j=1,3),i=1,3)
        write(*,'("Fe0:",3(/,3es28.16))') ((Fe0(i,j),j=1,3),i=1,3)
        write(*,'("Fp0:",3(/,3es28.16))') ((Fp0(i,j),j=1,3),i=1,3)
        write(*,'("cauchy0:",3(/,3es28.16))') ((cauchy0(i,j),j=1,3),i=1,3)
      end if
      Nsub    = Nsub + 1

    end if
  end do  !while(flag_continue)

  cauchy(:,:) = cauchy0(:,:)
  Fe_new(:,:) = Fe0(:,:)
  Fp_new(:,:) = Fp0(:,:)
  
  ! djl: dgamma_total added to totalize over all subincrements, then total 
  !      dgamma returned rather than final sub-increment's dgamma.
  dgamma = dgamma_total

END SUBROUTINE RK_SubIncrement

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>Computes the dislocation displacement gradient field according to the integral method
SUBROUTINE stroh_dislocation(props,vel,resolution,uij,Sb,Bb,Cmat)
  USE GlobalParams
  USE DMB
  use utilities, only: Identity, LinSpace, elbrak, Trapz, DetInv3x3
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
  INTEGER,                          INTENT(IN)    :: resolution
  REAL(KIND=8), DIMENSION(Nslip,Nchar), INTENT(IN) :: vel
! Outputs:
  REAL(KIND=8), DIMENSION(3,3,resolution,Nslip,Nchar), INTENT(OUT) :: uij
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar), INTENT(OUT) :: Sb, Bb ! stroh matrices S, B projected onto unit Burgers vector
  REAL(KIND=8), DIMENSION(3,3,3,3), INTENT(OUT) :: Cmat 
!-----------------------------------------------------------------------
  REAL(KIND=8), DIMENSION(3,Nslip)     :: slip_S, slip_M
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar) :: slip_V
  REAL(KIND=8), DIMENSION(3,resolution) :: stroh_M, stroh_N
  REAL(KIND=8), DIMENSION(3,3,resolution) :: NNinv, Sphi, Bphi
  REAL(KIND=8), DIMENSION(3,3) :: NN, MM, MN, NM, stroh_S, stroh_B
  REAL(KIND=8) :: det !, burger
  REAL(KIND=8) :: elastic_constants(3,3), Qrot(3,3)
  REAL(KIND=8) :: eye(3,3), rhobar0, Cdyn(3,3,3,3)
  REAL(KIND=8) :: phi(resolution) ! polar angle to be integrated over
  INTEGER :: i,j,k,l,kk,ll,oo
  !-- initialize output variables to 0:
  uij = 0.d0
  Sb = 0.d0
  Bb = 0.d0
  eye = Identity(3)
  ! ------------- assign properties from props list ----------------------------
  
  elastic_constants(:,:) = 0.d0
  elastic_constants(1,1) = props%C11
  elastic_constants(2,1) = props%C12
  elastic_constants(3,1) = props%C44
  
!~   burger = props%burger
  rhobar0 = props%rhobar0
  slip_M = props%slip_M
  slip_S = props%slip_S
  slip_V = props%slip_V
  
  ! don't rotate, keep in crystal coordinates for now (Qrot=eye)
  call DMB_Compute_Rotation((/0.d0,0.d0,0.d0/), Qrot)
  ! use room temperature and ambient pressure for now; MIGHT recompute for current conditions 
  ! in some later implementation, but that would require a lot of extra computation time for
  ! backstress computations, so need to check if it's worth it
  call DMB_form_elasticity(Cmat, Qrot, elastic_constants, 0.d0, 0.d0)
  ! now have tensor of SOEC in crystal axis stored in Cmat
  call LinSpace(0.d0,2.d0*pi,resolution,phi)

  do i=1,Nchar
    do j=1,Nslip
      stroh_B = 0.d0
      stroh_S = 0.d0
      Sphi=0.d0
      Bphi=0.d0
      ! compute shifted/dynamic elastic tensor
      do k=1,3
        do l=1,3
          do kk=1,3
           do ll=1,3
             Cdyn(ll,kk,l,k) = Cmat(ll,kk,l,k) - rhobar0*vel(j,i)*vel(j,i)*slip_V(ll,j,i)*eye(kk,l)*slip_V(k,j,i)
           end do !ll
          end do !kk
        end do !l
      end do !k
      ! compute vectors M(phi), N(phi) from m0 = -slip_V and n0 = slip_M
      do k=1,resolution
        do l=1,3
          stroh_M(l,k) = - slip_V(l,j,i)*cos(phi(k)) + slip_M(l,j)*sin(phi(k))
          stroh_N(l,k) = slip_M(l,j)*cos(phi(k)) + slip_V(l,j,i)*sin(phi(k))
        end do !l
        ! compute some building blocks (NN, NM, MN, MM, NNinv, Sphi, Bphi) for uij
        call elbrak(stroh_N(:,k),stroh_N(:,k),Cdyn,NN)
        call elbrak(stroh_M(:,k),stroh_M(:,k),Cdyn,MM)
        call elbrak(stroh_M(:,k),stroh_N(:,k),Cdyn,MN)
        call elbrak(stroh_N(:,k),stroh_M(:,k),Cdyn,NM)
        ! call Invert_Tensor(NN(:,:,k,j,i),3,NNinv(:,:,k,j,i)) ! lapack maybe overkill and too slow
        call DetInv3x3(NN,det,NNinv(:,:,k))
        do oo=1,3
          do ll=1,3
            do kk=1,3
              Sphi(kk,ll,k) = Sphi(kk,ll,k) - NNinv(kk,oo,k)*NM(oo,ll)
            end do !kk
          end do !ll
        end do !oo
        do ll=1,3
          do kk=1,3
            Bphi(kk,ll,k) = MM(kk,ll)
            do oo=1,3
              Bphi(kk,ll,k) = Bphi(kk,ll,k) + MN(kk,oo)*Sphi(oo,ll,k)
            end do !oo
          end do !kk
        end do !ll
      end do !resolution
      ! compute stroh matrices S, B by integrating building blocks Sphi, Bphi
      do ll=1,3
        do kk=1,3
          call Trapz(resolution, Sphi(kk,ll,:)/(4.d0*pi*pi), phi, stroh_S(kk,ll))
          call Trapz(resolution, Bphi(kk,ll,:)/(4.d0*pi*pi), phi, stroh_B(kk,ll))
        end do !kk
      end do !ll
      do kk=1,3
        do ll=1,3
          Sb(kk,j,i) = Sb(kk,j,i) + stroh_S(kk,ll)*slip_S(ll,j)
          Bb(kk,j,i) = Bb(kk,j,i) + stroh_B(kk,ll)*slip_S(ll,j)
        end do
      end do
      ! assemble everything into dislocation gradient field uij (note: uij is divided by burger here!)
      do k=1,resolution
        do kk=1,3
          do ll=1,3
            uij(kk,ll,k,j,i) = uij(kk,ll,k,j,i) - Sb(kk,j,i)*stroh_M(ll,k)
            do  oo=1,3
              uij(kk,ll,k,j,i) = uij(kk,ll,k,j,i) + stroh_N(ll,k)*(NNinv(kk,oo,k)*Bb(oo,j,i) - Sphi(kk,oo,k)*Sb(oo,j,i))
            end do !oo
          end do !ll
        end do !kk
      end do !resolution
    end do !Nslip
  end do !Nchar

END SUBROUTINE stroh_dislocation

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE Constitutive_DMB( &
           F_new, F_old, T_old, cauchy_old, Fp_old, Fe_old, field_variables, dt, props, elem_len, &
           cauchy, vel, vel_x, vel_y, vel_z, T_new, Fp_new, Fe_new, tau, dgamma, Nsub, i_cutback)
  USE GlobalParams
  USE DMB
  use EOS
  use utilities, only: DetInv3x3, operator(.det.)
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: F_old
  REAL(KIND=8), DIMENSION(1),             INTENT(IN)    :: T_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: cauchy_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fp_old
  REAL(KIND=8), DIMENSION(3,3),           INTENT(IN)    :: Fe_old
  TYPE(Type_DMB_FIELD_VARIABLES),   INTENT(IN)    :: field_variables
  REAL(KIND=8),                           INTENT(IN)    :: dt
  TYPE(Type_DMB_PROPS),             INTENT(IN)    :: props
  REAL(KIND=8),                           INTENT(IN)    :: elem_len
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: cauchy
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: vel
  REAL(KIND=8), DIMENSION(Nslip,Nchar),         INTENT(OUT)   :: vel_x, vel_y, vel_z
  REAL(KIND=8), DIMENSION(1),             INTENT(OUT)   :: T_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fp_new
  REAL(KIND=8), DIMENSION(3,3),           INTENT(OUT)   :: Fe_new
  REAL(KIND=8), DIMENSION(Nslip),         INTENT(OUT)   :: tau
  REAL(KIND=8), dimension(Nslip,Nchar),         intent(out)   :: dgamma
  
  INTEGER,                          INTENT(OUT)   :: Nsub
  INTEGER,                          INTENT(OUT)   :: i_cutback
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: third
  REAL(KIND=8), DIMENSION(Nslip,Nchar,2) :: varrho_ip
  REAL(KIND=8), DIMENSION(Nslip)   :: tau_b
  REAL(KIND=8), DIMENSION(3,3)     :: epsilon_i
  REAL(KIND=8) :: rhobar0, C11, C12, C44, burger, g0, bulk_c1, bulk_c2, J_old, J_new, plin, D
  REAL(KIND=8), DIMENSION(3,Nslip) :: slip_S, slip_Sg
  REAL(KIND=8), DIMENSION(3,Nslip) :: slip_M, slip_Mg
  REAL(KIND=8), DIMENSION(3,Nslip,nchar) :: slip_V, slip_Vg
  REAL(KIND=8), DIMENSION(3) :: euler_angle, v_int
  REAL(KIND=8), DIMENSION(3,3) :: Q, Fp_inv
  REAL(KIND=8), DIMENSION(3,3,3,3) :: Cmat
  REAL(KIND=8), DIMENSION(Nslip,Nchar) :: varrho_mobile_sum
  REAL(KIND=8), DIMENSION(3,3) :: Fbar_old, Fbar_new
  REAL(KIND=8), DIMENSION(4) :: pressure
  REAL(KIND=8), DIMENSION(4) :: free_energy
  REAL(KIND=8), DIMENSION(3) :: entropy
  REAL(KIND=8), DIMENSION(3) :: thermo_coeffs, temperature_change
  REAL(KIND=8) :: specific_heat, bulk_modulus, wave_speed, jdot, peos, p_diss
  REAL(KIND=8) :: pressure_old, T_ref
  REAL(KIND=8) :: d_plastic_work
  REAL(KIND=8), DIMENSION(3,Nslip,Nchar) :: vel0
  INTEGER :: i, ich
!-----------------------------------------------------------------------
  third = 1.d0 / 3.d0
  !==============================================================================!
  !       Step 0:  get properties, parameters assigned and defined               !
  !==============================================================================!

  ! ------------- assign properties from props list ----------------------------
  
  rhobar0     = props%rhobar0
  C11         = props%C11    
  C12         = props%C12    
  C44         = props%C44    
  
  burger      = props%burger    
  g0          = props%g0     
  bulk_c1     = props%bulk_c1
  bulk_c2     = props%bulk_c2
  euler_angle = props%euler_angle
  
  slip_S      = props%slip_S
  slip_M      = props%slip_M
  slip_V      = props%slip_V
  
  ! crystal-basis elasticity
  Q       = props%Qrot
  Cmat    = props%Cg
  slip_Sg = props%slip_Sg
  slip_Vg = props%slip_Vg
  
  T_ref   = props%T_ref
  
  pressure_old = -(cauchy_old(1,1)+cauchy_old(2,2)+cauchy_old(3,3)) / 3.d0

  call DMB_form_elasticity(Cmat, props%Qrot, props%elastic_constants, &
                            T_old(1)-T_ref, pressure_old)

  ! ------------- unpack field variables from list -----------------------------
  varrho_ip = field_variables%varrho_ip
  tau_b     = field_variables%tau_b
  epsilon_i = field_variables%epsilon_i
  
  ! get slip systems into global coordinates
  slip_Sg(:,:) = 0.d0
  slip_Mg(:,:) = 0.d0
  do i=1, Nslip
    slip_Sg(:,i) = Q(:,1)*slip_S(1,i) + Q(:,2)*slip_S(2,i) + Q(:,3)*slip_S(3,i)
    slip_Mg(:,i) = Q(:,1)*slip_M(1,i) + Q(:,2)*slip_M(2,i) + Q(:,3)*slip_M(3,i)
  end do  !i
  
  do ich=1, Nchar
    do i=1, Nslip
      slip_Vg(:,i,ich) = Q(:,1)*slip_V(1,i,ich) + Q(:,2)*slip_V(2,i,ich) + Q(:,3)*slip_V(3,i,ich)
    end do  !i
  end do  !ich
  

  ! debug constitutive
  ! pdb.set_trace()

  varrho_mobile_sum(:,:) = varrho_ip(:,:,1) + varrho_ip(:,:,2)

  !==============================================================================!
  !       Step 0:    Decompose into isochoric def. grad.                         !
  !==============================================================================!

  J_old = .det. F_old
  J_new = .det. F_new
    
  if (props%eos%eos_flag) then
      
    Fbar_old(:,:) = ( J_old**(-third) )*F_old(:,:)
    Fbar_new(:,:) = ( J_new**(-third) )*F_new(:,:)

  else

    Fbar_old(:,:) = F_old(:,:)
    Fbar_new(:,:) = F_new(:,:)

  end if

  !==============================================================================!
  !       Step 1:    Subincrementation, RK scheme to update Fp, Fe, stress       !
  !==============================================================================!
  
  call RK_SubIncrement(dt, T_old,  Fbar_old, Fbar_new, J_new, pressure_old, &
          Fe_old, Fp_old, cauchy_old, field_variables, slip_Sg, slip_Mg, Cmat, props, &
          cauchy, Fe_new, Fp_new, dgamma, Nsub, i_cutback)
  !write(*,'("cauchy:",3(/,3es28.16))') ((cauchy(i,j),j=1,3),i=1,3)
  !write(*,'("Fe_new:",3(/,3es28.16))') ((Fe_new(i,j),j=1,3),i=1,3)
  !write(*,'("Fp_new:",3(/,3es28.16))') ((Fp_new(i,j),j=1,3),i=1,3)
  !write(*,'("dgamma:",/,12es28.16)') (dgamma(j),j=1,Nslip)
  if(current_el==1.and.ncycle==debug_cycle) then
    write(*,'("Done with rk_subincrement()")')
    write(*,'("Nsub=",i2)') Nsub
    write(*,'("i_cutback=",i2)') i_cutback
  end if
  call Compute_Velocities(cauchy, tau_b, T_old, Fe_new, slip_Sg, slip_Mg, varrho_mobile_sum,  props, &
                          J_new, pressure_old, vel, tau)

  !==============================================================================!
  !       Step 2:    "reversible" pressure from equation of state                !
  !==============================================================================!

  if (props%eos%eos_flag) then
    plin     = -third*(cauchy(1,1)+cauchy(2,2)+cauchy(3,3))

    call DMB_Compute_EOS(props%eos, J_new, T_old(1), &
                         pressure, free_energy, entropy, thermo_coeffs)
    if(current_el==1.and.ncycle==debug_cycle) then
      write(*,'(a,4es28.16)') "pressure=",(pressure(i),i=1,4)
      write(*,'(a,4es28.16)') "free_energy=",(free_energy(i),i=1,4)
      write(*,'(a,3es28.16)') "entropy=",(entropy(i),i=1,3)
      write(*,'(a,3es28.16)') "thermo_coeffs=",(thermo_coeffs(i),i=1,3)
    end if
    peos = pressure(1) - props%eos%pressure0

    cauchy(1,1) = cauchy(1,1) - (peos - plin)
    cauchy(2,2) = cauchy(2,2) - (peos - plin)
    cauchy(3,3) = cauchy(3,3) - (peos - plin)

    if ( abs(thermo_coeffs(1)) < 1.d-15 ) then
       ! fall back to this if EOS parameters necessary to calculate it were not given
       thermo_coeffs(1) = props%specific_heat
    end if
    specific_heat = thermo_coeffs(1)
    bulk_modulus  = thermo_coeffs(3)
    wave_speed    = sqrt(bulk_modulus / rhobar0)
  
  else
    specific_heat = props%specific_heat  
    bulk_modulus = (C11 + 2.d0*C12) * third
    wave_speed    = sqrt(bulk_modulus / rhobar0)
    thermo_coeffs = (/specific_heat, 0.d0, bulk_modulus/) 
  end if
  
  !==============================================================================!
  !       Step 3: Include bulk viscosity                                         !
  !==============================================================================!  

  ! L    = sl.logm(nmp.dot(def_grad, sl.inv(F_old) ) ) / dt
  ! jdot = nmp.trace(L) 
  jdot   = (F_new(1,1) - F_old(1,1)) / dt ! approximation for 1D case
  
  ! pressure is positive in compression, negative in tension
  p_diss  = -bulk_c1*rhobar0*wave_speed*elem_len*jdot                     ! linear term
  p_diss = p_diss - sign(1.d0,jdot)*rhobar0*( bulk_c2*elem_len*jdot )**2  ! quadratic term 
  
  cauchy(1,1) = cauchy(1,1) - p_diss ! add dissipative pressure back into stress
  cauchy(2,2) = cauchy(2,2) - p_diss
  cauchy(3,3) = cauchy(3,3) - p_diss

  !==============================================================================!
  !       Step 4: Update temperature                                             !
  !==============================================================================!  

  ! shock heating
  call EOS_Shock_Heat(props%eos, thermo_coeffs, p_diss, J_new, J_old, T_old(1), &
       temperature_change)
  T_new = T_old + temperature_change(1)

  ! plastic work
  do ich=1,Nchar
    d_plastic_work = sum( (tau(:) - tau_b(:))*dgamma(:,ich) )
  end do

  ! T_new = T_old ! can turn off heating for debug ...

  T_new(1) = T_new(1) + d_plastic_work/(props%rhobar0*specific_heat)

  call DetInv3x3(Fp_new,D,Fp_inv)
  do ich=1,Nchar
    do i=1,Nslip
      v_int(1:3) = slip_Vg(1:3,i,ich) * vel(i,ich)
      vel0(:,i,ich)  = MATMUL(Fp_inv, v_int)
    end do  !i
  end do  !ich

  vel_x(:,:) = vel0(1,:,:) ! component of velocity in the global `x' direction
                                               ! assuming screw or edge character.
  vel_y(:,:) = vel0(2,:,:)
  vel_z(:,:) = vel0(3,:,:)
  
  ! if globals.stop == True:
  !     print error

END SUBROUTINE Constitutive_DMB

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE DMB_Update_State(state_old, dt, mass_array, xdot_bc, force_bc, & 
           all_props_dmb, cdt_state, ddc_state,                           &
           state_new, xdot_acc_bc)
  USE GlobalParams
  USE DMB
  USE CDT
  USE DDC
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_STATE),             INTENT(IN)    :: state_old
  REAL(KIND=8), INTENT(IN)                              :: dt, xdot_acc_bc
  REAL(KIND=8), DIMENSION(Nnode),         INTENT(IN)    :: mass_array
  TYPE(Type_DMB_VECTOR_BC),         INTENT(IN)    :: xdot_bc, force_bc
  TYPE(Type_DMB_PROPS), dimension(Nregion), INTENT(IN)    :: all_props_dmb
  TYPE(Type_CDT_STATE),             INTENT(IN)    :: cdt_state
  TYPE(Type_DDC_STATE),             INTENT(IN)    :: ddc_state
!-----------------------------------------------------------------------
!  Outputs:
  TYPE(Type_DMB_STATE),             INTENT(INOUT)   :: state_new
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(:,:,:,:), allocatable :: varrho_field
  REAL(KIND=8), DIMENSION(:,:), allocatable   :: tau_back_field
  REAL(KIND=8), DIMENSION(:,:,:), allocatable     :: epsilon_i_field
  REAL(KIND=8), DIMENSION(Nnode,3)     :: f_int, f_ext
  REAL(KIND=8), DIMENSION(3)           :: time_old, time
  REAL(KIND=8)                         :: dt_mid
  REAL(KIND=8), DIMENSION(Nel)         :: element_lengths
  REAL(KIND=8), DIMENSION(Nel,3)       :: dxdX1
  REAL(KIND=8), DIMENSION(2)           :: vec1
  REAL(KIND=8)                         :: elem_len
  REAL(KIND=8), DIMENSION(Nslip,Nchar, 2)    :: varrho_ip
  REAL(KIND=8), DIMENSION(Nslip)       :: tau_b
  REAL(KIND=8), DIMENSION(Nslip,Nchar) :: dgamma
  REAL(KIND=8), DIMENSION(3,3)         :: epsilon_i
  TYPE(Type_DMB_FIELD_VARIABLES) :: field_variables
  REAL(KIND=8), DIMENSION(3,3)         :: F_new, F_old, cauchy_old, Fp_old, Fe_old
  REAL(KIND=8), DIMENSION(1)           :: T_old, temp_new
!  REAL(KIND=8)                         :: P11
  REAL(KIND=8), DIMENSION(3,3)         :: sigma, Fp_new, Fe_new!, cauchy_new
  REAL(KIND=8), DIMENSION(Nslip,Nchar)       :: vel, vel_x, vel_y, vel_z
  REAL(KIND=8), DIMENSION(Nslip)       :: tau
  INTEGER                        :: Nsub
  INTEGER                        :: i_cutback, i, iregion, iel
!-----------------------------------------------------------------------

  allocate(varrho_field(Nel,Nslip,Nchar,2), tau_back_field(Nel,Nslip), epsilon_i_field(Nel,3,3))
!~   call state_new%allocatememory() ! will be called within DMB_State() below
  call field_variables%allocatememory()
  
  varrho_field    = cdt_state%varrho_field
  tau_back_field  = ddc_state%tau_b
  epsilon_i_field = ddc_state%epsilon_i

  !Nel = Nip ! note assumption here

  ! compute time step information
  time_old = state_old%time
  time     = (/time_old(1)+dt, time_old(1), time_old(2) /)
  dt_mid   = 0.5d0*(time(1) - time(3))
 
  ! compute element lengths and axial stretch field
  do iel = 1, Nel
    element_lengths(iel) = state_old%X_ref(iel+1) - state_old%X_ref(iel)
    dxdX1(iel,:)         = (state_old%x(iel+1,:) - state_old%x(iel,:)) / element_lengths(iel)
  end do  !iel

  ! initialize new state
  call DMB_State(time, state_new)      
  state_new%F(:,:,1) = dxdX1(:,:)
  state_new%X_ref    = state_old%X_ref

  !write(*,'(a)') "Before constitutive_DMB call"
  !do iel=1,Nel
  !  write(*,'(i2,/,"F(:,:,1):",3es28.16)') iel,(state_new%F(iel,i,1),i=1,3)
  !end do  !iel

  !call Print_DMB_State(state_new)
  
  ! array for computing element nodal forces
  vec1(:) = (/-1.d0, 1.d0/)

  f_int = 0.d0
  f_ext = 0.d0
  ! loop over elements
  !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(iregion, iel, elem_len, current_el, field_variables, &
  !$OMP   varrho_ip, tau_b, epsilon_i, F_new, F_old, T_old, cauchy_old, Fp_old, Fe_old, &
  !$OMP   sigma, vel, vel_x, vel_y, vel_z, temp_new, Fp_new, Fe_new, tau, dgamma, Nsub, i_cutback)
  do iregion = 1,nregion
    do iel = region_mask(iregion,1), region_mask(iregion,2)
      
    current_el = iel
    elem_len = element_lengths(iel)
  
    ! interpolate dis. density and backstress
    varrho_ip(:,:,:) = varrho_field(iel,:,:,:)
    tau_b(:)       = tau_back_field(iel,:)
    epsilon_i(:,:) = epsilon_i_field(iel,:,:)

    field_variables%varrho_ip = varrho_ip
    field_variables%tau_b     = tau_b
    field_variables%epsilon_i = epsilon_i
    
    ! get state variables at integration point
    F_new(:,:)      = state_new%F(iel,:,:)
    !write(*,*) "F_new, iel=",iel
    !write(*,*) F_new
    F_old(:,:)      = state_old%F(iel,:,:)
    T_old(:)        = state_old%T(iel,:)
    cauchy_old(:,:) = state_old%cauchy(iel,:,:)
    Fp_old(:,:)     = state_old%Fp(iel,:,:)
    Fe_old(:,:)     = state_old%Fe(iel,:,:)
    
    ! evaluate constitutive update
    call Constitutive_DMB( &
        F_new, F_old, T_old, cauchy_old, Fp_old, Fe_old, field_variables, dt, &
        all_props_dmb(iregion), elem_len,                                     &
        sigma, vel, vel_x, vel_y, vel_z, temp_new, Fp_new, Fe_new, tau, dgamma, Nsub, i_cutback)

!~     if(iel==-1) then
!~       write(*,*) "sigma"     ,sigma
!~       write(*,*) "vel"       ,vel
!~       write(*,*) "vel_x"     ,vel_x
!~       write(*,*) "temp_new"  ,temp_new
!~       write(*,*) "Fp_new"    ,Fp_new
!~       write(*,*) "Fe_new"    ,Fe_new
!~       write(*,*) "tau"       ,tau
!~       write(*,*) "Nsub"      ,Nsub
!~       write(*,*) "i_cutback" ,i_cutback
!~     end if

    ! note P1j = sigma1j for 1D problem
    f_int(iel,:)   = f_int(iel,:)   + vec1(1) * sigma(1,:)
    f_int(iel+1,:) = f_int(iel+1,:) + vec1(2) * sigma(1,:)

    state_new%vel(iel,:,:)       = vel(:,:)
    state_new%dis_acc(iel,:,:)   = (vel(:,:) - state_old%vel(iel,:,:))/dt
    state_new%vel_x(iel,:,:)     = vel_x(:,:)
    state_new%vel_y(iel,:,:)     = vel_y(:,:)
    state_new%vel_z(iel,:,:)     = vel_z(:,:)
    state_new%cauchy(iel,:,:)  = sigma(:,:)
    state_new%Fp(iel,:,:)      = Fp_new(:,:)
    state_new%Fe(iel,:,:)      = Fe_new(:,:)
    state_new%T(iel,:)         = temp_new(:)
    state_new%tau(iel,:)       = tau(:)
    state_new%gamma_acc(iel,:,:) = state_old%gamma_acc(iel,:,:) + abs(dgamma)
    state_new%Nsub(iel)        = Nsub
    state_new%i_cutback(iel)   = i_cutback
    end do  ! iel
  end do ! iregion
  !$OMP END PARALLEL DO
    
  ! DMB: FEM Step 4 - Update position, velocity

  ! 4b external force bc's
  if(force_bc%bc_active(1)) then
    do i=1,3
      if(force_bc%dof_active(1,i)) f_ext(1,i) = force_bc%value(1,i)
    end do  !i
  end if
  if(force_bc%bc_active(2)) then
    do i=1,3
      if(force_bc%dof_active(2,i)) f_ext(Nnode,i) = force_bc%value(2,i)
    end do  !i
  end if

  ! 4a update velocity using lumped mass --
  state_new%xddot(1:Nnode,1) = (f_ext(1:Nnode,1)-f_int(1:Nnode,1)) / mass_array(1:Nnode)
  state_new%xddot(1:Nnode,2) = (f_ext(1:Nnode,2)-f_int(1:Nnode,2)) / mass_array(1:Nnode)
  state_new%xddot(1:Nnode,3) = (f_ext(1:Nnode,3)-f_int(1:Nnode,3)) / mass_array(1:Nnode)
  state_new%xdot(1:Nnode,1:3)  = state_old%xdot(1:Nnode,1:3) + dt_mid*state_new%xddot(1:Nnode,1:3)
  
  ! 4b enforce vel bc's
  if(xdot_bc%bc_active(1)) then
    do i=1,3
      if(xdot_bc%dof_active(1,i)) state_new%xdot(1,i) = min(1.d0,time(1)/xdot_acc_bc)*xdot_bc%value(1,i)
    end do  !i
  end if
  if(xdot_bc%bc_active(2)) then
    do i=1,3
      if(xdot_bc%dof_active(2,i)) state_new%xdot(Nnode,i) = min(1.d0,time(1)/xdot_acc_bc)*xdot_bc%value(2,i)
    end do  !i
  end if
  
  ! 4c update current position
  state_new%x(:,:) = state_old%x(:,:) + dt * state_new%xdot(:,:)

END SUBROUTINE DMB_Update_State
end module DMB_routines

