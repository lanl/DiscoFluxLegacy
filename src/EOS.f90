module EOS
public
contains
!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  computes the change in temperature due to volumetric compression
SUBROUTINE EOS_Shock_Heat(props, thermo_coeffs, bulk_q, detFnew, detFold, temperature, &
           dTemp)
  
  USE GlobalParams
  USE DMB
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  TYPE(Type_DMB_EOS),               INTENT(IN)    :: props
  REAL(KIND=8), DIMENSION(3),             INTENT(IN)    :: thermo_coeffs
  REAL(KIND=8),                           INTENT(IN)    :: bulk_q
  REAL(KIND=8),                           INTENT(IN)    :: detFnew
  REAL(KIND=8),                           INTENT(IN)    :: detFold
  REAL(KIND=8),                           INTENT(IN)    :: temperature
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3),             INTENT(OUT)   :: dTemp
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8) :: rho0, Cv, te_coupling
!-----------------------------------------------------------------------

         
  rho0        = props%rho0
  Cv          = thermo_coeffs(1)  ! specific heat (const. volume)
  te_coupling = thermo_coeffs(2)  ! isentropic thermoelastic coupling coefficient, i.e., d^2(psi)/(dTdv) = (-dpdT)

  dTemp(2) = temperature / (rho0*Cv) * te_coupling * (detFnew - detFold) ! heating due to reversible thermoelastic coupling
  dTemp(3) = -bulk_q*(detFnew - detFold)/(rho0*Cv)                       ! heating due to dissipation from bulk viscosity
  dTemp(1) = dTemp(2) + dTemp(3)

  RETURN
END SUBROUTINE EOS_Shock_Heat
end module EOS
