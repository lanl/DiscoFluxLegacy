module UTILITIES
public
interface operator(.cross.)
  module procedure Cross
end interface
interface operator(.det.)
  module procedure Det3x3
end interface
interface operator(.otimes.)
  module procedure Outer_Product3
end interface
contains
!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Given vector F(X), coordinate vector X, and query location Xwant,
!>    return Fwant=F(Xwant) using linear interpolation.
SUBROUTINE Interp1D(Xwant, n, X, F, Fwant)

  IMPLICIT NONE
  REAL(KIND=8)                :: Xwant, Fwant
  INTEGER               :: n
  REAL(KIND=8), DIMENSION(n)  :: X, F
  INTENT (IN)  X, F, n, Xwant
  INTENT (OUT) Fwant

!~   REAL(KIND=8) :: b, m
  INTEGER :: i, ip1

  !If Xwant is outside the range of X, return the appropriate Fwant
  !  based on the first/last value.

  if(Xwant<=X(1)) then
    Fwant = F(1)
    return
  elseif(Xwant>=X(n)) then
    Fwant = F(n)
    return
  end if

  !Find interval i containing Xwant.

  call FindInterval(X, n, Xwant, i)

  if(i==0) then
    i = 1
  elseif(i>=n) then
    i = n - 1
  end if
  ip1 = i+1

!~   m = (F(ip1)-F(i)) / (X(ip1)-X(i))
!~   b = F(i) - m*X(i)

!~   Fwant = Xwant*m + b
  Fwant = F(i) + (Xwant - X(i))*(F(ip1)-F(i)) / (X(ip1)-X(i))
!  write(*,*) (F(ip1)-F(i)) / (X(ip1)-X(i)), &
!             (Fwant-F(i)) / (Xwant-X(i))

END SUBROUTINE Interp1D

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Given a query coordinate Xwant, find which interval in input X
!>    contains Xwant and return it as i
SUBROUTINE FindInterval(X, n, Xwant, i)

  IMPLICIT NONE
  INTEGER  :: n
  REAL(KIND=8), DIMENSION(n) :: X
  REAL(KIND=8) :: Xwant
  INTEGER  :: i
  INTENT (IN) X, n, Xwant
  INTENT (OUT) i

  INTEGER :: ilo, ihi, imid

  ilo = 0
  ihi = n+1

  do while (.true.)
    if(.not..true.) then
      return
    elseif(ihi-ilo>1) then
      imid = (ihi+ilo)/2
      if(X(n)>X(1).eqv.Xwant>X(imid)) then
        ilo = imid
      else
        ihi = imid
      end if
      cycle
    end if
    i = ilo
    exit
  end do

END SUBROUTINE FindInterval

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Generate a linearly-spaced vector from [x1,x2]
SUBROUTINE LinSpace(x1,x2,n,xout)

  IMPLICIT NONE
  REAL(KIND=8) :: x1, x2
  INTEGER :: n
  REAL(KIND=8), DIMENSION(N) :: xout
  INTENT (IN) x1, x2, n
  INTENT (OUT) xout
  
  REAL(KIND=8) ::dx
  INTEGER :: i
  
  dx = (x2 - x1) / (n-1)

  xout(1) = x1

  do i=2,n
    xout(i) = x1 + dx*(i-1)
  end do  !i

  xout(n) = x2

END SUBROUTINE LinSpace

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Integrate F(X) along X using trapezoidal rule
!>  Return Integral
SUBROUTINE Trapz(n, F, X, Integral)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL(KIND=8), INTENT(IN), DIMENSION(n)  :: F, X
  REAL(KIND=8), INTENT(OUT) :: Integral
 
  Integral = sum(0.5d0*(F(2:n)+F(1:n-1))*(x(2:n)-x(1:n-1)))

END SUBROUTINE Trapz

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Integrate F(X) along X using trapezoidal rule
!>  Return Integral, which is the cumulative integral at each interval
SUBROUTINE CumTrapz(n, F, X, Integral)

  IMPLICIT NONE
  INTEGER :: n
  REAL(KIND=8), DIMENSION(n)  :: F, X
  REAL(KIND=8), DIMENSION(n-1) :: Integral
  
  INTEGER :: i, ip1, im1
  REAL(KIND=8)  :: dx, fbar
  
  ! Compute each contribution
  
  do i=1,n-1
    ip1 = i+1
    dx = x(ip1) - x(i)
    fbar = 0.5d0*(f(i) + f(ip1))
    Integral(i) = fbar * dx
  end do  !i

  ! Integrate contributions

  do i=2,n-1
    im1 = i-1
    Integral(i) = Integral(i) + Integral(im1)
  end do  !i

END SUBROUTINE CumTrapz

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Return identity Identity(n,n)
pure function Identity(n)

  IMPLICIT NONE
  INTEGER, intent(in) :: n
  REAL(KIND=8), DIMENSION(n,n)  :: Identity
  INTEGER :: i

  Identity = 0.d0  
  do i=1,n
    Identity(i,i) = 1.d0
  end do  !i

end function Identity

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Return determinant for 3x3 REAL(KIND=8) array A
pure FUNCTION Det3x3(A)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3), INTENT(IN)  :: A
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8) :: Det3x3
!-----------------------------------------------------------------------

  Det3x3 = A(1,1)*A(2,2)*A(3,3)+A(1,2)*A(2,3)*A(3,1)+A(1,3)*A(2,1)*A(3,2) &
          -A(1,3)*A(2,2)*A(3,1)-A(1,1)*A(2,3)*A(3,2)-A(1,2)*A(2,1)*A(3,3)

END FUNCTION Det3x3

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Return determinant and inverse of 3x3 REAL(KIND=8) array A
SUBROUTINE DetInv3x3(A,D,AINV)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3,3), INTENT(IN)  :: A
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8),                 INTENT(OUT) :: D
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: AINV
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(3,3) :: CT !Transpose of cofactor matrix
  REAL(KIND=8) :: DI
!-----------------------------------------------------------------------

  D = A(1,1)*A(2,2)*A(3,3) &
     +A(1,2)*A(2,3)*A(3,1) &
     +A(1,3)*A(2,1)*A(3,2) &
     -A(1,3)*A(2,2)*A(3,1) &
     -A(1,1)*A(2,3)*A(3,2) &
     -A(1,2)*A(2,1)*A(3,3)
  if (D>=0.d0) then
    DI = 1.D0/max(1.d-99,D)
  else
    DI = -1.D0/max(1.d-99,abs(D))
  end if

  CT(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
  CT(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
  CT(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
  CT(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
  CT(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
  CT(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
  CT(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
  CT(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
  CT(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

  AINV(:,:) = CT(:,:) * DI

END SUBROUTINE DetInv3x3

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Return outer product AB of vectors a and b.
pure function Outer_Product3(a,b) result(AB)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3), INTENT(IN)  :: a, b
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3) :: AB
!-----------------------------------------------------------------------
!  Locals:
!  REAL(KIND=8), DIMENSION(3,3) :: CT !Transpose of cofactor matrix
!  REAL(KIND=8) :: DI
!-----------------------------------------------------------------------
  AB(:,1) = a(:)*b(1)
  AB(:,2) = a(:)*b(2)
  AB(:,3) = a(:)*b(3)

end function Outer_Product3
!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------
!> computes the cross product of two 3-dim vectors x and y
pure function Cross(a,b) result(axb)
  REAL(KIND=8), DIMENSION(3), INTENT(IN)  :: a,b
  REAL(KIND=8), DIMENSION(3)  :: axb
  
  axb(1) = a(2)*b(3) - a(3)*b(2)
  axb(2) = a(3)*b(1) - a(1)*b(3)
  axb(3) = a(1)*b(2) - a(2)*b(1)
 end function Cross
 
!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------
!> Compute the bracket (a,b) := a.Cmat.b, where Cmat is a tensor of 2nd order elastic constants.
SUBROUTINE ElBrak(a,b,Cmat,AB)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  REAL(KIND=8), DIMENSION(3), INTENT(IN)  :: a, b
  REAL(KIND=8), DIMENSION(3,3,3,3), INTENT(IN)  :: Cmat
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: AB
!-----------------------------------------------------------------------
  integer l,o,k,p
  AB(:,:) = 0.d0
  do p=1,3
    do o=1,3
      do l=1,3
        do k=1,3
          AB(l,o) = AB(l,o) + a(k)*Cmat(k,l,o,p)*b(p)
        end do
      end do
    end do
  end do

END SUBROUTINE ElBrak
  

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------
!>  Invert a NxN tensor.
SUBROUTINE Invert_Tensor(A,n,Ainv)
!-----------------------------------------------------------------------
  IMPLICIT NONE
!-----------------------------------------------------------------------
!  Inputs:
  INTEGER, INTENT(IN) :: n
  REAL(KIND=8), DIMENSION(n,n), INTENT(IN)  :: A
!-----------------------------------------------------------------------
!  Outputs:
  REAL(KIND=8), DIMENSION(n,n), INTENT(OUT) :: Ainv
!-----------------------------------------------------------------------
!  Locals:
  REAL(KIND=8), DIMENSION(n*n) :: work
  INTEGER, DIMENSION(n) :: ipiv
  INTEGER :: nn
  INTEGER :: info
!-----------------------------------------------------------------------
  Ainv = A
  nn=n*n

  call DGETRF(n,n,Ainv,n,ipiv,info)

  if(info/=0) then
    write(*,*) 'Invert_Tensor(): Lapack DGETRF() failed.'
    print*,A
    print*,"Ainv"
    print*,Ainv
    call Fatal('Bad matrix.')
  end if

  call DGETRI(n,Ainv,n,ipiv,work,nn,info)

  if(info/=0) then
    write(*,*) 'Invert_Tensor(): Lapack DGETRI() failed.'
    call Fatal('Bad matrix.')
  end if

END SUBROUTINE Invert_Tensor

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!~ !> compute the matrix exponential exp(A) of 3x3 matrix A
!~ !> this drop-in replacement for stdlib_linalg expm uses expokit's DGPADM routine under the hood
!~ function ExpM(A,order) result(expA)
!~ !-----------------------------------------------------------------------
!~   IMPLICIT NONE
!~ !-----------------------------------------------------------------------
!~ !  Inputs:
!~   REAL(KIND=8), DIMENSION(3,3), INTENT(IN)  :: A
!~   integer, intent(in), optional :: order
!~ !-----------------------------------------------------------------------
!~ !  Outputs:
!~   REAL(KIND=8), DIMENSION(3,3) :: expA
!~ !-----------------------------------------------------------------------
!~ !  Locals:
!~   INTEGER :: nonzero,k,j,l1,l2,ideg,ipiv(3),m,ldh,lwsp,iexph,ns,iflag
!~   REAL(KIND=8), DIMENSION(43) :: wsp
!~   REAL(KIND=8) :: rzero
!~   ideg=6; m=3; ldh=3; lwsp=43
!~   rzero = sqrt(tiny(0.d0))
!~   if (present(order)) ideg=order

!~ ! bypass (slow) matrix exponential DGPADM() when we already know the answer
!~   nonzero = count(abs(A)>rzero)
!~   if (0==nonzero) then
!~     expA = Identity(3)
!~   else if ((3==nonzero) .AND. (abs(A(1,1))>rzero) .AND. (abs(A(2,2))>rzero) .AND. (abs(A(3,3))>rzero)) then
!~     expA = Identity(3)
!~     do k=1,3
!~       expA(k,k) = exp(A(k,k))
!~     end do
!~   else
!~     call DGPADM(ideg,m,1.d0,A,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
    
!~     l1 = iexph - 1
!~     do k=1,3
!~       l2 = l1 + 3*(k-1)
!~       do j=1,3
!~         expA(j,k) = wsp(l2 + j)
!~       end do  !j
!~     end do  !k
!~   end if

!~ end function ExpM

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!> generate rotation matrix Mat that rotates by angle phi around unit vector vec
SUBROUTINE RotAround(vec,phi,Mat)
  IMPLICIT NONE
  ! Inputs
  REAL(KIND=8), INTENT(IN) :: vec(3), phi
  ! Outputs
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: Mat
  ! Locals
  REAL(KIND=8) :: vx(3,3), eye(3,3), s, c
  !-----------------------------------
  vx(3,3)=0.d0
  s = sin(phi)
  c = cos(phi)
  eye = Identity(3)
  vx(1,2) = vec(3)
  vx(2,1) = -vec(3)
  vx(1,3) = -vec(2)
  vx(3,1) = vec(2)
  vx(2,3) = vec(1)
  vx(3,2) = -vec(1)
  Mat = eye + s*vx + matmul(vx,vx)*(1.d0-c)

END SUBROUTINE RotAround

!> generate rotation matrix Mat that rotates unit vector old into unit vector new
SUBROUTINE RotInto(old,new,Mat)
  IMPLICIT NONE
  ! Inputs
  REAL(KIND=8), DIMENSION(3), INTENT(IN) :: old,new
  ! Outputs
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: Mat
  ! Locals
  REAL(KIND=8) :: v(3), vx(3,3), eye(3,3), s, c
  !-----------------------------------
  vx=0.d0
  v = old .cross. new
  eye = Identity(3)
  s = sqrt(dot_product(v,v))
  c = dot_product(old,new)
  vx(1,2) = -v(3)
  vx(2,1) = v(3)
  vx(1,3) = v(2)
  vx(3,1) = -v(2)
  vx(2,3) = -v(1)
  vx(3,2) = v(1)
  if (abs(s)<1.d-2) then
    Mat = eye + vx + matmul(vx,vx)/(1.d0+c)
  else
    Mat = eye + vx + matmul(vx,vx)*(1.d0-c)/s**2
  end if

END SUBROUTINE RotInto

!> generate rotation matrix that aligns unit vectors n with y and t with z, assuming n.t=0=y.z
SUBROUTINE RotAlign(n,t,y,z,Mat)
  USE GlobalParams
  IMPLICIT NONE
  ! Inputs
  REAL(KIND=8), DIMENSION(3), INTENT(IN) :: n,t,y,z
  ! Outputs
  REAL(KIND=8), DIMENSION(3,3), INTENT(OUT) :: Mat
  ! Locals
  REAL(KIND=8) :: rot1(3,3), newt(3)
  !-----------------------------------
  if (abs(dot_product(n,y)-1.d0)<1.d-15) then
    call RotAround(z,0.d0,rot1)
  elseif (abs(dot_product(n,y)+1.d0)<1.d-15) then
    call RotAround(z,pi,rot1)
  else
    call RotInto(n,y,rot1)
  end if
  newt = matmul(rot1,t)
  if (abs(dot_product(newt,z)-1.d0)<1.d-15) then
    call RotAround(y,0.d0,Mat)
  elseif (abs(dot_product(newt,z)+1.d0)<1.d-15) then
    call RotAround(y,pi,Mat)
  else
    call RotInto(newt,z,Mat)
  end if
  Mat = matmul(Mat,rot1)

END SUBROUTINE RotAlign

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

SUBROUTINE readcmdline(jobname)
  USE GlobalParams
!$   Use omp_lib
  IMPLICIT NONE
  ! Output
  CHARACTER(32), INTENT(OUT) :: jobname
  ! Locals
  CHARACTER(32) :: cmdlinearg, exe_name
  !----------------------------------------
  call get_command_argument(1, cmdlinearg)
  call get_command_argument(0, exe_name)
  if (len_trim(cmdlinearg) > 0) then
    jobname = cmdlinearg
    if ((cmdlinearg == '--version') .or. (cmdlinearg == '-v')) then
      print*,prog_version
      stop
    elseif ((cmdlinearg == '--help') .or. (cmdlinearg == '-h')) then
      print*,"USAGE: ",trim(exe_name)," [--version] [--help] <jobname>"
      print*,"       (needed to read 'input_parameters.<jobname>.dat')"
      print*,"       or ",trim(exe_name)," <inputfilename>"
      print*,""
      print*,"program version: ",prog_version
      print*,"simulation types: impact or shear"
      print*,"compiled with these defaults:"
      print*,"# of elements, Nel = ",Nel
      print*,"# of regions, Nregion = ",Nregion
      print*,"# of slip systems, Nslip = ",Nslip
      print*,"# of dislocation character angles, Nchar = ",Nchar
!$    print*, 'OpenMP parallelization enabled:'
!$    print*, 'using ',omp_get_max_threads(),' of ',omp_get_num_procs(),' processors'
!$    print*, 'type "export OMP_NUM_THREADS=n" before running this prog. to change'
      stop
    end if
  else
    print*,"usage: ",trim(exe_name)," [--version] [--help] <jobname> (needed to read 'input_parameters.<jobname>.dat')"
    stop 1
  end if
END SUBROUTINE readcmdline

!-----------------------------------------------------------------------
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
!-----------------------------------------------------------------------

!>  Termination routine.  Print the message and stop.
SUBROUTINE Fatal(message)

  IMPLICIT NONE
  CHARACTER*(*) message
!-----------------------------------------------------------------------
  write(*,10) message
10  format(                         &
      20('*'),/,                    &
      'Fatal error encountered:',/, &
      '----> ',a,/,                 &
      'Terminating',/,              &
      20('*'))

  STOP

END SUBROUTINE Fatal
end module UTILITIES
