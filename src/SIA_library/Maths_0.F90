module maths_0

!#########################################################################

!THIS MODULE IMPLEMENTS MATHEMATICAL FUNCTIONS

!#########################################################################

!IMPLEMENTATION IN VB6 BY RAINER FEISTEL
!IMPLEMENTATION IN FORTRAN BY D.G. WRIGHT
!FOR PUBLICATION IN OCEAN SCIENCE, AS DESCRIBED IN THE PAPERS

!FEISTEL, R., WRIGHT, D.G., JACKETT, D.R., MIYAGAWA, K., REISSMANN, J.H.,
!WAGNER, W., OVERHOFF, U., GUDER, C., FEISTEL, A., MARION, G.M.:
!NUMERICAL IMPLEMENTATION AND OCEANOGRAPHIC APPLICATION OF THE THERMODYNAMIC
!POTENTIALS OF WATER, VAPOUR, ICE, SEAWATER AND AIR. PART I: BACKGROUND AND EQUATIONS. 
!OCEAN SCIENCES, 2009, IN PREPARATION.

!WRIGHT, D.G., FEISTEL, R., JACKETT, D.R., MIYAGAWA, K., REISSMANN, J.H., 
!WAGNER, W., OVERHOFF, U., GUDER, C., FEISTEL, A., MARION, G.M.:
!NUMERICAL IMPLEMENTATION AND OCEANOGRAPHIC APPLICATION OF THE THERMODYNAMIC
!POTENTIALS OF WATER, VAPOUR, ICE, SEAWATER AND AIR. PART II: THE LIBRARY ROUTINES, 
!OCEAN SCIENCES., 2009, IN PREPARATION.

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULE
!     CONSTANTS_0, FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: get_cubicroots, matrix_solve 

contains

!==========================================================================
function get_cubicroots(r, s, t, x1, x2, x3)
!==========================================================================

!INPUT:  R,S,T ARE THE COEFFICIENTS OF THE POLYNOMIAL
!        X^3 + R * X^2 + S * X + T = 0
!OUTPUT: X1, X2, X3 ARE THE ROOTS OF THE POLYNOMIAL

!RETURNS:
!        GET_CUBICROOTS = 1:  X1 REAL SOLUTION, X2 REAL PART, X3 IMAG PART OF THE COMPLEX PAIR
!        GET_CUBICROOTS = 3:  X1, X2, X3 REAL SOLUTIONS

!REQUIRES THE FUNCTION ARCCOS(X)

implicit none

integer get_cubicroots
real*8 r, s, t, x1, x2, x3
real*8 p, q
real*8 a, phi, u, v

p = s - r ** 2 / 3d0
q = 2d0 * r ** 3 / 27d0 - r * s / 3d0 + t

a = (q / 2d0) ** 2 + (p / 3d0) ** 3
if(a >= 0d0) then
  u = -q / 2d0 + sqrt(a)
  if(u >= 0d0) then
    u = u ** (1d0 / 3d0)
  else
    u = -(-u) ** (1d0 / 3d0)
  end if
  v = -q / 2d0 - sqrt(a)
  if(v >= 0d0) then
    v = v ** (1d0 / 3d0)
  else
    v = -(-v) ** (1d0 / 3d0)
  end if
  if(a == 0d0) then
  !2 equal + 1 real solution
    get_cubicroots = 3
    x1 = u + v - r / 3d0
    x2 = -(u + v) / 2d0 - r / 3d0
    x3 = x2
  else
    get_cubicroots = 1
    !2 complex solutions + 1 real solution
    x1 = u + v - r / 3d0
    x2 = -(u + v) / 2d0 - r / 3d0    !real part
    x3 = (u - v) * 0.5d0 * sqrt(3d0)   !imag part
  end if
  return
end if

!3 REAL SOLUTIONS
get_cubicroots = 3
a = sqrt(-p ** 3 / 27d0)
phi = acos(-q / (2d0 * a))

if(phi == errorreturn) then
  get_cubicroots = 0
  return
end if

a = 2d0 * a ** (1d0 / 3d0)
x1 = a * cos(phi / 3d0) - r / 3d0
x2 = a * cos((phi + 2d0 * pi) / 3d0) - r / 3d0
x3 = a * cos((phi + 4d0 * pi) / 3d0) - r / 3d0

end function

!==========================================================================
function matrix_solve(a, b, x, n)
!==========================================================================

!SOLVES A SYSTEM OF LINEAR EQUATION BY MATRIX INVERSION

implicit none

integer matrix_solve, notok, i, j, n
real*8 a(n,n), b(n), x(n)

notok = matrix_invert(a, n)
if(notok /= 0) then  !singular
  matrix_solve = notok
  return
end if

do i = 1, n
  x(i) = 0d0
  do j = 1, n
    x(i) = x(i) + a(i, j) * b(j)
  enddo
enddo

matrix_solve = 0

end function

!==========================================================================
function matrix_invert(a, n)
!==========================================================================

!INVERTS A MATRIX IN PLACE BY GAUSS ELIMINATION WITHOUT PIVOTING

implicit none

integer matrix_invert, n, i, j, k
real*8 a(n,n)

do i = 1, n
  if(a(i, i) == 0d0) then  !matrix singular
    matrix_invert = i
    return
  end if
  a(i, i) = 1d0 / a(i, i)
  do j = 1, n
    if(j /= i) then
      a(j, i) = -a(j, i) * a(i, i)
      do k = 1, n
        if(k /= i) a(j, k) = a(j, k) + a(j, i) * a(i, k)
      enddo
    end if
  enddo
  do k = 1, n
    if(k /= i) a(i, k) = a(i, i) * a(i, k)
  enddo
enddo

matrix_invert = 0  !no error

end function

end module maths_0
