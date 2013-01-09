module sal_1

!#########################################################################

!THIS MODULE IMPLEMENTS THE SALINE PART OF THE GIBBS POTENTIAL OF SEAWATER
!IN THE FORM OF 7 COEFFICIENTS OF ITS POWER EXPANSION IN SALINITY,
!And THEIR FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE
!AND PRESSURE AS DEFINED IN IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION FOR THE THERMODYNAMIC PROPERTIES OF SEAWATER
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

!#########################################################################

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

!THIS MODULE REQUIRES MODULE
!     CONSTANTS_0, FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sal_g_term_si

real*8, public :: gsc(0:7,0:6,0:5), gsi(0:7,0:6,0:5)

integer, private, parameter :: maxs = 7, maxt = 6, maxp = 5

contains

!==========================================================================
function sal_g_term_si(term, drv_t, drv_p, t_si, p_si)
!==========================================================================

!RETURNS THE T-P DERIVATIVE OF THE I-TH TERM  (D/DT)^DRV_T (D/DP)^DRV_P GI(T, P)
!OF THE SALINITY EXPANSION OF THE SALINE GIBBS FUNCTION GS,

! GS(S,T,P) = G1(T,P)*S*LN(S) + SUM(I>1) GI(T,P)*S^(I/2)
! WHERE S = ABSOLUTE SALINITY IN KG/KG, I.E. THE MASS FRACTION OF SALT IN SEAWATER

! T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
! P_SI  = P = ABSOLUTE PRESSURE IN PA
! GS = SPECIFIC GIBBS ENERGY IN J/KG, RELATIVE TO PURE WATER

!CHECK VALUES WITH DEFAULT SETTINGS:
!SAL_G_TERM_SI(1,0,0,300,1E6)= 79427.9694846
!SAL_G_TERM_SI(1,1,0,300,1E6)= 264.759898282
!SAL_G_TERM_SI(1,0,1,300,1E6)= 0
!SAL_G_TERM_SI(1,2,0,300,1E6)= 0
!SAL_G_TERM_SI(1,1,1,300,1E6)= 0
!SAL_G_TERM_SI(1,0,2,300,1E6)= 0

!SAL_G_TERM_SI(2,0,0,300,1E6)= 301223.934546
!SAL_G_TERM_SI(2,1,0,300,1E6)= 1558.10730393
!SAL_G_TERM_SI(2,0,1,300,1E6)= -7.60101106780E-04
!SAL_G_TERM_SI(2,2,0,300,1E6)= 19.0397747694
!SAL_G_TERM_SI(2,1,1,300,1E6)= 9.35496617545E-07
!SAL_G_TERM_SI(2,0,2,300,1E6)= 1.28077803426E-12

!SAL_G_TERM_SI(3,0,0,300,1E6)= -345570.781498
!SAL_G_TERM_SI(3,1,0,300,1E6)= -1749.19911587
!SAL_G_TERM_SI(3,0,1,300,1E6)= 2.01608449881E-04
!SAL_G_TERM_SI(3,2,0,300,1E6)= -9.32088823620
!SAL_G_TERM_SI(3,1,1,300,1E6)= 2.17412289837E-08
!SAL_G_TERM_SI(3,0,2,300,1E6)= -4.78063827318E-13

!SAL_G_TERM_SI(4,0,0,300,1E6)= 1468073.64559
!SAL_G_TERM_SI(4,1,0,300,1E6)= 7741.24404962
!SAL_G_TERM_SI(4,0,1,300,1E6)= -4.33952339916E-04
!SAL_G_TERM_SI(4,2,0,300,1E6)= -6.13689642370
!SAL_G_TERM_SI(4,1,1,300,1E6)= -3.50876195080E-06
!SAL_G_TERM_SI(4,0,2,300,1E6)= -6.06204383305E-13

!SAL_G_TERM_SI(5,0,0,300,1E6)= -3776969.31546
!SAL_G_TERM_SI(5,1,0,300,1E6)= -15135.6522248
!SAL_G_TERM_SI(5,0,1,300,1E6)= 1.11272425476E-03
!SAL_G_TERM_SI(5,2,0,300,1E6)= 0
!SAL_G_TERM_SI(5,1,1,300,1E6)= 0
!SAL_G_TERM_SI(5,0,2,300,1E6)= 0

!SAL_G_TERM_SI(6,0,0,300,1E6)= 6151235.69234
!SAL_G_TERM_SI(6,1,0,300,1E6)= 14157.0509933
!SAL_G_TERM_SI(6,0,1,300,1E6)= 0
!SAL_G_TERM_SI(6,2,0,300,1E6)= 0
!SAL_G_TERM_SI(6,1,1,300,1E6)= 0
!SAL_G_TERM_SI(6,0,2,300,1E6)= 0

!SAL_G_TERM_SI(7,0,0,300,1E6)= -3734033.38866
!SAL_G_TERM_SI(7,1,0,300,1E6)= 0
!SAL_G_TERM_SI(7,0,1,300,1E6)= 0
!SAL_G_TERM_SI(7,2,0,300,1E6)= 0
!SAL_G_TERM_SI(7,1,1,300,1E6)= 0
!SAL_G_TERM_SI(7,0,2,300,1E6)= 0


!check values with Public Const IsExtension2010 = True:

!SAL_G_TERM_SI(1,0,0,300,1E6)= 79427.9694846
!SAL_G_TERM_SI(1,1,0,300,1E6)= 264.759898282
!SAL_G_TERM_SI(1,0,1,300,1E6)= 0
!SAL_G_TERM_SI(1,2,0,300,1E6)= 0
!SAL_G_TERM_SI(1,1,1,300,1E6)= 0
!SAL_G_TERM_SI(1,0,2,300,1e6)= 0

!SAL_G_TERM_SI(2,0,0,300,1E6)= 301223.934546274
!SAL_G_TERM_SI(2,1,0,300,1E6)= 1558.10730392553
!SAL_G_TERM_SI(2,0,1,300,1E6)= -7.60101106780E-04
!SAL_G_TERM_SI(2,2,0,300,1E6)= 19.0397747694
!SAL_G_TERM_SI(2,1,1,300,1E6)= 9.35496617545E-07
!SAL_G_TERM_SI(2,0,2,300,1E6)= 1.28077803426E-12

!SAL_G_TERM_SI(3,0,0,300,1E6)= -345570.781498
!SAL_G_TERM_SI(3,1,0,300,1E6)= -1749.19911587
!SAL_G_TERM_SI(3,0,1,300,1E6)= 2.01608449881E-04
!SAL_G_TERM_SI(3,2,0,300,1E6)= -9.3208882362
!SAL_G_TERM_SI(3,1,1,300,1E6)= 2.17412289837E-08
!SAL_G_TERM_SI(3,0,2,300,1E6)= -4.78063827318E-13

!SAL_G_TERM_SI(4,0,0,300,1E6)= 1468073.64559
!SAL_G_TERM_SI(4,1,0,300,1E6)= 7741.24404962
!SAL_G_TERM_SI(4,0,1,300,1E6)= -4.33952339916E-04
!SAL_G_TERM_SI(4,2,0,300,1E6)= -6.13689642370
!SAL_G_TERM_SI(4,1,1,300,1E6)= -3.50876195080E-06
!SAL_G_TERM_SI(4,0,2,300,1E6)= -6.06204383305E-13

!SAL_G_TERM_SI(5,0,0,300,1E6)= -3776969.31546
!SAL_G_TERM_SI(5,1,0,300,1E6)= -15135.6522248
!SAL_G_TERM_SI(5,0,1,300,1E6)= 1.11272425476E-03
!SAL_G_TERM_SI(5,2,0,300,1E6)= 0
!SAL_G_TERM_SI(5,1,1,300,1E6)= 0
!SAL_G_TERM_SI(5,0,2,300,1E6)= 0

!SAL_G_TERM_SI(6,0,0,300,1E6)= 6151235.69234
!SAL_G_TERM_SI(6,1,0,300,1E6)= 14157.0509933
!SAL_G_TERM_SI(6,0,1,300,1E6)= 0
!SAL_G_TERM_SI(6,2,0,300,1E6)= 0
!SAL_G_TERM_SI(6,1,1,300,1E6)= 0
!SAL_G_TERM_SI(6,0,2,300,1E6)= 0

!SAL_G_TERM_SI(7,0,0,300,1E6)= -3734033.38866
!SAL_G_TERM_SI(7,1,0,300,1E6)= 0
!SAL_G_TERM_SI(7,0,1,300,1E6)= 0
!SAL_G_TERM_SI(7,2,0,300,1E6)= 0
!SAL_G_TERM_SI(7,1,1,300,1E6)= 0
!SAL_G_TERM_SI(7,0,2,300,1E6)= 0

implicit none

integer term, drv_t, drv_p
real*8 t_si, p_si
real*8 sal_g_term_si

call init_iapws08

sal_g_term_si = errorreturn

if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (drv_t < 0) return
if (drv_p < 0) return

if (term < 1 .or. term > maxs) then
  sal_g_term_si = 0d0
  return
end if

if(check_limits == 1) then
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax &
    .or. p_si < sal_pmin .or. p_si > sal_pmax) return
endif

sal_g_term_si = poly_gtp(term, drv_t, drv_p, t_si, p_si)

end function

!==========================================================================
function poly_gtp(term, drv_t, drv_p, t_si, p_si)
!==========================================================================

!RETURNS THE T-P DERIVATIVE OF THE I-TH TERM  (D/DT)^DRV_T (D/DP)^DRV_P GI(T, P)
!OF THE SALINITY EXPANSION OF THE SALINE GIBBS FUNCTION GS,

!GS(S,T,P) = G1(T,P)*S*LN(S) + SUM(I>1) GI(T,P)*S^(I/2)

!S = ABSOLUTE SALINITY IN KG/KG, I.E. THE MASS FRACTION OF SALT IN SEAWATER
!T = ABSOLUTE TEMPERATURE IN K, ITS-90
!P = ABSOLUTE PRESSURE IN PA

implicit none
integer term, drv_t, drv_p
real*8 poly_gtp
real*8 tu, pu, t_si, p_si, y, z

tu = 40d0
pu = 100000000d0

!REDUCED INPUT VALUES OF T, P
y = (t_si - celsius_temperature_si) / tu
z = (p_si - sealevel_pressure_si) / pu

poly_gtp = poly_gyz(term, drv_t, drv_p, y, z) / (tu**drv_t * pu**drv_p)

end function

!==========================================================================
function poly_gyz(term, drv_y, drv_z, y, z)
!==========================================================================

!RETURNS THE VALUE OF THE POLYNOMIAL DERIVATIVE
!(D/DY)^DRV_Y (D/DZ)^DRV_Z SUM(J,K) GI(TERM,J,K)*Y^J*Z^K

implicit none
integer term, drv_y, drv_z
integer j, k, l, jmax, kmax
real*8 c, y, z
real*8 g, poly_gyz, yj, zk

g = 0d0
if (y == 0d0) then 
  jmax = drv_y 
else 
  jmax = maxt
endif

if (z == 0d0) then 
  kmax = drv_z 
else 
  kmax = maxp
endif

yj = 1d0
do j = drv_y, jmax   !loop over powers of y

  zk = 1d0
  do k = drv_z, kmax    !loop over powers of z

    if (gsi(term, j, k) /= 0d0) then
      c = gsi(term, j, k) * yj * zk

      do l = 1, drv_y            !factors from y-derivatives
        c = c * dfloat(j - l + 1)
      enddo

      do l = 1, drv_z          !factors from z-derivatives
        c = c * dfloat(k - l + 1)
      enddo

      g = g + c
    end if

    if (k < kmax) zk = zk * z
  enddo !k loop

  if (j < jmax) yj = yj * y
enddo

poly_gyz = g

end function

!==============================================================================
subroutine init_iapws08
!==============================================================================

!==============================================================================
!G(S,T,P)=(GSC(1,0,0) + GSC(1,1,0)*Y)*X^2*LN(X) + SUM GSC(I,J,K) * X^I * Y^J * Z^K 
!S = 40 PSS * X^2, T = 40 C * Y, P = 100 MPA * Z
!==============================================================================

implicit none

integer i, j, k
real*8 gsc100
real*8 su, ln_su, sqr_su, ups

gsc100 = 0.581281456626732d+4

if (gsc(1, 0, 0) == gsc100) return

!ERROR FLAG

!ENTHALPY=0 AND ENTROPY=0 STANDARD SEAWATER
!G(35,0,0) = 0:
!DG(35,0,0)/DT = 0:

!F08 COEFFICIENTS FOR G_S
gsc(1, 0, 0) = 0.581281456626732d+4
gsc(1, 1, 0) = 0.851226734946706d+3

gsc(2, 0, 0) = 0.141627648484197d+4
gsc(2, 1, 0) = 0.168072408311545d+3

gsc(3, 0, 0) =-0.243214662381794d+4
gsc(4, 0, 0) = 0.202580115603697d+4
gsc(5, 0, 0) =-0.109166841042967d+4
gsc(6, 0, 0) = 0.37460123787784d+3
gsc(7, 0, 0) =-0.485891069025409d+2
gsc(3, 1, 0) =-0.493407510141682d+3
gsc(4, 1, 0) = 0.543835333000098d+3
gsc(5, 1, 0) =-0.196028306689776d+3
gsc(6, 1, 0) = 0.367571622995805d+2
gsc(2, 2, 0) = 0.880031352997204d+3
gsc(3, 2, 0) =-0.430664675978042d+2
gsc(4, 2, 0) =-0.685572509204491d+2
gsc(2, 3, 0) =-0.225267649263401d+3
gsc(3, 3, 0) =-0.100227370861875d+2
gsc(4, 3, 0) = 0.493667694856254d+2
gsc(2, 4, 0) = 0.914260447751259d+2
gsc(3, 4, 0) = 0.875600661808945d0
gsc(4, 4, 0) =-0.171397577419788d+2
gsc(2, 5, 0) =-0.216603240875311d+2
gsc(4, 5, 0) = 0.249697009569508d+1
gsc(2, 6, 0) = 0.213016970847183d+1

gsc(2, 0, 1) =-0.331049154044839d+4
gsc(3, 0, 1) = 0.199459603073901d+3
gsc(4, 0, 1) =-0.547919133532887d+2
gsc(5, 0, 1) = 0.360284195611086d+2
gsc(2, 1, 1) = 0.729116529735046d+3
gsc(3, 1, 1) =-0.175292041186547d+3
gsc(4, 1, 1) =-0.226683558512829d+2
gsc(2, 2, 1) =-0.860764303783977d+3
gsc(3, 2, 1) = 0.383058066002476d+3
gsc(2, 3, 1) = 0.694244814133268d+3
gsc(3, 3, 1) =-0.460319931801257d+3
gsc(2, 4, 1) =-0.297728741987187d+3
gsc(3, 4, 1) = 0.234565187611355d+3

gsc(2, 0, 2) = 0.384794152978599d+3
gsc(3, 0, 2) =-0.522940909281335d+2
gsc(4, 0, 2) =-0.408193978912261d+1
gsc(2, 1, 2) =-0.343956902961561d+3
gsc(3, 1, 2) = 0.831923927801819d+2
gsc(2, 2, 2) = 0.337409530269367d+3
gsc(3, 2, 2) =-0.541917262517112d+2
gsc(2, 3, 2) =-0.204889641964903d+3
gsc(2, 4, 2) = 0.74726141138756d+2

gsc(2, 0, 3) =-0.965324320107458d+2
gsc(3, 0, 3) = 0.680444942726459d+2
gsc(4, 0, 3) =-0.301755111971161d+2
gsc(2, 1, 3) = 0.124687671116248d+3
gsc(3, 1, 3) =-0.29483064349429d+2
gsc(2, 2, 3) =-0.178314556207638d+3
gsc(3, 2, 3) = 0.256398487389914d+2
gsc(2, 3, 3) = 0.113561697840594d+3
gsc(2, 4, 3) =-0.364872919001588d+2

gsc(2, 0, 4) = 0.158408172766824d+2
gsc(3, 0, 4) =-0.341251932441282d+1
gsc(2, 1, 4) =-0.31656964386073d+2
gsc(2, 2, 4) = 0.442040358308d+2
gsc(2, 3, 4) =-0.111282734326413d+2

gsc(2, 0, 5) =-0.262480156590992d+1
gsc(2, 1, 5) = 0.704658803315449d+1
gsc(2, 2, 5) =-0.792001547211682d+1

ups = so_salinity_si / 35d0
su = 40d0 * ups

If(IsExtension2010) then
  !R. Feistel
  !Extended Equation of State for Seawater at Elevated Temperature and Salinity
  !Desalination, 2010, VOL. 250, 14-18.
  
  !to cover the measurements of
  !F. Millero and F. Huang
  !Ocean Sci. Discuss. 6 (2009) 153-169.
  !http://www.ocean-sci-discuss.net/6/153/2009/
  
  gsc(2, 1, 1) = gsc(2, 1, 1) - 0.347600838235511d2
  gsc(3, 1, 1) = gsc(3, 1, 1) + 0.935178208024272d2
  gsc(4, 1, 1) = gsc(4, 1, 1) - 0.603483495593212d2
  gsc(2, 2, 1) = gsc(2, 2, 1) + 0.228081199116236d2
  gsc(4, 2, 1) = gsc(4, 2, 1) - 0.242869748847311d2
  gsc(2, 3, 1) = gsc(2, 3, 1) - 0.560725782144008d2
  gsc(3, 3, 1) = gsc(3, 3, 1) - 0.143523729429211d2
  gsc(4, 3, 1) = gsc(4, 3, 1) + 0.926367388049097d2
  gsc(4, 4, 1) = gsc(4, 4, 1) - 0.416658900599273d2
  gsc(2, 5, 1) = gsc(2, 5, 1) + 0.645288813326254d2
  gsc(3, 5, 1) = gsc(3, 5, 1) - 0.403505133068118d2
  gsc(2, 6, 1) = gsc(2, 6, 1) - 4.32746069361075d0
  gsc(4, 6, 1) = gsc(4, 6, 1) + 2.05954716712622d0
End If

!=======================================================================
!MODIFIED COEFFICIENTS AS REQUIRED FOR COMPUTING EACH OF THE TERMS G1 TO G7
!OF THE SALINITY EXPANSION OF GS, STORED INTO THE GLOBAL ARRAY GSI()
!=======================================================================
ln_su = log(su)
sqr_su = sqrt(su)

do j = 0, maxt
  do k = 0, maxp
    gsi(1, j, k) = gsc(1, j, k) * 0.5d0 / su
    gsi(2, j, k) = (gsc(2, j, k) - 0.5d0 * gsc(1, j, k) * ln_su) / su
    do i = 3, maxs
      gsi(i, j, k) = gsc(i, j, k) / sqr_su ** i
    enddo
  enddo
enddo

end subroutine

end module sal_1
