module sal_2

!#########################################################################

!THIS MODULE IMPLEMENTS THE SALINE PART OF THE GIBBS FUNCTION OF SEAWATER,
!SAL_G_SI(DRV_S, DRV_T, DRV_P, SA_SI, T_SI, P_SI) AS WELL AS ITS 1ST AND 2ND
!PARTIAL DERIVATIVES WITH RESPECT TO SALINITY, TEMPERATURE AND PRESSURE,
!AS WELL AS THERMODYNAMIC PROPERTIES OF SEA SALT IN SEAWATER AS FUNCTIONS
!COMPUTED FROM THE SALINITY EXPANSION TERMS OF THE GIBBS POTENTIAL,
!AS DEFINED IN

!RELEASE ON THE IAPWS FORMULATION 2008 FOR THE THERMODYNAMIC PROPERTIES OF 
!SEAWATER.  THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
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

!THIS MODULE REQUIRES THE LIBRARY MODULE 
!     CONSTANTS_0,    FILE CONSTANTS_0.F90
!     SAL_1,          FILE SAL_1.F90 

!#########################################################################

use constants_0
use sal_1

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sal_act_coeff_si, sal_act_potential_si, sal_activity_w_si, &
          sal_chem_coeff_si, sal_chempot_h2o_si, sal_chempot_rel_si, &
          sal_dilution_si, sal_g_si, sal_mixenthalpy_si, sal_mixentropy_si, &
          sal_mixvolume_si, sal_molality_si, sal_osm_coeff_si, &
          sal_saltenthalpy_si, sal_saltentropy_si, sal_saltvolume_si

integer, private, parameter :: maxs = 7

contains

!==========================================================================
function sal_g_si(drv_s, drv_t, drv_p, sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE S-T-P DERIVATIVE (D/DS)^DRV_S (D/DT)^DRV_T (D/DP)^DRV_P GS(S,T,P)
!OF THE SALINE GIBBS FUNCTION GS,
!
! GS(S,T,P) = G1(T,P)*S*LN(S) + SUM(I>1) GI(T,P)*S^(I/2)
!
! SA_SI = S = ABSOLUTE SALINITY IN KG/KG, I.E. THE MASS FRACTION OF SALT IN SEAWATER
! T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
! P_SI  = P = ABSOLUTE PRESSURE IN PA
! SAL_G_SI = GS = SPECIFIC GIBBS ENERGY IN J/KG, RELATIVE TO PURE WATER

!CHECK VALUES:
!SAL_G_SI(0,0,0,0.035,300,1E6)= 127.033640309
!SAL_G_SI(1,0,0,0.035,300,1E6)= 77949.2100395
!SAL_G_SI(0,1,0,0.035,300,1E6)= 18.6360407073
!SAL_G_SI(0,0,1,0.035,300,1E6)= -2.55600080319E-05
!SAL_G_SI(2,0,0,0.035,300,1E6)= 2248200.54660
!SAL_G_SI(1,1,0,0.035,300,1E6)= 790.563810558
!SAL_G_SI(1,0,1,0.035,300,1E6)= -7.15686520588E-04
!SAL_G_SI(0,2,0,0.035,300,1E6)= 0.597842170749
!SAL_G_SI(0,1,1,0.035,300,1E6)= 2.85865076269E-08
!SAL_G_SI(0,0,2,0.035,300,1E6)= 4.09543164905E-14

implicit none

real*8 sal_g_si
integer drv_s, drv_t, drv_p, i
real*8 sa_si, t_si, p_si
real*8 sqr_s, gi, g

dimension gi(7)

sal_g_si = errorreturn

if (drv_s < 0 .or. drv_s > 3) return
if (drv_t < 0) return
if (drv_p < 0) return

if(check_limits == 1) then
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or.  &
  sa_si < sal_smin .or. sa_si > sal_smax .or.  &
  p_si < sal_pmin .or. p_si > sal_pmax) return
else
  if (sa_si < 0d0 .or. sa_si >= 1d0) return
  if (t_si <= 0d0) return
  if (p_si <= 0d0) return
endif

!COMPUTE THE LIMIT S = 0:
if (sa_si == 0d0) then

  if (drv_s == 0) then                                          !IF NO S-DERIVATIVE TAKEN, EACH TERM OF GS VANISHES
    sal_g_si = 0d0
    return
  end if
  if ((drv_p > 0 .or. drv_t > 1) .and. drv_s == 1) then   !IF THE LOG TERM IS GONE, THE 1ST S-DERIVATIVE
    sal_g_si = sal_g_term_si(2, drv_t, drv_p, t_si, p_si)        !RETURNS THE LINEAR S TERM
    return
  end if
  return                                                          !ALL OTHER CASES HAVE A SINGULARITY AT S = 0

end if

!COMPUTE THE CASES S > 0:
do i = 1, maxs
  gi(i) = sal_g_term_si(i, drv_t, drv_p, t_si, p_si)             !GET ALL COEFFICIENTS OF THE SALINITY EXPANSION
  if (gi(i) == errorreturn) return
enddo

g = 0d0
if (gi(1) /= 0d0) then                                          !TAKE THE REQUIRED S-DERIVATIVE OF THE LOG TERM
  select case (drv_s)
    case(0); g = gi(1) * sa_si * log(sa_si)
    case(1); g = gi(1) * (log(sa_si) + 1d0)
    case(2); g = gi(1) / sa_si
    case(3); g = -gi(1) / sa_si**2
  end select
end if

sqr_s = sqrt(sa_si)

do i = 2, maxs                                                    !ADD THE S-DERIVATIVES OF EACH ROOT(S) TERM
  if (gi(i) /= 0d0) then
    select case(drv_s)
      case(0); g = g + gi(i) * sqr_s**i
      case(1); g = g + gi(i)*0.5d0*dfloat(i)*sqr_s**(i - 2)
      case(2); g = g + gi(i)*0.5d0*dfloat(i)*(0.5d0*dfloat(i) - 1d0)  &
        *sqr_s**(i - 4)
      case(3); g = g + gi(i)*0.5d0*dfloat(i)*(0.5d0*dfloat(i) - 1d0)  &
        *(0.5d0*dfloat(i) - 2d0)*sqr_s**(i - 6)
    end select
  end if
enddo

sal_g_si = g

end function

!==========================================================================
function sal_act_coeff_si(sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE MEAN ACTIVITY COEFFICIENT, LN(GAMMA), OF SEAWATER

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_ACT_COEFF_SI(0.035,300,1E6) = -0.527003008913

implicit none

integer i
real*8 sal_act_coeff_si
real*8 sa_si, t_si, p_si
real*8 lng, sqr_s
real*8 gi

dimension gi(7)

sal_act_coeff_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (sa_si == 0d0) then
  sal_act_coeff_si = 0d0
  return
end if

!COMPUTE THE CASES S > 0:
do i = 1, maxs
  if (i /= 2) then
    gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !GET ALL COEFFICIENTS OF THE SALINITY EXPANSION
    if (gi(i) == errorreturn) return
  end if
enddo

sqr_s = sqrt(sa_si)
lng = 0d0

do i = 3, maxs
  lng = lng + gi(i)*(sa_si*(1d0 - 0.5d0*dfloat(i)) + 0.5d0*dfloat(i))  &
        *sqr_s**(i - 2)
enddo

sal_act_coeff_si = lng / gi(1) + log(1d0 - sa_si) - sa_si

end function

!==========================================================================
function sal_act_potential_si(sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE ACTIVITY POTENTIAL, PSI, OF SEAWATER

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_ACT_POTENTIAL_SI(0.035,300,1E6) = -0.429940465498

implicit none

integer i
real*8 sal_act_potential_si
real*8 sa_si, t_si, p_si
real*8 sqr_s
real*8 gi
real*8 psi

dimension gi(7)

sal_act_potential_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (sa_si == 0d0) then
  sal_act_potential_si = 0d0
  return
end if

!COMPUTE THE CASES S > 0:
do i = 1, maxs
  if (i /= 2) then
    gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !GET ALL COEFFICIENTS OF THE SALINITY EXPANSION
    if (gi(i) == errorreturn) return
  end if
enddo

sqr_s = sqrt(sa_si)
psi = 0d0

do i = 3, maxs
  psi = psi + gi(i)*sqr_s**(i - 2)
enddo

sal_act_potential_si = psi/gi(1) + log(1d0 - sa_si)

end function

!==========================================================================
function sal_activity_w_si(sa_si, t_si, p_si)
!==========================================================================

!FUNCTION RETURNS THE ACTIVITY OF WATER IN SEAWATER, A_W AS A FUNCTION OF 
!ABSOLUTE SALINITY IN KG/KG FOR SEAWATER WITH REFERENCE COMPOSITION, SEE 
!MILLERO ET AL., DEEP-SEA RESEARCH I, 55(2008) 50-72.


implicit none

real*8 sa_si, t_si, p_si, sal_activity_w_si
real*8 phi

!RETURNS   ACTIVITY OF WATER IN SEAWATER
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_ACTIVITY_W_SI(0.035,300,1E6) = 0.981388410188

if(sa_si == 0d0) then
   sal_activity_w_si = 0d0
   return
endif

phi = sal_osm_coeff_si(sa_si, t_si, p_si)

if (phi == errorreturn)then
   sal_activity_w_si = errorreturn
   return
endif

!Begin modified code
!dbg  sal_activity_w_si = exp(-phi*sa_si/(1d0-sa_si))
sal_activity_w_si = exp(-phi*sa_si/(1d0-sa_si)  &
                  * molar_mass_h2o_si/molar_mass_seasalt_si)
!End modified code

end function

!==========================================================================
function sal_dilution_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sa_si, t_si, p_si, sal_dilution_si
real*8 gi
real*8 d, sqr_s

dimension gi(7)

!RETURNS   DILUTION COEFFICIENT OF SEAWATER, D = S * (D2G/DS2)_T_P, IN J/KG
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_DILUTION_SI(0.035,300,1E6) = 78687.0191309

sal_dilution_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (sa_si == 0d0) then
  sal_dilution_si = sal_g_term_si(1, 0, 0, t_si, p_si)
end if

do i = 1, maxs
  if (i /= 2) then
    gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !GET ALL COEFFICIENTS OF THE SALINITY EXPANSION
    if (gi(i) == errorreturn) return
  end if
enddo

sqr_s = sqrt(sa_si)
d = gi(1)

do i = 3, maxs
  d = d + gi(i)*0.25d0*dfloat(i)*dfloat(i - 2)*sqr_s**(i - 2)
enddo

sal_dilution_si = d

end function

!==========================================================================
function sal_mixenthalpy_si(sa1_si, sa2_si, w1, t_si, p_si)
!==========================================================================

implicit none

real*8 sa1_si, sa2_si, t_si, p_si, sal_mixenthalpy_si
real*8 w1, w2
real*8 h, h1, h2, s, s1, s2, s12

!RETURNS   MIXING ENTHALPY IN J/KG

!SA1_SI    ABSOLUTE SALINITY OF COMPONENT 1, IN KG/KG
!SA2_SI    ABSOLUTE SALINITY OF COMPONENT 2, IN KG/KG
!W1        MASS FRACTION OF COMPONENT 1, W1 = M1/(M1+M2)
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_MIXENTHALPY_SI(0.01,0.035,0.6,300,1E6) = 16.1539911284

sal_mixenthalpy_si = errorreturn

if (sa1_si < 0d0 .or. sa1_si >= 1d0) return
if (sa2_si < 0d0 .or. sa2_si >= 1d0) return
if (w1 < 0d0 .or. w1 > 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

w2 = 1d0 - w1
s12 = w1 * sa1_si + w2 * sa2_si
h = sal_g_si(0, 0, 0, s12, t_si, p_si)
if (h == errorreturn) return
s = -sal_g_si(0, 1, 0, s12, t_si, p_si)
if (s == errorreturn) return
h = h + t_si * s

h1 = sal_g_si(0, 0, 0, sa1_si, t_si, p_si)
if (h1 == errorreturn) return
s1 = -sal_g_si(0, 1, 0, sa1_si, t_si, p_si)
if (s1 == errorreturn) return
h1 = h1 + t_si * s1

h2 = sal_g_si(0, 0, 0, sa2_si, t_si, p_si)
if (h2 == errorreturn) return
s2 = -sal_g_si(0, 1, 0, sa2_si, t_si, p_si)
if (s2 == errorreturn) return
h2 = h2 + t_si * s2

sal_mixenthalpy_si = h - w1*h1 - w2*h2

end function

!==========================================================================
function sal_mixentropy_si(sa1_si, sa2_si, w1, t_si, p_si)
!==========================================================================

implicit none

real*8 sa1_si, sa2_si, w1, w2, t_si, p_si, sal_mixentropy_si
real*8 s, s1, s2, s12

!RETURNS   MIXING ENTROPY OF SEAWATER IN J/KG

!SA1_SI    ABSOLUTE SALINITY OF COMPONENT 1, IN KG/KG
!SA2_SI    ABSOLUTE SALINITY OF COMPONENT 2, IN KG/KG
!W1        MASS FRACTION OF COMPONENT 1, W1 = M1/(M1+M2)
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_MIXENTROPY_SI(0.01,0.035,0.6,300,1E6) = 0.96682942261|7

sal_mixentropy_si = errorreturn

if (sa1_si < 0d0 .or. sa1_si >= 1d0)return
if (sa2_si < 0d0 .or. sa2_si >= 1d0)return
if (w1 < 0d0 .or. w1 > 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

w2 = 1d0 - w1
s12 = w1 * sa1_si + w2 * sa2_si
s = -sal_g_si(0, 1, 0, s12, t_si, p_si)
if (s == errorreturn) return

s1 = -sal_g_si(0, 1, 0, sa1_si, t_si, p_si)
if (s1 == errorreturn) return

s2 = -sal_g_si(0, 1, 0, sa2_si, t_si, p_si)
if (s2 == errorreturn) return

sal_mixentropy_si = s - w1*s1 - w2*s2

end function

!==========================================================================
function sal_mixvolume_si(sa1_si, sa2_si, w1, t_si, p_si)
!==========================================================================

implicit none

real*8 sal_mixvolume_si
real*8 sa1_si, sa2_si, w1, t_si, p_si
real*8 v, v1, v2, s12, w2

!RETURNS THE MIXING VOLUME OF SEAWATER, V_MIX, IN M3/KG

!SA1_SI:  SALINITY OF COMPONENT 1, IN KG/KG
!SA2_SI:  SALINITY OF COMPONENT 2, IN KG/KG
!W1:      MASS FRACTION OF COMPONENT 1, W1 = M1 / (M1 + M2)
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_MIXVOLUME_SI(0.01,0.035,0.6,300,1E6) = -5.94174956892E-08

sal_mixvolume_si = errorreturn

if (sa1_si < 0d0 .or. sa1_si >= 1d0) return
if (sa2_si < 0d0 .or. sa2_si >= 1d0) return
if (w1 < 0d0 .or. w1 > 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

w2 = 1d0 - w1
s12 = w1 * sa1_si + w2 * sa2_si
v = sal_g_si(0, 0, 1, s12, t_si, p_si)
if (v == errorreturn) return

v1 = sal_g_si(0, 0, 1, sa1_si, t_si, p_si)
if (v1 == errorreturn) return

v2 = sal_g_si(0, 0, 1, sa2_si, t_si, p_si)
if (v2 == errorreturn) return

sal_mixvolume_si = v - w1 * v1 - w2 * v2

end function

!==========================================================================
function sal_osm_coeff_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sal_osm_coeff_si
real*8 sa_si, t_si, p_si
real*8 phi, sqr_s, gi

dimension gi(7)

!RETURNS THE OSMOTIC COEFFICIENT OF SEAWATER, PHI

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_OSM_COEFF_SI(0.035,300,1E6) = 0.902937456585

sal_osm_coeff_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (sa_si == 0d0) then
  sal_osm_coeff_si = 1d0
  return
end if

!COMPUTE THE CASES S > 0:
do i = 1, maxs
  if (i /= 2) then
    gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !get all coefficients of the salinity expansion
    if (gi(i) == errorreturn) return
  end if
enddo

sqr_s = sqrt(sa_si)
phi = 0d0

do i = 3, maxs
  phi = phi + gi(i) * (0.5d0 * dfloat(i) - 1d0) * sqr_s**(i - 2)
enddo

!Begin code check
sal_osm_coeff_si = phi * (1d0 - sa_si) / gi(1) + 1d0 - sa_si
!rs = gas_constant_molar_si/molar_mass_seasalt_si
!sal_osm_coeff_si = (1d0 - sa_si) * (gi(1) + phi) / (rs * t_si)
!End code check

end function

!==========================================================================
function sal_chempot_h2o_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sal_chempot_h2o_si
real*8 sa_si, t_si, p_si
real*8 g, sqr_s, gi

dimension gi(7)

!RETURNS THE SALINE PART OF CHEMICAL POTENTIAL OF WATER IN SEAWATER, mu_WS, IN J/KG

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_CHEMPOT_H2O_SI(0.035,300,1E6) = -2601.18871107

sal_chempot_h2o_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

do i = 1, maxs
  if (i /= 2) then
    gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !get all coefficients of the salinity expansion
    if (gi(i) == errorreturn) return
  end if
enddo

sqr_s = sqrt(sa_si)
g = gi(1)

do i = 3, maxs
  g = g + gi(i) * (0.5d0 * dfloat(i) - 1d0) * sqr_s**(i - 2)
enddo

sal_chempot_h2o_si = -g * sa_si

end function

!==========================================================================
function sal_chempot_rel_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

real*8 sal_chempot_rel_si
real*8 sa_si, t_si, p_si

!RETURNS THE RELATIVE CHEMICAL POTENTIAL OF SEAWATER, mu, IN J/KG

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_CHEMPOT_REL_SI(0.035,300,1E6) = 77949.2100395

sal_chempot_rel_si = errorreturn

if (sa_si <= 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

sal_chempot_rel_si = sal_g_si(1, 0, 0, sa_si, t_si, p_si)

end function

!==========================================================================
function sal_chem_coeff_si(sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE CHEMICAL COEFFICIENT OF SEAWATER, DS, IN J/KG

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_CHEM_COEFF_SI(0.035,300,1E6) = 2754.04566958

real*8 sal_chem_coeff_si, d
real*8 sa_si, t_si, p_si

sal_chem_coeff_si = errorreturn

d = sal_dilution_si(sa_si, t_si, p_si)
if (d == errorreturn) return

sal_chem_coeff_si = sa_si * d

end function

!==========================================================================
function sal_saltenthalpy_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sal_saltenthalpy_si
real*8 sa_si, t_si, p_si
real*8 h, sqr_s, gi, gti

dimension gi(7), gti(7)

!RETURNS THE SPECIFIC ENTHALPY OF SEASALT IN SEAWATER, H_S, IN J/KG

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_SALTENTHALPY_SI(0.035,300,1E6) = -156107.959196

sal_saltenthalpy_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

do i = 2, maxs
  gi(i) = sal_g_term_si(i, 0, 0, t_si, p_si)  !get all coefficients of the salinity expansion
  if (gi(i) == errorreturn) return
  gti(i) = sal_g_term_si(i, 1, 0, t_si, p_si)  !get all t-derivatives of the salinity expansion
  if (gti(i) == errorreturn) return
enddo

sqr_s = sqrt(sa_si)
h = 0d0

do i = 2, maxs
  h = h + (gi(i) - t_si * gti(i)) * sqr_s**(i - 2)
enddo

sal_saltenthalpy_si = h

end function

!==========================================================================
function sal_saltentropy_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sal_saltentropy_si
real*8 sa_si, t_si, p_si
real*8 eta, sqr_s, gi

dimension gi(7)

!RETURNS THE SPECIFIC ENTROPY OF SEASALT IN SEAWATER, ETA_S, IN J/(KG K)

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_SALTENTROPY_SI(0.035,300,1E6) = -532.458305922

sal_saltentropy_si = errorreturn

if (sa_si <= 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

do i = 1, maxs
  gi(i) = sal_g_term_si(i, 1, 0, t_si, p_si)  !get all coefficients of the salinity expansion
  if (gi(i) == errorreturn) return
enddo

sqr_s = sqrt(sa_si)
eta = -gi(1) * log(sa_si)

do i = 2, maxs
  eta = eta - gi(i) * sqr_s**(i - 2)
enddo

sal_saltentropy_si = eta

end function

!==========================================================================
function sal_saltvolume_si(sa_si, t_si, p_si)
!==========================================================================

implicit none

integer i
real*8 sal_saltvolume_si
real*8 sa_si, t_si, p_si
real*8 v, sqr_s, gi

dimension gi(7)

!RETURNS THE SPECIFIC VOLUME OF SEASALT IN SEAWATER, V_S, IN M3/KG

!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE: 
!SAL_SALTVOLUME_SI(0.035,300,1e6) = -7.30285943768e-04

sal_saltvolume_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (t_si <= 0d0) return
if (p_si <= 0d0) return

if (sa_si == 0d0) then
  sal_saltvolume_si = sal_g_term_si(2, 0, 1, t_si, p_si)
  return
end if

do i = 2, maxs
  gi(i) = sal_g_term_si(i, 0, 1, t_si, p_si)  !GET ALL COEFFICIENTS OF THE SALINITY EXPANSION
  if (gi(i) == errorreturn) return
enddo

sqr_s = sqrt(sa_si)
v = gi(2)

do i = 3, maxs
  v = v + gi(i) * sqr_s**(i - 2)
enddo

sal_saltvolume_si = v

end function

!==========================================================================
function sal_molality_si(sa_si)
!==========================================================================

!FUNCTION RETURNS THE MOLALITY OF SEAWATER IN MOL/KG AS A FUNCTION OF ABSOLUTE 
!SALINITY IN KG/KG FOR SEAWATER WITH REFERENCE COMPOSITION, SEE MILLERO ET AL.,
!DEEP-SEA RESEARCH I, 55(2008) 50-72.


!SA_SI:   ABSOLUTE SALINITY IN KG/KG

!CHECK VALUE: 
!SAL_MOLALITY_SI(0.035) = 1.15493681893

real*8 sal_molality_si, sa_si

sal_molality_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return

sal_molality_si = sa_si / ((1d0 - sa_si) * Molar_mass_seasalt_si)

end function

end module sal_2
