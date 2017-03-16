module sea_3a

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS FUNCTION OF SEAWATER DEPENDING ON TEMPERATURE
!AND PRESSURE, AS WELL AS THEIR PARTIAL DERIVATIVES,
!AS DEFINED IN IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION 2008 FOR THE THERMODYNAMIC PROPERTIES OF SEAWATER
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

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     SAL_1,       FILE SAL_1.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_3A,      FILE FLU_3A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     FLU_1,       FILE FLU_1.F90

!#########################################################################

use constants_0
use sal_1
use sal_2
use flu_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: chk_iapws08_table8a, chk_iapws08_table8b, chk_iapws08_table8c, &
          sea_chempot_h2o_si, sea_chempot_rel_si, &
          sea_cp_si, sea_density_si, sea_enthalpy_si, sea_entropy_si, &
          sea_g_si, sea_g_contraction_t_si, sea_g_expansion_t_si, &
          sea_gibbs_energy_si, sea_internal_energy_si, &
          sea_kappa_s_si, sea_kappa_t_si, sea_lapserate_si, &
          sea_osm_coeff_si, sea_soundspeed_si, sea_temp_maxdensity_si

contains

!==========================================================================
function sea_g_si(drv_s, drv_t, drv_p, sa_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE GIBBS FUNCTION OF SEAWATER, IAPWS-08

!OUTPUT:
!SEA_G_SI AS THE S-T-P DERIVATIVE
!(D/DS)^DRV_S (D/DT)^DRV_T (D/DP)^DRV_P G(S,T,P)
!OF THE SPECIFIC GIBBS ENERGY OF SEAWATER, G(S,T,P), IN J/KG

!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_G_SI( 0, 0, 0, 0.035, 300, 1E5) =-5114.9961985|7
!SEA_G_SI( 1, 0, 0, 0.035, 300, 1E5) = 78593.7757371
!SEA_G_SI( 0, 1, 0, 0.035, 300, 1E5) =-374.452240540
!SEA_G_SI( 0, 0, 1, 0.035, 300, 1E5) = 9.77858615182E-04
!SEA_G_SI( 2, 0, 0, 0.035, 300, 1E5) = 2247550.41118
!SEA_G_SI( 1, 1, 0, 0.035, 300, 1E5) = 789.934255688
!SEA_G_SI( 1, 0, 1, 0.035, 300, 1E5) =-7.16682401265E-04
!SEA_G_SI( 0, 2, 0, 0.035, 300, 1E5) =-13.3358324655
!SEA_G_SI( 0, 1, 1, 0.035, 300, 1E5) = 3.04605539768E-07
!SEA_G_SI( 0, 0, 2, 0.035, 300, 1E5) =-4.10945807960E-13

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
integer drv_s, drv_t, drv_p
real*8 sea_g_si, sa_si, t_si, p_si
real*8 gw, gs

sea_g_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_s < 0 .or. drv_s > 3) return
if(drv_t < 0) return
if(drv_p < 0) return

!WATER PART
gw = 0d0
if(drv_s == 0) then
  gw = liq_g_si(drv_t, drv_p, t_si, p_si)
  if(gw == errorreturn) return
end if

!SALINE PART
gs = sal_g_si(drv_s, drv_t, drv_p, sa_si, t_si, p_si)
if(gs == errorreturn) return

sea_g_si = gw + gs

end function

!=========================================================================
function sea_chempot_h2o_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!W(S,T,P) = G - S*(DG/DS)_T_P  CHEM. POTENTIAL OF WATER IN SEAWATER IN J/KG
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_CHEMPOT_H2O_SI(0.035, 300, 1E5) = -7865.7783493|7

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_chempot_h2o_si, sa_si, t_si, p_si
real*8 g, mu_s

sea_chempot_h2o_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = liq_g_si(0, 0, t_si, p_si)
if(g == errorreturn) return

mu_s = sal_chempot_h2o_si(sa_si, t_si, p_si)
if(mu_s == errorreturn) return

sea_chempot_h2o_si = g + mu_s

end function

!=========================================================================
function sea_chempot_rel_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!mu(S,T,P) = (DG/DS)_T_P  RELATIVE CHEM. POTENTIAL OF SEAWATER IN J/KG
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_CHEMPOT_REL_SI(0.035, 300, 1E5) = 78593.7757371

implicit none
real*8 sea_chempot_rel_si, sa_si, t_si, p_si

sea_chempot_rel_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

sea_chempot_rel_si = sal_chempot_rel_si(sa_si, t_si, p_si)

end function

!=========================================================================
function sea_g_contraction_t_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!BETA(S,T,P) = (1/D)*(DD/DS)_T_P isothermal saline contraction of seawater in kg/kg

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_G_CONTRACTION_T_SI(0.035, 300, 1E5) = 0.732910044599

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_g_contraction_t_si, sa_si, t_si, p_si
real*8 g_p, gsp

sea_g_contraction_t_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

gsp = sea_g_si(1, 0, 1, sa_si, t_si, p_si)
if(gsp == errorreturn) return

sea_g_contraction_t_si = -gsp / g_p

end function

!=========================================================================
function sea_cp_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!CP(S,T,P) = T * (DETA/DT)_P  ISOBARIC HEAT CAPACITY OF SEAWATER IN J/(KG K),

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_CP_SI(0.035, 300, 1E5) = 4000.74973964

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_cp_si, sa_si, t_si, p_si
real*8 g_tt

sea_cp_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
if(g_tt == errorreturn) return

sea_cp_si = -t_si * g_tt

end function

!=========================================================================
function sea_density_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!D(S,T,P) = 1/(DG/DP)_S_T  DENSITY OF SEAWATER IN KG/M3

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_DENSITY_SI(0.035, 300, 1E5) = 1022.64272613

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_density_si, sa_si, t_si, p_si
real*8 g_p

sea_density_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

sea_density_si = 1d0 / g_p

end function

!=========================================================================
function sea_enthalpy_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!H(S,T,P) = G - T*(DG/DT)_S_P` ENTHALPY OF SEAWATER IN J/KG

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ENTHALPY_SI(0.035, 300, 1E5) = 107220.675963

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_enthalpy_si, sa_si, t_si, p_si
real*8 g, g_t

sea_enthalpy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = sea_g_si(0, 0, 0, sa_si, t_si, p_si)
if(g == errorreturn) return

g_t = sea_g_si(0, 1, 0, sa_si, t_si, p_si)
if(g_t == errorreturn) return

sea_enthalpy_si = g - t_si * g_t

end function

!==========================================================================
function sea_entropy_si(sa_si, t_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!ETA(S,T,P) = -(DG/DT)_S_P  SPECIFIC ENTROPY OF SEAWATER IN J/(KG K),

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ENTROPY_SI(0.035, 300, 1E5) = 374.452240540

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_entropy_si, sa_si, t_si, p_si
real*8 g_t

sea_entropy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_t = sea_g_si(0, 1, 0, sa_si, t_si, p_si)
if(g_t == errorreturn) return

sea_entropy_si = -g_t

end function

!=========================================================================
function sea_g_expansion_t_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!ALPHA(S,T,P) = -(1/D) * (DD/DT)_S_P  THERMAL EXPANSION OF SEAWATER IN 1/K

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_G_EXPANSION_T_SI(0.035, 300, 1E5) = 3.11502639583E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_g_expansion_t_si, sa_si, t_si, p_si
real*8 g_tp, g_p

sea_g_expansion_t_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
if(g_tp == errorreturn) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

sea_g_expansion_t_si = g_tp / g_p

end function

!=========================================================================
function sea_gibbs_energy_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!G(S,T,P) = G   GIBBS ENERGY OF SEAWATER IN J/KG

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_GIBBS_ENERGY_SI(0.035, 300, 1E5) = -5114.9961985|7

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_gibbs_energy_si, sa_si, t_si, p_si

sea_gibbs_energy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

sea_gibbs_energy_si = sea_g_si(0, 0, 0, sa_si, t_si, p_si)

end function

!=========================================================================
function sea_internal_energy_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!U(S,T,P) = G - T*(DG/DT)_S_P - P*(DG/DP)_S_T INTERNAL ENERGY OF SEAWATER IN J/KG

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_INTERNAL_ENERGY_SI(0.035, 300, 1E5) = 107122.890102

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_internal_energy_si, sa_si, t_si, p_si
real*8 g, g_t, g_p

sea_internal_energy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = sea_g_si(0, 0, 0, sa_si, t_si, p_si)
if(g == errorreturn) return

g_t = sea_g_si(0, 1, 0, sa_si, t_si, p_si)
if(g_t == errorreturn) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

sea_internal_energy_si = g - t_si * g_t - p_si * g_p

end function

!=========================================================================
function sea_kappa_s_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!KAPPA_S(S,T,P) = (1/D) * (DD/DP)_S_ETA  ISENTROPIC COMPRESSIBILITY OF SEAWATER IN 1/PA

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_KAPPA_S_SI(0.035, 300, 1E5) = 4.13135667732E-10

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_kappa_s_si, sa_si, t_si, p_si
real*8 g_p
real*8 g_tp, g_tt, g_pp

sea_kappa_s_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
if(g_tp == errorreturn) return

g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
if(g_tt == errorreturn) return
if(g_tt == 0d0) return

g_pp = sea_g_si(0, 0, 2, sa_si, t_si, p_si)
if(g_pp == errorreturn) return

sea_kappa_s_si = (g_tp * g_tp - g_tt * g_pp) / (g_p * g_tt)

end function

!=========================================================================
function sea_kappa_t_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!KAPPA_T(T,P) = (1/D) * (DD/DP)_S_T  ISOTHERMAL COMPRESSIBILITY OF SEAWATER IN 1/K

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_KAPPA_T_SI(0.035, 300, 1E5) = 4.20250741344E-10

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_kappa_t_si, sa_si, t_si, p_si
real*8 g_p, g_pp

sea_kappa_t_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

g_pp = sea_g_si(0, 0, 2, sa_si, t_si, p_si)
if(g_pp == errorreturn) return

sea_kappa_t_si = -g_pp / g_p

end function

!=========================================================================
function sea_lapserate_si(sa_si, t_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!GAMMA(S,T,P) = (DT/DP)_S_ETA  ADIABATIC LAPSE RATE OF SEAWATER IN K/PA

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_LAPSERATE_SI(0.035, 300, 1E5) = 2.28411342567E-08

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_lapserate_si, sa_si, t_si, p_si
real*8 g_tp, g_tt

sea_lapserate_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
if(g_tp == errorreturn) return

g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
if(g_tt == errorreturn) return
if(g_tt == 0d0) return

sea_lapserate_si = -g_tp / g_tt

end function

!=========================================================================
function sea_osm_coeff_si(sa_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!THE OSMOTIC COEFFICIENT OF SEAWATER, PHI

!INPUT:
!SA_SI:   ABSOLUTE SALINITY IN KG/KG
!T_SI:    ABSOLUTE TEMPERATURE IN K
!P_SI:    ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_OSM_COEFF_SI(0.035, 300, 1E5) = 0.902777495349

implicit none
real*8 sea_osm_coeff_si, sa_si, t_si, p_si

sea_osm_coeff_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

sea_osm_coeff_si = sal_osm_coeff_si(sa_si, t_si, p_si)

end function

!=========================================================================
function sea_soundspeed_si(sa_si, t_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!C(S,T,P) SQRT[(DP/DD)_S_ETA] SOUND SPEED OF SEAWATER IN M/S

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_SOUNDSPEED_SI(0.035, 300, 1E5) = 1538.47940766

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 sea_soundspeed_si, sa_si, t_si, p_si
real*8 g_p, c
real*8 g_tp, g_tt, g_pp

sea_soundspeed_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
if(g_p ==  errorreturn) return
if(g_p <= 0d0) return

g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
if(g_tp ==  errorreturn) return

g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
if(g_tt ==  errorreturn) return
if(g_tt ==  0d0) return

g_pp = sea_g_si(0, 0, 2, sa_si, t_si, p_si)
if(g_pp ==  errorreturn) return

c = g_tp * g_tp - g_tt * g_pp
if(c ==  0d0) return
c = g_tt / c
if(c < 0d0) return

sea_soundspeed_si = g_p * sqrt(c)

end function

!=========================================================================
function sea_temp_maxdensity_si(sa_si, p_si)
!=========================================================================

!RETURNS   T_MD(S,P) TEMPERATURE OF MAXIMUM DENSITY IN K, BY NEWTON ITERATION OF D2G/DTDP = 0

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_TEMP_MAXDENSITY_SI(0.01, 1E5) = 274.950121503

real*8 sea_temp_maxdensity_si, sa_si, p_si
real*8  g_tp, g_ttp, t, dt
real*8 h, t0, eps
integer iter, maxit

h = 0.000001d0      !FINITE DIFFERENCE IN T TO ESTIMATE D3G/DT2DP
t0 = 300d0          !INITIAL TEMPERATURE VALUE
eps = 0.000001d0    !TEMPERATURE TOLERANCE
maxit = 20          !MAX ITERATION STEPS

sea_temp_maxdensity_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(p_si <= 0d0) return

t = t0
do iter = 1, maxit

  g_tp = sea_g_si(0, 1, 1, sa_si, t, p_si)
  if(g_tp == errorreturn) return
  
  g_ttp = sea_g_si(0, 1, 1, sa_si, t + h, p_si)
  if(g_ttp == errorreturn) return
  g_ttp = (g_ttp - g_tp) / h                    !FINITE DIFFERENCE APPROXIMATION FOR 3RD DERIVATIVE
  if(g_ttp == 0d0) return
  
  !IMPROVE ESTIMATE
  dt = -g_tp / g_ttp
  t = t + dt
  if(t <= 0d0) return
  
  if(abs(dt) < eps)then
    sea_temp_maxdensity_si = t
    return
  end if

enddo

end function

!==========================================================================
subroutine chk_iapws08_table8a
!==========================================================================

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*)"Implementation of IAPWS-08 in Fortran"
write(*,*)"for Publication in Ocean Science, 2009"
write(*,*)"D. Wright, BIO, VERSION ", version
write(*,*)" "
write(*,*)"TABLE A5 OF OS2008 AND TABLE 8A OF IAPWS-08"
write(*,*)"PROPERTIES AT S = 0.03516504, T = 273.15 K, P = 101325 Pa"
write(*,*)" "

write(*,*)'QUANTITY  WATER PART           SALINE PART          SEAWATER VALUE      UNITS'
write(*,*)' '

write(*,*)'g         1.01342742e2        -1.0134274172939e2    4.e-9               J/KG'
write(*,*)'dg/ds     0.0                  6.39974067312299e4   6.39974067312299e4  J/KG'
write(*,*)'dg/dt     1.47643376e-1       -1.47643376346e-1    -6.e-11              J/KG/K'
write(*,*)'dg/dp     1.00015693912169e-3 -2.74957224268433e-5  9.7266121669485e-4  M3/KG'
write(*,*)'d2g/dsdp  0.0                 -7.59615411515309e-4 -7.59615411515309e-4 M3/KG'
write(*,*)'d2g/dt2  -1.54473542320e1      8.5286115117592e-1  -1.459449308080e1    J/KG/K2'
write(*,*)'d2g/dtdp -6.777003179e-8       1.1928678741396e-7   5.1516755627e-8     M3/KG/K'
write(*,*)'d2g/dp2  -5.08928894643e-13    5.8153517233288e-14 -4.507753774102e-13  M3/KG/Pa'
write(*,*)'h         6.10139535e1        -6.10139534804e1      2.e-8               J/KG'
write(*,*)'f         1.8399e-3           -9.8556737654491e1   -9.855489778e1       J/KG'
write(*,*)'u        -4.03269484e1        -5.82279494055e1     -9.85548978e1        J/KG'
write(*,*)'sigma    -1.47643376e-1        1.47643376346e-1     6.e-11              J/KG/K'
write(*,*)'rho       9.9984308550433e2   -9d99                 1.02810719995401e3  KG/M3'
write(*,*)'cp        4.21944480846e3     -2.3295902344370e2    3.98648578502e3     J/KG/K'
write(*,*)'c         1.40238253109e3     -9d99                 1.44900246362e3     M/S'
write(*,*)'muw       1.01342742e2        -2.35181410932936e3  -2.250471368e3       J/KG'
write(*,*)' '

s = 0.03516504d0
t = 273.15d0
p = 101325d0

zg(1) = 1.01342742d2  	
zdgds(1) = 0d0
zdgdt(1) = 1.47643376d-1  
zdgdp(1) = 1.00015693912169d-3
zd2gdsdp(1) = 0d0
zd2gdt2(1) = -1.54473542320d1
zd2gdtdp(1) = -6.777003179d-8
zd2gdp2(1) = -5.08928894643d-13	
zh(1) = 6.10139535d1
zf(1) = 1.8399d-3
zu(1) = -4.03269484d1	
zsigma(1) = -1.47643376d-1
zrho(1) = 9.9984308550433d2
zcp(1) = 4.21944480846d3
zc(1) = 1.40238253109d3	
zmuw(1) = 1.01342742d2

zg(2) = -1.0134274172939d2      
zdgds(2) = 6.39974067312299d4 
zdgdt(2) = -1.47643376346d-1
zdgdp(2) = -2.74957224268433d-5	
zd2gdsdp(2) = -7.59615411515309d-4
zd2gdt2(2) = 8.5286115117592d-1
zd2gdtdp(2) = 1.1928678741396d-7	
zd2gdp2(2) = 5.8153517233288d-14
zh(2) = -6.10139534804d1
zf(2) = -9.8556737654491d1	
zu(2) = -5.82279494055d1
zsigma(2) = 1.47643376346d-1
zrho(2) = -9.d99	       
zcp(2) = -2.3295902344370d2  
zc(2) = -9.d99
zmuw(2) = -2.35181410932936d3

zg(3) = 4.d-9
zdgds(3) = 6.39974067312299d4
zdgdt(3) = -6.d-11          
zdgdp(3) = 9.7266121669485d-4
zd2gdsdp(3) = -7.59615411515309d-4 
zd2gdt2(3) = -1.459449308080d1   
zd2gdtdp(3) = 5.1516755627d-8    
zd2gdp2(3) = -4.507753774102d-13
zh(3) = 2.d-8             
zf(3) = -9.855489778d1   
zu(3) = -9.85548978d1   
zsigma(3) = 6.d-11         
zrho(3) = 1.02810719995401d3
zcp(3) = 3.98648578502d3  
zc(3) = 1.44900246362d3 
zmuw(3) = -2.250471368d3 
             
icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A <1 OR >3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol ==0 .or. icol >= 4) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif

write(*,*) " "
write(*,*) " "
write(*,*) " COMPARISON WITH FUNCTION VALUES AS GIVEN IN TABLE A5, COLUMN", icol, &
           " OF OS2008"
write(*,*) " PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
 &       UNITS"
write(*,*) " "

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q

go to 100

end subroutine

!==========================================================================
subroutine chk_iapws08_table8b
!==========================================================================

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*)"Implementation of IAPWS-08 in Fortran"
write(*,*)"for Publication in Ocean Science, 2009"
write(*,*)"D. Wright, BIO, VERSION ", version
write(*,*)" "
write(*,*)"TABLE A6 OF OS2008 AND TABLE 8B OF IAPWS-08"
write(*,*)"PROPERTIES AT S = 0.1, T = 353 K, P = 101325 Pa"
write(*,*)" "

write(*,*)'QUANTITY  WATER PART          SALINE PART          SEAWATER VALUE      UNITS'
write(*,*)' '

write(*,*)'g        -4.46114968996e4      1.5087174003705e4   -2.95243228959e4    J/KG'
write(*,*)'dg/ds     0.0                  2.51957275851413e5   2.51957275851413e5 J/KG'
write(*,*)'dg/dt    -1.0737599318875e3    1.56230907404291e2  -9.175290244832e2   J/KG/K'
write(*,*)'dg/dp     1.02892955635611e-3 -5.7922728577126e-5   9.7100682777898e-4 M3/KG'
write(*,*)'d2g/dsdp  0.0                 -3.0595780244234e-4  -3.0595780244234e-4 M3/KG'
write(*,*)'d2g/dt2  -1.1888500004755e1    1.27922649315507e0  -1.0609273511600e1  J/KG/K2'
write(*,*)'d2g/dtdp  6.59051552339e-7     8.0306159575153e-7   1.462113148091e-6  M3/KG/K'
write(*,*)'d2g/dp2  -4.746728193611e-13   2.1308615424374e-13 -2.61586665117e-13  M3/KG/Pa'
write(*,*)'h         3.344257590567e5    -4.006233631001e4     2.94363422747e5    J/KG'
write(*,*)'f        -4.47157531869e4      1.5093043024178e4   -2.96227101627e4    J/KG'
write(*,*)'u         3.343215027694e5    -4.0056467289536e4    2.942650354799e5   J/KG'
write(*,*)'sigma     1.0737599318875e3   -1.56230907404291e2   9.175290244832e2   J/KG/K'
write(*,*)'rho       9.7188383191308e2   -9d99                 1.02985887574790e3 KG/M3'
write(*,*)'cp        4.1966405016784e3   -4.51566952083741e2   3.745073549595e3   J/KG/K'
write(*,*)'c         1.5544629665347e3   -9d99                 3.9612783529e3     M/S'
write(*,*)'muw      -4.46114968996e4     -1.01085535814360e4  -5.47200504810e4    J/KG'
write(*,*)' '

s = 0.1d0
t = 353d0
p = 101325d0

zg(1) = -4.46114968996d4
zdgds(1) = 0d0
zdgdt(1) = -1.0737599318875d3
zdgdp(1) = 1.02892955635611d-3
zd2gdsdp(1) = 0d0
zd2gdt2(1) = -1.1888500004755d1
zd2gdtdp(1) = 6.59051552339d-7
zd2gdp2(1) = -4.746728193611d-13
zh(1) = 3.344257590567d5
zf(1) = -4.47157531869d4
zu(1) = 3.343215027694d5
zsigma(1) = 1.0737599318875d3
zrho(1) = 9.7188383191308d2
zcp(1) = 4.1966405016784d3
zc(1) = 1.5544629665347d3
zmuw(1) = -4.46114968996d4

zg(2) = 1.5087174003705d4
zdgds(2) = 2.51957275851413d5
zdgdt(2) = 1.56230907404291d2
zdgdp(2) = -5.7922728577126d-5
zd2gdsdp(2) = -3.0595780244234d-4
zd2gdt2(2) = 1.27922649315507d0
zd2gdtdp(2) = 8.0306159575153d-7
zd2gdp2(2) = 2.1308615424374d-13
zh(2) = -4.006233631001d4
zf(2) = 1.5093043024178d4
zu(2) = -4.0056467289536d4
zsigma(2) = -1.56230907404291d2
zrho(2) = -9d99	           
zcp(2) = -4.51566952083741d2
zc(2) = -9d99	           
zmuw(2) = -1.01085535814360d4

zg(3) = -2.95243228959d4 
zdgds(3) = 2.51957275851413d5
zdgdt(3) = -9.175290244832d2
zdgdp(3) = 9.7100682777898d-4
zd2gdsdp(3) = -3.0595780244234d-4
zd2gdt2(3) = -1.0609273511600d1
zd2gdtdp(3) = 1.462113148091d-6
zd2gdp2(3) = -2.61586665117d-13
zh(3) = 2.94363422747d5
zf(3) = -2.96227101627d4
zu(3) = 2.942650354799d5
zsigma(3) = 9.175290244832d2
zrho(3) = 1.02985887574790d3
zcp(3) = 3.745073549595d3
zc(3) = 3.9612783529d3
zmuw(3) = -5.47200504810d4
             
icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A VALUE < 1 or > 3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol < 1 .or. icol > 3) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif

write(*,*) " "
write(*,*) " "
write(*,*) " FUNCTION VALUES AS GIVEN IN TABLE A6, COLUMN", icol, &
            " OF OS2008"
write(*,*) " PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
&       UNITS"
write(*,*) " "

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q
 
go to 100

end subroutine

!==========================================================================
subroutine chk_iapws08_table8c
!==========================================================================

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) " "
write(*,*)"Implementation of IAPWS-08 in Fortran"
write(*,*)"for Publication in Ocean Science, 2009"
write(*,*)"D. Wright, BIO, VERSION ", version
write(*,*)" "
write(*,*)"TABLE A7 OF OS2008 AND TABLE 8C OF IAPWS-08"
write(*,*)"PROPERTIES S = 0.1, T = 353 K, P = 10**8 Pa"
write(*,*) " "

write(*,*)'QUANTITY  WATER PART           SALINE PART          SEAWATER VALUE      UNITS'
write(*,*)' '

write(*,*)'g         9.773038621954e+4  -2.60093050730637e+3   9.512945571223e+4   J/KG'
write(*,*)'dg/ds     0.0                -5.4586158064880e+3   -5.4586158064880e+3  J/KG'
write(*,*)'dg/dt     8.5146650206e       7.5404568488117e+0    1.60551218694e+1    J/KG/K'
write(*,*)'dg/dp     9.5668332915351e-4 -2.2912384179113e-5    9.3377094497440e-4  M3/KG'
write(*,*)'d2g/dsdp  0.0                -6.4075761854575e-4   -6.40757618545748e-4 M3/KG'
write(*,*)'d2g/dt2  -1.429698733876e+1   4.88076973942251e-1  -1.380891036482e+1   J/KG/K2'
write(*,*)'d2g/dtdp  1.99079570803e-7    4.6628441224121e-8    2.45708012027e-7    M3/KG/K'
write(*,*)'d2g/dp2  -3.715308894234e-13  3.57345735845327e-14 -3.35796315839e-13   M3/KG/Pa'
write(*,*)'h         9.54046054692e+4   -4.6606062955593e+3    9.07439991736e+4    J/KG'
write(*,*)'f         2.06205330419e+3   -3.0969208939506e+2    1.75236121479e+3    J/KG'
write(*,*)'u        -2.637274462e+2     -2.3693678776480e+3   -2.6330953238e+3     J/KG' 
write(*,*)'sigma    -8.5146650206e+0    -7.5404568488117      -1.6055121869e+1     J/KG/K'
write(*,*)'rho       1.04527796139692e3 -9d99                  1.0709264465574e+3  KG/M3'
write(*,*)'cp        3.905222091582e+3  -1.33318225432326e+2   3.77190386615e+3    J/KG/K'
write(*,*)'c         1.575422398486e+3  -9d99                  1.621989976499e+3   M/S'
write(*,*)'muw       9.77303862195e+4   -2.4089780641266e+3    9.532140815541e+4   J/KG'
write(*,*)' '

s = 0.03516504d0
t = 273.15d0
p = 1d8

zg(1) = 9.773038621954d+4
zdgds(1) = 0.0	     
zdgdt(1) = 8.5146650206d+0
zdgdp(1) = 9.5668332915351d-4
zd2gdsdp(1) = 0.0	      
zd2gdt2(1) = -1.429698733876d+1
zd2gdtdp(1) = 1.99079570803d-7
zd2gdp2(1) = -3.715308894234d-13
zh(1) = 9.54046054692d+4
zf(1) = 2.06205330419d+3
zu(1) = -2.637274462d+2
zsigma(1) = -8.5146650206d+0
zrho(1) = 1.04527796139692d+3
zcp(1) = 3.905222091582d+3
zc(1) = 1.575422398486d+3
zmuw(1) = 9.77303862195d+4

zg(2) = -2.60093050730637d+3
zdgds(2) = -5.4586158064880d+3
zdgdt(2) = 7.5404568488117d+0
zdgdp(2) = -2.2912384179113d-5
zd2gdsdp(2) = -6.4075761854575d-4
zd2gdt2(2) = 4.88076973942251d-1
zd2gdtdp(2) = 4.6628441224121d-8
zd2gdp2(2) = 3.57345735845327d-14
zh(2) = -4.6606062955593d+3
zf(2) = -3.0969208939506d+2
zu(2) = -2.3693678776480d+3
zsigma(2) = -7.5404568488117d+0
zrho(2) = -9d99
zcp(2) = -1.33318225432326d+2
zc(2) = -9d99
zmuw(2) = -2.4089780641266d+3

zg(3) = 9.512945571223d+4
zdgds(3) = -5.4586158064880d+3
zdgdt(3) = 1.60551218694d+1
zdgdp(3) = 9.3377094497440d-4
zd2gdsdp(3) = -6.40757618545748d-4
zd2gdt2(3) = -1.380891036482d+1
zd2gdtdp(3) = 2.45708012027d-7
zd2gdp2(3) = -3.35796315839d-13
zh(3) = 9.07439991736d+4
zf(3) = 1.75236121479d+3
zu(3) = -2.6330953238d+3
zsigma(3) = -1.6055121869d+1
zrho(3) = 1.0709264465574d+3
zcp(3) = 3.77190386615d+3
zc(3) = 1.621989976499d+3
zmuw(3) = 9.532140815541d+4

icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A VALUE < 1 or > 3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol < 1 .or. icol > 3) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif


write(*,*) " "
write(*,*) " "
write(*,*) " FUNCTION VALUES AS GIVEN IN TABLE A7, COLUMN", icol, &
           " OF OS2008"
write(*,*) " PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
&       UNITS"
write(*,*) " "

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q
 
go to 100

end subroutine


end module sea_3a
