module sea_vap_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIUM PROPERTIES OF VAPOUR WITH SEAWATER
!COMPUTED FROM IAPWS-95 AND IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

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
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     MATHS_0,      FILE MATHS_0.F90
!     FLU_1,        FILE FLU_1.F90
!     SAL_1,        FILE SAL_1.F90
!     SAL_2,        FILE SAL_2.F90
!     FLU_3A,       FILE FLU_3A.F90
!     FLU_3B,       FILE FLU_3B.F90
!     SEA_3A,       FILE SEA_3A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,    FILE CONVERT_0.F90
!     FLU_2,        FILE FLU_2.F90

!#########################################################################

use constants_0
use maths_0
use flu_1
use sal_1
use sal_2
use flu_3a
use flu_3b
use sea_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_vap_boilingtemperature_si, sea_vap_brinefraction_seavap_si, &
          sea_vap_brinesalinity_si, sea_vap_cp_seavap_si, &
          sea_vap_density_sea_si, sea_vap_density_seavap_si, &
          sea_vap_density_vap_si, sea_vap_enthalpy_evap_si, &
          sea_vap_enthalpy_sea_si, sea_vap_enthalpy_seavap_si, &
          sea_vap_enthalpy_vap_si, sea_vap_entropy_sea_si, &
          sea_vap_entropy_seavap_si, sea_vap_entropy_vap_si, &
          sea_vap_expansion_seavap_si, sea_vap_g_si, &
          sea_vap_kappa_t_seavap_si, sea_vap_pressure_si, &
          sea_vap_salinity_si, sea_vap_temperature_si, &
          sea_vap_vapourpressure_si, sea_vap_volume_evap_si, &
          set_it_ctrl_sea_vap, set_sea_vap_eq_at_t_p, &
          set_sea_vap_eq_at_s_p, set_sea_vap_eq_at_s_t

integer, private :: ctrl_mode_sea, ctrl_mode_temperature, ctrl_mode_liquid, &
         ctrl_mode_vapour, ctrl_loop_maximum, ctrl_initialized
real*8,  private :: ctrl_init_s_sea, ctrl_init_t, ctrl_init_d_liq, ctrl_init_d_vap, &
         ctrl_eps_exit_s, ctrl_eps_exit_t, ctrl_eps_exit_p
real*8,  private :: equi_sea_vap_t, equi_sea_vap_s, equi_sea_vap_p
real*8,  private :: equi_sea_vap_d_vap, equi_sea_vap_d_liq
integer, private :: equi_sea_vap_done

!IF97 COEFFICIENTS FROM:
!IF97 coefficients from:
!Revised Release on the IAPWS Industrial Formulation 1997
!for the Thermodynamic Properties of Water and Steam
!The International Association for the Properties of Water and Steam
!Lucerne, Switzerland, August 2007

!COEFFICIENTS OF REGION 1
integer, private, save :: i1i(34), j1i(34)
real*8,  private, save :: n1i(34)

!COEFFICIENTS OF REGION 2, IDEAL PART
integer, private, save :: j0i(9)
real*8,  private, save :: n0i(9)

!COEFFICIENTS OF REGION 2, RESIDUAL PART
integer, private, save :: iri(43), jri(43)
real*8,  private, save :: nri(43)

contains

!==========================================================================
function sea_vap_boilingtemperature_si(sa_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE TEMPERATURE IN K AT WHICH SEAWATER WITH
!ABSOLUTE SALINITY SA_SI IN KG/KG BOILS AT THE ABSOLUTE PRESSURE P_SI IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_BOILINGTEMPERATURE_SI(0.035, 640) = 274.042416829

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_boilingtemperature_si, sa_si, p_si

sea_vap_boilingtemperature_si = errorreturn

if(set_sea_vap_eq_at_s_p(sa_si, p_si) == errorreturn) return

sea_vap_boilingtemperature_si = equi_sea_vap_t

end function

!==========================================================================
function sea_vap_brinesalinity_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE SALINITY SA_SI IN KG/KG WITH WHICH SEAWATER BOILS AT THE
!ABSOLUTE TEMPERATURE IN K AND THE ABSOLUTE PRESSURE P_SI IN PA

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE
!AND OF THOSE IN FLU_3A

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_BRINESALINITY_SI(274, 640) = 2.9439629829|4E-02


implicit none
real*8 sea_vap_brinesalinity_si, t_si, p_si

sea_vap_brinesalinity_si = errorreturn

if(set_sea_vap_eq_at_t_p(t_si, p_si) == errorreturn) return

sea_vap_brinesalinity_si = equi_sea_vap_s

end function

!=========================================================================
function sea_vap_cp_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!CP(SSV,T,P) = -T * (D2G/DT2)_SSV_P ISOBARIC HEAT CAPACITY OF SEA VAPOUR IN J/(KG K)

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_CP_SEAVAP_SI(0.035, 274, 600) = 756270.431593

!NOTE: THE ACCURACY OF THIS FUNCTion DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_cp_seavap_si, ssv_si, t_si, p_si
real*8 g_tt

sea_vap_cp_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tt = sea_vap_g_si(0, 2, 0, ssv_si, t_si, p_si)
if(g_tt == errorreturn) return

sea_vap_cp_seavap_si = -t_si * g_tt

end function

!==========================================================================
function sea_vap_density_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY IN KG/M3 OF SEAWATER AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE,
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN fLU_3A

!check values with default settings:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_DENSITY_SEA_SI = 1027.87349556

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_DENSITY_SEA_SI = 1027.87626132

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_DENSITY_SEA_SI = 1023.42713047

implicit none
real*8 sea_vap_density_sea_si

sea_vap_density_sea_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_density_sea_si = sea_density_si(equi_sea_vap_s, equi_sea_vap_t, equi_sea_vap_p)

end function

!=========================================================================
function sea_vap_density_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!D(SSV,T,P) = 1/(DG/DP)_SSV_T  DENSITY OF SEA VAPOUR IN KG/M3

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_DENSITY_SEAVAP_SI(0.035, 274, 600) = 7.27092786882E-03

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_density_seavap_si, ssv_si, t_si, p_si
real*8 g_p

sea_vap_density_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_vap_g_si(0, 0, 1, ssv_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

sea_vap_density_seavap_si = 1d0 / g_p

end function

!==========================================================================
function sea_vap_density_vap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY IN KG/M3 OF VAPOUR AT EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_DENSITY_VAP_SI = 5.06324890264E-03

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_DENSITY_VAP_SI = 5.04855547811E-03

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_DENSITY_VAP_SI = 5.06403699513E-03
 
implicit none
real*8 sea_vap_density_vap_si

sea_vap_density_vap_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_density_vap_si = vap_density_si(equi_sea_vap_t, equi_sea_vap_p)

end function

!==========================================================================
function sea_vap_enthalpy_evap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE EVAPORATION ENTHALPY IN J/KG OF SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_ENTHALPY_EVAP_SI = 2498295.32187

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_ENTHALPY_EVAP_SI = 2498395.40101

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_ENTHALPY_EVAP_SI = 2498551.19875

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_enthalpy_evap_si
real*8 s, t, p, d

sea_vap_enthalpy_evap_si = errorreturn

if(equi_sea_vap_done /= isok) return

s = equi_sea_vap_s
t = equi_sea_vap_t
p = equi_sea_vap_p

d = sal_dilution_si(s, t, p)
if(d == errorreturn) return

sea_vap_enthalpy_evap_si = t * d * sea_vap_s_si(1, 0, t, p)

end function

!==========================================================================
function sea_vap_enthalpy_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTHALPY IN KG/M3 OF SEAWATER AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_ENTHALPY_SEA_SI = 3465.11896144

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_ENTHALPY_SEA_SI = 3295.96629299

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_ENTHALPY_SEA_SI = 3405.93353730

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQ_AT_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_enthalpy_sea_si

sea_vap_enthalpy_sea_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_enthalpy_sea_si = sea_enthalpy_si(equi_sea_vap_s, equi_sea_vap_t, equi_sea_vap_p)

end function

!=========================================================================
function sea_vap_enthalpy_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!H(SSV,T,P) = GSV - T * (DG/DT)_SSV_P  ENTHALPY OF SEA VAPOUR IN J/KG

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_ENTHALPY_SEAVAP_SI(0.035, 274, 600) = 1661118.41089

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_enthalpy_seavap_si, ssv_si, t_si, p_si
real*8 g, g_t

sea_vap_enthalpy_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = sea_vap_g_si(0, 0, 0, ssv_si, t_si, p_si)
if(g == errorreturn) return
g_t = sea_vap_g_si(0, 1, 0, ssv_si, t_si, p_si)
if(g_t == errorreturn) return

sea_vap_enthalpy_seavap_si = g - t_si * g_t

end function

!==========================================================================
function sea_vap_enthalpy_vap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTHALPY IN KG/M3 OF VAPOUR AT EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_ENTHALPY_VAP_SI = 2502546.89358

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_ENTHALPY_VAP_SI = 2502469.07187

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_ENTHALPY_VAP_SI = 2502466.96633

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real *8 sea_vap_enthalpy_vap_si

sea_vap_enthalpy_vap_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_enthalpy_vap_si = vap_enthalpy_si(equi_sea_vap_t, equi_sea_vap_p)

end function

!==========================================================================
function sea_vap_entropy_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTROPY IN J/(KG K) OF SEAWATER AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_ENTROPY_SEA_SI = 13.0616891215

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_ENTROPY_SEA_SI = 12.4443983378

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_ENTROPY_SEA_SI = 14.0256815112

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE,
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_entropy_sea_si

sea_vap_entropy_sea_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_entropy_sea_si = sea_entropy_si(equi_sea_vap_s, equi_sea_vap_t, equi_sea_vap_p)

end function

!=========================================================================
function sea_vap_entropy_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!ETA(SSV,T,P) = - (DG/DT)_SSV_P  ENTROPY OF SEA VAPOUR IN J/(KG K)

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_ENTROPY_SEAVAP_SI(0.035, 274, 600) = 6072.50817709

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_entropy_seavap_si, ssv_si, t_si, p_si
real*8 g_t

sea_vap_entropy_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_t = sea_vap_g_si(0, 1, 0, ssv_si, t_si, p_si)
if(g_t == errorreturn) return

sea_vap_entropy_seavap_si = -g_t

end function

!==========================================================================
function sea_vap_entropy_vap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTROPY IN J/(KG K) OF VAPOUR AT EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_ENTROPY_VAP_SI = 9140.56256065

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_ENTROPY_VAP_SI = 9141.68990452

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_ENTROPY_VAP_SI = 9140.27087793

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_entropy_vap_si

sea_vap_entropy_vap_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_entropy_vap_si = vap_entropy_si(equi_sea_vap_t, equi_sea_vap_p)

end function

!=========================================================================
function sea_vap_expansion_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!ALPHA(SSV,T,P) = (D2G/DTDP)_SSV/(DG/DP)_SSV_T  THERMAL EXPANSION OF SEA VAPOUR IN 1/K,

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K,
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_EXPANSION_SEAVAP_SI(0.035, 274, 600) = 0.458863421347

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_expansion_seavap_si, ssv_si, t_si, p_si
real*8 g_p, g_tp

sea_vap_expansion_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_vap_g_si(0, 0, 1, ssv_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_tp = sea_vap_g_si(0, 1, 1, ssv_si, t_si, p_si)
if(g_tp == errorreturn) return

sea_vap_expansion_seavap_si = g_tp / g_p

end function

!==========================================================================
function sea_vap_g_si(drv_s, drv_t, drv_p, ssv_si, t_si, p_si)
!==========================================================================

!RETURNS THE SSV-T-P DERIVATIVE (D/DSSV)^DRV_S (D/DT)^DRV_T (D/DP)^DRV_P G(SSV,T,P)
!OF THE GIBBS FUNCTION G OF SEA VAPOUR, AS A FUNCTION OF SEA-VAPOUR SALINITY IN KG/KG,
!ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA
!
!G(S,T,P) = (1-W) * GV(T,P) + W * GSW(S(T,P),T,P)
!
!W = MASS FRACTION OF BRINE, W(SSV, T, P) = SSV/S(T,P)
!S = BRINE SALINITY
!
!SSV_SI = ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG, I.E. THE MASS FRACTION OF SALT IN SEAWATER + VAPOUR
!T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
!P_SI  = P = ABSOLUTE PRESSURE IN PA
!SAL_G_SI  = SPECIFIC GIBBS ENERGY IN J/KG


!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_VAP_G_SI( 0, 0, 0, 0.035, 274, 610) =-2748.82963245
!SEA_VAP_G_SI( 1, 0, 0, 0.035, 274, 610) = 151028.257424
!SEA_VAP_G_SI( 0, 1, 0, 0.035, 274, 610) =-6072.50817709
!SEA_VAP_G_SI( 0, 0, 1, 0.035, 274, 610) = 137.534028399
!SEA_VAP_G_SI( 2, 0, 0, 0.035, 274, 610) = 0
!SEA_VAP_G_SI( 1, 1, 0, 0.035, 274, 610) = 14965.0677011
!SEA_VAP_G_SI( 1, 0, 1, 0.035, 274, 610) =-321.591932572
!SEA_VAP_G_SI( 0, 2, 0, 0.035, 274, 610) =-2760.11106421
!SEA_VAP_G_SI( 0, 1, 1, 0.035, 274, 610) = 63.1093348229
!SEA_VAP_G_SI( 0, 0, 2, 0.035, 274, 610) =-1.65027885871

implicit none
integer drv_s, drv_t, drv_p
real*8 sea_vap_g_si, ssv_si, t_si, p_si
real*8 g, s, w, d
real*8 gb, gv
real*8 gb_p, gv_p
real*8 gb_t, gv_t
real*8 gb_pp, gv_pp, gm_pp
real*8 gb_tp, gv_tp, gm_tp
real*8 gb_tt, gv_tt, gm_tt
real*8 gb_s, gb_sp, gb_st
real*8 s_p, s_t

sea_vap_g_si = errorreturn
d = errorreturn
g = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_s < 0) return
if(drv_t < 0 .or. drv_t > 2) return
if(drv_p < 0 .or. drv_p > 2) return

if(drv_s > 1) then    !G IS LINEAR IN SSV
  sea_vap_g_si = 0d0
  return
end if

if(set_sea_vap_eq_at_t_p(t_si, p_si) == errorreturn) return

s = sea_vap_salinity_si()
if(s == errorreturn) return
if(s <= 0d0 .or. s >= 1d0) return
if(s < ssv_si) return

w = ssv_si / s               !LIQUID MASS FRACTION OF SEA VAPOUR

if(drv_t + drv_p > 1d0) then
  d = sal_dilution_si(s, t_si, p_si)
end if

select case (drv_s)

  case (0)

    select case (drv_t)

      case (0)     !D/DT

        select case (drv_p)

          case (0)
            gb = sea_g_si(0, 0, 0, s, t_si, p_si)
            if(gb == errorreturn) return
            gv = vap_g_si(0, 0, t_si, p_si)
            if(gv == errorreturn) return
            g = (1d0 - w) * gv + w * gb                      !G

          case (1)
            gb_p = sea_g_si(0, 0, 1, s, t_si, p_si)
            if(gb_p == errorreturn) return
            gv_p = vap_g_si(0, 1, t_si, p_si)
            if(gv_p == errorreturn) return
            g = (1d0 - w) * gv_p + w * gb_p                  !V = G_P

          case (2)
            s_p = sea_vap_s_si(0, 1, t_si, p_si)
            if(s_p == errorreturn) return
            gm_pp = -d * s_p ** 2 / s                      !LATENT DERIVATIVE
            gb_pp = sea_g_si(0, 0, 2, s, t_si, p_si)
            if(gb_pp == errorreturn) return
            gv_pp = vap_g_si(0, 2, t_si, p_si)
            if(gv_pp == errorreturn) return
            g = (1d0 - w) * gv_pp + w * (gb_pp + gm_pp)      !G_PP

          case default; return
        end select

      case (1)     !D/DT

        select case (drv_p)

          case (0)
            gb_t = sea_g_si(0, 1, 0, s, t_si, p_si)
            if(gb_t == errorreturn) return
            gv_t = vap_g_si(1, 0, t_si, p_si)
            if(gv_t == errorreturn) return
            g = (1d0 - w) * gv_t + w * gb_t                  !-ETA = G_T

          case (1)
            s_t = sea_vap_s_si(1, 0, t_si, p_si)
            if(s_t == errorreturn) return
            s_p = sea_vap_s_si(0, 1, t_si, p_si)
            if(s_p == errorreturn) return
            gm_tp = -d * s_t * s_p / s                     !LATENT DERIVATIVE
            gb_tp = sea_g_si(0, 1, 1, s, t_si, p_si)
            if(gb_tp == errorreturn) return
            gv_tp = vap_g_si(1, 1, t_si, p_si)
            if(gv_tp == errorreturn) return
            g = (1d0 - w) * gv_tp + w * (gb_tp + gm_tp)                 !G_TP

          case default; return
        end select

      case (2)     !D2/DT2

        select case (drv_p)

          case (0)
            s_t = sea_vap_s_si(1, 0, t_si, p_si)
            if(s_t == errorreturn) return
            gm_tt = -d * s_t ** 2 / s                      !LATENT DERIVATIVE
            gb_tt = sea_g_si(0, 2, 0, s, t_si, p_si)
            if(gb_tt == errorreturn) return
            gv_tt = vap_g_si(2, 0, t_si, p_si)
            if(gv_tt == errorreturn) return
            g = (1d0 - w) * gv_tt + w * (gb_tt + gm_tt)                 !G_TT

          case default; return
        end select

      case default; return
    end select

  case (1)       !D/DS
  
    select case (drv_t)
    
      case (0)
        select case (drv_p)

          case (0)
             g = sal_g_si(1, 0, 0, s, t_si, p_si)                     !G_S

          case (1)
             gb_sp = sal_g_si(1, 0, 1, s, t_si, p_si)                 !G_SP
             if(gb_sp == errorreturn) return
             gb_s = sal_g_si(1, 0, 0, s, t_si, p_si)
             if(gb_s == errorreturn) return
             s_p = sea_vap_s_si(0, 1, t_si, p_si)
             if(s_p == errorreturn) return
             g = gb_sp + gb_s * s_p

          case default; return
        end select

      case (1)   !D2/DSDT

        select case (drv_p)

          case (0)
             gb_st = sal_g_si(1, 1, 0, s, t_si, p_si)                 !G_ST
             if(gb_st == errorreturn) return
             gb_s = sal_g_si(1, 0, 0, s, t_si, p_si)
             if(gb_s == errorreturn) return
             s_t = sea_vap_s_si(1, 0, t_si, p_si)
             if(s_t == errorreturn) return
             g = gb_st + gb_s * s_t

          case default; return
        end select

       case default; return
    end select
end select

sea_vap_g_si = g

end function

!=========================================================================
function sea_vap_kappa_t_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!KAPPA_T(SSV,T,P) = - (D2G/DP2)_SSV_T/(DG/DP)_SSV_T  ISOTHERMAL COMPRESSIBILITY OF SEA VAPOUR IN 1/PA

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_KAPPA_T_SEAVAP_SI(0.035, 274, 600) = 1.19990585451E-02

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_kappa_t_seavap_si, ssv_si, t_si, p_si
real*8 g_p, g_pp

sea_vap_kappa_t_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_vap_g_si(0, 0, 1, ssv_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_pp = sea_vap_g_si(0, 0, 2, ssv_si, t_si, p_si)
if(g_pp == errorreturn) return

sea_vap_kappa_t_seavap_si = -g_pp / g_p

end function

!==========================================================================
function sea_vap_pressure_si()
!==========================================================================

!THIS FUNCTION RETURNS THE PRESSURE IN PA OF VAPOUR AT EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_PRESSURE_SI = 640

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_PRESSURE_SI = 638.044692615

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_PRESSURE_SI = 640

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_pressure_si

sea_vap_pressure_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_pressure_si = equi_sea_vap_p

end function

!=========================================================================
function sea_vap_s_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS ABSOLUTE BRINE SALINITY S IN KG/KG, AND ITS PARTIAL DERIVATIVES
!  (D/DT)^DRV_T (D/DP)^DRV_P S(T,P)
!AS A FUNCTION OF ABSOLUTE TEMPERATURE, T_SI, IN K, AND ABSOLUTE PRESSURE, P_SI, IN PA

implicit none
integer drv_t, drv_p
real*8 sea_vap_s_si, t_si, p_si
real*8 d  !DILUTION COEFF
real*8 s  !BRINE SALINITY

real*8 gb_p, gv_p, gb_sp
real*8 gb_t, gv_t, gb_st

sea_vap_s_si = errorreturn

d = errorreturn

if(drv_t < 0 .or. drv_t > 1) return
if(drv_p < 0 .or. drv_p > 1) return

if(t_si < 0d0) return
if(p_si < 0d0) return

!THIS MUST BE CALLED OUTSIDE BEFORE:
!IF(SET_SEA_ICE_EQUILIBRIUM_TO_T_P(T_SI, P_SI) == ERRORRETURN) RETURN

!IF NOT, EXIT HERE
if(equi_sea_vap_done /= isok) return

s = sea_vap_salinity_si()
if(s == errorreturn) return
if(s < 0d0 .or. s >= 1d0) return

if(drv_t > 0 .or. drv_p > 0) then
  d = sal_dilution_si(s, t_si, p_si)
  if(d == errorreturn) return
  if(d <= 0d0) return
end if

select case (drv_t)

  case (0)

    select case (drv_p)

      case (0)
        sea_vap_s_si = s

      case (1)
        gb_p = sea_g_si(0, 0, 1, s, t_si, p_si)
        if(gb_p == errorreturn) return
        gb_sp = sea_g_si(1, 0, 1, s, t_si, p_si)
        if(gb_sp == errorreturn) return
        gv_p = vap_g_si(0, 1, t_si, p_si)
        if(gv_p == errorreturn) return
        sea_vap_s_si = (gb_p - gv_p - s * gb_sp) / d
    end select

  case (1)

    select case (drv_p)

      case (0)
        gb_t = sea_g_si(0, 1, 0, s, t_si, p_si)
        if(gb_t == errorreturn) return
        gb_st = sea_g_si(1, 1, 0, s, t_si, p_si)
        if(gb_st == errorreturn) return
        gv_t = vap_g_si(1, 0, t_si, p_si)
        if(gv_t == errorreturn) return
        sea_vap_s_si = (gb_t - gv_t - s * gb_st) / d
    end select

end select

end function

!==========================================================================
function sea_vap_salinity_si()
!==========================================================================

!THIS FUNCTION RETURNS THE BRINE SALINITY IN KG/KG OF SEAWATER AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_SALINITY_SI = 0.035

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_SALINITY_SI = 0.035

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_SALINITY_SI = 2.9439629829|4E-02

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_salinity_si

sea_vap_salinity_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_salinity_si = equi_sea_vap_s

end function

!==========================================================================
function sea_vap_temperature_si()
!==========================================================================

!THIS FUNCTION RETURNS THE TEMPERATURE IN K OF VAPOUR AT EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_TEMPERATURE_SI = 274.042416829

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_TEMPERATURE_SI = 274

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_TEMPERATURE_SI = 274

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_temperature_si

sea_vap_temperature_si = errorreturn

if(equi_sea_vap_done /= isok) return

sea_vap_temperature_si = equi_sea_vap_t

end function

!==========================================================================
function sea_vap_vapourpressure_si(sa_si, t_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE VAPOUR PRESSURE IN PA AT THE
!ABSOLUTE TEMPERATURE T_SI IN K AND THE AMBIENT SEAWATER SALINITY SA_SI IN KG/KG

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_VAP_VAPOURPRESSURE_SI(0.035, 274) = 638.044692615

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_vap_vapourpressure_si, sa_si, t_si

sea_vap_vapourpressure_si = errorreturn

if(set_sea_vap_eq_at_s_t(sa_si, t_si) == errorreturn) return

sea_vap_vapourpressure_si = equi_sea_vap_p

end function

!==========================================================================
function sea_vap_volume_evap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE EVAPORATION VOLUME IN M3/KG OF SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_VAP_EQUILIBRIUM_TO_S_P,
!SET_SEA_VAP_EQUILIBRIUM_TO_S_T OR SET_SEA_VAP_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_VAP_EQ_AT_S_P 0.035, 640
!SEA_VAP_VOLUME_EVAP_SI = 197.500648110

!SET_SEA_VAP_EQ_AT_S_T 0.035, 274
!SEA_VAP_VOLUME_EVAP_SI = 198.075461154

!SET_SEA_VAP_EQ_AT_T_P 274, 640
!SEA_VAP_VOLUME_EVAP_SI = 197.469911653

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
!AND, IN THE CASE OF SET_SEA_VAP_EQUILIBRIUM_TO_T_P, OF THOSE IN FLU_3A

implicit none
real*8 sea_vap_volume_evap_si
real*8 s, t, p, d

sea_vap_volume_evap_si = errorreturn

if(equi_sea_vap_done /= isok) return

s = equi_sea_vap_s
t = equi_sea_vap_t
p = equi_sea_vap_p

d = sal_dilution_si(s, t, p)
if(d == errorreturn) return

sea_vap_volume_evap_si = -d * sea_vap_s_si(0, 1, t, p)

end function

!==========================================================================
function set_sea_vap_eq_at_s_p(sa_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-VAPOUR EQUILIBRIUM AT GIVEN SALINITY AND PRESSURE
!BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_vap_eq_at_s_p
real*8 sa_si, p_si
real*8 t, d_liq, d_vap, eps, dvt, dlt, tt

tt = tp_temperature_si
dlt = tp_density_liq_iapws95_si   !TRIPLE POINT DENSITY OF LIQUID WATER
dvt = tp_density_vap_iapws95_si   !TRIPLE POINT DENSITY OF WATER VAPOUR

if(equi_sea_vap_done == isok .and. &
  sa_si == equi_sea_vap_s .and. p_si == equi_sea_vap_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_vap_eq_at_s_p = isok
  return
end if

call clear_sea_vap_state !WASTE ANY PREVIOUS STATE

set_sea_vap_eq_at_s_p = errorreturn

if(sa_si <= 0d0) return
if(sa_si >= 1d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_vap

!SET INITIAL TEMPERATURE GUESS
select case (ctrl_mode_temperature)
  case (0);  t = aux_boilingtemperature_si(sa_si, p_si)
  case (-1); t = tt
  case (1);  t = ctrl_init_t
  case default; return
end select
if(t <= 0d0) return
if(t == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  d_liq = aux_liq_density_if97_si(t, p_si)
  case (-1); d_liq = dlt
  case (1);  d_liq = ctrl_init_d_liq
  case default; return
end select
if(d_liq <= 0d0) return
if(d_liq == errorreturn) return

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0);  d_vap = aux_density_ideal_si(t, p_si)
  case (-1); d_vap = dvt
  case (1);  d_vap = ctrl_init_d_vap
  case default; return
end select
if(d_vap <= 0d0) return
if(d_vap == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
   maxit = 100
elseif(ctrl_loop_maximum == (-1))then
   call set_sea_vap_state (sa_si, t, p_si, d_liq, d_vap)
   set_sea_vap_eq_at_s_p = isok
   return
elseif(ctrl_loop_maximum > 0)then
   maxit = ctrl_loop_maximum
else
   return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_t
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_vap_iteration_at_s_p(sa_si, p_si, maxit, eps, d_liq, d_vap, t) == errorreturn) then
  return
end if

call set_sea_vap_state (sa_si, t, p_si, d_liq, d_vap)

set_sea_vap_eq_at_s_p = isok

end function

!==========================================================================
function set_sea_vap_eq_at_s_t(sa_si, t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-VAPOUR EQUILIBRIUM AT GIVEN TEMPERATURE AND SALINITY
!BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_vap_eq_at_s_t
real*8 sa_si, t_si
real*8 eps
real*8 d_liq, d_vap, p_si, dvt, dlt

dlt = tp_density_liq_iapws95_si   !TRIPLE POINT DENSITY OF LIQUID WATER
dvt = tp_density_vap_iapws95_si   !TRIPLE POINT DENSITY OF WATER VAPOUR

if(equi_sea_vap_done == isok .and. &
  sa_si == equi_sea_vap_s .and. t_si == equi_sea_vap_t) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_vap_eq_at_s_t = isok
  return
end if

call clear_sea_vap_state !WASTE ANY PREVIOUS STATE

set_sea_vap_eq_at_s_t = errorreturn

if(sa_si <= 0d0) return
if(sa_si >= 1d0) return
if(t_si <= 0d0) return

call init_it_ctrl_sea_vap

if(ctrl_mode_liquid == 0 .or. ctrl_mode_vapour == 0) then
  p_si = aux_vapourpressure_wp80_si(sa_si, t_si)
end if

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  d_liq = aux_liq_density_if97_si(t_si, p_si)
  case (-1); d_liq = dlt
  case (1);  d_liq = ctrl_init_d_liq
  case default; return
end select
if(d_liq <= 0d0) return
if(d_liq == errorreturn) return

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0);  d_vap = aux_density_ideal_si(t_si, p_si)
  case (-1); d_vap = dvt
  case (1);  d_vap = ctrl_init_d_vap
  case default; return
end select
if(d_vap <= 0d0) return
if(d_vap == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  p_si = d_vap ** 2 * flu_f_si(0, 1, t_si, d_vap)
  call set_sea_vap_state (sa_si, t_si, p_si, d_liq, d_vap)
  set_sea_vap_eq_at_s_t = isok
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_vap_iteration_at_s_t(sa_si, t_si, maxit, eps, d_liq, d_vap, p_si) == errorreturn) then
  return
end if

call set_sea_vap_state (sa_si, t_si, p_si, d_liq, d_vap)

set_sea_vap_eq_at_s_t = isok

end function

!==========================================================================
function set_sea_vap_eq_at_t_p(t_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE EQUILIBRIUM AT GIVEN ABSOLUTE TEMPERATURE 
!T_SI IN K AND ABSOLUTE PRESSURE P_SI IN PA BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_vap_eq_at_t_p, t_si, p_si
real*8 d_liq, d_vap, s_norm
real*8 s, eps

s_norm = so_salinity_si  !KG/KG OF KCL NORMALISED SEAWATER

if(equi_sea_vap_done == isok .and. &
  t_si == equi_sea_vap_t .and. p_si == equi_sea_vap_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_vap_eq_at_t_p = isok
  return
end if

call clear_sea_vap_state !WASTE ANY PREVIOUS STATE

set_sea_vap_eq_at_t_p = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_vap

!SET INITIAL BRINE SALINITY GUESS
select case (ctrl_mode_sea)
  case (0);  s = aux_brinesalinity_si(t_si, p_si)
  case (-1); s = s_norm
  case (1);  s = ctrl_init_s_sea
  case default; return
end select

if(s <= 0d0 .or. s >= 1d0) return
if(s == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  d_liq = liq_density_si(t_si, p_si)
  d_vap = vap_density_si(t_si, p_si)
  call set_sea_vap_state (s, t_si, p_si, d_liq, d_vap)
  set_sea_vap_eq_at_t_p = isok
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_s
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_vap_iteration_at_t_p(t_si, p_si, maxit, eps, s) == errorreturn) then
  return
end if

d_liq = liq_density_si(t_si, p_si) !ONLY FOR COMPLETENESS OF THE STATE DESCRIPTOR
d_vap = vap_density_si(t_si, p_si)
call set_sea_vap_state (s, t_si, p_si, d_liq, d_vap)

set_sea_vap_eq_at_t_p = isok

end function

!==========================================================================
function sea_vap_iteration_at_s_p(sa_si, p_si, maxit, eps, d_liq_si, d_vap_si, t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN SALINITY, SA_SI, AND PRESSURE, P_SI,
!FROM INITIAL GUESSES FOR THE TEMPERATURE, T_SI, THE LIQUID WATER DENSITY, D_LIQ_SI,
!AND THE VAPOUR DENSITY, D_VAP_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN TEMPERATURE.

!OUTPUT:  
!         SEA_VAP_ITERATION_AT_S_P = ISOK IF SUCCESSFULLY DONE
!         SEA_VAP_ITERATION_AT_S_P = ERRORRETURN IS RETURNED IF
!         - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!         - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!         - DENSITY OR TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!         - S-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
! D_LIQ_SI: DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)
! D_VAP_SI: DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!     T_SI: ABSOLUTE TEMPERATURE IN K OF THE SEAWATER-VAPOUR EQUILIBRIUM

!INPUT: 
!       SA_SI: ABSOLUTE SALINITY IN KG/KG
!        P_SI: ABSOLUTE PRESSURE IN PA
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF TEMPERATURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE
!    D_LIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)
!    D_VAP_SI: INITIAL GUESS FOR DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!        T_SI: INITIAL GUESS FOR THE FREEZING TEMPERATURE


implicit none
integer it, maxit
real*8 sea_vap_iteration_at_s_p, sa_si, p_si, eps, d_liq_si, d_vap_si, t_si
real*8 dl, ddl
real*8 dv, ddv
real*8 t, dt
real*8 mus, p
real*8 gl, gv
real*8 gs_t, gs_st
real*8 fl, fl_t, fl_td
real*8 fl_d, fl_dd
real*8 fv, fv_t, fv_td
real*8 fv_d, fv_dd

real*8 a(3, 3), b(3), x(3)

sea_vap_iteration_at_s_p = errorreturn

if((d_liq_si <= 0d0) .or. &
  (d_vap_si <= 0d0) .or. &
  (t_si <= 0d0)) then
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  t_si = errorreturn
  return
endif

if(check_limits == 1) then
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) then
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    t_si = errorreturn
    return
  endif
else
  if((p_si <= 0d0) .or. &
    (sa_si <= 0d0 .or. sa_si >= 1d0)) then
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    t_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
dv = d_vap_si
dl = d_liq_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF FLUID WATER FOR 3D NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fv = flu_f_si(0, 0, t, dv)
  if(fv == errorreturn) go to 999

  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fv_t = flu_f_si(1, 0, t, dv)
  if(fv_t == errorreturn) go to 999

  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t, dv)
  if(fv_d == errorreturn) go to 999

  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999
  fv_td = flu_f_si(1, 1, t, dv)
  if(fv_td == errorreturn) go to 999

  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t, dv)
  if(fv_dd == errorreturn) go to 999

  p = dv ** 2 * fv_d
  gl = fl + dl * fl_d
  gv = fv + dv * fv_d

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(sa_si, t, p) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  gs_t = sal_g_si(0, 1, 0, sa_si, t, p)
  if(gs_t == errorreturn) go to 999
  gs_st = sal_g_si(1, 1, 0, sa_si, t, p)
  if(gs_st == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = fl_t + dl * fl_td + gs_t - sa_si * gs_st - fv_t - dv * fv_td
  a(1, 2) = 2d0 * fl_d + dl * fl_dd
  a(1, 3) = 2d0 * fv_d + dv * fv_dd

  a(2, 1) = dl ** 2 * fl_td
  a(2, 2) = dl * (2d0 * fl_d + dl * fl_dd)
  a(2, 3) = 0d0

  a(3, 1) = dv ** 2 * fv_td
  a(3, 2) = 0d0
  a(3, 3) = dv * (2d0 * fv_d + dv * fv_dd)

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = gv - gl - mus
  b(2) = p_si - dl ** 2 * fl_d
  b(3) = p_si - dv ** 2 * fv_d

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 3) /= 0) go to 999 !MATRIX SINGULAR
  dt = x(1)
  ddl = x(2)
  ddv = x(3)

  !UPDATE DENSITIES AND TEMPERATURE
  t = t + dt
  if(t <= 0d0) go to 999
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  dv = dv + ddv
  if(dv <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(dt) < eps) .or. &
     (eps < 0d0 .and. abs(dt) < -eps * abs(t))) then
    if(t <= 0d0) go to 999
    sea_vap_iteration_at_s_p = isok
    t_si = t
    d_liq_si = dl
    d_vap_si = dv
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_vap_iteration_at_s_p = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_vap_iteration_at_s_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_vap_iteration_at_s_t(sa_si, t_si, maxit, eps, d_liq_si, d_vap_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN SALINITY, SA_SI, AND TEMPERATURE, T_SI,
!FROM AN INITIAL GUESS FOR THE LIQUID & VAPOUR DENSITY, D_LIQ_SI, D_VAP_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT:  
!          SEA_VAP_ITERATION_AT_S_T = ISOK IF SUCCESSFULLY DONE
!          SEA_VAP_ITERATION_AT_S_T = ERRORRETURN IS RETURNED IF
!          - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!          - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!          - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
! D_LIQ_SI: DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)
! D_VAP_SI: DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!     P_SI: ABSOLUTE PRESSURE IN PA IN THE SEAWATER-ICE EQUILIBRIUM

!INPUT: 
!      SA_SI: ABSOLUTE SALINITY IN KG/KG
!       T_SI: ABSOLUTE TEMPERATURE IN K
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF PRESSURE
!             EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE PRESSURE TOLERANCE
!   D_LIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)
!   D_VAP_SI: INITIAL GUESS FOR DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)


implicit none
integer it, maxit
real*8 sea_vap_iteration_at_s_t, sa_si, t_si, eps, d_liq_si, d_vap_si, p_si
real*8 dl, ddl
real*8 dv, ddv
real*8 mus, p, p_old
real*8 gl, gv
real*8 gs_p, gs_sp
real*8 fl, fl_d, fl_dd
real*8 fv, fv_d, fv_dd

real*8 a(2, 2), b(2), x(2)

sea_vap_iteration_at_s_t = errorreturn

if((d_liq_si <= 0d0) .or. &
  (d_vap_si <= 0d0)) then
  t_si = errorreturn
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax .or. &
  t_si < sal_tmin .or. t_si > sal_tmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
else
  if((t_si <= 0d0) .or. &
    (sa_si < 0d0 .or. sa_si > 1d0)) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

dv = d_vap_si
dl = d_liq_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF FLUID WATER FOR 3D NEWTON ITERATION
  fl = flu_f_si(0, 0, t_si, dl)
  if(fl == errorreturn) go to 999
  fv = flu_f_si(0, 0, t_si, dv)
  if(fv == errorreturn) go to 999

  fl_d = flu_f_si(0, 1, t_si, dl)
  if(fl_d == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t_si, dv)
  if(fv_d == errorreturn) go to 999

  fl_dd = flu_f_si(0, 2, t_si, dl)
  if(fl_dd == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t_si, dv)
  if(fv_dd == errorreturn) go to 999

  p_old = p
  p = dv ** 2 * fv_d
  gl = fl + dl * fl_d
  gv = fv + dv * fv_d

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(sa_si, t_si, p) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  gs_p = sal_g_si(0, 0, 1, sa_si, t_si, p)
  if(gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, sa_si, t_si, p)
  if(gs_sp == errorreturn) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(p - p_old) < eps) .or. &
      (eps < 0d0 .and. abs(p - p_old) < -eps * abs(p))) then
      if(p <= 0d0) go to 999
      sea_vap_iteration_at_s_t = isok
      p_si = p
      d_liq_si = dl
      d_vap_si = dv
      go to 999
    end if
  end if

  if(it == maxit) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = (2d0 * fl_d + dl * fl_dd) * (1d0 + (gs_p - sa_si * gs_sp) * dl)
  a(1, 2) = -2d0 * fv_d - dv * fv_dd

  a(2, 1) = dl * (2d0 * fl_d + dl * fl_dd)
  a(2, 2) = -dv * (2d0 * fv_d + dv * fv_dd)

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = gv - gl - mus
  b(2) = dv ** 2 * fv_d - dl ** 2 * fl_d

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR
  ddl = x(1)
  ddv = x(2)

  !UPDATE DENSITIES
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  dv = dv + ddv
  if(dv <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_vap_iteration_at_s_t = errorreturn
    return
  endif
  !SAL_LIMITS
  if(p_si < sal_pmin .or. p_si > sal_pmax) then
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_vap_iteration_at_s_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_vap_iteration_at_t_p(t_si, p_si, maxit, eps, s_sea_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN TEMPERATURE, T_SI, AND PRESSURE, P_SI,
!FROM AN INITIAL GUESS FOR THE BRINE SALINITY, S_SEA_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN BRINE SALINITY.

!OUTPUT: 
!        SEA_VAP_ITERATION_AT_T_P = ISOK IF SUCCESSFULLY DONE
!        SEA_VAP_ITERATION_AT_T_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS FUNCTION HAS RETURNED AN ERROR
!        - SALINITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
! S_SEA_SI: BRINE SALINITY IN KG/KG IN THE SEAWATER-ICE EQUILIBRIUM

!INPUT: 
!       T_SI: ABSOLUTE TEMPERATURE IN K
!       P_SI: ABSOLUTE PRESSURE IN PA
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF SALINITY
!             EPS > 0: ABSOLUTE SALINITY TOLERANCE IN KG/KG
!             EPS < 0: RELATIVE SALINITY TOLERANCE
!   S_SEA_SI: INITIAL GUESS OF BRINE SALINITY IN KG/KG

implicit none
integer it, maxit
real*8 sea_vap_iteration_at_t_p, t_si, p_si, eps, s_sea_si
real*8 s, ds
real*8 mus, sgss
real*8 gl, gv

sea_vap_iteration_at_t_p = errorreturn

if(s_sea_si <= 0d0 .or. s_sea_si >= 1d0) then
  s_sea_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    s_sea_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) then
    s_sea_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0 .or. &
    t_si <= 0d0) then
    s_sea_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

!GET THE CHEMICAL POTENTIAL OF FLUID WATER (THIS IMPLIES ITERATION IN FLU_3A))
gl = liq_gibbs_energy_si(t_si, p_si)
if(gl == errorreturn) go to 999

gv = vap_gibbs_energy_si(t_si, p_si)
if(gv == errorreturn) go to 999

s = s_sea_si

do it = 1, maxit

  !DERIVATIVES OF THE GIBBS FUNCTIONS FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(s, t_si, p_si) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  sgss = sal_dilution_si(s, t_si, p_si)  !S* D2G/DS2
  if(sgss == errorreturn) go to 999
  if(sgss <= 0) go to 999  !PHYSICALLY, THIS IS IMPOSSIBLE

  !BRINE SALINITY ITERATION STEP
  ds = (gl + mus - gv) / sgss

  !UPDATE BRINE SALINITY
  s = s + ds
  if(s <= 0d0 .or. s >= 1d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(ds) < eps) .or. &
     (eps < 0 .and. abs(ds) < -eps * s)) then
    if(s >= 1d0) go to 999
    sea_vap_iteration_at_t_p = isok
    s_sea_si = s
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !SAL_LIMITS
  if(s_sea_si < sal_smin .or. s_sea_si > sal_smax) then
    s_sea_si = errorreturn
    sea_vap_iteration_at_t_p = errorreturn
    return
  endif
endif

end function

!=========================================================================
function sea_vap_brinefraction_seavap_si(ssv_si, t_si, p_si)
!==========================================================================

!OUTPUT:
!W(SSV,T,P) = SSV/S(T,P)  MASS FRACTION OF BRINE IN SEA VAPOUR (0 < W < 1)

!INPUT:
!SSV_SI    ABSOLUTE SEA-VAPOUR SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A AND ON THOSE OF THIS MODULE

!CHECK VALUE:
!SEA_VAP_BRINEFRACTION_SEAVAP_SI(0.0035,274,640)  =   0.11888736442|5E+00 

implicit none
real*8 sea_vap_brinefraction_seavap_si
real*8 ssv_si, t_si, p_si
real*8 s

sea_vap_brinefraction_seavap_si = errorreturn

if(ssv_si < 0d0 .or. ssv_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

s = sea_vap_brinesalinity_si(t_si, p_si)
if(s == errorreturn) return
if(s <= ssv_si) return

sea_vap_brinefraction_seavap_si = ssv_si / s

end function

!==========================================================================
function aux_boilingtemperature_si(sa_si, p_si)
!==========================================================================

!COMPUTES A BOILING TEMPERATURE ESTIMATE IN K FOR SEAWATER-VAPOUR EQUILIBRIUM AT GIVEN
!BRINE SALINITY IN KG/KG AND ABSOLUTE PRESSURE IN PA FROM CLAUSIUS-CLAPEYRON
!AND RAOULT LAWS

implicit none
real*8 aux_boilingtemperature_si, sa_si, p_si
real*8 tt, pt, a1, a2, alpha
real*8 p, q
real*8 tr, ps

tt = tp_temperature_si
pt = tp_pressure_iapws95_si

a1 = -19.8731005709116d0
a2 = -3.08975437352998d0

alpha = -0.57d0       !RAOULT COEFFICIENT, (P(S)-P(0))/P(S) = ALPHA * S,  S IN KG/KG

aux_boilingtemperature_si = errorreturn
if(p_si <= 0d0) return
if(sa_si < 0d0) return
if(sa_si >= 1d0) return

!PRESSURE IS CORRECTED FOR SALINITY BY RAOULT!S LAW
ps = p_si / (1d0 - alpha * sa_si)

!BOILING TEMPERATURE OF H2O ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_P (TT < T < 640 K WITH RMS = 0.01 IN LN(P/PT))

! CLAUSIUS-CLAPEYRON TYPE EQUATION:
! LN(P/PT) = A1 * (TT/T-1) + A2 * (TT/T-1)^2

p = a1 / a2
q = -log(ps / pt) / a2

tr = -0.5d0 * p + sqrt(0.25d0 * p ** 2 - q)

aux_boilingtemperature_si = tt / (tr + 1d0)

end function

!==========================================================================
function aux_brinesalinity_si(t_si, p_si)
!==========================================================================

!COMPUTES A BRINE SALINITY ESTIMATE IN KG/KG FOR SEAWATER-VAPOUR EQUILIBRIUM AT GIVEN
!ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA FROM WEISS & PRICE 1980 AND
!AND RAOULT LAWS

implicit none
real*8 aux_brinesalinity_si, t_si, p_si
real*8 p0, alpha

alpha = -0.57d0       !RAOULT COEFFICIENT, (P(S)-P(0))/P(S) = ALPHA * S,  S IN KG/KG

aux_brinesalinity_si = errorreturn
if(p_si <= 0d0) return
if(t_si <= 0d0) return

p0 = aux_vapourpressure_wp80_si(0d0, t_si)
aux_brinesalinity_si = (1d0 - p0/p_si)/alpha

end function

!==========================================================================
function aux_vapourpressure_wp80_si(sa_si, t_si)
!==========================================================================

!COMPUTES A VAPOUR PRESSURE ESTIMATE IN PA FOR SEAWATER-VAPOUR EQUILIBRIUM AT GIVEN
!BRINE SALINITY IN KG/KG AND ABSOLUTE TEMPERATURE IN K FROM WEISS & PRICE (1980)

implicit none
real*8 aux_vapourpressure_wp80_si, sa_si, t_si
real*8 spsu
real*8 tabs48_k
real*8 psw

!CRUDE CONVERSION IS SUFFICIENT FOR A GUESS
spsu = 1000d0 * sa_si
tabs48_k = t_si

!WEISS, R.F. AND B.A. PRICE, 1980. NITROUS OXIDE SOLUBILITY IN WATER AND SEAWATER. MARINE CHEM., 8, 347-359.
!SATURATION WATER VAPOUR PRESSURE OF SEAWATER OVER THE TEMPERATURE RANGE 273 TO 313 K AND THE
!SALINITY RANGE 0 TO 40:

psw = exp(24.4543d0 - 67.4509d0*(100d0/tabs48_k) &
      - 4.8489d0*log(tabs48_k/100d0) - 0.000544d0*spsu)

!PSW IS THE WATER VAPOUR PRESSURE (IN ATM), TABS48_K IS THE TEMPERATURE (IN K), 
!AND SPSS IS THE SALINITY ON THE PRACTICAL SALINITY SCALE.

aux_vapourpressure_wp80_si = psw*101325d0

end function

!==========================================================================
function aux_density_ideal_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF IDEAL-GAS VAPOUR AS A FUNCTION OF TEMPERATURE
!AND PRESSURE

!OUTPUT:  
!AUX_DENSITY_IDEAL_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA


implicit none
real*8 aux_density_ideal_si, t_si, p_si, r

r = gas_constant_h2o_si         !SPECIFIC GAS CONSTANT OF WATER IN J/(KG K)

aux_density_ideal_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

aux_density_ideal_si = p_si / (r * t_si)

end function

!=========================================================================
function aux_liq_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 1 (LIQUID) IN IAPWS-IF97

!OUTPUT:  
!AUX_LIQ_DENSITY_IF97_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 aux_liq_density_if97_si
real*8 t_si, p_si
real*8 d

aux_liq_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = aux_liq_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

aux_liq_density_if97_si = 1d0 / d

end function

!=========================================================================
function aux_vap_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF WATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 2 (VAPOUR) IN IAPWS-IF97

!OUTPUT:  
!AUX_VAP_DENSITY_IF97_SI: DENSITY OF LIQUID WATER IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 aux_vap_density_if97_si
real*8 t_si, p_si
real*8 d

aux_vap_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = aux_vap_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

aux_vap_density_if97_si = 1d0 / d

end function

!=========================================================================
function aux_liq_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 1 (LIQUID) IN IAPWS-IF97
!
!OUTPUT: 
!AUX_LIQ_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!DRV_P: ORDER OF THE PRESSURE DERIVATIVE (0-2)
! T_SI:  ABSOLUTE TEMPERATURE IN K
! P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 aux_liq_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 1386d0
pu = 16530000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_SI = 461.51805

aux_liq_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_1

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0); g = gamma_1(0, 0, tau, psi)
        if(g == errorreturn) return
        aux_liq_g_if97_si = rt * g

      case (1); g_p = gamma_1(0, 1, tau, psi)
        if(g_p == errorreturn) return
        aux_liq_g_if97_si = rt * g_p / pu

      case (2); g_pp = gamma_1(0, 2, tau, psi)
        if(g_pp == errorreturn) return
        aux_liq_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0)
        g = gamma_1(0, 0, tau, psi)
        if(g == errorreturn) return
        g_t = gamma_1(1, 0, tau, psi)
        if(g_t == errorreturn) return
        aux_liq_g_if97_si = r * (g - tau * g_t)

      case (1) 
        g_p = gamma_1(0, 1, tau, psi)
        if(g_p == errorreturn) return
        g_tp = gamma_1(1, 1, tau, psi)
        if(g_tp == errorreturn) return
        aux_liq_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
        g_tt = gamma_1(2, 0, tau, psi)
        if(g_tt == errorreturn) return
        aux_liq_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function
 
!=========================================================================
function aux_vap_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 2 (VAPOUR) IN IAPWS-IF97

!OUTPUT: 
!AUX_VAP_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!DRV_P: ORDER OF THE PRESSURE DERIVATIVE (0-2)
! T_SI:  ABSOLUTE TEMPERATURE IN K
! P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 aux_vap_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 540d0
pu = 1000000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_SI = 461.51805

aux_vap_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_2

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0) 
        g = gamma_2(0, 0, tau, psi)
        if(g == errorreturn) return
        aux_vap_g_if97_si = rt * g

      case (1) 
        g_p = gamma_2(0, 1, tau, psi)
        if(g_p == errorreturn) return
        aux_vap_g_if97_si = rt * g_p / pu

      case (2) 
        g_pp = gamma_2(0, 2, tau, psi)
        if(g_pp == errorreturn) return
        aux_vap_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0) 
        g = gamma_2(0, 0, tau, psi)
        if(g == errorreturn) return
        g_t = gamma_2(1, 0, tau, psi)
        if(g_t == errorreturn) return
        aux_vap_g_if97_si = r * (g - tau * g_t)

      case (1) 
        g_p = gamma_2(0, 1, tau, psi)
        if(g_p == errorreturn) return
        g_tp = gamma_2(1, 1, tau, psi)
        if(g_tp == errorreturn) return
        aux_vap_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
        g_tt = gamma_2(2, 0, tau, psi)
        if(g_tt == errorreturn) return
        aux_vap_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function
 
!=========================================================================
function gamma_0(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_0 AS GIVEN IN TABLE 13 OF IF-97

implicit none

integer i, k, drv_t, drv_p
real*8 tau, psi
real*8 gamma_0
real*8 g, pwrt

gamma_0 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

g = 0d0
if(drv_t == 0) then
  if(drv_p == 0) then
    g = log(psi)
  else
    g = 1d0 / psi
    do k = 2, drv_p
      g = dfloat(1 - k) * g / psi
    enddo
  end if
end if

if(drv_p == 0) then
  do i = 1, 9
    pwrt = tau ** (j0i(i) - drv_t)
    do k = 0, drv_t - 1
      pwrt = pwrt * dfloat(j0i(i) - k)
    enddo
    g = g + n0i(i) * pwrt
  enddo
end if

gamma_0 = g

end function
 
!=========================================================================
function gamma_1(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN TABLE 4 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 tau, psi
real*8 gamma_1
real*8 g, pp, tt
real*8 pwrt, pwrp

gamma_1 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

pp = 7.1d0 - psi
tt = tau - 1.222d0

g = 0d0
do i = 1, 34

  if(tt == 0d0) then
    if(j1i(i) == drv_t) then
      pwrt = 1d0
    elseif(j1i(i) > drv_t) then
      pwrt = 0d0
    else
      return
    endif
  else
    pwrt = tt ** (j1i(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(j1i(i) - k)
  enddo

  if(pp == 0d0) then
    if(i1i(i) == drv_p) then
      pwrp = 1d0
    elseif(i1i(i) > drv_p) then
      pwrp = 0d0
    else
      return
    endif
  else
    pwrp = pp ** (i1i(i) - drv_p)
  end if
  do k = 0, drv_p - 1
    pwrp = -pwrp * dfloat(i1i(i) - k)
  enddo

  g = g + n1i(i) * pwrp * pwrt

enddo

gamma_1 = g

end function
 
!=========================================================================
function gamma_2(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN EQ. 15 OF IF-97

implicit none

integer drv_t, drv_p
real*8 tau, psi
real*8 gamma_2
real*8 g0, gr

call init_if97_2

gamma_2 = errorreturn

g0 = gamma_0(drv_t, drv_p, tau, psi)
if(g0 == errorreturn) return

gr = gamma_r(drv_t, drv_p, tau, psi)
if(gr == errorreturn) return

gamma_2 = g0 + gr

end function

!=========================================================================
function gamma_r(drv_t, drv_p, tau, psi)
!=========================================================================
                         
!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_R AS GIVEN IN TABLE 14 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 gamma_r
real*8 tau, psi
real*8 g, tt
real*8 pwrt, pwrp

gamma_r = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

tt = tau - 0.5d0

g = 0d0
do i = 1, 43

  if(tt == 0d0) then
    if(jri(i) == drv_t) then
      pwrt = 1d0
    elseif(jri(i) > drv_t) then
      pwrt = 0d0
    else
      return
    endif
  else
    pwrt = tt ** (jri(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(jri(i) - k)
  enddo

  pwrp = psi ** (iri(i) - drv_p)
  do k = 0, drv_p - 1
    pwrp = pwrp * dfloat(iri(i) - k)
  enddo

  g = g + nri(i) * pwrp * pwrt

enddo

gamma_r = g

end function

!==========================================================================
subroutine set_sea_vap_state(sa_si, t_si, p_si, d_liq_si, d_vap_si)
!==========================================================================

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none
real*8 sa_si, t_si, p_si, d_liq_si, d_vap_si

equi_sea_vap_done = isok

equi_sea_vap_s = sa_si
equi_sea_vap_t = t_si
equi_sea_vap_p = p_si
equi_sea_vap_d_vap = d_vap_si   !DENSITY OF WATER VAPOUR
equi_sea_vap_d_liq = d_liq_si   !DENSITY OF LIQUID WATER

return
end subroutine

!==========================================================================
subroutine clear_sea_vap_state
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

equi_sea_vap_done = 0

return
end subroutine

!==========================================================================
subroutine init_it_ctrl_sea_vap
!==========================================================================

implicit none
real*8 s_norm, tt, dlt, dvt

s_norm = so_salinity_si          !KG/KG OF KCL NORMALISED SEAWATER
tt = tp_temperature_si
dlt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dvt = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR THE ITERATION
ctrl_mode_liquid = 0
ctrl_mode_vapour = 0
ctrl_mode_sea = 0
ctrl_mode_temperature = 0
ctrl_loop_maximum = 100
ctrl_init_s_sea = s_norm
ctrl_init_t = tt
ctrl_init_d_liq = dlt
ctrl_init_d_vap = dvt

ctrl_eps_exit_s = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_t = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_p = -0.0000001d0 !RELATIVE, 0.1 PPM

return
end subroutine

!==========================================================================
subroutine set_it_ctrl_sea_vap(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!SEAWATER-VAPOUR EQUILIBRIUM

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VALUES TO COMPUTE EQUILIBRIUM

!'INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY TO START =
!                                   LIQUID DENSITY CORRELATION FUNCTION, D = AUX_LIQ_DENSITY_IF97_SI(T)
!'INIT_LIQ_DENS'       -1           USE TRIPLE POINT LIQUID DENSITY TO START
!'INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!'INIT_VAP_DENS'        0           USE DEFAULT VAPOUR DENSITY TO START =
!                                   VAPOUR DENSITY CORRELATION FUNCTION, D = AUX_DENSITY_IDEAL_SI(T)
!'INIT_VAP_DENS'       -1           USE TRIPLE POINT VAPOUR DENSITY TO START
!'INIT_VAP_DENS'        D > 0       USE VALUE D AS VAPOUR DENSITY TO START

!'INIT_BRINE_SA'        0           USE DEFAULT BRINE SALINITY TO START ( = AUX_BRINESALINITY_SI(T, P))
!'INIT_BRINE_SA'       -1           USE NORMAL SALINITY TO START
!'INIT_BRINE_SA'        S > 0       USE VALUE S AS BRINE SALINITY TO START

!'INIT_TEMP    '        0           USE DEFAULT TEMPERATURE TO START ( = AUX_BOILINGTEMPERATURE_SI(S, P))
!'INIT_TEMP    '       -1           USE 0 C TEMPERATURE TO START
!'INIT_TEMP    '        T > 0       USE VALUE T AS TEMPERATURE TO START

!'TOL_BRINE_SA '        0           USE DEFAULT EXIT TOLERANCE FOR BRINE SALINITY (0.1 PPM)
!'TOL_BRINE_SA '        EPS         USE EPS AS EXIT TOLERANCE FOR BRINE SALINITY (EPS < 0 MEANS RELATIVE ERROR)

!'TOL_TEMP     '        0           USE DEFAULT EXIT TOLERANCE FOR TEMPERATURE (0.1 PPM)
!'TOL_TEMP     '        EPS         USE EPS AS EXIT TOLERANCE FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!'TOL_PRESS    '        0           USE DEFAULT EXIT TOLERANCE FOR PRESSURE (0.1 PPM)
!'TOL_PRESS    '        EPS         USE EPS AS EXIT TOLERANCE FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

implicit none
character*13 key
real*8 value

call init_it_ctrl_sea_vap

call clear_sea_vap_state

select case (trim(key))

  case ('it_steps')                !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100   !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')           !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_liquid = 0      !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_liquid = nint(value)
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ('init_vap_dens')           !START VAPOUR DENSITY
    if(value == 0d0)then
      ctrl_mode_vapour = 0      !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_vapour = nint(value)
    else
      ctrl_mode_vapour = 1
      ctrl_init_d_vap = value
    endif

  case ('init_brine_sa')           !START BRINE SALINITY
    if(value == 0d0)then
      ctrl_mode_sea = 0         !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_sea = nint(value)
    else
      ctrl_mode_sea = 1
      ctrl_init_s_sea = value
    endif

  case ('init_temp')               !START TEMPERATURE
    if(value == 0d0)then
      ctrl_mode_temperature = 0 !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_temperature = nint(value)
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif

  case ('tol_brine_sa')            !TOLERANCE FOR BRINE SALINITY
    if(value == 0d0)then
      ctrl_eps_exit_s = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_s = value
    endif

  case ('tol_temp')                !TOLERANCE FOR TEMPERATURE
    if(value == 0d0)then
      ctrl_eps_exit_t = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_t = value
    endif

  case ('tol_press')               !TOLERANCE FOR TEMPERATURE
    if(value == 0d0)then
      ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select

return
end subroutine
 
!=========================================================================
subroutine init_if97_1
!=========================================================================

implicit none
integer i

if(n1i(1) == 0.14632971213167d0) return

!TABLE-2. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE DIMENSIONLESS GIBBS FREE ENERGY
!FOR REGION 1, EQ. (7)
!I II JI NI I II JI NI
i = 1; i1i(i)  = 0;  j1i(i) = -2;  n1i(i) = 0.14632971213167d0
i = 2; i1i(i)  = 0;  j1i(i) = -1;  n1i(i) = -0.84548187169114d0
i = 3; i1i(i)  = 0;  j1i(i) = 0;   n1i(i) = -3.756360367204d0
i = 4; i1i(i)  = 0;  j1i(i) = 1;   n1i(i) = 3.3855169168385d0
i = 5; i1i(i)  = 0;  j1i(i) = 2;   n1i(i) = -0.95791963387872d0
i = 6; i1i(i)  = 0;  j1i(i) = 3;   n1i(i) = 0.15772038513228d0
i = 7; i1i(i)  = 0;  j1i(i) = 4;   n1i(i) = -0.016616417199501d0
i = 8; i1i(i)  = 0;  j1i(i) = 5;   n1i(i) = 8.1214629983568d-04
i = 9; i1i(i)  = 1;  j1i(i) = -9;  n1i(i) = 2.8319080123804d-04
i = 10; i1i(i) = 1;  j1i(i) = -7;  n1i(i) = -6.0706301565874d-04
i = 11; i1i(i) = 1;  j1i(i) = -1;  n1i(i) = -0.018990068218419d0
i = 12; i1i(i) = 1;  j1i(i) = 0;   n1i(i) = -0.032529748770505d0
i = 13; i1i(i) = 1;  j1i(i) = 1;   n1i(i) = -0.021841717175414d0
i = 14; i1i(i) = 1;  j1i(i) = 3;   n1i(i) = -5.283835796993d-05  
i = 15; i1i(i) = 2;  j1i(i) = -3;  n1i(i) = -4.7184321073267d-04
i = 16; i1i(i) = 2;  j1i(i) = 0;   n1i(i) = -3.0001780793026d-04
i = 17; i1i(i) = 2;  j1i(i) = 1;   n1i(i) = 4.7661393906987d-05
i = 18; i1i(i) = 2;  j1i(i) = 3;   n1i(i) = -4.4141845330846d-06
i = 19; i1i(i) = 2;  j1i(i) = 17;  n1i(i) = -7.2694996297594d-16
i = 20; i1i(i) = 3;  j1i(i) = -4;  n1i(i) = -3.1679644845054d-05
i = 21; i1i(i) = 3;  j1i(i) = 0;   n1i(i) = -2.8270797985312d-06
i = 22; i1i(i) = 3;  j1i(i) = 6;   n1i(i) = -8.5205128120103d-10
i = 23; i1i(i) = 4;  j1i(i) = -5;  n1i(i) = -2.2425281908d-06
i = 24; i1i(i) = 4;  j1i(i) = -2;  n1i(i) = -6.5171222895601d-07
i = 25; i1i(i) = 4;  j1i(i) = 10;  n1i(i) = -1.4341729937924d-13
i = 26; i1i(i) = 5;  j1i(i) = -8;  n1i(i) = -4.0516996860117d-07
i = 27; i1i(i) = 8;  j1i(i) = -11; n1i(i) = -1.2734301741641d-09
i = 28; i1i(i) = 8;  j1i(i) = -6;  n1i(i) = -1.7424871230634d-10
i = 29; i1i(i) = 21; j1i(i) = -29; n1i(i) = -6.8762131295531d-19
i = 30; i1i(i) = 23; j1i(i) = -31; n1i(i) = 1.4478307828521d-20
i = 31; i1i(i) = 29; j1i(i) = -38; n1i(i) = 2.6335781662795d-23
i = 32; i1i(i) = 30; j1i(i) = -39; n1i(i) = -1.1947622640071d-23
i = 33; i1i(i) = 31; j1i(i) = -40; n1i(i) = 1.8228094581404d-24
i = 34; i1i(i) = 32; j1i(i) = -41; n1i(i) = -9.3537087292458d-26

return
end subroutine

!=========================================================================
subroutine init_if97_2
!=========================================================================

implicit none
integer i

if(n0i(1) == -9.6927686500217d0) return

!TABLE 10. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE IDEAL-GAS PART GAMMA_0 OF THE
!DIMENSIONLESS GIBBS FREE ENERGY FOR REGION 2, EQ. (16)
!I JI NI I JI NI
i = 1; j0i(i) = 0;  n0i(i) = -9.6927686500217d0
i = 2; j0i(i) = 1;  n0i(i) = 10.086655968018d0
i = 3; j0i(i) = -5; n0i(i) = -0.005608791128302d0
i = 4; j0i(i) = -4; n0i(i) = 0.071452738081455d0
i = 5; j0i(i) = -3; n0i(i) = -0.40710498223928d0
i = 6; j0i(i) = -2; n0i(i) = 1.4240819171444d0
i = 7; j0i(i) = -1; n0i(i) = -4.383951131945d0
i = 8; j0i(i) = 2;  n0i(i) = -0.28408632460772d0
i = 9; j0i(i) = 3;  n0i(i) = 0.021268463753307d0
!IF(EQ. (16) IS INCORPORATED INTO EQ. (18), INSTEAD OF THE-VALUES FOR N0I(1) AND N0I(1)
!GIVEN ABOVE, THE FOLLOWING VALUES
!N0I(1) = -0.96937268393049D1
!N0I(2) = 0.10087275970006D2
!SHOULD BE USED

!TABLE 11. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE
!RESIDUAL PART GAMMA_R OF THE DIMENSIONLESS GIBBS FREE ENERGY FOR
!REGION 2, EQ. (17)
!I II JI NI
i = 1;  iri(i) = 1;  jri(i) = 0;  nri(i) = -1.7731742473213d-03
i = 2;  iri(i) = 1;  jri(i) = 1;  nri(i) = -0.017834862292358d0
i = 3;  iri(i) = 1;  jri(i) = 2;  nri(i) = -0.045996013696365d0
i = 4;  iri(i) = 1;  jri(i) = 3;  nri(i) = -0.057581259083432d0
i = 5;  iri(i) = 1;  jri(i) = 6;  nri(i) = -0.05032527872793d0
i = 6;  iri(i) = 2;  jri(i) = 1;  nri(i) = -3.3032641670203d-05
i = 7;  iri(i) = 2;  jri(i) = 2;  nri(i) = -1.8948987516315d-04
i = 8;  iri(i) = 2;  jri(i) = 4;  nri(i) = -3.9392777243355d-03
i = 9;  iri(i) = 2;  jri(i) = 7;  nri(i) = -0.043797295650573d0
i = 10; iri(i) = 2;  jri(i) = 36; nri(i) = -2.6674547914087d-05
i = 11; iri(i) = 3;  jri(i) = 0;  nri(i) = 2.0481737692309d-08
i = 12; iri(i) = 3;  jri(i) = 1;  nri(i) = 4.3870667284435d-07
i = 13; iri(i) = 3;  jri(i) = 3;  nri(i) = -3.227767723857d-05
i = 14; iri(i) = 3;  jri(i) = 6;  nri(i) = -1.5033924542148d-03
i = 15; iri(i) = 3;  jri(i) = 35; nri(i) = -0.040668253562649d0
i = 16; iri(i) = 4;  jri(i) = 1;  nri(i) = -7.8847309559367d-10
i = 17; iri(i) = 4;  jri(i) = 2;  nri(i) = 1.2790717852285d-08
i = 18; iri(i) = 4;  jri(i) = 3;  nri(i) = 4.8225372718507d-07
i = 19; iri(i) = 5;  jri(i) = 7;  nri(i) = 2.2922076337661d-06
i = 20; iri(i) = 6;  jri(i) = 3;  nri(i) = -1.6714766451061d-11
i = 21; iri(i) = 6;  jri(i) = 16; nri(i) = -2.1171472321355d-03
i = 22; iri(i) = 6;  jri(i) = 35; nri(i) = -23.895741934104d0
i = 23; iri(i) = 7;  jri(i) = 0;  nri(i) = -5.905956432427d-18
i = 24; iri(i) = 7;  jri(i) = 11; nri(i) = -1.2621808899101d-06
i = 25; iri(i) = 7;  jri(i) = 25; nri(i) = -0.038946842435739d0
i = 26; iri(i) = 8;  jri(i) = 8;  nri(i) = 1.1256211360459d-11
i = 27; iri(i) = 8;  jri(i) = 36; nri(i) = -0.082311340897998d0
i = 28; iri(i) = 9;  jri(i) = 13; nri(i) = 1.9809712802088d-08
i = 29; iri(i) = 10; jri(i) = 4;  nri(i) = 1.0406965210174d-19
i = 30; iri(i) = 10; jri(i) = 10; nri(i) = -1.0234747095929d-13
i = 31; iri(i) = 10; jri(i) = 14; nri(i) = -1.0018179379511d-09
i = 32; iri(i) = 16; jri(i) = 29; nri(i) = -8.0882908646985d-11
i = 33; iri(i) = 16; jri(i) = 50; nri(i) = 0.10693031879409d0
i = 34; iri(i) = 18; jri(i) = 57; nri(i) = -0.33662250574171d0
i = 35; iri(i) = 20; jri(i) = 20; nri(i) = 8.9185845355421d-25
i = 36; iri(i) = 20; jri(i) = 35; nri(i) = 3.0629316876232d-13
i = 37; iri(i) = 20; jri(i) = 48; nri(i) = -4.2002467698208d-06
i = 38; iri(i) = 21; jri(i) = 21; nri(i) = -5.9056029685639d-26
i = 39; iri(i) = 22; jri(i) = 53; nri(i) = 3.7826947613457d-06
i = 40; iri(i) = 23; jri(i) = 39; nri(i) = -1.2768608934681d-15
i = 41; iri(i) = 24; jri(i) = 26; nri(i) = 7.3087610595061d-29
i = 42; iri(i) = 24; jri(i) = 40; nri(i) = 5.5414715350778d-17
i = 43; iri(i) = 24; jri(i) = 58; nri(i) = -9.436970724121d-07

return
end subroutine

end module sea_vap_4
