module air_3b

!#########################################################################

!THIS MODULE IMPLEMENTS THERMODYNAMIC PROPERTIES OF HUMID AIR AS FUNCTIONS
!OF AIR FRACTION, TEMPERATURE AND PRESSURE, COMPUTED NUMERICALLY FROM
!THE GIBBS FUNCTION OF HUMID AIR.

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., !WRIGHT, D.G., HERRMANN, S.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     AIR_1,       FILE AIR_1.F90
!     AIR_2,       FILE AIR_2.F90
!     AIR_3A,      FILE AIR_3A.F90

!NESTED USE ASSOCIATIONS
!     MATHS_0,     FILE MATHS_0.F90
!     FLU_1,       FILE FLU_1.F90

!#########################################################################

use constants_0
use convert_0
use air_1
use air_2
use air_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: air_g_chempot_vap_si, air_g_compressibilityfactor_si, &
          air_g_contraction_si, air_g_cp_si, air_g_cv_si, air_g_density_si, &
          air_g_enthalpy_si, air_g_entropy_si, air_g_expansion_si, &
          air_g_gibbs_energy_si, air_g_internal_energy_si, &
          air_g_kappa_s_si, air_g_kappa_t_si, air_g_lapserate_si, &
          air_g_soundspeed_si, chk_lemmon_etal_2000

contains

!=========================================================================
function air_g_chempot_vap_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE CHEMICAL POTENTIAL OF VAPOUR IN HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT: 
!CHEMICAL POTENTIAL IN J/KG

!INPUT:  
!A_SI:  MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!AIR_G_CHEMPOT_VAP_SI(0.9,300,1E5) = 193678.583379

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!AIR_G_CHEMPOT_VAP_SI(0.9,300,1E5) = 193705.688576

real*8 air_g_chempot_vap_si, a_si, t_si, p_si
real*8  g, g_a

air_g_chempot_vap_si = errorreturn

if(a_si < 0 .or. a_si > 1) return
if(t_si <= 0) return
if(p_si <= 0) return

!G
g = air_g_si(0, 0, 0, a_si, t_si, p_si)
if(g == errorreturn ) return

!DG/DA
g_a = air_g_si(1, 0, 0, a_si, t_si, p_si)
if(g == errorreturn) return

air_g_chempot_vap_si = g - a_si * g_a  !g - A(dg/dA)

end function

!=========================================================================
function air_g_compressibilityfactor_si(a_si, t_si, p_si)  !dbg
!=========================================================================

!THIS FUNCTION COMPUTES THE COMPRESSIBILITY FACTOR OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!Z = P * V/(RAV * T)  COMPRESSIBILITY FACTOR OF HUMID AIR

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_COMPRESSIBILITYFACTOR_SI(0.9, 300, 1E5) = 0.99782665440423

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_COMPRESSIBILITYFACTOR_SI(0.9, 300, 1E5) = 0.997825670991

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_compressibilityfactor_si, a_si, t_si, p_si
real*8 d, rav

air_g_compressibilityfactor_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return
if(d <= 0d0) return

rav = air_molar_mass_si(a_si)
if(rav == errorreturn) return
if(rav <= 0d0) return
rav = gas_constant_molar_L2000 / rav

air_g_compressibilityfactor_si = p_si / (d * rav * t_si)

end function

!=========================================================================
function air_g_contraction_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE CONTRACTION COEFFICIENT OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:   
!BETA = -(1/V) * (DV/DA)_T_P  CONTRACTION COEFFICIENT OF HUMID AIR IN KG/KG,
!       I.E. THE RELATIVE DENSITY INCREASE IF SOME VAPOUR IS REPLACED BY AIR

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_CONTRACTION_SI(0.9, 300, 1E5) = 0.529966294386

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_CONTRACTION_SI(0.9, 300, 1E5) = 0.530280645260

real*8 air_g_contraction_si, a_si, t_si,p_si
real*8 d, g_ap

air_g_contraction_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

g_ap = air_g_si(1, 0, 1, a_si, t_si, p_si)
if(g_ap == errorreturn) return

air_g_contraction_si = -d * g_ap

end function

!=========================================================================
function air_g_cp_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ISOBARIC HEAT CAPACITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!CP = T * (DETA/DT)_P  ISOBARIC HEAT CAPACITY OF HUMID AIR IN J/(KG K)

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_CP_SI(0.9, 300, 1E5) = 1246.47586191

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_CP_SI(0.9, 300, 1E5) = 1246.34991644

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_cp_si, a_si, t_si, p_si
real*8 d

air_g_cp_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_cp_si = air_f_cp_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_cv_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ISOCHORIC HEAT CAPACITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!CV = T * (D2F/DT2)_D  ISOCHORIC HEAT CAPACITY OF HUMID AIR IN J/(KG K)

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_CV_SI(0.9, 300, 1E5) = 920.676202551

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_CV_SI(0.9, 300, 1E5) = 920.600781012

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_cv_si, a_si, t_si, p_si
real*8 d

air_g_cv_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_cv_si = air_f_cv_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_density_si(a_si, t_si, p_si)
!=========================================================================

!THIS SIMPLY WRAPS THE DENSITY CALL AIR_DENSITY_SI() OF AIR_3A
!FOR A UNIFORM APPEARANCE AIR_G_... IF A, T, P ARE THE ARGUMENTS

!OUTPUT:
!RHO = 1/(DG/DT)_A_T  DENSITY OF HUMID AIR IN KG/M3

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_DENSITY_SI(0.9, 300, 1E5) = 1.09686619009

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_DENSITY_SI(0.9, 300, 1E5) = 1.09708772444

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_density_si, a_si, t_si, p_si

air_g_density_si = air_density_si(a_si, t_si, p_si)

end function

!=========================================================================
function air_g_enthalpy_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ENTHALPY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!H = G - T*(DG/DT)_A_P ENTHALPY OF HUMID AIR IN J/KG

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_ENTHALPY_SI(0.9, 300, 1E5) = 277935.361122

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_ENTHALPY_SI(0.9, 300, 1E5) = 277928.954795

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_enthalpy_si, a_si, t_si, p_si
real*8 d

air_g_enthalpy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_enthalpy_si = air_f_enthalpy_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_entropy_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ENTROPY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!ETA = - (DG/DT)_A_P ENTROPY OF HUMID AIR IN J/(KG K)

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_ENTROPY_SI(0.9, 300, 1E5) = 911.203092171

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_ENTROPY_SI(0.9, 300, 1E5) = 911.170080461

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_entropy_si, a_si, t_si, p_si
real*8 d

air_g_entropy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_entropy_si = air_f_entropy_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_expansion_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!ALPHA = -(1/V) * (DV/DT)_A_P  THERMAL EXPANSION COEFFICIENT OF HUMID AIR IN 1/K

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_EXPANSION_SI(0.9, 300, 1E5) = 3.45696232117E-03

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_EXPANSION_SI(0.9, 300, 1E5) = 3.45704654420E-03

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_expansion_si, a_si, t_si, p_si
real*8 d

air_g_expansion_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_expansion_si = air_f_expansion_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_gibbs_energy_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE GIBBS ENERGY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!G = F + D*(DF/DD)_A_T   GIBBS ENERGY OF HUMID AIR IN J/KG

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_GIBBS_ENERGY_SI(0.9, 300, 1E5) = 4574.43347102

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_GIBBS_ENERGY_SI(0.9, 300, 1E5) = 4577.93065689

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_gibbs_energy_si, a_si, t_si, p_si
real*8 d

air_g_gibbs_energy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_gibbs_energy_si = air_f_gibbs_energy_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_internal_energy_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE INTERNAL ENERGY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!U(T,P) = G - T*(DG/DT)_A_P -P*(DG/DP)_A_T INTERNAL ENERGY OF HUMID AIR IN J/KG

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_INTERNAL_ENERGY_SI(0.9, 300, 1E5) = 186766.537702

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_INTERNAL_ENERGY_SI(0.9, 300, 1E5) = 186778.541048

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_internal_energy_si, a_si, t_si, p_si
real*8 d

air_g_internal_energy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_internal_energy_si = air_f_internal_energy_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_kappa_s_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ISENTROPIC COMPRESSIBILITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!KAPPA_S = (1/V) * (DV/DP)_A_ETA  ISENTROPIC COMPRESSIBILITY OF HUMID AIR IN 1/PA

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_KAPPA_S_SI(0.9, 300, 1E5) = 7.41018952703E-06

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_KAPPA_S_SI(0.9, 300, 1E5) = 7.41034505449E-06

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_kappa_s_si, a_si, t_si, p_si
real*8 d

air_g_kappa_s_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_kappa_s_si = air_f_kappa_s_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_kappa_t_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE ISOTHERMAL COMPRESSIBILITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!KAPPA_T = (1/V) * (DV/DP)_A_T  ISOTHERMAL COMPRESSIBILITY OF HUMID AIR IN 1/PA

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_KAPPA_T_SI(0.9, 300, 1E5) = 1.00324330660E-05

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_KAPPA_T_SI(0.9, 300, 1E5) = 1.00324517749E-05

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_kappa_t_si, a_si, t_si, p_si
real*8 d

air_g_kappa_t_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_kappa_t_si = air_f_kappa_t_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_lapserate_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE "DRY"-ADIABATIC LAPSE RATE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!GAMMA = (DT/DP)_A_ETA  "DRY"-ADIABATIC LAPSE RATE OF HUMID AIR IN K/PA

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_LAPSERATE_SI(0.9, 300, 1E5) = 7.58539809057E-04

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_LAPSERATE_SI(0.9, 300, 1E5) = 7.58481752251E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_lapserate_si, a_si, t_si, p_si
real*8 d

air_g_lapserate_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_lapserate_si = air_f_lapserate_si(a_si, t_si, d)

end function

!=========================================================================
function air_g_soundspeed_si(a_si, t_si, p_si)
!=========================================================================

!THIS FUNCTION COMPUTES THE SOUND SPEED IN HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND PRESSURE

!OUTPUT:
!C = SQRT[(DP/dD)_A_ETA] SOUND SPEED IN HUMID AIR IN M/S

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_G_SOUNDSPEED_SI(0.9, 300, 1E5) = 350.758752725

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_G_SOUNDSPEED_SI(0.9, 300, 1E5) = 350.719656182

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A

real*8 air_g_soundspeed_si, a_si, t_si, p_si
real*8 d

air_g_soundspeed_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = air_density_si(a_si, t_si, p_si)
if(d == errorreturn) return

air_g_soundspeed_si = air_f_soundspeed_si(a_si, t_si, d)

end function

!==========================================================================
subroutine chk_lemmon_etal_2000(mode)
!==========================================================================

!COMPARISON OF COMPUTED SELECTED VALUES FROM THE PAPER
!LEMMON, E.W., JACOBSEN, R.T., PENONCELLO, S.G., FRIEND, D.G.
!THERMODYNAMIC PROPERTIES OF AIR AND MIXTURES OF NITROGEN, ARGON AND OXYGEN
!FROM 60 TO 2000 K AT PRESSURES TO 2000 MPA.
!J. PHYS. CHEM. REF. DATA 29(2000)331-362.

real*8 ma, d, t, x, p
real*8 temp, dens, pres, inte, enth, entr, cv, cp, cs
integer mode, i
!character*130 row
character*90 fmt
character*9 fmt_d
character*4 fmt_t,fmt_p,fmt_i,fmt_enth,fmt_entr,fmt_cv,fmt_cp,fmt_cs

ma = molar_mass_air_L2000    !MOLAR MASS OF AIR IN KG/MOL USED BY LEMMON ET AL 2000

!TEMPERATURE DENSITY  PRESSURE INT.ENERGY ENTHALPY ENTROPY    CV         CP          SOUNDSPEED
!K           MOL/DM3  MPA      J/MOL      J/MOL    J/(MOL K)  J/(MOL K)  J/(MOL K)   M/S
!
!140         0.087718 0.101325 2873.2     4028.3   176.60     20.81      29.38       236.4
!270         0.045164 0.101325 5578.9     7822.4   195.78     20.76      29.13       329.6
!400         0.030461 0.101325 8294.3     11621.0  207.26     21.04      29.38       400.5
!2000        0.006092 0.101325 48610.0    65242.0  259.62     27.90      36.21       863.5
!140         19.810   10.0     -329.46    175.34   118.30     25.82      72.88       418.1
!270         4.6064   10.0     4911.3     7082.2   155.25     21.61      35.28       349.7
!400         2.9202   10.0     7923.0     11347.0  168.19     21.38      31.50       425.6
!2000        0.59094  10.0     48600.0    65522.0  221.44     27.93      36.25       878.6
!270         47.327   2000.0   4354.8     46614.0  96.232     37.64      45.74       2899.8
!400         45.208   2000.0   8076.2     52316.0  113.52     34.45      42.27       2822.9
!2000        32.893   2000.0   53433.0    114240.0 176.45     31.86      38.21       2472.1

!mode = 1: check the Helmholtz function
!mode = 2: check the Gibbs function
                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
if (mode == 1) then
  write(*,*) "Implementation of the Dry-Air HELMHOLTZ function in Fortran"
else
  write(*,*) "Implementation of the Dry-Air GIBBS function in Fortran"
end if
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) ' '
write(*,*) "Selected function values as given in Lemmon et al., "
write(*,*) "J. Phys. Chem. Ref. Data 29(2000)331-362:"
write(*,*) "Note that these values are determined with the original coefficients"
write(*,*) "Top row gives published values, bottom row computed now"
write(*,*) " "

call dry_init_Lemmon2000  !set original adjustable constants

x = errorreturn

write(*,*) "Temp       Density  Pr    IntEnergy Enthalpy  Entropy    Cv        Cp    SndSpd"
write(*,*) "K          mol/dm3  MPa     J/mol    J/mol   J/(mol K) J/(mol K) J/(mol K) m/s"

Do i = 1, 11
  write(*,*) ' '
  select case (i)
    case (1)
      write(*,*) "   140   0.087718 0.101325   2873.2   4028.3   176.60    20.81    29.38    236.4"
      t =  140d0
      d = 0.087718d0
      p = 0.101325d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (2)
      write(*,*) "   270   0.045164 0.101325   5578.9   7822.4   195.78    20.76    29.13    329.6"
      t = 270d0
      d = 0.045164d0
      p = 0.101325d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (3)
      write(*,*) "   400   0.030461 0.101325   8294.3  11621.0   207.26    21.04    29.38    400.5"
      t = 400d0
      d = 0.030461d0
      p = 0.101325d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (4)
      write(*,*) "  2000   0.006092 0.101325  48610.0  65242.0   259.62    27.90    36.21    863.5"
      t = 2000d0
      d = 0.006092d0
      p = 0.101325d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (5)
      write(*,*) "   140     19.810     10.0  -329.46   175.34   118.30    25.82    72.88    418.1"
      t = 140d0
      d = 19.810d0
      p = 10.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (6)
      write(*,*) "   270     4.6064     10.0   4911.3   7082.2   155.25    21.61    35.28    349.7"
      t = 270d0
      d = 4.6064d0
      p = 10.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (7)
      write(*,*) "   400     2.9202     10.0   7923.0  11347.0   168.19    21.38    31.50    425.6"
      t = 400d0
      d = 2.9202d0
      p = 10.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (8)
      write(*,*) "  2000    0.59094     10.0  48600.0  65522.0   221.44    27.93    36.25    878.6"
      t = 2000d0
      d = 0.59094d0
      p = 10.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (9)
      write(*,*) "   270     47.327   2000.0   4354.8  46614.0   96.232    37.64    45.74   2899.8"
      t = 270d0
      d = 47.327d0
      p = 2000.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (10)
      write(*,*) "   400     45.208   2000.0   8076.2  52316.0   113.52    34.45    42.27   2822.9"
      t = 400d0
      d = 45.208d0
      p = 2000.0d0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000d0)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000d0)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000d0)
      else
        temp = t
           x = 0.001 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
    case (11)
      write(*,*) "  2000     32.893   2000.0  53433.0 114240.0   176.45    31.86    38.21   2472.1"
      t = 2000d0
      d = 32.893d0
      p = 2000.0D0
      if (mode == 1) then
        temp = t
        dens = d
           x = 0.000001d0 * air_f_pressure_si(1d0, t, d * ma * 1000d0)
        pres = x
           x = ma * air_f_internal_energy_si(1d0, t, d * ma * 1000d0)
        inte = x
           x = ma * air_f_enthalpy_si(1d0, t, d * ma * 1000d0)
        enth = x
           x = ma * air_f_entropy_si(1d0, t, d * ma * 1000d0)
        entr = x
        cv = ma * air_f_cv_si(1d0, t, d * ma * 1000)
        cp = ma * air_f_cp_si(1d0, t, d * ma * 1000)
        cs = air_f_soundspeed_si(1d0, t, d * ma * 1000)
      else
        temp = t
           x = 0.001d0 * air_g_density_si(1d0, t, p * 1000000d0) / ma
        dens = x
        pres = p
           x = ma * air_g_internal_energy_si(1d0, t, p * 1000000d0)
        inte = x
           x = ma * air_g_enthalpy_si(1d0, t, p * 1000000d0)
        enth = x
           x = ma * air_g_entropy_si(1d0, t, p * 1000000d0)
        entr = x
        cv = ma * air_g_cv_si(1d0, t, p * 1000000d0)
        cp = ma * air_g_cp_si(1d0, t, p * 1000000d0)
        cs = air_g_soundspeed_si(1d0, t, p * 1000000d0)
      endif
  end select

fmt_t = "F8.0"
fmt_d = "F9.6"
if (pres < 1d0) then
   fmt_p = "F8.6"
else 
   fmt_p = "F8.1"
endif
if (dens < 0.1d0) then
  fmt_d = "F9.6"
elseif (dens < 1d0) then
  fmt_d = "1X,F8.5"
elseif (dens < 10d0) then
  fmt_d = "2X,F7.4"
else
  fmt_d = "3X,F6.3"
end if
if (abs(inte) < 1000d0) then
   fmt_i = "F8.2"
else
   fmt_i = "F8.1"
endif
if (abs(enth) < 1000d0) then 
   fmt_enth = "F8.2"
else 
   fmt_enth = "F8.1"
endif
if (x < 100d0) then 
   fmt_entr = "F8.3"
else 
   fmt_entr = "F8.2"
endif
fmt_cv = "F8.2"
fmt_cp = "F8.2"
fmt_cs = "F8.1"
fmt = "("//fmt_t//",1X,"//fmt_d//",1X,"//fmt_p//",1X,"//fmt_i//",1X,"//fmt_enth//",1X,"// &
fmt_entr//",1X,"//fmt_cv//",1X,"//fmt_cp//",1X,"//fmt_cs//")"
!write(*,*)fmt
write(*,fmt)temp, dens, pres, inte, enth, entr, cv, cp, cs
    
enddo

call dry_init_clear  !restore sea-air adjustable coefficients

write(*,*) ' '
write(*,*) "Note that for this comparison, the originally published adjustable coefficients"
write(*,*) "in Lemmon et al. (2000) were used"
write(*,*) " "
write(*,*) "Note: The molar mass of air used in SIA Version 1.1 results in minor offsets"
write(*,*) "Note: For this comparison, check_limits =", check_limits
write(*,*) ' '
write(*,*) ' '

return

!V1.1
!=================================
!PRINT CHK_LEMMON_ETAL_2000(0):

!CHECK FOR PUBLICATION IN OCEAN SCIENCE, 2009
!R. FEISTEL, IOW, VERSION 12 APR 2009
!COMPILED ON 13.04.2009 09:45:57
!
!SELECTED FUNCTION VALUES AS GIVEN IN LEMMON ET AL.,
!J. PHYS. CHEM. REF. DATA 29(2000)331-362:
!
!              TEMPERATURE   DENSITY   PRESSURE  INT.ENERGY  ENTHALPY  ENTROPY     CV          CP            SOUNDSPEED
!              K             MOL/DM3   MPA       J/MOL       J/MOL     J/(MOL K)   J/(MOL K)   J/(MOL K)     M/S
!TABLE VALUE:  140           0.087718  0.101325  2873.2      4028.3    176.60      20.81       29.38         236.4
!COMPUTED:     140           0.087718  0.101324  2873.2      4028.3    176.6       20.8        29.38         236.4
!TABLE VALUE:  270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.6
!COMPUTED:     270           0.045164  0.101324  5578.9      7822.4    195.78      20.76       29.13         329.5
!TABLE VALUE:  400           0.030461  0.101325  8294.3      11621.0   207.26      21.04       29.38         400.5
!COMPUTED:     400           0.030461  0.101326  8294.3      11620.8   207.26      21.04       29.38         400.4
!TABLE VALUE:  2000          0.006092  0.101325  48610.0     65242.0   259.62      27.90       36.21         863.5
!COMPUTED:     2000          0.006092  0.101322  48610.4     65242.4   259.62      27.9        36.21         863.4
!TABLE VALUE:  140           19.810    10.0      -329.46     175.34    118.30      25.82       72.88         418.1
!COMPUTED:     140           19.81     10        -329.5      175.32    118.3       25.82       72.88         418.1
!TABLE VALUE:  270           4.6064    10.0      4911.3      7082.2    155.25      21.61       35.28         349.7
!COMPUTED:     270           4.6064    10        4911.3      7082.2    155.25      21.61       35.28         349.7
!TABLE VALUE:  400           2.9202    10.0      7923.0      11347.0   168.19      21.38       31.50         425.6
!COMPUTED:     400           2.9202    10        7923        11347.4   168.19      21.38       31.5          425.5
!TABLE VALUE:  2000          0.59094   10.0      48600.0     65522.0   221.44      27.93       36.25         878.6
!COMPUTED:     2000          0.59094   10        48599.8     65522     221.44      27.93       36.25         878.5
!TABLE VALUE:  270           47.327    2000.0    4354.8      46614.0   96.232      37.64       45.74         2899.8
!COMPUTED:     270           47.327    2000      4354.9      46615     96.232      37.64       45.74         2899.5
!TABLE VALUE:  400           45.208    2000.0    8076.2      52316.0   113.52      34.45       42.27         2822.9
!COMPUTED:     400           45.208    2000      8076.1      52315.6   113.52      34.45       42.27         2822.5
!TABLE VALUE:  2000          32.893    2000.0    53433.0     114240.0  176.45      31.86       38.21         2472.1
!COMPUTED:     2000          32.893    2000      53432.6     114236.2  176.45      31.86       38.21         2471.8
!
!NOTE: FOR THIS COMPARISON, THE ORIGINALLY PUBLISHED ADJUSTABLE COEFFICIENTS WERE USED

!NOTE: THE MOLAR MASS OF AIR USED IN sia vERSION 1.1 RESULTS IN MINOR OFFSETS

!NOTE: fOR THIS COMPARISON, CHECK_LIMITS = 0 WAS INITILIALIZED 


!V1.1
!=================================
!PRINT CHK_LEMMON_ETAL_2000(1):

!CHECK FOR PUBLICATION IN OCEAN SCIENCE, 2009
!R. FEISTEL, IOW, VERSION 12 APR 2009
!COMPILED ON 13.04.2009 09:46:04
!
!SELECTED FUNCTION VALUES AS GIVEN IN LEMMON ET AL.,
!J. PHYS. CHEM. REF. DATA 29(2000)331-362:
!
!              TEMPERATURE   DENSITY   PRESSURE  INT.ENERGY  ENTHALPY  ENTROPY     CV          CP            SOUNDSPEED
!              K             MOL/DM3   MPA       J/MOL       J/MOL     J/(MOL K)   J/(MOL K)   J/(MOL K)     M/S
!TABLE VALUE:  140           0.087718  0.101325  2873.2      4028.3    176.60      20.81       29.38         236.4
!COMPUTED:     140           0.087718  0.101325  2873.2      4028.3    176.6       20.81       29.38         236.4
!TABLE VALUE:  270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.6
!COMPUTED:     270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.5
!TABLE VALUE:  400           0.030461  0.101325  8294.3      11621.0   207.26      21.04       29.38         400.5
!COMPUTED:     400           0.030461  0.101325  8294.3      11620.8   207.26      21.04       29.38         400.4
!TABLE VALUE:  2000          0.006092  0.101325  48610.0     65242.0   259.62      27.90       36.21         863.5
!COMPUTED:     2000          0.006092  0.101325  48610.4     65242.4   259.62      27.9        36.21         863.4
!TABLE VALUE:  140           19.810    10.0      -329.46     175.34    118.30      25.82       72.88         418.1
!COMPUTED:     140           19.81     10        -329.46     175.34    118.3       25.82       72.88         418.1
!TABLE VALUE:  270           4.6064    10.0      4911.3      7082.2    155.25      21.61       35.28         349.7
!COMPUTED:     270           4.6064    10        4911.3      7082.2    155.25      21.61       35.28         349.7
!TABLE VALUE:  400           2.9202    10.0      7923.0      11347.0   168.19      21.38       31.50         425.6
!COMPUTED:     400           2.9202    10        7923        11347.4   168.19      21.38       31.5          425.5
!TABLE VALUE:  2000          0.59094   10.0      48600.0     65522.0   221.44      27.93       36.25         878.6
!COMPUTED:     2000          0.59094   10        48599.8     65522     221.44      27.93       36.25         878.5
!TABLE VALUE:  270           47.327    2000.0    4354.8      46614.0   96.232      37.64       45.74         2899.8
!COMPUTED:     270           47.327    2000      4354.8      46614.2   96.232      37.64       45.74         2899.5
!TABLE VALUE:  400           45.208    2000.0    8076.2      52316.0   113.52      34.45       42.27         2822.9
!COMPUTED:     400           45.208    2000      8076.2      52316     113.52      34.45       42.27         2822.5
!TABLE VALUE:  2000          32.893    2000.0    53433.0     114240.0  176.45      31.86       38.21         2472.1
!COMPUTED:     2000          32.893    2000      53432.6     114235.9  176.45      31.86       38.21         2471.8
!
!NOTE: FOR THIS COMPARISON, THE ORIGINALLY PUBLISHED ADJUSTABLE COEFFICIENTS WERE USED

!NOTE: THE MOLAR MASS OF AIR USED IN SIA VERSION 1.1 RESULTS IN MINOR OFFSETS

!NOTE: FOR THIS COMPARISON, CHECK_LIMITS = 0 WAS INITILIALIZED 


!V1.0
!=================================
!PRINT CHK_LEMMON_ETAL_2000(0):

!CHECK FOR PUBLICATION IN OCEAN SCIENCE, 2009
!R. FEISTEL, IOW, VERSION 12 APR 2009
!COMPILED ON 13.04.2009 09:45:57
!
!SELECTED FUNCTION VALUES AS GIVEN IN LEMMON ET AL.,
!J. PHYS. CHEM. REF. DATA 29(2000)331-362:
!
!              TEMPERATURE   DENSITY   PRESSURE  INT.ENERGY  ENTHALPY  ENTROPY     CV          CP            SOUNDSPEED
!              K             MOL/DM3   MPA       J/MOL       J/MOL     J/(MOL K)   J/(MOL K)   J/(MOL K)     M/S
!TABLE VALUE:  140           0.087718  0.101325  2873.2      4028.3    176.60      20.81       29.38         236.4
!COMPUTED:     140           0.087718  0.101324  2873.2      4028.3    176.6       20.81       29.38         236.4
!TABLE VALUE:  270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.6
!COMPUTED:     270           0.045164  0.101324  5578.9      7822.4    195.78      20.76       29.13         329.5
!TABLE VALUE:  400           0.030461  0.101325  8294.3      11621.0   207.26      21.04       29.38         400.5
!COMPUTED:     400           0.030461  0.101326  8294.3      11620.8   207.26      21.04       29.38         400.4
!TABLE VALUE:  2000          0.006092  0.101325  48610.0     65242.0   259.62      27.90       36.21         863.5
!COMPUTED:     2000          0.006092  0.101322  48610.4     65242.4   259.62      27.9        36.21         863.4
!TABLE VALUE:  140           19.810    10.0      -329.46     175.34    118.30      25.82       72.88         418.1
!COMPUTED:     140           19.81     10        -329.5      175.32    118.3       25.82       72.88         418.1
!TABLE VALUE:  270           4.6064    10.0      4911.3      7082.2    155.25      21.61       35.28         349.7
!COMPUTED:     270           4.6064    10        4911.3      7082.2    155.25      21.61       35.28         349.6
!TABLE VALUE:  400           2.9202    10.0      7923.0      11347.0   168.19      21.38       31.50         425.6
!COMPUTED:     400           2.9202    10        7923        11347.4   168.19      21.38       31.5          425.5
!TABLE VALUE:  2000          0.59094   10.0      48600.0     65522.0   221.44      27.93       36.25         878.6
!COMPUTED:     2000          0.59094   10        48599.8     65522     221.44      27.93       36.25         878.5
!TABLE VALUE:  270           47.327    2000.0    4354.8      46614.0   96.232      37.64       45.74         2899.8
!COMPUTED:     270           47.327    2000      4354.9      46615     96.232      37.64       45.74         2899.5
!TABLE VALUE:  400           45.208    2000.0    8076.2      52316.0   113.52      34.45       42.27         2822.9
!COMPUTED:     400           45.208    2000      8076.1      52315.6   113.52      34.45       42.27         2822.5
!TABLE VALUE:  2000          32.893    2000.0    53433.0     114240.0  176.45      31.86       38.21         2472.1
!COMPUTED:     2000          32.893    2000      53432.6     114236.2  176.45      31.86       38.21         2471.8
!
!NOTE: FOR THIS COMPARISON, THE ORIGINALLY PUBLISHED ADJUSTABLE COEFFICIENTS WERE USED

!NOTE: FOR THIS COMPARISON, CHECK_LIMITS_VALUE = 0 WAS INITILIALIZED

!V1.0
!=================================
!PRINT CHK_LEMMON_ETAL_2000(1):

!CHECK FOR PUBLICATION IN OCEAN SCIENCE, 2009
!R. FEISTEL, IOW, VERSION 12 APR 2009
!COMPILED ON 13.04.2009 09:46:04
!
!SELECTED FUNCTION VALUES AS GIVEN IN LEMMON ET AL.,
!J. PHYS. CHEM. REF. DATA 29(2000)331-362:
!
!              TEMPERATURE   DENSITY   PRESSURE  INT.ENERGY  ENTHALPY  ENTROPY     CV          CP            SOUNDSPEED
!              K             MOL/DM3   MPA       J/MOL       J/MOL     J/(MOL K)   J/(MOL K)   J/(MOL K)     M/S
!TABLE VALUE:  140           0.087718  0.101325  2873.2      4028.3    176.60      20.81       29.38         236.4
!COMPUTED:     140           0.087718  0.101325  2873.2      4028.3    176.6       20.8        29.38         236.4
!TABLE VALUE:  270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.6
!COMPUTED:     270           0.045164  0.101325  5578.9      7822.4    195.78      20.76       29.13         329.5
!TABLE VALUE:  400           0.030461  0.101325  8294.3      11621.0   207.26      21.04       29.38         400.5
!COMPUTED:     400           0.030461  0.101325  8294.3      11620.8   207.26      21.04       29.38         400.4
!TABLE VALUE:  2000          0.006092  0.101325  48610.0     65242.0   259.62      27.90       36.21         863.5
!COMPUTED:     2000          0.006092  0.101325  48610.4     65242.4   259.62      27.9        36.21         863.4
!TABLE VALUE:  140           19.810    10.0      -329.46     175.34    118.30      25.82       72.88         418.1
!COMPUTED:     140           19.81     10        -329.46     175.34    118.3       25.82       72.88         418.1
!TABLE VALUE:  270           4.6064    10.0      4911.3      7082.2    155.25      21.61       35.28         349.7
!COMPUTED:     270           4.6064    10        4911.3      7082.2    155.25      21.61       35.28         349.7
!TABLE VALUE:  400           2.9202    10.0      7923.0      11347.0   168.19      21.38       31.50         425.6
!COMPUTED:     400           2.9202    10        7923        11347.4   168.19      21.38       31.5          425.5
!TABLE VALUE:  2000          0.59094   10.0      48600.0     65522.0   221.44      27.93       36.25         878.6
!COMPUTED:     2000          0.59094   10        48599.8     65522     221.44      27.93       36.25         878.5
!TABLE VALUE:  270           47.327    2000.0    4354.8      46614.0   96.232      37.64       45.74         2899.8
!COMPUTED:     270           47.327    2000      4354.8      46614.2   96.232      37.64       45.74         2899.5
!TABLE VALUE:  400           45.208    2000.0    8076.2      52316.0   113.52      34.45       42.27         2822.9
!COMPUTED:     400           45.208    2000      8076.2      52316     113.52      34.45       42.27         2822.5
!TABLE VALUE:  2000          32.893    2000.0    53433.0     114240.0  176.45      31.86       38.21         2472.1
!COMPUTED:     2000          32.893    2000      53432.6     114235.9  176.45      31.86       38.21         2471.8
!
!NOTE: fOR THIS COMPARISON, THE ORIGINALLY PUBLISHED ADJUSTABLE COEFFICIENTS WERE USED

!NOTE: FOR THIS COMPARISON, CHECK_LIMITS_VALUE = 0 WAS INITILIALIZED

end subroutine

end module air_3b
