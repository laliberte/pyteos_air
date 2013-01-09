module liq_ice_air_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE EQUILIBRIUM PROPERTIES OF WET ICE AIR, 
!I.E., BETWEEN LIQUID WATER, ICE AND HUMID AIR

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., WRIGHT, D.G., HERRMANN, S.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     CONVERT_0,    FILE CONVERT_0.F90
!     MATHS_0,      FILE MATHS_0.F90
!     AIR_1,        FILE AIR_1.F90
!     FLU_1,        FILE FLU_1.F90
!     ICE_1,        FILE ICE_1.F90
!     AIR_2,        FILE AIR_2.F90
!     FLU_2,        FILE FLU_2.F90
!     ICE_2,        FILE ICE_2.F90
!     AIR_3B,       FILE AIR_3B.F90
!     ICE_LIQ_4,    FILE ICE_LIQ_4.F90

!NESTED USE ASSOCIATIONS
!     AIR_3A,       FILE AIR_3A.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use air_1
use flu_1
use ice_1
use air_2
use flu_2
use ice_2
use air_3b
use ice_liq_4

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: liq_ice_air_airfraction_si, liq_ice_air_density_si, &
          liq_ice_air_dryairfraction_si, liq_ice_air_enthalpy_si, &
          liq_ice_air_entropy_si, liq_ice_air_ifl_si, &
          liq_ice_air_iml_si, liq_ice_air_liquidfraction_si, &
          liq_ice_air_pressure_si, liq_ice_air_solidfraction_si, &
          liq_ice_air_temperature_si, liq_ice_air_vapourfraction_si, &
          set_it_ctrl_liq_ice_air, set_liq_ice_air_eq_at_a, &
          set_liq_ice_air_eq_at_p, set_liq_ice_air_eq_at_t, &
          set_liq_ice_air_eq_at_wa_eta_wt, &
          set_liq_ice_air_eq_at_wa_wl_wi

!CONTROL PARAMETERS OF THE LIQ-ICE-AIR ITERATION
integer ctrl_initialized

integer ctrl_mode_humidity
integer ctrl_mode_temperature
integer ctrl_mode_pressure
integer ctrl_mode_liquid
integer ctrl_mode_air
integer ctrl_loop_maximum

real*8 ctrl_init_a             !initial dry-air fraction of humid air
real*8 ctrl_init_t             !initial temperature
real*8 ctrl_init_p             !initial pressure
real*8 ctrl_init_d_liq         !initial liquid density
real*8 ctrl_init_d_air         !initial humid-air density

real*8 ctrl_eps_exit_t         !temperature exit tolerance
real*8 ctrl_eps_exit_p         !pressure exit tolerance

!PROPERTIES OF THE CURRENT WATER-ICE-HUMID-AIR EQUILIBRIUM STATE
integer equi_liq_ice_air_done  !flag = isok if equilibrium was successfully computed

real*8 equi_liq_ice_air_t      !equilibrium temperature
real*8 equi_liq_ice_air_p      !equilibrium pressure
real*8 equi_liq_ice_air_a      !equilibrium dry-air fraction in humid air
real*8 equi_liq_ice_air_d_liq  !equilibrium density of liquid water at a, t, p
real*8 equi_liq_ice_air_d_air  !equilibrium density of humid air at a, t, p

!IF ONLY A-T-P WAS FIXED BY EQUILIBRIUM SETTING, THE FOLLOWING VALUES MAY REMAIN UNSPECIFIED
real*8 equi_liq_ice_air_wa     !equilibrium dry-air fraction
real*8 equi_liq_ice_air_wt    !equilibrium liquid fraction of the condensed part
real*8 equi_liq_ice_air_eta    !equilibrium entropy

contains

!==========================================================================
function liq_ice_air_airfraction_si()
!==========================================================================

!THIS FUNCTION RETURNS THE FRACTION A IN KG/KG OF DRY AIR IN HUMID AIR
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_A, ETC

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.99

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.961015437341

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.99436473982|6

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.996582564510

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.2

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.99

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.961024307544

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.994366063923

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.996583352944

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_AIRFRACTION_SI = 0.2

real*8 liq_ice_air_airfraction_si

liq_ice_air_airfraction_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

liq_ice_air_airfraction_si = equi_liq_ice_air_a

end function

!==========================================================================
function liq_ice_air_density_si()
!==========================================================================

!FUNCTION COMPUTES THE DENSITY OF WET ICE AIR AFTER PREVIOUSLY CALLING
!SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_WL_WI,... ETC

!Note: the accuracy of this function depends on the iteration settings in this module

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_DENSITY_SI = 7.74649959202

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_DENSITY_SI = 474.974398772

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_DENSITY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_DENSITY_SI = 7.74757979404

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_DENSITY_SI = 474.974398769

real*8 liq_ice_air_density_si
real*8 t, p
real*8 dh, dl, di
real*8 wl, wi

liq_ice_air_density_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

p = liq_ice_air_pressure_si()
if(p == errorreturn) return
t = liq_ice_air_temperature_si()
if(t == errorreturn) return
dh = equi_liq_ice_air_d_air
if(dh == errorreturn) return
dl = equi_liq_ice_air_d_liq
if(dl == errorreturn) return

if(equi_liq_ice_air_wt == errorreturn) return
if(equi_liq_ice_air_wa == errorreturn) return
if(equi_liq_ice_air_a == errorreturn) return
if(equi_liq_ice_air_a == 0d0) then !air-free triple point is not supported here
  return
end if
wl = equi_liq_ice_air_wt * (1d0 - equi_liq_ice_air_wa / equi_liq_ice_air_a)
wi = (1d0 - equi_liq_ice_air_wt) * (1d0 - equi_liq_ice_air_wa / equi_liq_ice_air_a)
di = ice_density_si(t, p)
if(di == errorreturn) return

liq_ice_air_density_si = (1d0 - wi - wl) * dh + wl * dl + wi * di

end function

!==========================================================================
function liq_ice_air_dryairfraction_si()
!==========================================================================

!THIS FUNCTION RETURNS THE FRACTION A OF DRY AIR IN WET ICE AIR IN KG/KG
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_A, ETC

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0 AND V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_DRYAIRFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_DRYAIRFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_DRYAIRFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_DRYAIRFRACTION_SI = 0.99
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_DRYAIRFRACTION_SI = 0.1

real*8 liq_ice_air_dryairfraction_si

liq_ice_air_dryairfraction_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

liq_ice_air_dryairfraction_si = equi_liq_ice_air_wa

end function

!==========================================================================
function liq_ice_air_enthalpy_si()
!==========================================================================

!FUNCTION COMPUTES THE ENTHALPY OF WET ICE AIR AFTER PREVIOUSLY CALLING
!SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_WL_WI,... ETC

!Note: the accuracy of this function depends on the iteration settings in this module

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_ENTHALPY_SI = 7358.21843392

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_ENTHALPY_SI = 900361.142094

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_ENTHALPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_ENTHALPY_SI = 7356.12943724

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_ENTHALPY_SI = 900361.135280

real*8 liq_ice_air_enthalpy_si

real*8 a, t, p
real*8 dh, dl, wet
real*8 wa, wh, wl, wi
real*8 hh, hl, hi

liq_ice_air_enthalpy_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

if(equi_liq_ice_air_wa == errorreturn .or. &
   equi_liq_ice_air_wt == errorreturn) return

a = equi_liq_ice_air_a
t = equi_liq_ice_air_t
p = equi_liq_ice_air_p
dh = equi_liq_ice_air_d_air
dl = equi_liq_ice_air_d_liq
wa = equi_liq_ice_air_wa
wet = equi_liq_ice_air_wt

if(a == 0d0) return                       !THE AIR-FREE TRIPLE POINT CANNOT BE COMPUTED HERE

wh = wa / a                               !WH IS THE HUMID-AIR FRACTION, 1-WH IS THE CONDENSED FRACTION
if(wh < 0d0 .or. wh > 1d0) return
wl = wet * (1d0 - wh)                     !THE LIQUID FRACTION
if(wl < 0d0 .or. wl > 1d0) return
wi = (1d0 - wet) * (1d0 - wh)             !THE ICE FRACTION
if(wi < 0d0 .or. wi > 1d0) return

hh = air_f_enthalpy_si(a, t, dh)          !ENTHALPY OF HUMID AIR
if(hh == errorreturn) return
hl = flu_enthalpy_si(t, dl)               !ENTHALPY OF LIQUID WATER
if(hl == errorreturn) return
hi = ice_enthalpy_si(t, p)                !ENTHALPY OF ICE
if(hi == errorreturn) return

liq_ice_air_enthalpy_si = wh * hh + wl * hl + wi * hi

end function

!==========================================================================
function liq_ice_air_entropy_si()
!==========================================================================

!FUNCTION COMPUTES THE ENTROPY OF WET ICE AIR AFTER PREVIOUSLY CALLING
!SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_WL_WI,... ETC


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999e+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_ENTROPY_SI = 0

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_ENTROPY_SI = 3496.20368786

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999e+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_ENTROPY_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_ENTROPY_SI = 0

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_ENTROPY_SI = 3496.16306903

real*8 liq_ice_air_entropy_si

liq_ice_air_entropy_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

liq_ice_air_entropy_si = equi_liq_ice_air_eta

end function

!==========================================================================
function liq_ice_air_ifl_si(wa_si, eta_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISENTROPIC FREEZING LEVEL (IFL) IN PA
!FROM THE DRY-AIR FRACTION, WA_SI IN KG/KG, AND THE ENTROPY, ETA_SI IN J/(KG K)
!OF WET AIR

!Note: the accuracy of this function depends on the iteration settings in this module

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_ICE_AIR_IFL_SI(.99, 100) = 83240.7163330

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_ICE_AIR_IFL_SI(.99, 100) = 83234.7314360

real*8 liq_ice_air_ifl_si, wa_si, eta_si

liq_ice_air_ifl_si = errorreturn

if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si, eta_si, 1d0) == errorreturn) return

liq_ice_air_ifl_si = equi_liq_ice_air_p

end function

!==========================================================================
function liq_ice_air_iml_si(wa_si, eta_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISENTROPIC MELTING LEVEL (IML) IN PA
!FROM THE DRY-AIR FRACTION, WA_SI IN KG/KG, AND THE ENTROPY, ETA_SI IN J/(KG K)
!OF ICE AIR

!Note: the accuracy of this function depends on the iteration settings in this module

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_ICE_AIR_IML_SI(.99, 100) = 81612.0101048

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_ICE_AIR_IML_SI(.99, 100) = 81605.5557728

real*8 liq_ice_air_iml_si, wa_si, eta_si

liq_ice_air_iml_si = errorreturn

if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si, eta_si, 0d0) == errorreturn) return

liq_ice_air_iml_si = equi_liq_ice_air_p

end function

!==========================================================================
function liq_ice_air_liquidfraction_si()
!==========================================================================

!THIS FUNCTION RETURNS THE FRACTION OF LIQUID WATER IN WET ICE AIR IN KG/KG
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_ETA_WT, ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 3.30256857013E-03

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 3.30296152581E-03

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_LIQUIDFRACTION_SI = 0.2

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 liq_ice_air_liquidfraction_si

liq_ice_air_liquidfraction_si = errorreturn

if(equi_liq_ice_air_done /= isok) return
if(equi_liq_ice_air_wa == errorreturn .or. &
  equi_liq_ice_air_wt == errorreturn) return
if(equi_liq_ice_air_a == 0d0) return

liq_ice_air_liquidfraction_si = equi_liq_ice_air_wt &
  * (1d0 - equi_liq_ice_air_wa / equi_liq_ice_air_a)

end function

!==========================================================================
function liq_ice_air_pressure_si()
!==========================================================================

!THIS FUNCTION RETURNS THE PRESSURE IN PA OF WATER + ICE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_A, ETC

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_PRESSURE_SI = 38347.8983147

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_PRESSURE_SI = 10000

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_PRESSURE_SI = 67931.6010|764

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_PRESSURE_SI = 112016.572992

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_PRESSURE_SI = 706.839968365

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_PRESSURE_SI = 38338.9622424

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_PRESSURE_SI = 10000

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_PRESSURE_SI = 67931.6010|764

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_PRESSURE_SI = 112016.075795

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_PRESSURE_SI = 706.817425301

real*8 liq_ice_air_pressure_si

liq_ice_air_pressure_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

liq_ice_air_pressure_si = equi_liq_ice_air_p

end function

!==========================================================================
function liq_ice_air_solidfraction_si()
!==========================================================================

!THIS FUNCTION RETURNS THE FRACTION OF ICE IN WET ICE AIR IN KG/KG
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_ETA_WT, ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT((0.99, 0, 0.5)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 3.30256857013E-03
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 0.3

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 9.99999999E+98
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT((0.99, 0, 0.5)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 3.30296152581E-03
!
!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_SOLIDFRACTION_SI = 0.3

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 liq_ice_air_solidfraction_si

liq_ice_air_solidfraction_si = errorreturn

if(equi_liq_ice_air_done /= isok) return
if(equi_liq_ice_air_wa == errorreturn .or. &
  equi_liq_ice_air_wt == errorreturn) return
if(equi_liq_ice_air_a == 0d0) return

liq_ice_air_solidfraction_si = (1d0 - equi_liq_ice_air_wt) &
  * (1d0 - equi_liq_ice_air_wa / equi_liq_ice_air_a)

end function

!==========================================================================
function liq_ice_air_temperature_si()
!==========================================================================

!THIS FUNCTION RETURNS THE TEMPERATURE IN K OF WATER + ICE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_A, ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.157197424

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.159302793

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.155

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.151724933

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.159992931

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.157198087

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.159302793

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.155

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.151724970

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_TEMPERATURE_SI = 273.159992933

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 liq_ice_air_temperature_si

liq_ice_air_temperature_si = errorreturn

if(equi_liq_ice_air_done /= isok) return

liq_ice_air_temperature_si = equi_liq_ice_air_t

end function

!==========================================================================
function liq_ice_air_vapourfraction_si()
!==========================================================================

!THIS FUNCTION RETURNS THE FRACTION OF VAPOUR IN WET ICE AIR IN KG/KG
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR,
!SET BY A PREVIOUS CALL OF SET_LIQ_ICE_AIR_EQUILIBRIUM_TO_WA_ETA_WT, ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 3.39486285975E-03

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 0.4

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_LIQ_ICE_AIR_EQ_AT_A(0.99)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_P(1E4)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_T(273.155)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 9.99999999E+98

!SET_LIQ_ICE_AIR_EQ_AT_WA_ETA_WT(0.99, 0, 0.5)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 3.39407694837E-03

!SET_LIQ_ICE_AIR_EQ_AT_WA_WL_WI(.1, .2, .3)
!LIQ_ICE_AIR_VAPOURFRACTION_SI = 0.4

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 liq_ice_air_vapourfraction_si

liq_ice_air_vapourfraction_si = errorreturn

if(equi_liq_ice_air_done /= isok) return
if(equi_liq_ice_air_wa == errorreturn .or. &
  equi_liq_ice_air_wt == errorreturn) return
if(equi_liq_ice_air_a == 0d0) return

liq_ice_air_vapourfraction_si = equi_liq_ice_air_wa * (1d0 / equi_liq_ice_air_a - 1d0)

end function

!==========================================================================
function set_liq_ice_air_eq_at_a(a_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE WATER-ICE-HUMID-AIR EQUILIBRIUM AT GIVEN DRY-AIR FRACTION
!IN HUMID AIR BY NEWTON ITERATION

real*8 set_liq_ice_air_eq_at_a, a_si
real*8 tt, pt, dt_liq, dt_vap

real*8 t, p, dl, dh, eps
real*8 wa, wet, eta
integer maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !triple point density of liquid water
dt_vap = tp_density_vap_iapws95_si  !triple point density of water vapour

if(equi_liq_ice_air_done == isok .and. &
  a_si == equi_liq_ice_air_a) then
 !the requested state has already been computed earlier
 set_liq_ice_air_eq_at_a = isok
 return
end if

call clear_liq_ice_air_state !waste any previous state

set_liq_ice_air_eq_at_a = errorreturn

if(a_si < 0d0) return
if(a_si >= 1d0) return

!AT THIS EQUILIBRIUM, THE MASS FRACTIONS OF DRY AIR, WA, LIQUID, WET,
!AND ENTROPY ETA ARE STILL AMBIGUOUS AND MUST BE SPECIFIED BY ADDITIONAL CONSTRAINTS
wa = errorreturn
wet = errorreturn
eta = errorreturn

if(a_si == 0d0) then  !pure water: set iapws-95 triple point
  call set_liq_ice_air_state(a_si, tt, pt, dt_vap, dt_liq, wa, wet, eta)
  set_liq_ice_air_eq_at_a = isok
  return
end if

call init_it_ctrl_liq_ice_air

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0);  t = tt
  case(-1); t = tt
  case(1);  t = ctrl_init_t
  case default; return
end select
if(t <= 0d0) return
if(t == errorreturn) return

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case(0);  p = aux_pressure_si(a_si)
  case(-1); p = pt
  case(1);  p = ctrl_init_p
  case default; return
end select
if(p <= 0d0) return
if(p == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case(0);  dl = dt_liq
  case(-1); dl = dt_liq
  case(1);  dl = ctrl_init_d_liq
  case default; return
end select
if(dl <= 0d0) return
if(dl == errorreturn) return

!SET INITIAL HUMID-AIR DENSITY GUESS
select case(ctrl_mode_air)
  case(0);  dh = air_g_density_si(a_si, t, p)
  case(-1); dh = dt_vap
  case(1);  dh = ctrl_init_d_air
  case default; return
end select
if(dh <= 0d0) return
if(dh == errorreturn) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    call set_liq_ice_air_state(a_si, t, p, dh, dl, wa, wet, eta)
    set_liq_ice_air_eq_at_a = isok
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(liq_ice_air_iter_at_a(a_si, maxit, eps, t, p, dh, dl) == errorreturn) then
  return
end if

call set_liq_ice_air_state(a_si, t, p, dh, dl, wa, wet, eta)

set_liq_ice_air_eq_at_a = isok

end function

!==========================================================================
function set_liq_ice_air_eq_at_p(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE WATER-ICE-HUMID-AIR EQUILIBRIUM AT GIVEN PRESSURE
! BY NEWTON ITERATION

real*8 set_liq_ice_air_eq_at_p, p_si
real*8 tt, pt, dt_liq, dt_vap

real*8 t, a, dl, dh, eps
real*8 wa, wet, eta
integer maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !triple point density of liquid water
dt_vap = tp_density_vap_iapws95_si  !triple point density of water vapour

if(equi_liq_ice_air_done == isok .and. &
  p_si == equi_liq_ice_air_p) then
  !the requested state has already been computed earlier
  set_liq_ice_air_eq_at_p = isok
  return
end if

call clear_liq_ice_air_state !waste any previous state

set_liq_ice_air_eq_at_p = errorreturn

if(p_si < pt) return

!AT THIS EQUILIBRIUM, THE MASS FRACTIONS OF DRY AIR, WA, LIQUID, WET,
!AND ENTROPY ETA ARE STILL AMBIGUOUS AND MUST BE SPECIFIED BY ADDITIONAL CONSTRAINTS
wa = errorreturn
wet = errorreturn
eta = errorreturn

if(p_si == pt) then  !pure water: set iapws-95 triple point
  call set_liq_ice_air_state(0d0, tt, pt, dt_vap, dt_liq, wa, wet, eta)
  set_liq_ice_air_eq_at_p = isok
  return
end if

call init_it_ctrl_liq_ice_air

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0);  t = tt
  case(-1); t = tt
  case(1);  t = ctrl_init_t
  case default; return
end select
if(t <= 0d0) return
if(t == errorreturn) return

a = air_massfraction_air_si(1d0 - pt / p_si)
if(a == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case(0);  dl = dt_liq
  case(-1); dl = dt_liq
  case(1);  dl = ctrl_init_d_liq
  case default; return
end select
if(dl <= 0d0) return
if(dl == errorreturn) return

!SET INITIAL HUMID-AIR DENSITY GUESS
select case(ctrl_mode_air)
  case(0);  dh = air_g_density_si(a, t, p_si)
  case(-1); dh = dt_vap
  case(1);  dh = ctrl_init_d_air
  case default; return
end select
if(dh <= 0d0) return
if(dh == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0) then
  maxit = 100
elseif(ctrl_loop_maximum == (-1)) then
  call set_liq_ice_air_state(a, t, p_si, dh, dl, wa, wet, eta)
  set_liq_ice_air_eq_at_p = isok
  return
elseif(ctrl_loop_maximum > 0) then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(liq_ice_air_iter_at_p(p_si, maxit, eps, a, t, dh, dl) == errorreturn) then
  return
end if

call set_liq_ice_air_state(a, t, p_si, dh, dl, wa, wet, eta)

set_liq_ice_air_eq_at_p = isok

end function

!==========================================================================
function set_liq_ice_air_eq_at_t(t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE WATER-ICE-HUMID-AIR EQUILIBRIUM AT GIVEN TEMPERATURE
!BY NEWTON ITERATION

real*8 set_liq_ice_air_eq_at_t, t_si
real*8 tt, pt, dt_liq, dt_vap
real*8 p, a, dl, dh, eps
real*8 wa, wet, eta
integer maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !triple point density of liquid water
dt_vap = tp_density_vap_iapws95_si  !triple point density of water vapour

if(equi_liq_ice_air_done == isok .and. &
  t_si == equi_liq_ice_air_t) then
  !the requested state has already been computed earlier
  set_liq_ice_air_eq_at_t = isok
  return
end if

call clear_liq_ice_air_state !waste any previous state

set_liq_ice_air_eq_at_t = errorreturn

if(t_si > tt) return

!AT THIS EQUILIBRIUM, THE MASS FRACTIONS OF DRY AIR, WA, LIQUID, WET,
!AND ENTROPY ETA ARE STILL AMBIGUOUS AND MUST BE SPECIFIED BY ADDITIONAL CONSTRAINTS
wa = errorreturn
wet = errorreturn
eta = errorreturn

if(t_si == tt) then  !pure water: set iapws-95 triple point
  call set_liq_ice_air_state(0d0, tt, pt, dt_vap, dt_liq, wa, wet, eta)
  set_liq_ice_air_eq_at_t = isok
  return
end if

call init_it_ctrl_liq_ice_air

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case(0);  p = ice_liq_meltingpressure_si(t_si)
  case(-1); p = pt
  case(1);  p = ctrl_init_p
  case default; return
end select
if(p <= 0d0) return
if(p == errorreturn) return

a = air_massfraction_air_si(1d0 - pt / p)
if(a == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case(0);  dl = dt_liq
  case(-1); dl = dt_liq
  case(1);  dl = ctrl_init_d_liq
  case default; return
end select
if(dl <= 0d0) return
if(dl == errorreturn) return

!SET INITIAL HUMID-AIR DENSITY GUESS
select case(ctrl_mode_air)
  case(0);  dh = air_g_density_si(a, t_si, p)
  case(-1); dh = dt_vap
  case(1);  dh = ctrl_init_d_air
  case default; return
end select
if(dh <= 0d0) return
if(dh == errorreturn) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    call set_liq_ice_air_state(a, t_si, p, dh, dl, wa, wet, eta)
    set_liq_ice_air_eq_at_t = isok
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(liq_ice_air_iter_at_t(t_si, maxit, eps, a, p, dh, dl) == errorreturn) then
  return
end if

call set_liq_ice_air_state(a, t_si, p, dh, dl, wa, wet, eta)

set_liq_ice_air_eq_at_t = isok

end function

!==========================================================================
function set_liq_ice_air_eq_at_wa_eta_wt(wa_si, eta_si, wet_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE WATER-ICE-HUMID-AIR EQUILIBRIUM AT GIVEN AIR FRACTION WA_SI,
!ENTROPY ETA_SI, AND THE LIQUID FRACTION, WET_SI, OF THE CONDENSED PART BY NEWTON ITERATION.
!IN PARTICULAR, AT WET_SI = 0, THE ISENTROPIC MELTING LEVEL IML IS COMPUTED,
!AND AT WET_SI = 1, THE ISENTROPIC FREEZING LEVEL IFL

real*8 set_liq_ice_air_eq_at_wa_eta_wt, wa_si, eta_si, wet_si
real*8 tt, pt, dt_liq, dt_vap
real*8 a, t, p, dl, dh, eps
integer maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dt_vap = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR

if(equi_liq_ice_air_done == isok .and. &
  wa_si == equi_liq_ice_air_wa .and. &
  eta_si == equi_liq_ice_air_eta .and. &
  wet_si == equi_liq_ice_air_wt) then
  !the requested state has already been computed earlier
  set_liq_ice_air_eq_at_wa_eta_wt = isok
  return
end if

call clear_liq_ice_air_state !waste any previous state

set_liq_ice_air_eq_at_wa_eta_wt = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return

if(wa_si == 0d0) then               !PURE WATER: SET IAPWS-95 TRIPLE POINT
  call set_liq_ice_air_state(0d0, tt, pt, dt_vap, dt_liq, wa_si, wet_si, eta_si)
  set_liq_ice_air_eq_at_wa_eta_wt = isok
  return
end if

call init_it_ctrl_liq_ice_air

!SET INITIAL MASS FRACTIONS GUESS
select case(ctrl_mode_humidity)
  case(0);  a = aux_airfraction_si(wa_si)
  case(-1); a = 0.999
  case(1);  a = ctrl_init_a
  case default; return
end select
if(a < 0d0 .or. a > 1d0) return
if(a == errorreturn) return

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0);  t = tt
  case(-1); t = tt
  case(1);  t = ctrl_init_t
  case default; return
end select
if(t <= 0d0) return
if(t == errorreturn) return

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case(0);  p = aux_pressure_si(a)
  case(-1); p = pt
  case(1);  p = ctrl_init_p
  case default; return
end select
if(p <= 0d0) return
if(p == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case(0);  dl = dt_liq
  case(-1); dl = dt_liq
  case(1);  dl = ctrl_init_d_liq
  case default; return
end select
if(dl <= 0d0) return
if(dl == errorreturn) return

!SET INITIAL HUMID-AIR DENSITY GUESS
select case(ctrl_mode_air)
  case(0);  dh = air_g_density_si(a, t, p)
  case(-1); dh = dt_vap
  case(1);  dh = ctrl_init_d_air
  case default; return
end select
if(dh <= 0d0) return
if(dh == errorreturn) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    call set_liq_ice_air_state(a, t, p, dh, dl, wa_si, wet_si, eta_si)
    set_liq_ice_air_eq_at_wa_eta_wt = isok
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(liq_ice_air_iter_at_wa_eta_wt(wa_si, eta_si, wet_si, maxit,  &
  eps, a, t, p, dh, dl) == errorreturn) then
  return
end if

call set_liq_ice_air_state(a, t, p, dh, dl, wa_si, wet_si, eta_si)

set_liq_ice_air_eq_at_wa_eta_wt = isok

end function

!==========================================================================
function set_liq_ice_air_eq_at_wa_wl_wi(wa_si, wl_si, wi_si)
!==========================================================================

!THIS FUNCTION COMPUTES BY NEWTON INTERATION THE WATER-ICE-HUMID-AIR EQUILIBRIUM
!AT GIVEN AIR FRACTION WA_SI (0 <= WA_SI < 1D0) OF AIR IN KG/KG,
!GIVEN LIQUID FRACTION WL_SI (0 <= WL_SI < 1D0) IN KG/KG,
!AND THE ICE FRACTION WI_SI (0 <= WI_SI < 1D0) IN KG/KG.
!NOTE THAT WA + WL + WI < 1 MUST HOLD FOR A VAPOUR FRACTION WV > 0

real*8 set_liq_ice_air_eq_at_wa_wl_wi, wa_si, wl_si, wi_si
real*8 t, p, a, s, wet
real*8 sh, si, sl
real*8 dh, dl

set_liq_ice_air_eq_at_wa_wl_wi = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(wl_si < 0d0 .or. wl_si > 1d0) return
if(wi_si < 0d0 .or. wi_si > 1d0) return
if(wa_si + wl_si + wi_si >= 1d0) return
if(wl_si + wi_si == 0d0) return

a = wa_si / (1d0 - wl_si - wi_si)  !AIR FRACTION IN HUMID AIR
wet = wl_si / (wl_si + wi_si)      !LIQUID FRACTION OF THE CONDENSED PART

!ITERATIVELY COMPUTE EQUILIBRIUM TEMPERATURE AND PRESSURE
if(set_liq_ice_air_eq_at_a(a) == errorreturn) return

p = liq_ice_air_pressure_si()
if(p == errorreturn) return
t = liq_ice_air_temperature_si()
if(t == errorreturn) return

!COMPUTE ENTROPIES OF THE COMPONENTS
!HUMID AIR:
dh = equi_liq_ice_air_d_air  != air_g_density_si(a, t, p)
if(dh == errorreturn) return
sh = air_f_entropy_si(a, t, dh)
if(sh == errorreturn) return

!ICE:
si = ice_entropy_si(t, p)
if(si == errorreturn) return

!LIQUID WATER
dl = equi_liq_ice_air_d_liq  !liq_density_si(t, p)
if(dl == errorreturn) return
sl = flu_entropy_si(t, dl)
if(sl == errorreturn) return

!TOTAL ENTROPY
s = (1d0 - wl_si - wi_si) * sh + wl_si * sl + wi_si * si

!STORE THE COMPUTED EQUILIBRIUM VALUES
call set_liq_ice_air_state(a, t, p, dh, dl, wa_si, wet, s)

set_liq_ice_air_eq_at_wa_wl_wi = isok

end function

!==========================================================================
function liq_ice_air_iter_at_a(a_si, maxit, eps, t_si, &
                                    p_si, d_air_si, d_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE WATER-ICE-HUMID-AIR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN AIR FRACTION, A_SI,
!FROM AN INITIAL GUESSES FOR THE TEMPERATURE, T_SI, PRESSURE, P_SI,
!THE HUMID-AIR DENSITY, DAIR_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT: 
!        LIQ_ICE_AIR_ITER_AT_A = ISOK IF SUCCESSFULLY DONE
!        LIQ_ICE_AIR_ITER_AT_A = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, TEMPERATURE OR PRESSURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - A-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE OR WATER
!  T_SI: TEMPERATURE IN K OF THE WATER-ICE-AIR EQUILIBRIUM
!  P_SI: PRESSURE IN PA OF THE WATER-ICE-AIR EQUILIBRIUM
!  DAIR_SI: DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID WATER AT (T_SI, P_SI)

!INPUT: 
!        A_SI: AIR FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!              EPS < 0: RELATIVE PRESSURE TOLERANCE
!        T_SI: INITIAL GUESS FOR ABSOLUTE TEMPERATURE IN K
!        P_SI: INITIAL GUESS FOR ABSOLUTE PRESSURE IN PA
!        DAIR_SI: INITIAL GUESS FOR DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!        DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID AT (T_SI, P_SI)

real*8 liq_ice_air_iter_at_a, a_si, eps, t_si, &
       p_si, d_air_si, d_liq_si
real*8 dh, dl, t, p
real*8 ddh, ddl, dt, dp
real*8 gi, gi_t, gi_p
real*8 fh, fh_a, fh_t, fh_td
real*8 fh_d, fh_dd, fh_at, fh_ad
real*8 fl, fl_t, fl_td
real*8 fl_d, fl_dd
real*8 a(4, 4), b(4), x(4)

integer it, maxit

liq_ice_air_iter_at_a = errorreturn

if(a_si <= 0d0 .or. a_si >= 1d0 .or. &
  d_air_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  t_si <= 0d0 .or. &
  p_si <= 0d0) then
  t_si = errorreturn
  p_si = errorreturn
  d_air_si = errorreturn
  d_liq_si = errorreturn
  return
endif

check_limits = check_limits - 1

t = t_si
p = p_si
dl = d_liq_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF HUMID AIR FOR 4D-NEWTON ITERATION
  fh = air_f_si(0, 0, 0, a_si, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a_si, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a_si, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a_si, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a_si, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a_si, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a_si, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a_si, t, dh)
  if(fh_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p)
  if(gi_t == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t, p)
  if(gi_p == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = gi_t - fl_t
  a(1, 2) = gi_p - 1d0 / dl
  a(1, 3) = p / dl**2 - fl_d
  a(1, 4) = 0d0
  a(2, 1) = -dl * fl_td
  a(2, 2) = 1d0 / dl
  a(2, 3) = -fl_d - dl * fl_dd - p / dl**2
  a(2, 4) = 0d0
  a(3, 1) = gi_t - fh_t + a_si * fh_at
  a(3, 2) = gi_p - 1d0 / dh
  a(3, 3) = 0d0
  a(3, 4) = p / dh**2 + a_si * fh_ad - fh_d
  a(4, 1) = -dh * fh_td
  a(4, 2) = 1d0 / dh
  a(4, 3) = 0d0
  a(4, 4) = -p / dh**2 - fh_d - dh * fh_dd

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + p / dl - gi
  b(2) = dl * fl_d - p / dl
  b(3) = -gi + fh + p / dh - a_si * fh_a
  b(4) = dh * fh_d - p / dh
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  dt = x(1)
  dp = x(2)
  ddl = x(3)
  ddh = x(4)
  
  !UPDATE UNKNOWNS
  dh = dh + ddh
  if(dh <= 0d0) go to 999
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  t = t + dt
  if(t <= 0d0) go to 999
  p = p + dp
  if(p <= 0d0) go to 999
  
  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(dp) < eps) .or. &
     (eps < 0d0 .and. abs(dp) < -eps * abs(p))) then
    liq_ice_air_iter_at_a = isok
    t_si = t
    p_si = p
    d_air_si = dh
    d_liq_si = dl
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
    p_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    liq_ice_air_iter_at_a = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    liq_ice_air_iter_at_a = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    liq_ice_air_iter_at_a = errorreturn
    return
  endif
endif

end function

!==========================================================================
function liq_ice_air_iter_at_p(p_si, maxit, eps, a_si, t_si, &
         d_air_si, d_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE WATER-ICE-HUMID-AIR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN PRESSURE, P_SI,
!FROM AN INITIAL GUESSES FOR THE TEMPERATURE, T_SI, AIR FRACTION, A_SI,
!THE HUMID-AIR DENSITY, DAIR_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN AIR FRACTION.

!OUTPUT: 
!        LIQ_ICE_AIR_ITER_AT_P = ISOK IF SUCCESSFULLY DONE
!        LIQ_ICE_AIR_ITER_AT_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, TEMPERATURE OR AIR FRACTION HAVE TAKEN AN INVALID VALUE DURING THE ITERATION
!        - A-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE OR WATER
!  T_SI: TEMPERATURE IN K OF THE WATER-ICE-AIR EQUILIBRIUM
!  A_SI: AIR FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!  DAIR_SI: DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID WATER AT (T_SI, P_SI)

!INPUT: 
!        P_SI: ABSOLUTE PRESSURE IN PA
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF AIR FRACTION
!              EPS > 0: AIR FRACTION TOLERANCE IN KG/KG
!              EPS < 0: RELATIVE AIR FRACTION TOLERANCE
!        A_SI: AIR FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!        T_SI: INITIAL GUESS FOR ABSOLUTE TEMPERATURE IN K
!        DAIR_SI: INITIAL GUESS FOR DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!        DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID AT (T_SI, P_SI)

real*8 liq_ice_air_iter_at_p, p_si, eps, a_si, t_si, &
         d_air_si, d_liq_si
real*8 dh, dl, t
real*8 gi, gi_t, a
real*8 fh, fh_a, fh_t, fh_td
real*8 fh_d, fh_dd, fh_at, fh_ad
real*8 fl, fl_t, fl_td, fh_aa
real*8 fl_d, fl_dd
real*8 aa(4, 4), b(4), x(4)

integer it, maxit

liq_ice_air_iter_at_p = errorreturn

If(d_air_si <= 0) return
If(d_liq_si <= 0) return
If(t_si <= 0) return
If(p_si <= 0) return
If(a_si <= 0) return
If(a_si >= 1) return

if(d_air_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  t_si <= 0d0 .or. &
  a_si <= 0d0 .or. a_si >= 1d0) then
  t_si = errorreturn
  a_si = errorreturn
  d_air_si = errorreturn
  d_liq_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0) then
    t_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
a = a_si
dl = d_liq_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF HUMID AIR FOR 4D-NEWTON ITERATION
  fh = air_f_si(0, 0, 0, a, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_aa = air_f_si(2, 0, 0, a, t, dh)
  if(fh_aa == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a, t, dh)
  if(fh_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p_si)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p_si)
  if(gi_t == errorreturn) go to 999

  !COEFFICIENT MATRIX
  aa(1, 1) = gi_t - fl_t
  aa(1, 2) = p_si / dl**2 - fl_d
  aa(1, 3) = 0d0
  aa(1, 4) = 0d0
  
  aa(2, 1) = -dl * fl_td
  aa(2, 2) = -fl_d - dl * fl_dd - p_si / dl**2
  aa(2, 3) = 0d0
  aa(2, 4) = 0d0
  
  aa(3, 1) = gi_t - fh_t + a * fh_at
  aa(3, 2) = 0d0
  aa(3, 3) = p_si / dh**2 + a * fh_ad - fh_d
  aa(3, 4) = a * fh_aa
  
  aa(4, 1) = -dh * fh_td
  aa(4, 2) = 0d0
  aa(4, 3) = -p_si / dh**2 - fh_d - dh * fh_dd
  aa(4, 4) = -dh * fh_ad

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + p_si / dl - gi
  b(2) = dl * fl_d - p_si / dl
  b(3) = -gi + fh + p_si / dh - a * fh_a
  b(4) = dh * fh_d - p_si / dh

  !SOLVE EQUATIONS
  if(matrix_solve(aa, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  
  !UPDATE UNKNOWNS
  t = t + x(1)
  if(t <= 0d0) go to 999
  dl = dl + x(2)
  if(dl <= 0d0) go to 999
  dh = dh + x(3)
  if(dh <= 0d0) go to 999
  a = a + x(4)
  if(a <= 0d0 .or. a >= 1d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(4)) < eps) .or. &
    (eps < 0d0 .and. abs(x(4)) < -eps * a)) then
    liq_ice_air_iter_at_p = isok
    a_si = a
    t_si = t
    d_air_si = dh
    d_liq_si = dl
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
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_p = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    t_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_p = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  a_si < 0d0 .or. a_si > 1d0 .or. & 
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    t_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function liq_ice_air_iter_at_t(t_si, maxit, eps, a_si, p_si, &
                                    d_air_si, d_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE WATER-ICE-HUMID-AIR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN TEMPERATURE, T_SI,
!FROM AN INITIAL GUESSES FOR THE PRESSURE, P_SI, AIR FRACTION, A_SI,
!THE HUMID-AIR DENSITY, DAIR_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN AIR FRACTION.

!OUTPUT: 
!        LIQ_ICE_AIR_ITER_AT_T = ISOK IF SUCCESSFULLY DONE
!        LIQ_ICE_AIR_ITER_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, TEMPERATURE OR AIR FRACTION HAVE TAKEN AN INVALID VALUE DURING THE ITERATION
!        - A-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE OR WATER
!  P_SI: PRESSURE IN PA OF THE WATER-ICE-AIR EQUILIBRIUM
!  A_SI: AIR FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!  DAIR_SI: DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID WATER AT (T_SI, P_SI)

!INPUT: 
!        T_SI: ABSOLUTE TEMPERATURE IN K
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!              EPS < 0: RELATIVE PRESSURE TOLERANCE
!        A_SI: AIR FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!        P_SI: INITIAL GUESS FOR ABSOLUTE PRESSURE IN PA
!        DAIR_SI: INITIAL GUESS FOR DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!        DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID AT (T_SI, P_SI)

real*8 liq_ice_air_iter_at_t, t_si, eps, a_si, p_si, d_air_si, d_liq_si
real*8 dh, dl, p
real*8 gi, gi_p, a
real*8 fh, fh_a
real*8 fh_d, fh_dd, fh_ad
real*8 fl, fh_aa
real*8 fl_d, fl_dd
real*8 aa(4, 4), b(4), x(4)

integer it, maxit

liq_ice_air_iter_at_t = errorreturn

if((d_air_si <= 0d0) .or. &
  (d_liq_si <= 0d0) .or. &
  (p_si <= 0d0) .or. &
  (a_si <= 0d0) .or. &
  (a_si >= 1d0)) then
  p_si = errorreturn
  a_si = errorreturn
  d_air_si = errorreturn
  d_liq_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
else
  if(t_si <= 0d0) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

p = p_si
a = a_si
dl = d_liq_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF HUMID AIR FOR 4D-NEWTON ITERATION
  fh = air_f_si(0, 0, 0, a, t_si, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a, t_si, dh)
  if(fh_a == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a, t_si, dh)
  if(fh_d == errorreturn) go to 999
  fh_aa = air_f_si(2, 0, 0, a, t_si, dh)
  if(fh_aa == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a, t_si, dh)
  if(fh_ad == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a, t_si, dh)
  if(fh_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t_si, dl)
  if(fl == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t_si, dl)
  if(fl_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t_si, dl)
  if(fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t_si, p)
  if(gi == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t_si, p)
  if(gi_p == errorreturn) go to 999

  !COEFFICIENT MATRIX
  aa(1, 1) = gi_p - 1d0 / dl
  aa(1, 2) = p / dl**2 - fl_d
  aa(1, 3) = 0d0
  aa(1, 4) = 0d0
  
  aa(2, 1) = 1d0 / dl
  aa(2, 2) = -fl_d - dl * fl_dd - p / dl**2
  aa(2, 3) = 0d0
  aa(2, 4) = 0d0
  
  aa(3, 1) = gi_p - 1d0 / dh
  aa(3, 2) = 0d0
  aa(3, 3) = p / dh**2 + a * fh_ad - fh_d
  aa(3, 4) = a * fh_aa
  
  aa(4, 1) = 1d0 / dh
  aa(4, 2) = 0d0
  aa(4, 3) = -p / dh**2 - fh_d - dh * fh_dd
  aa(4, 4) = -dh * fh_ad

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + p / dl - gi
  b(2) = dl * fl_d - p / dl
  b(3) = -gi + fh + p / dh - a * fh_a
  b(4) = dh * fh_d - p / dh

  !SOLVE EQUATIONS
  if(matrix_solve(aa, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  
  !UPDATE UNKNOWNS
  p = p + x(1)
  if(p <= 0d0) go to 999
  dl = dl + x(2)
  if(dl <= 0d0) go to 999
  dh = dh + x(3)
  if(dh <= 0d0) go to 999
  a = a + x(4)
  if(a <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(1)) < eps) .or. &
     (eps < 0d0 .and. abs(x(1)) < -eps * p)) then
    liq_ice_air_iter_at_t = isok
    a_si = a
    p_si = p
    d_air_si = dh
    d_liq_si = dl
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax .or. &
  p_si <= 0d0) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_t = errorreturn
    return
  endif
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_t = errorreturn
    return
  endif
  !AIR_LIMITS
  if(a_si < 0d0 .or. a_si > 1d0 .or. &
  p_si <= 0d0 .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    p_si = errorreturn
    a_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function liq_ice_air_iter_at_wa_eta_wt(wa_si, eta_si, wet_si, &
                                             maxit, eps, a_si, t_si, p_si, &
                                             d_air_si, d_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE WATER-ICE-HUMID-AIR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN AIR FRACTION, WA_SI,
!ENTROPY, ETA_SI, AND THE WETNESS, WET_SI, OF THE CONDENSATE, FROM AN INITIAL GUESS
!FOR THE AIR FRACTION A_SI OF THE HUMID-AIR PART, TEMPERATURE, T_SI, PRESSURE, P_SI,
!THE HUMID-AIR DENSITY, DAIR_SI, AND THE LIQUID DENSITY, DLIQ_SI.
!WET_SI = 0 MEANS THE CONDENSED PART IS ICE, WET_SI = 1 MEANS THE CONDENSED PART IS LIQUID.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT: 
!        LIQ_ICE_AIR_ITER_AT_WA_ETA_WT = ISOK IF SUCCESSFULLY DONE
!        LIQ_ICE_AIR_ITER_AT_WA_ETA_WT = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, TEMPERATURE, PRESSURE OR AIR FRACTION HAVE TAKEN AN INVALID VALUE DURING THE ITERATION
!        - A-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE OR WATER
!     A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!     T_SI: TEMPERATURE IN K OF THE WATER-ICE-AIR EQUILIBRIUM
!     P_SI: PRESSURE IN PA OF THE WATER-ICE-AIR EQUILIBRIUM
!  DAIR_SI: DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID WATER AT (T_SI, P_SI)

!INPUT: 
!       WA_SI: AIR FRACTION OF DRY AIR IN WET ICE AIR IN KG/KG
!      ETA_SI: SPECIFIC ENTROPY OF WET ICE AIR IN J/(KG K)
!      WET_SI: LIQUID FRACTION OF WET ICE IN KG/KG, (1-WET_SI) IS THE ICE FRACTION
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!              EPS < 0: RELATIVE PRESSURE TOLERANCE
!        A_SI: INITIAL GUESS FOR THE DRY-AIR FRACTION IN HUMID AIR IN KG/KG
!        T_SI: INITIAL GUESS FOR ABSOLUTE TEMPERATURE IN K
!        P_SI: INITIAL GUESS FOR ABSOLUTE PRESSURE IN PA
!     DAIR_SI: INITIAL GUESS FOR DENSITY OF HUMID AIR AT (A_SI, T_SI, P_SI)
!     DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID AT (T_SI, P_SI)

real*8 liq_ice_air_iter_at_wa_eta_wt, wa_si, eta_si, wet_si, &
       eps, a_si, t_si, p_si, d_air_si, d_liq_si
real*8 dh, dl, t, p, a
real*8 gi, gi_t, gi_p, gi_tt, gi_tp
real*8 fh, fh_a, fh_t, fh_td
real*8 fh_d, fh_dd, fh_at, fh_ad
real*8 fh_aa, fh_tt
real*8 fl, fl_t, fl_tt, fl_td
real*8 fl_d, fl_dd

!THE 5 UNKNOWNS ARE, INDEX 1-5: T, P, DL, DH, A
real*8 aa(5, 5), b(5), x(5)

integer it, maxit

liq_ice_air_iter_at_wa_eta_wt = errorreturn

If(d_air_si <= 0) return
If(d_liq_si <= 0) return
If(a_si < 0 .or. a_si >= 1) return
If(t_si <= 0) return
If(p_si <= 0) return
If(wa_si < 0 .or. wa_si >= a_si) return
If(wet_si < 0 .or. wet_si > 1) return

if(d_air_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  a_si < 0d0 .or. a_si >= 1d0 .or. &
  wa_si < 0d0 .or. wa_si > 1d0 .or.  &
  wet_si < 0d0 .or. wet_si > 1d0 .or. &
  t_si <= 0d0 .or. &
  p_si <= 0d0) then
  a_si = errorreturn
  t_si = errorreturn
  p_si = errorreturn
  d_air_si = errorreturn
  d_liq_si = errorreturn
  return
endif

check_limits = check_limits - 1

a = a_si
t = t_si
p = p_si
dl = d_liq_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF HUMID AIR FOR 5D-NEWTON ITERATION
  fh = air_f_si(0, 0, 0, a, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_aa = air_f_si(2, 0, 0, a, t, dh)
  if(fh_aa == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_tt = air_f_si(0, 2, 0, a, t, dh)
  if(fh_tt == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a, t, dh)
  if(fh_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fl_tt = flu_f_si(2, 0, t, dl)
  if(fl_tt == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p)
  if(gi_t == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t, p)
  if(gi_p == errorreturn) go to 999
  gi_tt = ice_g_si(2, 0, t, p)
  if(gi_tt == errorreturn) go to 999
  gi_tp = ice_g_si(1, 1, t, p)
  if(gi_tp == errorreturn) go to 999

  !COEFFICIENT MATRIX
  aa(1, 1) = gi_t - fl_t
  aa(1, 2) = gi_p - 1d0 / dl
  aa(1, 3) = p / dl**2 - fl_d
  aa(1, 4) = 0d0
  aa(1, 5) = 0d0
  
  aa(2, 1) = -dl * fl_td
  aa(2, 2) = 1d0 / dl
  aa(2, 3) = -fl_d - dl * fl_dd - p / dl**2
  aa(2, 4) = 0d0
  aa(2, 5) = 0d0

  aa(3, 1) = gi_t - fh_t + a * fh_at
  aa(3, 2) = gi_p - 1d0 / dh
  aa(3, 3) = 0d0
  aa(3, 4) = p / dh**2 + a * fh_ad - fh_d
  aa(3, 5) = a * fh_aa

  aa(4, 1) = -dh * fh_td
  aa(4, 2) = 1d0 / dh
  aa(4, 3) = 0d0
  aa(4, 4) = -p / dh**2 - fh_d - dh * fh_dd
  aa(4, 5) = -dh * fh_ad

  aa(5, 1) = (a - wa_si) * (wet_si * fl_tt + (1d0 - wet_si) * gi_tt) + wa_si * fh_tt
  aa(5, 2) = (a - wa_si) * (1d0 - wet_si) * gi_tp
  aa(5, 3) = (a - wa_si) * wet_si * fl_td
  aa(5, 4) = wa_si * fh_td
  aa(5, 5) = wet_si * fl_t + (1d0 - wet_si) * gi_t + wa_si * fh_at

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + p / dl - gi
  b(2) = dl * fl_d - p / dl
  b(3) = -gi + fh + p / dh - a * fh_a
  b(4) = dh * fh_d - p / dh
  b(5) = (wa_si - a) * (wet_si * fl_t + (1d0 - wet_si) * gi_t) - wa_si * fh_t - a * eta_si
  
  !SOLVE EQUATIONS
  if(matrix_solve(aa, b, x, 5) /= 0) go to 999 !MATRIX SINGULAR
  
  !UPDATE UNKNOWNS
  t = t + x(1)
  if(t <= 0d0) go to 999
  p = p + x(2)
  if(p <= 0d0) go to 999
  dl = dl + x(3)
  if(dl <= 0d0) go to 999
  dh = dh + x(4)
  if(dh <= 0d0) go to 999
  a = a + x(5)
  if(a < 0d0 .or. a >= 1d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(2)) < eps) .or. &
     (eps < 0d0 .and. abs(x(2)) < -eps * abs(p))) then
    liq_ice_air_iter_at_wa_eta_wt = isok
    a_si = a
    t_si = t
    p_si = p
    d_air_si = dh
    d_liq_si = dl
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    a_si = errorreturn
    t_si = errorreturn
    p_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_wa_eta_wt = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) then
    a_si = errorreturn
    t_si = errorreturn
    p_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_wa_eta_wt = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    a_si = errorreturn
    t_si = errorreturn
    p_si = errorreturn
    d_air_si = errorreturn
    d_liq_si = errorreturn
    liq_ice_air_iter_at_wa_eta_wt = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_airfraction_si(wa_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED DRY-AIR FRACTION A IN KG/KG OF HUMID AIR
!AT GIVEN DRY-AIR FRACTION WA IN KG/KG OF WET ICE AIR

!A = WA /(WA + WV)                   !DRY-AIR FRACTION OF HUMID AIR, WV = VAPOUR FRACTION
!WL = WET / (1D0 - WA / A)           !LIQUID FRACTION OF WET ICE AIR
!WI = (1D0 - WET) / (1D0 - WA / A)   !ICE FRACTION OF WET ICE AIR
!WA + WV + WL + WI = 1               !MASS BALANCE: AIR + VAPOUR + LIQUID + ICE = 100%

real*8 aux_airfraction_si, wa_si
real*8 wv

aux_airfraction_si = errorreturn

!WI+WL=WV, ROUGH GUESS FOR THE VAPOUR FRACTION WV, 50% OF THE WATER MAY BE VAPOUR
wv = 0.5d0 * (1d0 - wa_si)
if(wa_si + wv <= 0d0) return

aux_airfraction_si = wa_si / (wa_si + wv)

end function

!==========================================================================
function aux_pressure_si(a_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED TRIPLE-POINT PRESSURE IN PA AT GIVEN AIR FRACTION A_SI
!IN KG/KG OF THE HUMID-AIR PART

real*8 aux_pressure_si, a_si
real*8 xv
real*8 pt

pt = tp_pressure_iapws95_si

aux_pressure_si = errorreturn

!ASSUME THE PARTIAL VAPOUR PRESSURE = TRIPLE-POINT PRESSURE
xv = air_molfraction_vap_si(a_si)
if(xv == errorreturn) return
if(xv <= 0d0 .or. xv > 1d0) return

aux_pressure_si = pt / xv

end function

!==========================================================================
subroutine set_it_ctrl_liq_ice_air(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!WET ICE AIR PROPERTIES

!KEY              VALUE
!IT_STEPS         0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS         N > 0    SET ITERATION NUMBER TO N
!IT_STEPS        -1           DO NOT ITERATE, USE INITIAL VALUES

!INIT_AIR         0           USE DEFAULT AIR DENSITY TO START ( = RHO(CTRL_INIT_P, CTRL_INIT_T))
!INIT_AIR         D > 0    USE VALUE D AS AIR DENSITY TO START

!INIT_LIQ         0           USE DEFAULT LIQUID DENSITY TO START ( = RHO(CTRL_INIT_P, CTRL_INIT_T))
!INIT_LIQ         D > 0    USE VALUE D AS LIQUID DENSITY TO START

!INIT_HUM         0           USE DEFAULT AIR FRACTION TO START ( = AUX_AIRFRACTION_SI(WA_SI))
!INIT_HUM         A > 0    USE VALUE A AS AIR FRACTION TO START

!INIT_TEMP        0           USE DEFAULT TEMPERATURE TO START ( = TRIPLE POINT TEMPERATURE)
!INIT_TEMP        T > 0    USE VALUE T AS TEMPERATURE TO START

!INIT_PRESS       0           USE DEFAULT PRESSURE TO START ( = AUX_PRESSURE_SI(WA_SI))
!INIT_PRESS       P > 0    USE VALUE P AS PRESSURE TO START

!TOL_TEMP         0           USE DEFAULT EXIT ACCURACY FOR TEMPERATURE (0.1 PPM)
!TOL_TEMP         EPS         USE EPS AS EXIT ACCURACY FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_PRESS        0           USE DEFAULT EXIT ACCURACY FOR PRESSURE (0.1 PPM)
!TOL_PRESS        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

character*13 key
real*8 value

call init_it_ctrl_liq_ice_air

call clear_liq_ice_air_state

select case (trim(key))

  case ("it_steps");               !ITERATION STEPS
    if(value == 0d0) then
      ctrl_loop_maximum = 100   !DEFAULT = 100
    elseif(value < 0d0) then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ("init_air");               !START AIR DENSITY
    if(value == 0d0) then
      ctrl_mode_air = 0         !DEFAULT = RHO(CTRL_INIT_P, CTRL_INIT_T)
    else
      ctrl_mode_air = 1
      ctrl_init_d_air = value
    endif
    
  case ("init_liq");               !START LIQUID WATER DENSITY
    if(value == 0d0) then
      ctrl_mode_liquid = 0      !DEFAULT = RHO(CTRL_INIT_P, CTRL_INIT_T)
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ("init_hum");               !START AIR FRACTION
    if(value == 0d0) then
      ctrl_mode_humidity = 0    !DEFAULT = AUX FUNCTION
    else
      ctrl_mode_humidity = 1
      ctrl_init_a = value
    endif
    
  case ("init_temp");              !START TEMPERATURE
    if(value == 0d0) then
      ctrl_mode_temperature = 0 !DEFAULT = TRIPLE POINT TEMPERATURE
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif

  case ("init_press");             !START PRESSURE
    if(value == 0d0) then
      ctrl_mode_pressure = 0    !DEFAULT = AUX FUNCTION
    else
      ctrl_mode_pressure = 1
      ctrl_init_p = value
    endif

  case ("tol_temp");               !REQUIRED TEMPERATURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_t = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_t = value
    endif

  case ("tol_press");              !REQUIRED TOTAL PRESSURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select


end subroutine

!==========================================================================
subroutine clear_liq_ice_air_state()
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_liq_ice_air_done = 0

end subroutine

!==========================================================================
subroutine set_liq_ice_air_state(a_si, t_si, p_si, d_air_si, d_liq_si, &
                                 wa_si, wet_si, eta_si)
!==========================================================================
                                  
!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

real*8 a_si, t_si, p_si, d_air_si, d_liq_si, wa_si, wet_si, eta_si

equi_liq_ice_air_done = isok

!INTENSIVE PROPERTIES, WELL DEFINED AT THE TRIPLE POINT
equi_liq_ice_air_a = a_si          !TRIPLE-POINT AIR FRACTION
equi_liq_ice_air_t = t_si          !TRIPLE-POINT TEMPERATURE
equi_liq_ice_air_p = p_si          !TRIPLE-POINT PRESSURE
equi_liq_ice_air_d_liq = d_liq_si   !DENSITY OF PURE WATER
equi_liq_ice_air_d_air = d_air_si   !DENSITY OF HUMID AIR

!"EXTENSIVE PROPERTIES", DEPEND ON ADDITIONAL CONSTRAINTS
!AND MAY BE UNDEFINED (I.E. = ERRORRETURN) IN CERTAIN COMPUTED STATES

equi_liq_ice_air_wa = wa_si        !DRY-AIR MASS FRACTION IN KG/KG OF WET ICE AIR
equi_liq_ice_air_wt = wet_si      !LIQUID-WATER MASS FRACTION IN KG/KG OF THE CONDENSED PART
equi_liq_ice_air_eta = eta_si      !TOTAL SPECIFIC ENTROPY OF THE PARCEL

end subroutine

!==========================================================================
subroutine init_it_ctrl_liq_ice_air()
!==========================================================================

real*8 tt, pt, dlt, dvt

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dlt = tp_density_liq_iapws95_si    !TRIPLE POINT DENSITY OF LIQUID WATER
dvt = tp_density_vap_iapws95_si    !TRIPLE POINT DENSITY OF WATER VAPOUR

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR ITERATION
ctrl_mode_humidity = 0
ctrl_mode_liquid = 0
ctrl_mode_air = 0
ctrl_mode_temperature = 0
ctrl_mode_pressure = 0
ctrl_loop_maximum = 100

ctrl_init_a = 0.999d0
ctrl_init_t = tt
ctrl_init_p = pt
ctrl_init_d_liq = dlt
ctrl_init_d_air = dvt

ctrl_eps_exit_t = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_p = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

end module liq_ice_air_4
