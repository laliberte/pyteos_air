module ice_air_4a

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILIBRIUM PROPERTIES OF ICE IH
!WITH HUMID AIR, COMMONLY REGARDED AS "SATURATED AIR" OR "HUMIDITY 100%".
!THE AIR PROPERTIES COMPUTED HERE REFER TO SATURATED AIR BELOW THE FROST POINT.

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
!     ICE_1,        FILE ICE_1.F90
!     AIR_2,        FILE AIR_2.F90
!     ICE_2,        FILE ICE_2.F90
!     AIR_3A,       FILE AIR_3A.F90
!     AIR_3B,       FILE AIR_3B.F90
!     ICE_VAP_4,    FILE ICE_VAP_4.F90

!NESTED USE ASSOCIATIONS
!     FLU_1,        FILE FLU_1.F90
!     FLU_2,        FILE FLU_2.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use air_1
use ice_1
use air_2
use ice_2
use air_3a
use air_3b
use ice_vap_4

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: ice_air_a_from_rh_cct_si, ice_air_a_from_rh_wmo_si, &
          ice_air_condensationpressure_si, ice_air_density_air_si, &
          ice_air_density_ice_si, ice_air_density_vap_si, &
          ice_air_enthalpy_subl_si, ice_air_frostpoint_si, &
          ice_air_icl_si, ice_air_ict_si, ice_air_massfraction_air_si, &
          ice_air_pressure_si, ice_air_rh_cct_from_a_si, &
          ice_air_rh_wmo_from_a_si, ice_air_sublimationpressure_si, &
          ice_air_temperature_si, &
          set_ice_air_eq_at_a_eta, set_ice_air_eq_at_a_p, &
          set_ice_air_eq_at_a_t, set_ice_air_eq_at_t_p, &
          set_it_ctrl_ice_air

!CONTROL PARAMETERS OF THE VAPOUR PRESSURE ITERATION
integer ctrl_initialized

!SELECT INITIAL GUESSING METHOD
integer ctrl_mode_humidity      !HUMIDTY
integer ctrl_mode_temperature   !TEMPERATURE
integer ctrl_mode_pressure      !PRESSURE

!STORE CUSTOMIZED INITIAL VALUES HERE
integer ctrl_loop_maximum
real*8 ctrl_init_t
real*8 ctrl_init_a
real*8 ctrl_init_p

!SET TOLERANCE OF THE ITERATION LOOP FOR
real*8 ctrl_eps_exit_p_vap   !VAPOUR PARTIAL PRESSURE
real*8 ctrl_eps_exit_temp    !TEMPERATURE
real*8 ctrl_eps_exit_press   !PRESSURE

!PROPERTIES OF THE CURRENT LIQUID-VAPOUR EQUILIBRIUM STATE
integer equi_ice_air_done

!STORAGE OF THE ITERATION RESULTS = THE EQUILIBRIUM STATE (SATURATED AIR)
real*8 equi_ice_air_a        !DRY AIR MASS FRACTION OF HUMID AIR
real*8 equi_ice_air_t        !TEMPERATURE
real*8 equi_ice_air_p        !PRESSURE
real*8 equi_ice_air_eta      !HUMID AIR ENTROPY
real*8 equi_ice_air_d_air    !HUMID AIR DENSITY

contains

!==========================================================================
function ice_air_a_from_rh_wmo_si(rh_si, t_si, p_si)
!==========================================================================

!RETURNS THE AIR MASS FRACTION FROM RELATIVE HUMIDITYWITH RESPECT TO ICE, TEMPERATURE AND PRESSURE

!OUTPUT: 
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG

!INPUT:  
!RH_SI RELATIVE HUMIDITY, WMO DEFINITION, 0 < RH < 1
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_A_FROM_RH_WMO_SI(0.8,270,1E5) = 0.997645144450

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_A_FROM_RH_WMO_SI(0.8,270,1E5) = 0.997645698908

real*8 ice_air_a_from_rh_wmo_si, rh_si, t_si, p_si
real*8 a_sat

ice_air_a_from_rh_wmo_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return
if(rh_si < 0d0 .or. rh_si > 1d0) return

if(rh_si == 0d0) then !DRY AIR
  ice_air_a_from_rh_wmo_si = 1d0
  return
end if

a_sat = ice_air_massfraction_air_si(t_si, p_si)
if(a_sat == errorreturn) return
if(a_sat <= 0d0 .or. a_sat > 1d0) return

!NOTE: RH_WMO = (1-A)*A_SAT/(1-A_SAT)*A (ALL MASS FRACTIONS)
ice_air_a_from_rh_wmo_si = 1d0 / (1d0 + rh_si * (1d0 / a_sat - 1d0))

end function

!==========================================================================
function ice_air_a_from_rh_cct_si( rh_si, t_si, p_si)
!==========================================================================

!RETURNS THE AIR MASS FRACTION FROM RELATIVE HUMIDITY WITH RESPECT TO ICE, TEMPERATURE AND PRESSURE

!OUTPUT: 
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG

!INPUT:  
!RH_SI RELATIVE HUMIDITY, CCT/IUPAC DEFINITION, 0 < RH < 1
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_A_FROM_RH_CCT_SI(0.8,270,1E5) = 0.997647370806

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_A_FROM_RH_CCT_SI(0.8,270,1E5) = 0.997647924743

real*8 ice_air_a_from_rh_cct_si, rh_si, t_si, p_si
real*8 a_sat

ice_air_a_from_rh_cct_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return
if(rh_si < 0d0 .or. rh_si > 1d0) return

if(rh_si == 0d0) then !DRY AIR
  ice_air_a_from_rh_cct_si = 1d0
  return
end if

a_sat = ice_air_massfraction_air_si(t_si, p_si)  !SATURATED AIR MASS FRACTION
if(a_sat == errorreturn) return
if(a_sat <= 0d0 .or. a_sat > 1d0) return
 
a_sat = air_molfraction_vap_si(a_sat)   !SATURATED AIR MOLE FRACTION

ice_air_a_from_rh_cct_si = air_massfraction_air_si(1d0 - rh_si * a_sat)

end function

!==========================================================================
function ice_air_massfraction_air_si(t_si, p_si)
!==========================================================================

!RETURNS THE AIR MASS FRACTION IN KG/KG OF SATURATED AIR WITH RESPECT TO ICE

!OUTPUT: 
!A_SI  MASS FRACTION OF DRY AIR IN SATURATED HUMID AIR, IN KG/KG

!INPUT:  
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_MASSFRACTION_AIR_SI(270,1E5) = 0.997058162464

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_MASSFRACTION_AIR_SI(270,1E5) = 0.997058854720

real*8 ice_air_massfraction_air_si, t_si, p_si

ice_air_massfraction_air_si = errorreturn

if(set_ice_air_eq_at_t_p(t_si, p_si) == errorreturn) return

ice_air_massfraction_air_si = equi_ice_air_a

end function

!==========================================================================
function ice_air_condensationpressure_si( a_si, t_si)
!==========================================================================

!RETURNS THE CONDENSATION PRESSURE OF ICE FROM HUMID AIR

!OUTPUT: 
!P_SI  ABSOLUTE PRESSURE IN PA

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!T_SI  ABSOLUTE TEMPERATURE IN K

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_CONDENSATIONPRESSURE_SI(0.997,270) = 98057.5687060

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_CONDENSATIONPRESSURE_SI(0.997,270) = 98034.4511233

real*8 ice_air_condensationpressure_si, a_si, t_si

ice_air_condensationpressure_si = errorreturn

if(set_ice_air_eq_at_a_t(a_si, t_si) == errorreturn) return

ice_air_condensationpressure_si = equi_ice_air_p

end function

!==========================================================================
function ice_air_sublimationpressure_si(t_si, p_si)
!==========================================================================

!RETURNS THE PARTIAL PRESSURE OF VAPOUR IN HUMID AIR AT EQUILIBRIUM WITH ICE,
!FROM TEMPERATURE AND TOTAL PRESSURE

!INPUT:  T_SI  ABSOLUTE TEMPERATURE IN K
!        P_SI  ABSOLUTE (TOTAL) PRESSURE IN PA
!OUTPUT: SUBLIMATION PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!ICE_AIR_SUBLIMATIONPRESSURE_SI(270,1E5) = 472.041436186

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!ICE_AIR_SUBLIMATIONPRESSURE_SI(270,1E5) = 472.041823975

real*8 ice_air_sublimationpressure_si, t_si, p_si
real*8 a_sat, xv

ice_air_sublimationpressure_si = errorreturn

if(t_si <= 0) return
if(p_si <= 0) return

!GET THE AIR FRACTION AT SATURATION WITH RESPECT TO ICE
a_sat = ice_air_massfraction_air_si(t_si, p_si)
if(a_sat == errorreturn) return
if(a_sat <= 0 .or. a_sat > 1) return

!CONVERT AIR MASS FRACTION TO VAPOUR MOLE FRACTION
xv = air_molfraction_vap_si(a_sat)
if(xv == errorreturn) return
if(xv <= 0 .or. xv > 1) return

!COMPUTE THE PARTIAL PRESSURE
ice_air_sublimationpressure_si = xv * p_si

end function

!==========================================================================
function ice_air_density_air_si()
!==========================================================================

!RETURNS THE HUMID AIR DENSITY IN KG/M3 OF SATURATED AIR WITH RESPECT TO ICE
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_density_air_si = 1.28849644315

!set_ice_air_eq_at_a_t .997, 270
!ice_air_density_air_si = 1.26340826697

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_density_air_si = 1.28734106365

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_density_air_si = 0.949136077891

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_density_air_si = 1.28880078014

!set_ice_air_eq_at_a_t .997, 270
!ice_air_density_air_si = 1.26340801028

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_density_air_si = 1.28763121402

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_density_air_si = 0.949325026119

real*8 ice_air_density_air_si

ice_air_density_air_si = errorreturn

if(equi_ice_air_done /= isok) return

ice_air_density_air_si = equi_ice_air_d_air

end function

!==========================================================================
function ice_air_density_ice_si()
!==========================================================================

!RETURNS THE ICE DENSITY IN KG/M3 OF SATURATED AIR WITH RESPECT TO ICE
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p(270, 1E5)
!ice_air_density_ice_si = 917.181167192

!set_ice_air_eq_at_a_t(.997, 270)
!ice_air_density_ice_si = 917.180958346

!set_ice_air_eq_at_a_p(.997, 1E5)
!ice_air_density_ice_si = 917.147466108

!set_ice_air_eq_at_a_eta(.997, 100)
!ice_air_density_ice_si = 917.682093944

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p(270, 1E5)
!ice_air_density_ice_si = 917.181167192

!set_ice_air_eq_at_a_t(.997, 270)
!ice_air_density_ice_si = 917.180955861

!set_ice_air_eq_at_a_p(.997, 1E5)
!ice_air_density_ice_si = 917.147060527

!set_ice_air_eq_at_a_eta(.997, 100)
!ice_air_density_ice_si = 917.681749114


real*8 ice_air_density_ice_si

ice_air_density_ice_si = errorreturn

if(equi_ice_air_done /= isok) return

ice_air_density_ice_si = ice_density_si(equi_ice_air_t, equi_ice_air_p)

end function

!==========================================================================
function ice_air_density_vap_si()
!==========================================================================

!RETURNS THE VAPOUR DENSITY IN KG/M3 OF SATURATED AIR WITH RESPECT TO ICE
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_density_vap_si = 3.79054720202E-03

!set_ice_air_eq_at_a_t .997, 270
!ice_air_density_vap_si = 3.79022480092E-03

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_density_vap_si = 3.86202319096E-03

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_density_vap_si = 2.84740823367E-03

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_density_vap_si = 3.79055033080E-03

!set_ice_air_eq_at_a_t .997, 270
!ice_air_density_vap_si = 3.79022403085E-03

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_density_vap_si = 3.86289364206E-03

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_density_vap_si = 2.84797507836E-03


real*8 ice_air_density_vap_si

ice_air_density_vap_si = errorreturn

if(equi_ice_air_done /= isok) return

ice_air_density_vap_si = equi_ice_air_d_air * (1d0 - equi_ice_air_a)

end function

!==========================================================================
function ice_air_enthalpy_subl_si()
!==========================================================================

!RETURNS THE SUBLIMATION HEAT IN J/KG OF SATURATED AIR BELOW THE FROST POINT
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_enthalpy_subl_si = 2833359.68775

!set_ice_air_eq_at_a_t .997, 270
!ice_air_enthalpy_subl_si = 2833386.63215

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_enthalpy_subl_si = 2833297.67344

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_enthalpy_subl_si = 2834613.25618

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_enthalpy_subl_si = 2833359.27614

!set_ice_air_eq_at_a_t .997, 270
!ice_air_enthalpy_subl_si = 2833386.54980

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_enthalpy_subl_si = 2833296.51317

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_enthalpy_subl_si = 2834612.42351

real*8 ice_air_enthalpy_subl_si
real*8 a, t, dh
real*8 hh, hh_a, hi
real*8 fh_a, fh_d
real*8 fh_at, fh_ad
real*8 fh_td, fh_dd
real*8 p_a, p_d
real*8 h_a, h_d

ice_air_enthalpy_subl_si = errorreturn

if(equi_ice_air_done /= isok) return

!GET THE EQUILIBRIUM PROPERTIES
a = equi_ice_air_a
if(a == errorreturn) return
if(a < 0d0 .or. a > 1d0) return
t = equi_ice_air_t
if(t == errorreturn) return
if(t <= 0d0) return
dh = equi_ice_air_d_air
if(dh == errorreturn) return
if(dh <= 0d0) return

!ENTHALPY OF HUMID AIR
hh = air_f_enthalpy_si(a, t, dh)
if(hh == errorreturn) return

!DERIVATIVE HH_A = (DH/DA)_T_P OF ENTHALPY OF HUMID AIR
fh_a = air_f_si(1, 0, 0, a, t, dh)
if(fh_a == errorreturn) return
fh_d = air_f_si(0, 0, 1, a, t, dh)
if(fh_d == errorreturn) return
fh_at = air_f_si(1, 1, 0, a, t, dh)
if(fh_at == errorreturn) return
fh_ad = air_f_si(1, 0, 1, a, t, dh)
if(fh_ad == errorreturn) return
fh_td = air_f_si(0, 1, 1, a, t, dh)
if(fh_td == errorreturn) return
fh_dd = air_f_si(0, 0, 2, a, t, dh)
if(fh_dd == errorreturn) return

p_a = dh * fh_ad                            !DP/DA AT CONST T,D
p_d = dh * (2d0 * fh_d + dh * fh_dd)        !DP/DD AT CONST A,T
if(p_d == 0d0) return
h_a = fh_a + dh * fh_ad - t * fh_at         !DH/DA AT CONST T,D
h_d = 2d0 * fh_d + dh * fh_dd - t * fh_td   !DH/DD AT CONST A,T

hh_a = h_a - h_d * p_a / p_d                !DH/DA AT CONST T,P

!ENTHALPY OF ICE
hi = ice_enthalpy_si(t, equi_ice_air_p)
if(hi == errorreturn) return

!COMPUTE THE LATENT HEAT
ice_air_enthalpy_subl_si = hh - a * hh_a - hi

end function

!==========================================================================
function ice_air_frostpoint_si(a_si, p_si)
!==========================================================================

!RETURNS THE FROSTPOINT TEMPERATURE OF HUMID AIR

!OUTPUT: 
!T_SI  ABSOLUTE FROSTPOINT TEMPERATURE IN K

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_FROSTPOINT_SI(0.997,1E5) = 270.232024746

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_FROSTPOINT_SI(0.997,1E5) = 270.234816126

real*8 ice_air_frostpoint_si, a_si, p_si

ice_air_frostpoint_si = errorreturn

if(set_ice_air_eq_at_a_p(a_si, p_si) == errorreturn) return

ice_air_frostpoint_si = equi_ice_air_t

end function

!==========================================================================
function ice_air_icl_si(a_si, t_si, p_si)
!==========================================================================

!RETURNS THE ISENTROPIC CONDENSATION LEVEL ICL OF ICE FROM HUMID AIR

!OUTPUT: 
!ICL_SI  ABSOLUTE ISENTROPIC CONDENSATION PRESSURE IN PA

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_ICL_SI(0.997,300,1E5) = 64985.6857989

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_ICL_SI(0.997,300,1E5) = 64988.3931838

real*8 ice_air_icl_si, a_si, t_si, p_si
real*8 eta

ice_air_icl_si = errorreturn

eta = air_g_entropy_si(a_si, t_si, p_si)
if(eta == errorreturn) return

if(set_ice_air_eq_at_a_eta(a_si, eta) == errorreturn) return

ice_air_icl_si = equi_ice_air_p

end function

!==========================================================================
function ice_air_ict_si(a_si, t_si, p_si)
!==========================================================================

!RETURNS THE ISENTROPIC CONDENSATION TEMPERATURE ICT OF ICE FROM HUMID AIR

!OUTPUT: 
!ICT_SI  ABSOLUTE ISENTROPIC CONDENSATION TEMPERATURE IN K

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_ICT_SI(0.997,300,1E5) = 265.221832079

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_ICT_SI(0.997,300,1E5) = 265.224998411

real*8 ice_air_ict_si, a_si, t_si, p_si
real*8 eta

ice_air_ict_si = errorreturn

eta = air_g_entropy_si(a_si, t_si, p_si)
if(eta == errorreturn) return

if(set_ice_air_eq_at_a_eta(a_si, eta) == errorreturn) return

ice_air_ict_si = equi_ice_air_t

end function

!==========================================================================
function ice_air_pressure_si()
!==========================================================================

!RETURNS THE PRESSURE IN PA OF SATURATED AIR WITH RESPECT TO ICE
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_pressure_si = 100000

!set_ice_air_eq_at_a_t .997, 270
!ice_air_pressure_si = 98057.5687060

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_pressure_si = 100000

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_pressure_si = 72723.467705

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_pressure_si = 100000

!set_ice_air_eq_at_a_t .997, 270
!ice_air_pressure_si = 98034.4511233

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_pressure_si = 100000

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_pressure_si = 72721.4579415

real*8 ice_air_pressure_si

ice_air_pressure_si = errorreturn

if(equi_ice_air_done /= isok) return

ice_air_pressure_si = equi_ice_air_p

end function

!==========================================================================
function ice_air_rh_cct_from_a_si(a_si, t_si, p_si)
!==========================================================================

!RETURNS THE RELATIVE HUMIDITY WITH RESPECT TO ICE FROM AIR MASS FRACTION, TEMPERATURE AND PRESSURE

!OUTPUT: 
!RELATIVE HUMIDITY, CCT/IUPAC DEFINITION

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_RH_CCT_FROM_A_SI(0.998,270,1E5) = 0.680235676087

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_RH_CCT_FROM_A_SI(0.998,270,1E5) = 0.680395740553

real*8 ice_air_rh_cct_from_a_si, a_si, t_si, p_si
real*8 a_sat

ice_air_rh_cct_from_a_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return
if(a_si <= 0d0 .or. a_si > 1d0) return

if(a_si == 1d0) then !DRY AIR
  ice_air_rh_cct_from_a_si = 0d0
  return
end if

a_sat = ice_air_massfraction_air_si(t_si, p_si)
if(a_sat == errorreturn) return
if(a_sat <= 0d0 .or. a_sat > a_si) return

if(a_sat == 1d0) then !DRY AIR
  ice_air_rh_cct_from_a_si = 0d0
  return
end if

ice_air_rh_cct_from_a_si = air_molfraction_vap_si(a_si) &
                         / air_molfraction_vap_si(a_sat)

end function

!==========================================================================
function ice_air_rh_wmo_from_a_si(a_si, t_si, p_si)
!==========================================================================

!RETURNS THE RELATIVE HUMIDITY WITH RESPECT TO ICE FROM AIR MASS FRACTION, TEMPERATURE AND PRESSURE

!OUTPUT: 
!RELATIVE HUMIDITY, WMO DEFINITION

!INPUT:  
!A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!T_SI  ABSOLUTE TEMPERATURE IN K
!P_SI  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ICE_AIR_RH_WMO_FROM_A_SI(0.998,270,1E5) = 0.679205607135

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ICE_AIR_RH_WMO_FROM_A_SI(0.998,270,1E5) = 0.679365943331

real*8 ice_air_rh_wmo_from_a_si, a_si, t_si, p_si
real*8 a_sat

ice_air_rh_wmo_from_a_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return
if(a_si <= 0d0 .or. a_si > 1d0) return

if(a_si == 1d0) then !DRY AIR
  ice_air_rh_wmo_from_a_si = 0
  return
end if

a_sat = ice_air_massfraction_air_si(t_si, p_si)
if(a_sat == errorreturn) return
if(a_sat <= 0d0 .or. a_sat > a_si) return

if(a_sat == 1d0) then !DRY AIR
  ice_air_rh_wmo_from_a_si = 0d0
  return
end if

ice_air_rh_wmo_from_a_si = (1d0 / a_si - 1d0) / (1d0 / a_sat - 1d0)

end function

!==========================================================================
function ice_air_temperature_si()
!==========================================================================

!RETURNS THE TEMPERATURE IN K OF SATURATED AIR WITH RESPECT TO ICE
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_ICE_AIR_EQUILIBRIUM_TO_A_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_temperature_si = 270

!set_ice_air_eq_at_a_t .997, 270
!ice_air_temperature_si = 270

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_temperature_si = 270.232024746

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_temperature_si = 266.511951178

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!set_ice_air_eq_at_t_p 270, 1E5
!ice_air_temperature_si = 270

!set_ice_air_eq_at_a_t .997, 270
!ice_air_temperature_si = 270

!set_ice_air_eq_at_a_p .997, 1E5
!ice_air_temperature_si = 270.234816126

!set_ice_air_eq_at_a_eta .997, 100
!ice_air_temperature_si = 266.514349350

real*8 ice_air_temperature_si

ice_air_temperature_si = errorreturn

if(equi_ice_air_done /= isok) return

ice_air_temperature_si = equi_ice_air_t

end function

!==========================================================================
function set_ice_air_eq_at_a_eta(a_si, eta_si)
!==========================================================================
                                             
!THIS FUNCTION COMPUTES THE SATURATION EQUILIBRIUM STATE BETWEEN HUMID AIR AND
!ICE, IN PARTICULAR, THE ISENTROPIC CONDENSATION TEMPERATURE AND PRESSURE.

!NOTE: THERE IS NO CHECK MADE WHETHER THE ICE PHASE IS OVERHEATED, I.E. THE
!RESULT MAY APPLY TO METASTABLE ICE IF THE TEMPERATURE IS ABOVE THE MELTING POINT.

!INPUT:  A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!        ETA_SI  SPECIFIC ENTROPY OF HUMID AIR IN J/(KG K)

real*8 set_ice_air_eq_at_a_eta, a_si, eta_si
real*8 p_si, t_si, dh, eps

integer maxit

if(equi_ice_air_done == isok .and. &
  a_si == equi_ice_air_a .and. &
  eta_si == equi_ice_air_eta) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_ice_air_eq_at_a_eta = isok
  return
end if

call clear_ice_air_state !WASTE ANY PREVIOUS STATE

set_ice_air_eq_at_a_eta = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return

call init_it_ctrl_ice_air

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0); t_si = aux_ict_si(a_si, eta_si)
  case(1); t_si = ctrl_init_t
  case default; return
end select
if(t_si == errorreturn) return
if(t_si <= 0d0) return

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case(0); p_si = aux_pressure_si(a_si, t_si)
  case(1); p_si = ctrl_init_p
  case default; return
end select
if(p_si == errorreturn) return
if(p_si <= 0d0) return

!SET INITIAL HUMID AIR DENSITY
dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    set_ice_air_eq_at_a_eta = isok
    call set_ice_air_state(a_si, t_si, p_si, eta_si, dh)
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_air_iteration_at_a_eta(a_si, eta_si, maxit, eps, t_si, p_si, dh) &
  == errorreturn) then
  return
end if

call set_ice_air_state(a_si, t_si, p_si, eta_si, dh)

set_ice_air_eq_at_a_eta = isok

end function

!==========================================================================
function set_ice_air_eq_at_a_p(a_si, p_si)
!==========================================================================
                                             
!THIS FUNCTION COMPUTES THE SATURATION EQUILIBRIUM STATE BETWEEN HUMID AIR AND
!ICE, IN PARTICULAR, THE FROSTPOINT TEMPERATURE.

!NOTE: THERE IS NO CHECK MADE WHETHER THE ICE PHASE IS SUPERHEATED, I.E. THE
!RESULT MAY APPLY TO METASTABLE ICE IF THE TEMPERATURE IS ABOVE THE MELTING POINT.

!INPUT:  A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!        P_SI  ABSOLUTE PRESSURE IN PA

real*8 set_ice_air_eq_at_a_p, a_si, p_si
real*8 t_si, dh, eta, eps

integer maxit

if(equi_ice_air_done == isok .and. &
  a_si == equi_ice_air_a .and. &
  p_si == equi_ice_air_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_ice_air_eq_at_a_p = isok
  return
end if

call clear_ice_air_state !WASTE ANY PREVIOUS STATE

set_ice_air_eq_at_a_p = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(p_si <= 0d0) return

call init_it_ctrl_ice_air

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0); t_si = aux_temperature_si(a_si, p_si)
  case(1); t_si = ctrl_init_t
  case default; return
end select
if(t_si == errorreturn) return
if(t_si <= 0d0) return

!SET INITIAL HUMID AIR DENSITY
dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum  == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    set_ice_air_eq_at_a_p = isok
    eta = air_f_entropy_si(a_si, t_si, dh)
    call set_ice_air_state(a_si, t_si, p_si, eta, dh)
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_air_iteration_at_a_p(a_si, p_si, maxit, eps, t_si, dh) == errorreturn) then
  return
end if

eta = air_f_entropy_si(a_si, t_si, dh)
call set_ice_air_state(a_si, t_si, p_si, eta, dh)

set_ice_air_eq_at_a_p = isok

end function

!==========================================================================
function set_ice_air_eq_at_a_t(a_si, t_si)
!==========================================================================
                                             
!THIS FUNCTION COMPUTES THE SATURATION EQUILIBRIUM STATE BETWEEN HUMID AIR AND
!ICE, IN PARTICULAR, THE FROSTPOINT CONDENSATION PRESSURE.

!NOTE: THERE IS NO CHECK MADE WHETHER THE ICE PHASE IS SUPERHEATED, I.E. THE
!RESULT MAY APPLY TO METASTABLE ICE IF THE TEMPERATURE IS ABOVE THE MELTING POINT.

!INPUT:  A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!        T_SI  ABSOLUTE TEMPERATURE IN K

real*8 set_ice_air_eq_at_a_t, a_si, t_si
real*8 p_si, dh, eta, eps

integer maxit

if(equi_ice_air_done == isok .and. &
  a_si == equi_ice_air_a .and. &
  t_si == equi_ice_air_t) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMP~TED EARLIER
  set_ice_air_eq_at_a_t = isok
  return
end if

call clear_ice_air_state !WASTE ANY PREVIOUS STATE

set_ice_air_eq_at_a_t = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si <= 0d0) return

call init_it_ctrl_ice_air

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case(0); p_si = aux_temperature_si(a_si, t_si)
  case(1); p_si = ctrl_init_p
  case default; return
end select
if(p_si == errorreturn) return
if(p_si <= 0d0) return

!SET INITIAL HUMID AIR DENSITY
dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    set_ice_air_eq_at_a_t = isok
    eta = air_f_entropy_si(a_si, t_si, dh)
    call set_ice_air_state(a_si, t_si, p_si, eta, dh)
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_air_iteration_at_a_t(a_si, t_si, maxit, eps, p_si, dh) == errorreturn) then
  return
end if

eta = air_f_entropy_si(a_si, t_si, dh)
call set_ice_air_state(a_si, t_si, p_si, eta, dh)

set_ice_air_eq_at_a_t = isok

end function

!==========================================================================
function set_ice_air_eq_at_t_p( t_si, p_si)
!==========================================================================
                                             
!THIS FUNCTION COMPUTES THE SATURATION EQUILIBRIUM STATE BETWEEN HUMID AIR AND
!ICE, IN PARTICULAR, THE SATURATED SPECIFIC HUMIDITY OF AIR.

!NOTE: THERE IS NO CHECK MADE WHETHER THE ICE PHASE IS SUPERHEATED, I.E. THE
!RESULT MAY APPLY TO METASTABLE ICE IF THE TEMPERATURE IS ABOVE THE MELTING POINT.

!INPUT:  T_SI  ABSOLUTE TEMPERATURE IN K
!        P_SI  ABSOLUTE PRESSURE IN PA

real*8 set_ice_air_eq_at_t_p, t_si, p_si
real*8 a_si, dh, eta, eps

integer maxit

if(equi_ice_air_done == isok .and. &
  t_si == equi_ice_air_t .and. &
  p_si == equi_ice_air_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_ice_air_eq_at_t_p = isok
  return
end if

call clear_ice_air_state !WASTE ANY PREVIOUS STATE

set_ice_air_eq_at_t_p = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_it_ctrl_ice_air

!SET INITIAL AIR FRACTION GUESS
select case(ctrl_mode_humidity)
  case(0)  
    a_si = aux_airfraction_si(t_si, p_si)
    if(a_si == errorreturn) return
    if(a_si < 0d0 .or. a_si > 1d0) return
  case(1)  
    a_si = ctrl_init_a
  case default; return
end select
if(a_si == errorreturn) return
if(a_si < 0d0 .or. a_si > 1d0) return

dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    set_ice_air_eq_at_t_p = isok
    eta = air_f_entropy_si(a_si, t_si, dh)
    call set_ice_air_state(a_si, t_si, p_si, eta, dh)
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_air_iteration_at_t_p(t_si, p_si, maxit, eps, a_si, dh) &
  == errorreturn) then
  return
end if

eta = air_f_entropy_si(a_si, t_si, dh)
call set_ice_air_state(a_si, t_si, p_si, eta, dh)

set_ice_air_eq_at_t_p = isok

end function

!==========================================================================
function ice_air_iteration_at_t_p( t_si, p_si, maxit, eps, a_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURE,
!TEMPERATURE AND CHEMICAL POTENTIAL OF WATER OF The TWO PHASES AT GIVEN PRESSURE, P_SI,
!AND  TEMPERATURE, T_SI, FROM INITIAL GUESSES FOR THE HUMID-AIR DENSITY, D_AIR_SI,
!AND THE AIR FRACTION, A_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PARTIAL PRESSURE

!OUTPUT: 
!        ICE_AIR_ITERATION_AT_T_P = ISOK IF SUCCESSFUL
!        ICE_AIR_ITERATION_AT_T_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO ICE_G_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - AIR FRACTION IS NOT BETWEEN 0 AND 1 DURING THE ITERATION
!     A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
! D_AIR_SI: HUMID-AIR DENSITY AT THE SATURATION EQUILIBRIUM IN KG/M3

!INPUT: 
!       P_SI: ABSOLUTE PRESSURE IN PA
!       T_SI: ABSOLUTE TEMPERATURE IN K
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       A_SI: INITIAL GUESS OF THE MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY IN KG/M3

real*8 ice_air_iteration_at_t_p, t_si, p_si, eps, a_si, d_air_si

real*8 ah, fh_a, fh_aa, fh_ad
real*8 gi
real*8 fh, fh_d, fh_dd
real*8 dh
real*8 pv, pv_old

real*8 a(2, 2), b(2), x(2)

integer it, maxit

ice_air_iteration_at_t_p = errorreturn
pv = errorreturn

!CHECK FOR INVALID INPUT VALUES
if(a_si < 0d0 .or. a_si > 1d0 .or. &
  d_air_si <= 0d0 .or. &
  t_si <= 0 .or. &
  p_si <= 0d0) then
  a_si = errorreturn
  d_air_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

ah = a_si
dh = d_air_si

!CHEMICAL POTENTIAL OF ICE
gi = ice_g_si(0, 0, t_si, p_si)
if(gi == errorreturn) go to 999

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 2D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, ah, t_si, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, ah, t_si, dh)
  if(fh_a == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, ah, t_si, dh)
  if(fh_d == errorreturn) go to 999
  fh_aa = air_f_si(2, 0, 0, ah, t_si, dh)
  if(fh_aa == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, ah, t_si, dh)
  if(fh_ad == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, ah, t_si, dh)
  if(fh_dd == errorreturn) go to 999
  
  !VAPOUR PARTIAL PRESSURE FOR ERROR CHECKING
  pv_old = pv
  pv = air_molfraction_vap_si(ah)
  if(pv == errorreturn) go to 999
  pv = p_si * pv
  
  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or. &
       (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      ice_air_iteration_at_t_p = isok
      d_air_si = dh
      a_si = ah
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = -ah * fh_aa
  a(1, 2) = fh_d - ah * fh_ad - p_si / dh**2
  a(2, 1) = dh * fh_ad
  a(2, 2) = 2d0 * fh_d + dh * fh_dd

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gi - fh - p_si / dh + ah * fh_a
  b(2) = p_si / dh - dh * fh_d
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR

  !UPDATE AIR FRACTION & DENSITY
  ah = ah + x(1)
  if(ah < 0d0 .or. ah > 1d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !AIR_LIMITS
  if(a_si < 0d0 .or. a_si > 1d0 .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_t_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_air_iteration_at_a_eta(a_si, eta_si, maxit, eps, t_si, p_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURE,
!TEMPERATURE AND CHEMICAL POTENTIAL OF WATER OF THE TWO PHASES AT GIVEN ENTROPY, ETA_SI,
!AND AIR FRACTION, A_SI, FROM INITIAL GUESSES FOR THE HUMID-AIR DENSITY, D_AIR_SI,
!THE TEMPERATURE, T_SI, AND THE PRESSURE, P_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN THE PRESSURE

!OUTPUT: 
!        ICE_AIR_ITERATION_AT_A_ETA = ISOK IF SUCCESSFUL
!        ICE_AIR_ITERATION_AT_A_ETA = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO ICE_G_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - PRESSURE, DENSITY OR TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
!     P_SI: ABSOLUTE PRESSURE IN PA
! D_AIR_SI: HUMID-AIR DENSITY AT THE SATURATION EQUILIBRIUM IN KG/M3

!INPUT: 
!     ETA_SI: SPECIFIC ENTROPY IN J/(KG K)
!       A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!       P_SI: INITIAL GUESS OF ABSOLUTE PRESSURE IN PA
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY IN KG/M3

real*8 ice_air_iteration_at_a_eta, a_si, eta_si, eps, t_si, p_si, d_air_si

real*8 fh, fh_a, fh_d
real*8 fh_ad, fh_dd
real*8 fh_t, fh_at, fh_td, fh_tt
real*8 gi, gi_t, gi_p
real*8 dh, p, t

real*8 a(3, 3), b(3), x(3)

integer it, maxit

ice_air_iteration_at_a_eta = errorreturn

!CHECK FOR INVALID INPUT VALUES
if(a_si < 0d0 .or. a_si >= 1d0 .or. &
  d_air_si <= 0d0 .or. &
  t_si <= 0d0 .or. &
  p_si <= 0d0) then
  t_si = errorreturn
  p_si = errorreturn
  d_air_si = errorreturn
  return
endif

check_limits = check_limits - 1

t = t_si
p = p_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 3D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, a_si, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a_si, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a_si, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a_si, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a_si, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a_si, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_tt = air_f_si(0, 2, 0, a_si, t, dh)
  if(fh_tt == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a_si, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a_si, t, dh)
  if(fh_dd == errorreturn) go to 999

  !ICE:
  gi = ice_g_si(0, 0, t, p)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p)
  if(gi_t == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t, p)
  if(gi_p == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = 1d0 / dh - gi_p
  a(1, 2) = fh_d - a_si * fh_ad - p / dh**2
  a(1, 3) = fh_t - a_si * fh_at - gi_t
  a(2, 1) = -1d0 / dh
  a(2, 2) = 2d0 * fh_d + dh * fh_dd
  a(2, 3) = dh * fh_td
  a(3, 1) = 0d0
  a(3, 2) = fh_td
  a(3, 3) = fh_tt

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gi - fh + a_si * fh_a - p / dh
  b(2) = p / dh - dh * fh_d
  b(3) = -eta_si - fh_t
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 3) /= 0) go to 999 !MATRIX SINGULAR

  !UPDATE PRESSURE, HUMID-AIR DENSITY AND TEMPERATURE
  p = p + x(1)
  if(p <= 0d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999
  t = t + x(3)
  if(t <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(1)) < eps) .or. &
     (eps < 0d0 .and. abs(x(1)) < -eps * p)) then
    ice_air_iteration_at_a_eta = isok
    d_air_si = dh
    t_si = t
    p_si = p
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    p_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_eta = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    t_si = errorreturn
    p_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_eta = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_air_iteration_at_a_t(a_si, t_si, maxit, eps, p_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURE,
!TEMPERATURE AND CHEMICAL POTENTIAL OF WATER OF THE TWO PHASES AT GIVEN TEMPERATURE, T_SI,
!AND AIR FRACTION, A_SI, FROM INITIAL GUESSES FOR THE HUMID-AIR DENSITY, D_AIR_SI,
!AND THE PRESSURE, P_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN THE PRESSURE

!OUTPUT: 
!        ICE_AIR_ITERATION_AT_A_T = ISOK IF SUCCESSFUL
!        ICE_AIR_ITERATION_AT_A_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO ICE_G_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - PRESSURE OR TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     P_SI: ABSOLUTE PRESSURE IN PA
! D_AIR_SI: HUMID-AIR DENSITY AT THE SATURATION EQUILIBRIUM IN KG/M3

!INPUT: 
!       T_SI: ABSOLUTE TEMPERATURE IN K
!       A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       P_SI: INITIAL GUESS OF ABSOLUTE PRESSURE IN PA
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY IN KG/M3


real*8 ice_air_iteration_at_a_t, a_si, t_si, eps, p_si, d_air_si
real*8  fh, fh_a, fh_d
real*8 fh_ad, fh_dd
real*8 gi, gi_p
real*8 dh, p

real*8  a(2, 2), b(2), x(2)

integer it, maxit

ice_air_iteration_at_a_t = errorreturn

!CHECK FOR INVALID INPUT VALUES
if(d_air_si <= 0d0 .or. &
  p_si <= 0d0) then
  p_si = errorreturn
  d_air_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    p_si = errorreturn
    d_air_si = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  a_si < 0d0 .or. a_si > 1d0) then
    p_si = errorreturn
    d_air_si = errorreturn
    return
  endif
else
  if(a_si < 0d0 .or. a_si >= 1d0 .or. &
    t_si <= 0d0) then
    p_si = errorreturn
    d_air_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

p = p_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 2D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, a_si, t_si, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a_si, t_si, dh)
  if(fh_a == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a_si, t_si, dh)
  if(fh_d == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a_si, t_si, dh)
  if(fh_ad == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a_si, t_si, dh)
  if(fh_dd == errorreturn) go to 999

  !ICE:
  gi = ice_g_si(0, 0, t_si, p)
  if(gi == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t_si, p)
  if(gi_p == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = 1d0 / dh - gi_p
  a(1, 2) = fh_d - a_si * fh_ad - p / dh**2
  a(2, 1) = -1d0 / dh
  a(2, 2) = 2d0 * fh_d + dh * fh_dd

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gi - fh + a_si * fh_a - p / dh
  b(2) = p / dh - dh * fh_d
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR

  !UPDATE PRESSURE & DENSITY
  p = p + x(1)
  if(p <= 0d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(1)) < eps) .or. &
     (eps < 0d0 .and. abs(x(1)) < -eps * p)) then
    ice_air_iteration_at_a_t = isok
    d_air_si = dh
    p_si = p
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    p_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_t = errorreturn
    return
  endif
  !AIR_LIMITS
  if(d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    p_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_air_iteration_at_a_p(a_si, p_si, maxit, eps, t_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURE, 
!TEMPERATURE AND CHEMICAL POTENTIAL OF WATER OF THE TWO PHASES AT GIVEN PRESSURE, P_SI,
!AND AIR FRACTION, A_SI, FROM INITIAL GUESSES FOR THE THE HUMID-AIR DENSITY, D_AIR_SI,
!AND THE TEMPERATURE, T_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN THE PRESSURE

!OUTPUT: 
!        ICE_AIR_ITERATION_AT_A_P = ISOK IF SUCCESSFUL
!        ICE_AIR_ITERATION_AT_A_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO ICE_G_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - DENSITY OR TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
! D_AIR_SI: HUMID-AIR DENSITY AT THE SATURATION EQUILIBRIUM IN KG/M3

!INPUT: 
!       P_SI: ABSOLUTE PRESSURE IN PA
!       A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY IN KG/M3


real*8 ice_air_iteration_at_a_p, a_si, p_si, eps, t_si, d_air_si

real*8 fh, fh_a, fh_t, fh_d
real*8 fh_ad, fh_at, fh_td, fh_dd
real*8 gi, gi_t
real*8 dh, t

real*8 a(2, 2), b(2), x(2)

integer it, maxit

ice_air_iteration_at_a_p = errorreturn

!CHECK FOR INVALID INPUT VALUES
if(a_si < 0d0 .or. a_si >= 1d0 .or. &
  d_air_si <= 0d0 .or. &
  t_si <= 0d0) then
  t_si = errorreturn
  d_air_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    d_air_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0) then
    t_si = errorreturn
    d_air_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
dh = d_air_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 2D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, a_si, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a_si, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a_si, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a_si, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a_si, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a_si, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a_si, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a_si, t, dh)
  if(fh_dd == errorreturn) go to 999

  !ICE:
  gi = ice_g_si(0, 0, t, p_si)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p_si)
  if(gi_t == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = fh_t - a_si * fh_at - gi_t
  a(1, 2) = fh_d - a_si * fh_ad - p_si / dh**2
  a(2, 1) = dh * fh_td
  a(2, 2) = 2d0 * fh_d + dh * fh_dd

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gi - fh + a_si * fh_a - p_si / dh
  b(2) = p_si / dh - dh * fh_d
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR

  !UPDATE TEMPERATURE & DENSITY
  t = t + x(1)
  if(t <= 0d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(1)) < eps) .or. &
     (eps < 0d0 .and. abs(x(1)) < -eps * t)) then
    ice_air_iteration_at_a_p = isok
    d_air_si = dh
    t_si = t
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_air_si <= flu_dmin .or. d_air_si > flu_dmax) then
    t_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_p = errorreturn
    return
  endif
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    t_si = errorreturn
    d_air_si = errorreturn
    ice_air_iteration_at_a_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_airfraction_si(t_si, p_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE AIR MASS FRACTION OF SATURATED
!HUMID AIR AT GIVEN ABSOLUTE TEMPERATURE, T_SI, IN K AND ABSOLUTE
!PRESSURE, P_SI, IN PA

real*8 aux_airfraction_si, t_si, p_si
real*8 dh, dv, pv

aux_airfraction_si = errorreturn
if(t_si <= 0d0) return
if(p_si <= 0d0) return

!AS AN ESTIMATE, ADD THE DENSITIES OF DRY AIR AND VAPOUR

!VAPOUR PRESSURE
pv = ice_vap_sublimationpressure_si(t_si)
if(pv == errorreturn) return
if(pv > p_si) return

if(pv == p_si) then
  aux_airfraction_si = 0d0
  return
end if

!VAPOUR DENSITY
dv = ice_vap_density_vap_si()
if(dv == errorreturn) return
if(dv <= 0d0) return

!DRY AIR DENSITY AT REMAINING PARTIAL PRESSURE
dh = air_density_si(1d0, t_si, p_si - pv)
if(dh == errorreturn) return
if(dh <= 0d0) return

!AIR FRACTION OF THE IDEAL MIXTURE
aux_airfraction_si = dh / (dh + dv)
           
end function

!==========================================================================
function aux_ict_si( a_si, eta_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE ISENTROPIC CONDENSATION TEMPERATURE
!IN K OF HUMID AIR AT GIVEN AIR MASS FRACTION, A_SI, IN KG/KG AND ENTROPY
!ETA_SI, IN J/(KG K)

real*8 aux_ict_si, a_si, eta_si
real*8 tt, pt, ra, rw, cpa, cpv, l

real*8 etat, rav, xv
real*8 numer, denom

tt = tp_temperature_si
pt = tp_pressure_iapws95_si

ra = gas_constant_air_si  !SPECIFIC GAS CONSTANT OF AIR IN J/(KG K)
rw = gas_constant_h2o_si  !SPECIFIC GAS CONSTANT OF H2O IN J/(KG K)

!PROPERTIES AT THE TRIPLE POINT
cpa = 1003.69d0  !HEAT CAPACITY OF AIR, IN J/(KG K)
cpv = 1884.352d0 !HEAT CAPACITY OF VAPOUR, IN J/(KG K)
l = tp_enthalpy_vap_si - tp_enthalpy_ice_si    !LATENT HEAT

aux_ict_si = errorreturn
if(a_si < 0d0 .or. a_si >= 1d0) return

!TO ESTIMATE T FROM ETA AND A, USE THE CLAUSIUS-CLAPEYRON EQUATION PVAP(T)
!FOR THE PARTIAL SUBLIMATION PRESSURE TO ELIMINATE PVAP FROM THE IDEAL-GAS
!EQUATION ETA(A,T,P) = CONST, AND SOLVE IT FOR T

etat = air_g_entropy_si(a_si, tt, pt)  !ENTROPY OF HUMID AIR AT THE TRIPLE POINT
if(etat == errorreturn) return
rav = a_si * ra + (1d0 - a_si) * rw    !GAS CONSTANT OF HUMID AIR
xv = air_molfraction_vap_si(a_si)
if(xv == errorreturn) return
if(xv <= 0d0) return

denom = a_si * (cpa - ra * l / (rw * tt)) + (1d0 - a_si) * (cpv - l / tt)
if(denom == 0d0) return

numer = eta_si - etat - rav * log(xv)

aux_ict_si = tt * exp(numer / denom)

end function

!==========================================================================
function aux_temperature_si(a_si, p_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE TEMPERATURE IN K OF SATURATED
!HUMID AIR AT GIVEN AIR MASS FRACTION, A_SI, IN KG/KG AND ABSOLUTE
!PRESSURE, P_SI, IN PA

real*8 aux_temperature_si, a_si, p_si
real*8 pv, t

aux_temperature_si = errorreturn
if(a_si < 0d0 .or. a_si >= 1d0) return
if(p_si <= 0d0) return

!AS AN ESTIMATE, IGNORE THE AIR-VAPOUR INTERACTION

!PARTIAL VAPOUR PRESSURE
pv = p_si * air_molfraction_vap_si(a_si)
!RELATED SATURATION TEMPERATURE
t = ice_vap_sublimationtemp_si(pv)
if(t == errorreturn) return
if(t <= 0d0) return

aux_temperature_si = t

end function

!==========================================================================
function aux_pressure_si(a_si, t_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE PRESSURE IN PA OF SATURATED
!HUMID AIR AT GIVEN AIR MASS FRACTION, A_SI, IN KG/KG AND ABSOLUTE
!TEMPERATURE, T_SI, IN K

real*8 aux_pressure_si, a_si, t_si
real*8 pv, pa, d

aux_pressure_si = errorreturn
if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si <= 0d0) return

!AS AN ESTIMATE, IGNORE THE AIR-VAPOUR INTERACTION

!PARTIAL VAPOUR PRESSURE
pv = ice_vap_sublimationpressure_si(t_si)
if(pv == errorreturn) return
if(pv <= 0d0) return
d = ice_vap_density_vap_si()
if(d == errorreturn) return
if(d <= 0d0) return  !VAPOUR DENSITY

!PARTIAL AIR PRESSURE
d = a_si / (1d0 - a_si) * d !AIR DENSITY
pa = air_f_pressure_si(1d0, t_si, d)
if(pa == errorreturn) return
if(pa <= 0d0) return

aux_pressure_si = pa + pv

end function

!==========================================================================
subroutine set_ice_air_state(a, t, p, eta, dh)

real*8 a, t, p, eta, dh

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR
equi_ice_air_done = isok

equi_ice_air_a = a
equi_ice_air_t = t
equi_ice_air_p = p
equi_ice_air_eta = eta
equi_ice_air_d_air = dh

end subroutine

!==========================================================================
subroutine clear_ice_air_state()

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_ice_air_done = 0

end subroutine

!==========================================================================
subroutine init_it_ctrl_ice_air()

real*8 t_triple, p_triple

!TRIPLE POINT VALUES
t_triple = tp_temperature_si
p_triple = tp_pressure_iapws95_si

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR THE ITERATION
ctrl_mode_pressure = 0
ctrl_mode_temperature = 0
ctrl_mode_humidity = 0

ctrl_loop_maximum = 100
ctrl_init_t = t_triple
ctrl_init_p = p_triple
ctrl_init_a = 1d0 - 0.001d0

ctrl_eps_exit_p_vap = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_temp = -0.0000001d0  !RELATIVE, 0.1 PPM
ctrl_eps_exit_press = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

!==========================================================================
subroutine set_it_ctrl_ice_air(key, value)

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!PROPERTIES OF HUMID AIR SATURATED WITH RESPECT TO ICE

!KEY              VALUE
!IT_STEPS         0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS         N > 0    SET ITERATION NUMBER TO N
!IT_STEPS        -1           DO NOT ITERATE, USE INITIAL VALUES AS THE SATURATION STATE

!INIT_AIR         0           USE DEFAULT AIR FRACTION TO START ( = AUX_AIRFRACTION_SI(T, P))
!INIT_AIR         A > 0    USE VALUE A AS AIR FRACTION TO START

!INIT_TEMP        0           USE DEFAULT TEMPERATURE TO START ( = AUX_TEMPERATURE_SI(A, P))
!INIT_TEMP        T > 0    USE VALUE T AS TEMPERATURE TO START

!INIT_PRESS       0           USE DEFAULT PRESSURE TO START ( = AUX_PRESSURE_SI(A, T))
!INIT_PRESS       P > 0    USE VALUE P AS PRESSURE TO START

!TOL_VAP_PRESS    0           USE DEFAULT EXIT ACCURACY FOR VAPOUR PRESSURE (0.1 PPM)
!TOL_VAP_PRESS    EPS         USE EPS AS EXIT ACCURACY FOR VAPOUR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_TEMP         0           USE DEFAULT EXIT ACCURACY FOR TEMPERATURE (0.1 PPM)
!TOL_TEMP         EPS         USE EPS AS EXIT ACCURACY FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_PRESS        0           USE DEFAULT EXIT ACCURACY FOR PRESSURE (0.1 PPM)
!TOL_PRESS        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

real*8 value
character*13 key

call init_it_ctrl_ice_air

call clear_ice_air_state

select case (trim(key))

  case ("it_steps");   !ITERATION STEPS
    if(value == 0d0) then
      ctrl_loop_maximum = 100  !DEFAULT = 100
    elseif(value < 0d0) then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ("init_air");   !START AIR FRACTION
    if(value == 0d0) then
      ctrl_mode_humidity = 0    !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_humidity = 1
      ctrl_init_a = value
    endif
    
  case ("init_temp");   !START TEMPERATURE
    if(value == 0d0) then
      ctrl_mode_temperature = 0    !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif
    
  case ("init_press");   !START PRESSURE
    if(value == 0d0) then
      ctrl_mode_pressure = 0    !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_pressure = 1
      ctrl_init_p = value
    endif

  case ("tol_vap_press");      !REQUIRED VAPOUR PRESSURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_p_vap = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p_vap = value
    endif

  case ("tol_temp");      !REQUIRED TEMPERATURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_temp = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_temp = value
    endif

  case ("tol_press");      !REQUIRED TOTAL PRESSURE TOLERANCE
    if(value == 0d0)then
      ctrl_eps_exit_press = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_press = value
    endif

  case default; write(*,*) key,' not supported'

end select

end subroutine

end module ice_air_4a

