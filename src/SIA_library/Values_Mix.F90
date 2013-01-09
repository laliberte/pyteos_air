module values_mix

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE MULTI-PHASE ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90
!     AIR_1,         FILE AIR_1.F90
!     AIR_2,         FILE AIR_2.F90
!     FLU_1,         FILE FLU_1.F90
!     FLU_2,         FILE FLU_2.F90
!     FLU_3A,        FILE FLU_3A.F90
!     FLU_3B,        FILE FLU_3B.F90
!     ICE_1,         FILE ICE_1.F90
!     ICE_2,         FILE ICE_2.F90
!     SAL_1,         FILE ICE_1.F90
!     SAL_2,         FILE ICE_2.F90
!     SEA_3A,        FILE SEA_3A.F90
!     SEA_3B,        FILE SEA_3B.F90
!     SEA_3C,        FILE SEA_3C.F90
!     SEA_3D,        FILE SEA_3D.F90
!     ICE_AIR_4A,    FILE ICE_AIR_4A.F90
!     ICE_AIR_4B,    FILE ICE_AIR_4B.F90
!     ICE_AIR_4C,    FILE ICE_AIR_4C.F90
!     ICE_LIQ_4,     FILE ICE_LIQ_4.F90
!     ICE_VAP_4,     FILE ICE_VAP_4.F90
!     LIQ_AIR_4A,    FILE LIQ_AIR_4A.F90
!     LIQ_AIR_4B,    FILE LIQ_AIR_4B.F90
!     LIQ_AIR_4C,    FILE LIQ_AIR_4B.F90
!     LIQ_VAP_4,     FILE LIQ_VAP_4.F90
!     LIQ_ICE_AIR_4, FILE LIQ_ICE_AIR_4.F90
!     SEA_AIR_4,     FILE SEA_AIR_4.F90
!     SEA_ICE_4,     FILE SEA_ICE_4.F90
!     SEA_ICE_VAP_4, FILE SEA_ICE_VAP_4.F90
!     SEA_LIQ_4,     FILE SEA_LIQ_4.F90
!     SEA_VAP_4,     FILE SEA_VAP_4.F90
!     ICE_FLU_5,     FILE FIT_ICE_FLU_5.F90

use constants_0
use maths_0
use convert_0
use air_1
use flu_1
use ice_1
use sal_1
use air_2
use flu_2
use ice_2
use sal_2
use flu_3a
use flu_3b
use sea_3a
use sea_3b
use sea_3c
use sea_3d
use ice_liq_4
use ice_vap_4
use liq_vap_4
use sea_air_4
use sea_ice_4
use sea_liq_4
use sea_vap_4
use liq_ice_air_4
use sea_ice_vap_4
use ice_air_4a
use ice_air_4b
use ice_air_4c
use liq_air_4a
use liq_air_4b
use liq_air_4c
use ice_flu_5

implicit none 
character*16, private :: version = '8 Jun 2010'
character (len = *), parameter :: fmt = "(A,T38,1X,E19.12,2X,A)"

contains

!*******************************************************************
subroutine mix_chk_values
!*******************************************************************

implicit none 

real*8 temp

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF EQUILIBRIUM STATES'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

if(molar_mass_air_L2000 == molar_mass_air_si)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.1 (V1.1)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ELSEIF) USES THIS VALUE OF THE MOLAR MASS.
!TO USES THE ORIGINAL VALUE USED BY LEMMON ET AL. (2000) YOU MUST EDIT CONSTANTS_0 TO 
!CORRESPOND TO VERSION 1.0 (V1.0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4A'

write(*,fmt) 'ice_air_massfraction_air(270,1d5) =', ice_air_massfraction_air_si(270d0,1d5), & 
' 0.997058854720'
write(*,fmt) 'ice_air_frostpoint(0.997,1d5) =', ice_air_frostpoint_si(0.997d0,1d5), & 
' 270.234816126'
write(*,fmt) 'ice_air_ict(0.997,300,1d5) =', ice_air_ict_si(0.997d0,300d0,1d5), & 
' 265.224998411'
write(*,fmt) 'ice_air_icl(0.997,300,1d5) =', ice_air_icl_si(0.997d0,300d0,1d5), & 
' 64988.3931838'
write(*,fmt) 'ice_air_condenspressure(0.997,270) =', ice_air_condensationpressure_si(0.997d0,270d0), & 
' 98034.4511233'
write(*,fmt) 'ice_air_sublimationpress(270,1E5) =', ice_air_sublimationpressure_si(270d0,1d5), & 
' 472.041823975'
write(*,fmt) 'ice_air_rh_cct_from_a(0.998,270,1d5) =', ice_air_rh_cct_from_a_si(0.998d0,270d0,1d5), & 
' 0.680395740553'
write(*,fmt) 'ice_air_rh_wmo_from_a(0.998,270,1d5) =', ice_air_rh_wmo_from_a_si(0.998d0,270d0,1d5), & 
' 0.679365943331'
write(*,fmt) 'ice_air_a_from_rh_wmo(0.8,270,1d5) =', ice_air_a_from_rh_wmo_si(0.8d0,270d0,1d5), & 
' 0.997645698908'
write(*,fmt) 'ice_air_a_from_rh_cct(0.8,270,1d5) =', ice_air_a_from_rh_cct_si(0.8d0,270d0,1d5), & 
' 0.997647924743'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833359.27614'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833386.54980'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833296.51317'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2834612.42351'

write(*,"(A)")' '
temp =set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.28880078014'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.26340801028'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.28763121402'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 0.949325026119'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.79055033080E-03'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.79022403085E-03'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.86289364206E-03'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 2.84797507836E-03'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 98034.4511233'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 72721.4579415'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270'

write(*,"(A)")' '
temp =set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270.234816126'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 266.514349350'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.181167192'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.180955861'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.147060527'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.681749114'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4B'

write(*,fmt)'ice_air_g(0,0,0,0.5,270,1d5) =', ice_air_g_si( 0, 0, 0, 0.5d0, 270d0, 1d5), &
'-2595.57166634'
write(*,fmt)'ice_air_g(1,0,0,0.5,270,1d5) =', ice_air_g_si( 1, 0, 0, 0.5d0, 270d0, 1d5), &
' 2382.35592988'
write(*,fmt)'ice_air_g(0,1,0,0.5,270,1d5) =', ice_air_g_si( 0, 1, 0, 0.5d0, 270d0, 1d5), &
' 610.264515318'
write(*,fmt)'ice_air_g(0,0,1,0.5,270,1d5) =', ice_air_g_si( 0, 0, 1, 0.5d0, 270d0, 1d5), &
' 0.389645501224'
write(*,fmt)'ice_air_g(2,0,0,0.5,270,1d5) =', ice_air_g_si( 2, 0, 0, 0.5d0, 270d0, 1d5), &
' 0'
write(*,fmt)'ice_air_g(1,1,0,0.5,270,1d5) =', ice_air_g_si( 1, 1, 0, 0.5d0, 270d0, 1d5), &
'-1269.41767949'
write(*,fmt)'ice_air_g(1,0,1,0.5,270,1d5) =', ice_air_g_si( 1, 0, 1, 0.5d0, 270d0, 1d5), &
' 0.777110408175'
write(*,fmt)'ice_air_g(0,2,0,0.5,270,1d5) =', ice_air_g_si( 0, 2, 0, 0.5d0, 270d0, 1d5), &
'-7.00810930740'
write(*,fmt)'ice_air_g(0,1,1,0.5,270,1d5) =', ice_air_g_si( 0, 1, 1, 0.5d0, 270d0, 1d5), &
' 1.60095965101E-03'
write(*,fmt)'ice_air_g(0,0,2,0.5,270,1d5) =', ice_air_g_si( 0, 0, 2, 0.5d0, 270d0, 1d5), &
'-3.91178603885E-06'

write(*,fmt)' '
write(*,fmt)'ice_air_solidfrac(0.5,270,1d5) =', ice_air_solidfraction_si(0.5d0, 270d0, 1d5), & 
' 0.498525089434'
write(*,fmt)'ice_air_vaporfrac(0.5,270,1d5) =', ice_air_vapourfraction_si(0.5d0, 270d0, 1d5), & 
' 1.47491056602E-03'
write(*,fmt)'ice_air_g_cp(0.5,270,1d5) =', ice_air_g_cp_si(0.5d0, 270d0, 1d5), & 
' 1892.18951300'
write(*,fmt)'ice_air_g_density(0.5,270,1d5) =', ice_air_g_density_si(0.5d0, 270d0, 1d5), & 
' 2.56643538000'
write(*,fmt)'ice_air_g_enthalpy(0.5,270,1d5) =', ice_air_g_enthalpy_si(0.5d0, 270d0, 1d5), & 
'-167366.990802'
write(*,fmt)'ice_air_g_entropy(0.5,270,1d5) =', ice_air_g_entropy_si(0.5d0, 270d0, 1d5), & 
'-610.264515318'
write(*,fmt)'ice_air_g_exp(0.5,270,1d5) =', ice_air_g_expansion_si(0.5d0, 270d0, 1d5), & 
' 4.10875949031E-03'
write(*,fmt)'ice_air_g_kap_t(0.5,270,1d5) =', ice_air_g_kappa_t_si(0.5d0, 270d0, 1d5), & 
' 1.00393460891E-05'
write(*,fmt)'ice_air_g_lapserate(0.5,270,1d5) =', ice_air_g_lapserate_si(0.5d0, 270d0, 1d5), & 
' 2.28443875629E-04'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4C'

write(*,fmt)'ice_air_h(0,0,0,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 0, 0.5d0, -600d0, 1d5), &
'-164588.106002'
write(*,fmt)'ice_air_h(1,0,0,0.5,-600,1d5) =', ice_air_h_si( 1, 0, 0, 0.5d0, -600d0, 1d5), &
' 543.016638396'
write(*,fmt)'ice_air_h(0,1,0,0.5,-600,1d5) =', ice_air_h_si( 0, 1, 0, 0.5d0, -600d0, 1d5), &
' 271.449994437'
write(*,fmt)'ice_air_h(0,0,1,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 1, 0.5d0, -600d0, 1d5), &
' 0.391981878510'
write(*,fmt)'ice_air_h(2,0,0,0.5,-600,1d5) =', ice_air_h_si( 2, 0, 0, 0.5d0, -600d0, 1d5), &
' 224958.525864'
write(*,fmt)'ice_air_h(1,1,0,0.5,-600,1d5) =', ice_air_h_si( 1, 1, 0, 0.5d0, -600d0, 1d5), &
'-177.457078495'
write(*,fmt)'ice_air_h(1,0,1,0.5,-600,1d5) =', ice_air_h_si( 1, 0, 1, 0.5d0, -600d0, 1d5), &
' 0.781782661019'
write(*,fmt)'ice_air_h(0,2,0,0.5,-600,1d5) =', ice_air_h_si( 0, 2, 0, 0.5d0, -600d0, 1d5), &
' 0.139985868894'
write(*,fmt)'ice_air_h(0,1,1,0.5,-600,1d5) =', ice_air_h_si( 0, 1, 1, 0.5d0, -600d0, 1d5), &
' 2.26912930199E-04'
write(*,fmt)'ice_air_h(0,0,2,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 2, 0.5d0, -600d0, 1d5), &
'-3.56976697603E-06'

write(*,"(A)")' '
write(*,fmt)'ice_air_pottemp(0.9,230,5d4,1d5) =', ice_air_pottemp_si(0.9d0, 230d0, 5d4, 1d5), & 
' 266.105208871'
write(*,fmt)'ice_air_potdens(0.9,230,5d4,1d5) =', ice_air_potdensity_si(0.9d0, 230d0, 5d4, 1d5), & 
' 1.45048110422'
write(*,fmt)'ice_air_potenth(0.9,230,5d4,1d5) =', ice_air_potenthalpy_si(0.9d0, 230d0, 5d4, 1d5), & 
'-35781.2564451'
write(*,fmt)'ice_air_h_temp(0.9,-100,1d5) =', ice_air_h_temperature_si(0.9d0, -100d0, 1d5), & 
' 270.383680119'
write(*,fmt)'ice_air_h_lapsrate(0.9,-100,1d5) =', ice_air_h_lapserate_si(0.9d0, -100d0, 1d5), & 
' 4.42457786755E-04'
write(*,fmt)'ice_air_h_cp(0.9,-100,1d5) =', ice_air_h_cp_si(0.9d0, -100d0, 1d5), & 
' 1766.52051488'
write(*,fmt)'ice_air_h_kappa_s(0.9,-100,1d5) =', ice_air_h_kappa_s_si(0.9d0, -100d0, 1d5), & 
' 8.23031581047E-06'
write(*,fmt)'ice_air_h_density(0.9,-100,1d5) =', ice_air_h_density_si(0.9d0, -100d0, 1d5), & 
' 1.42531895993'


else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.0 (V1.0)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ENDIF) WILL ONLY BE ACTIVATED IF 
!YOU CHANGE TO VERSION 1.0 (V1.0) IN CONSTANTS_0.F90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4A'

write(*,fmt) 'ice_air_massfraction_air(270,1d5) =', ice_air_massfraction_air_si(270d0,1d5), & 
' 0.997058162464'
write(*,fmt) 'ice_air_frostpoint(0.997,1d5) =', ice_air_frostpoint_si(0.997d0,1d5), & 
' 270.232024746'
write(*,fmt) 'ice_air_ict(0.997,300,1d5) =', ice_air_ict_si(0.997d0,300d0,1d5), & 
' 265.221832079'
write(*,fmt) 'ice_air_icl(0.997,300,1d5) =', ice_air_icl_si(0.997d0,300d0,1d5), & 
' 64985.6857989'
write(*,fmt) 'ice_air_condenspressure(0.997,270) =', ice_air_condensationpressure_si(0.997d0,270d0), & 
' 98057.5687060'
write(*,fmt) 'ice_air_sublimationpress(270,1E5) =', ice_air_sublimationpressure_si(270d0,1d5), & 
' 472.041436186'
write(*,fmt) 'ice_air_rh_cct_from_a(0.998,270,1d5) =', ice_air_rh_cct_from_a_si(0.998d0,270d0,1d5), & 
' 0.680235676087'
write(*,fmt) 'ice_air_rh_wmo_from_a(0.998,270,1d5) =', ice_air_rh_wmo_from_a_si(0.998d0,270d0,1d5), & 
' 0.679205607135'
write(*,fmt) 'ice_air_a_from_rh_wmo(0.8,270,1d5) =', ice_air_a_from_rh_wmo_si(0.8d0,270d0,1d5), & 
' 0.997645144450'
write(*,fmt) 'ice_air_a_from_rh_cct(0.8,270,1d5) =', ice_air_a_from_rh_cct_si(0.8d0,270d0,1d5), & 
' 0.997647370806'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833359.68775'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833386.63215'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2833297.67344'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_enthalpy_subl =', ice_air_enthalpy_subl_si(),&
' 2834613.25618'

write(*,"(A)")' '
temp =set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.28849644315'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.26340826697'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 1.28734106365'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_air =', ice_air_density_air_si(),&
' 0.949136077891'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.79054720202d-3'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.79022480092d-3'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 3.86202319096d-3'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_vap =', ice_air_density_vap_si(),&
' 2.84740823367d-3'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 98057.5687060'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_pressure =', ice_air_pressure_si(),&
' 72723.4677050'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270'

write(*,"(A)")' '
temp =set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 270.232024746'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_temp =', ice_air_temperature_si(),&
' 266.511951178'

write(*,"(A)")' '
temp = set_ice_air_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_t_p(270d0,1d5):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.181167192'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_t(.997d0, 270d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_t(.997d0,270d0):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.180958346'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_p(.997d0, 1d5)
write(*,"(A)")'set_ice_air_equilibrium_at_a_p(.997d0,1d5):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.147466108'

write(*,"(A)")' '
temp = set_ice_air_eq_at_a_eta(.997d0, 100d0)
write(*,"(A)")'set_ice_air_equilibrium_at_a_eta(.997d0,100d0):'
write(*,fmt)'ice_air_density_ice =', ice_air_density_ice_si(),&
' 917.682093944'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4B'

write(*,fmt)'ice_air_g(0,0,0,0.5,270,1d5) =', ice_air_g_si( 0, 0, 0, 0.5d0, 270d0, 1d5), &
'-2595.73785825'
write(*,fmt)'ice_air_g(1,0,0,0.5,270,1d5) =', ice_air_g_si( 1, 0, 0, 0.5d0, 270d0, 1d5), &
' 2382.02354606'
write(*,fmt)'ice_air_g(0,1,0,0.5,270,1d5) =', ice_air_g_si( 0, 1, 0, 0.5d0, 270d0, 1d5), &
' 610.261631357'
write(*,fmt)'ice_air_g(0,0,1,0.5,270,1d5) =', ice_air_g_si( 0, 0, 1, 0.5d0, 270d0, 1d5), &
' 0.389737675174'
write(*,fmt)'ice_air_g(1,1,0,0.5,270,1d5) =', ice_air_g_si( 1, 1, 0, 0.5d0, 270d0, 1d5), &
'-1269.42344741'
write(*,fmt)'ice_air_g(1,0,1,0.5,270,1d5) =', ice_air_g_si( 1, 0, 1, 0.5d0, 270d0, 1d5), &
' 0.777294756075'
write(*,fmt)'ice_air_g(0,2,0,0.5,270,1d5) =', ice_air_g_si( 0, 2, 0, 0.5d0, 270d0, 1d5), &
'-7.00885819160'
write(*,fmt)'ice_air_g(0,1,1,0.5,270,1d5) =', ice_air_g_si( 0, 1, 1, 0.5d0, 270d0, 1d5), &
' 1.60133877466d-3'
write(*,fmt)'ice_air_g(0,0,2,0.5,270,1d5) =', ice_air_g_si( 0, 0, 2, 0.5d0, 270d0, 1d5), &
'-3.91271268873d-6'

write(*,fmt)' '
write(*,fmt)'ice_air_solidfrac(0.5,270,1d5) =', ice_air_solidfraction_si(0.5d0, 270d0, 1d5), & 
' 0.498524741260'
write(*,fmt)'ice_air_vaporfrac(0.5,270,1d5) =', ice_air_vapourfraction_si(0.5d0, 270d0, 1d5), & 
' 1.47525873976d-3'
write(*,fmt)'ice_air_g_cp(0.5,270,1d5) =', ice_air_g_cp_si(0.5d0, 270d0, 1d5), & 
' 1892.39171173'
write(*,fmt)'ice_air_g_density(0.5,270,1d5) =', ice_air_g_density_si(0.5d0, 270d0, 1d5), & 
' 2.56582841152'
write(*,fmt)'ice_air_g_enthalpy(0.5,270,1d5) =', ice_air_g_enthalpy_si(0.5d0, 270d0, 1d5), & 
'-167366.378325'
write(*,fmt)'ice_air_g_entropy(0.5,270,1d5) =', ice_air_g_entropy_si(0.5d0, 270d0, 1d5), & 
'-610.261631357'
write(*,fmt)'ice_air_g_exp(0.5,270,1d5) =', ice_air_g_expansion_si(0.5d0, 270d0, 1d5), & 
' 4.10876052449d-3'
write(*,fmt)'ice_air_g_kap_t(0.5,270,1d5) =', ice_air_g_kappa_t_si(0.5d0, 270d0, 1d5), & 
' 1.00393493829d-5'
write(*,fmt)'ice_air_g_lapserate(0.5,270,1d5) =', ice_air_g_lapserate_si(0.5d0, 270d0, 1d5), & 
' 2.28473558871d-4'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_AIR_4C'

write(*,fmt)'ice_air_h(0,0,0,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 0, 0.5d0, -600d0, 1d5), &
'-164588.277175'
write(*,fmt)'ice_air_h(1,0,0,0.5,-600,1d5) =', ice_air_h_si( 1, 0, 0, 0.5d0, -600d0, 1d5), &
' 543.382787268'
write(*,fmt)'ice_air_h(0,1,0,0.5,-600,1d5) =', ice_air_h_si( 0, 1, 0, 0.5d0, -600d0, 1d5), &
' 271.449435548'
write(*,fmt)'ice_air_h(0,0,1,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 1, 0.5d0, -600d0, 1d5), &
' 0.392073699416'
write(*,fmt)'ice_air_h(2,0,0,0.5,-600,1d5) =', ice_air_h_si( 2, 0, 0, 0.5d0, -600d0, 1d5), &
' 224938.714580'
write(*,fmt)'ice_air_h(1,1,0,0.5,-600,1d5) =', ice_air_h_si( 1, 1, 0, 0.5d0, -600d0, 1d5), &
'-177.440249223'
write(*,fmt)'ice_air_h(1,0,1,0.5,-600,1d5) =', ice_air_h_si( 1, 0, 1, 0.5d0, -600d0, 1d5), &
' 0.781966303025'
write(*,fmt)'ice_air_h(0,2,0,0.5,-600,1d5) =', ice_air_h_si( 0, 2, 0, 0.5d0, -600d0, 1d5), &
' 0.139971645623'
write(*,fmt)'ice_air_h(0,1,1,0.5,-600,1d5) =', ice_air_h_si( 0, 1, 1, 0.5d0, -600d0, 1d5), &
' 2.26942461987E-04'
write(*,fmt)'ice_air_h(0,0,2,0.5,-600,1d5) =', ice_air_h_si( 0, 0, 2, 0.5d0, -600d0, 1d5), &
'-3.57055653968E-06'

write(*,"(A)")' '
write(*,fmt)'ice_air_pottemp(0.9,230,5d4,1d5) =', ice_air_pottemp_si(0.9d0, 230d0, 5d4, 1d5), & 
' 266.106323310'
write(*,fmt)'ice_air_potdens(0.9,230,5d4,1d5) =', ice_air_potdensity_si(0.9d0, 230d0, 5d4, 1d5), & 
' 1.45013106603'
write(*,fmt)'ice_air_potenth(0.9,230,5d4,1d5) =', ice_air_potenthalpy_si(0.9d0, 230d0, 5d4, 1d5), & 
'-35779.7348519'
write(*,fmt)'ice_air_h_temp(0.9,-100,1d5) =', ice_air_h_temperature_si(0.9d0, -100d0, 1d5), & 
' 270.382806934'
write(*,fmt)'ice_air_h_lapsrate(0.9,-100,1d5) =', ice_air_h_lapserate_si(0.9d0, -100d0, 1d5), & 
' 4.42476802729d-4'
write(*,fmt)'ice_air_h_cp(0.9,-100,1d5) =', ice_air_h_cp_si(0.9d0, -100d0, 1d5), & 
' 1766.84429536'
write(*,fmt)'ice_air_h_kappa_s(0.9,-100,1d5) =', ice_air_h_kappa_s_si(0.9d0, -100d0, 1d5), & 
' 8.23024116522d-6'
write(*,fmt)'ice_air_h_density(0.9,-100,1d5) =', ice_air_h_density_si(0.9d0, -100d0, 1d5), & 
' 1.42498657691'

endif


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_LIQ_4'

temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_chempot =',ice_liq_chempot_si(),&
' 9972.8817138|1'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_chempot =',ice_liq_chempot_si(),&
' 38870.060519|2'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_density_ice =',ice_liq_density_ice_si(),&
' 917.896690831'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_density_ice =',ice_liq_density_ice_si(),&
' 921.359428514'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_density_liq =',ice_liq_density_liq_si(),&
' 1004.79353660'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_density_liq =',ice_liq_density_liq_si(),&
' 1019.05568894'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_enthalpy_ice =',ice_liq_enthalpy_ice_si(),&
'-324602.983812'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_enthalpy_ice =',ice_liq_enthalpy_ice_si(),&
'-299055.938629'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_enthalpy_liq =',ice_liq_enthalpy_liq_si(),&
' 6945.9270048|3'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_enthalpy_liq =',ice_liq_enthalpy_liq_si(),&
' 26110.7481094'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_enthalpy_melt =',ice_liq_enthalpy_melt_si(),&
' 331548.910817'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_enthalpy_melt =',ice_liq_enthalpy_melt_si(),&
' 325166.686739'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_entropy_ice =',ice_liq_entropy_ice_si(),&
'-1228.24464139'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_entropy_ice =',ice_liq_entropy_ice_si(),&
'-1251.57777462'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_entropy_liq =',ice_liq_entropy_liq_si(),&
'-11.112101272|3'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_entropy_liq =',ice_liq_entropy_liq_si(),&
'-47.2567126291'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_entropy_melt =',ice_liq_entropy_melt_si(),&
' 1217.13254011'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_entropy_melt =',ice_liq_entropy_melt_si(),&
' 1204.32106199'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_pressure_liq =',ice_liq_pressure_liq_si(),&
' 10000000.0070'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_pressure_liq =',ice_liq_pressure_liq_si(),&
' 39313338.8825'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_temp =',ice_liq_temperature_si(),&
' 272.401648869'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_temp =',ice_liq_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_p(1d7)
write(*,"(A)")'set_ice_liq_equilibrium_at_p(1d7): '
write(*,fmt)'ice_liq_volume_melt =',ice_liq_volume_melt_si(),&
'-9.42178903288d-5'

write(*,"(A)")' '
temp = set_ice_liq_eq_at_t(270d0)
write(*,"(A)")'set_ice_liq_equilibrium_at_t(270): '
write(*,fmt)'ice_liq_volume_melt =',ice_liq_volume_melt_si(),&
'-1.04052121182d-4'

write(*,fmt)'ice_liq_meltingpressure(270) =',ice_liq_meltingpressure_si(270d0),&
' 39313338.8825'
write(*,fmt)'ice_liq_meltingtemp(1d7) =',ice_liq_meltingtemperature_si(1d7),&
' 272.401648869'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_VAP_4'

write(*,fmt)'ice_vap_sublimationpressure(270) =',ice_vap_sublimationpressure_si(270d0),&
' 470.059067981'
write(*,fmt)'ice_vap_sublimationtemp(100) =',ice_vap_sublimationtemp_si(100d0),&
' 252.817910215'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_chempot =', ice_vap_chempot_si(),&
'-3895.26747392'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_chempot =', ice_vap_chempot_si(),&
'-26421.2820403'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_density_ice =', ice_vap_density_ice_si(),&
' 917.170465733'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_density_ice =', ice_vap_density_ice_si(),&
' 919.600269745'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_density_vap =', ice_vap_density_vap_si(),&
' 3.77406140772d-3'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_density_vap =', ice_vap_density_vap_si(),&
' 8.57185487853d-4'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_enthalpy_ice =', ice_vap_enthalpy_ice_si(),&
'-340033.434649'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_enthalpy_ice =', ice_vap_enthalpy_ice_si(),&
'-374576.247867'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_enthalpy_vap =', ice_vap_enthalpy_vap_si(),&
' 2495132.21977'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_enthalpy_vap =', ice_vap_enthalpy_vap_si(),&
' 2463525.19629'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_entropy_ice =', ice_vap_entropy_ice_si(),&
'-1244.95617472'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_entropy_ice =', ice_vap_entropy_ice_si(),&
'-1377.09771247'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_entropy_vap =', ice_vap_entropy_vap_si(),&
' 9255.65736018'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_entropy_vap =', ice_vap_entropy_vap_si(),&
' 9848.77406912'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_pressure_vap =', ice_vap_pressure_vap_si(),&
' 470.059067981'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_pressure_vap =', ice_vap_pressure_vap_si(),&
' 100'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_temp =', ice_vap_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_temp =', ice_vap_temperature_si(),&
' 252.817910215'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_volume_subl =', ice_vap_volume_subl_si(),&
' 264.965451558'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_volume_subl =', ice_vap_volume_subl_si(),&
' 1166.60755699'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_entropy_subl =', ice_vap_entropy_subl_si(),&
' 10500.6135349'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_entropy_subl =', ice_vap_entropy_subl_si(),&
' 11225.8717816'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'ice_vap_enthalpy_subl =', ice_vap_enthalpy_subl_si(),&
' 2835165.65442'

write(*,"(A)")' '
temp = set_ice_vap_eq_at_p(100d0)
write(*,"(A)")'set_ice_vap_equilibrium_at_p(100): '
write(*,fmt)'ice_vap_enthalpy_subl =', ice_vap_enthalpy_subl_si(),&
' 2838101.44416'



if(molar_mass_air_L2000 == molar_mass_air_si)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.1 (V1.1)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ELSEIF) USES THIS VALUE OF THE MOLAR MASS.
!TO USES THE ORIGINAL VALUE USED BY LEMMON ET AL. (2000) YOU MUST EDIT CONSTANTS_0 TO 
!CORRESPOND TO VERSION 1.0 (V1.0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4A'

write(*,fmt)'liq_air_massfraction_air(300,1d5) =', liq_air_massfraction_air_si(300d0, 1d5), & 
' 0.977605797727'
write(*,fmt)'liq_air_dewpoint(0.99,1d5) =', liq_air_dewpoint_si(0.99d0, 1d5), & 
' 287.078299795'
write(*,fmt)'liq_air_ict(0.99,300,1d5) =', liq_air_ict_si(0.99d0, 300d0, 1d5), & 
' 284.200207629'
write(*,fmt)'liq_air_icl(0.99,300,1d5) =', liq_air_icl_si(0.99d0, 300d0, 1d5), & 
' 82723.6047631'
write(*,fmt)'liq_air_condenspressure(0.9,300) =', liq_air_condensationpressure_si(0.9d0, 300d0), & 
' 23381.2332935'

write(*,"(A)")' '
write(*,fmt)'liq_air_rh_cct_from_a(0.99,300,1d5) =', liq_air_rh_cct_from_a_si(0.99d0, 300d0, 1d5), & 
' 0.449887886959'
write(*,fmt)'liq_air_rh_wmo_from_a(0.99,300,1d5) =', liq_air_rh_wmo_from_a_si(0.99d0, 300d0, 1d5), & 
' 0.440953686019'
write(*,fmt)'liq_air_a_from_rh_wmo(0.8,300,1d5) =', liq_air_a_from_rh_wmo_si(0.8d0, 300d0, 1d5), & 
' 0.982004037135'
write(*,fmt)'liq_air_a_from_rh_cct(0.8,300,1d5) =', liq_air_a_from_rh_cct_si(0.8d0, 300d0, 1d5), & 
' 0.982133277948'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2458121.74961'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2465656.38630'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2433111.29696'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2434585.53919'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 100'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 145.863545194'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
'-41.9991507402'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 296.711483507'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 290.107386673'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 287.078299795'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 300'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 300'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 121546.373652'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 100000'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 223057.741750'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 100000'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.45154665083'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.20675806022'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 2.57657653270'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.14614215827'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 1.45154665083E-02'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 0.0120675806022'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99,300):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 2.57657653270E-02'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 2.56669393257E-02'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 998.794738784'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 999.256685197'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 996.611581662'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 996.556340389'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4B'

write(*,fmt)'liq_air_g(0,0,0,0.5,300,1d5) =', liq_air_g_si( 0, 0, 0, 0.5d0, 300d0, 1d5), &
'-5396.77820137'
write(*,fmt)'liq_air_g(1,0,0,0.5,300,1d5) =', liq_air_g_si( 1, 0, 0, 0.5d0, 300d0, 1d5), &
'-263.45549120|3'
write(*,fmt)'liq_air_g(0,1,0,0.5,300,1d5) =', liq_air_g_si( 0, 1, 0, 0.5d0, 300d0, 1d5), &
'-343.783393872'
write(*,fmt)'liq_air_g(0,0,1,0.5,300,1d5) =', liq_air_g_si( 0, 0, 1, 0.5d0, 300d0, 1d5), &
' 0.446729465555'
write(*,fmt)'liq_air_g(2,0,0,0.5,300,1d5) =', liq_air_g_si( 2, 0, 0, 0.5d0, 300d0, 1d5), &
' 0'
write(*,fmt)'liq_air_g(1,1,0,0.5,300,1d5) =', liq_air_g_si( 1, 1, 0, 0.5d0, 300d0, 1d5), &
' 98.5580798842'
write(*,fmt)'liq_air_g(1,0,1,0.5,300,1d5) =', liq_air_g_si( 1, 0, 1, 0.5d0, 300d0, 1d5), &
' 0.891452019991'
write(*,fmt)'liq_air_g(0,2,0,0.5,300,1d5) =', liq_air_g_si( 0, 2, 0, 0.5d0, 300d0, 1d5), &
'-14.0995955397'
write(*,fmt)'liq_air_g(0,1,1,0.5,300,1d5) =', liq_air_g_si( 0, 1, 1, 0.5d0, 300d0, 1d5), &
' 2.43183979422E-03'
write(*,fmt)'liq_air_g(0,0,2,0.5,300,1d5) =', liq_air_g_si( 0, 0, 2, 0.5d0, 300d0, 1d5), &
'-4.62360294023E-06'

write(*,fmt)'liq_air_g_cp(0.5,300,1d5) =', liq_air_g_cp_si(0.5d0, 300d0, 1d5), &
' 4229.87866191'
write(*,fmt)'liq_air_g_density(0.5,300,1d5) =', liq_air_g_density_si(0.5d0, 300d0, 1d5), &
' 2.23849125053'
write(*,fmt)'liq_air_g_enth(0.5,300,1d5) =', liq_air_g_enthalpy_si(0.5d0, 300d0, 1d5), &
' 97738.2399604'
write(*,fmt)'liq_air_g_entropy(0.5,300,1d5) =', liq_air_g_entropy_si(0.5d0, 300d0, 1d5), &
' 343.783393872'
write(*,fmt)'liq_air_g_exp(0.5,300, 1d5) =', liq_air_g_expansion_si(0.5d0, 300d0, 1d5), &
' 5.44365210207E-03'
write(*,fmt)'liq_air_g_kap_t(0.5,300,1d5) =', liq_air_g_kappa_t_si(0.5d0, 300d0, 1d5), &
' 1.03498947276E-05'
write(*,fmt)'liq_air_g_lapserate(0.5,300,1d5) =', liq_air_g_lapserate_si(0.5d0, 300d0, 1d5), &
' 1.72475854884E-04'

write(*,"(A)")' '
write(*,fmt)'liq_air_liquidfrac(0.5,300,1d5) =', liq_air_liquidfraction_si(0.5d0, 300d0, 1d5), &
' 0.488546404734'
write(*,fmt)'liq_air_vapourfrac(0.5,300,1d5) =', liq_air_vapourfraction_si(0.5d0, 300d0, 1d5), &
' 1.14535952655E-02'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4C'

write(*,fmt)'liq_air_h(0,0,0,0.5,100,1d5) =', liq_air_h_si( 0, 0, 0, 0.5d0, 100d0, 1d5), &
' 26898.5215492'
write(*,fmt)'liq_air_h(1,0,0,0.5,100,1d5) =', liq_air_h_si( 1, 0, 0, 0.5d0, 100d0, 1d5), &
'-1681.7936611|3'
write(*,fmt)'liq_air_h(0,1,0,0.5,100,1d5) =', liq_air_h_si( 0, 1, 0, 0.5d0, 100d0, 1d5), &
' 280.393544899'
write(*,fmt)'liq_air_h(0,0,1,0.5,100,1d5) =', liq_air_h_si( 0, 0, 1, 0.5d0, 100d0, 1d5), &
' 0.406872930019'
write(*,fmt)'liq_air_h(2,0,0,0.5,100,1d5) =', liq_air_h_si( 2, 0, 0, 0.5d0, 100d0, 1d5), &
' 35.768970891|5'
write(*,fmt)'liq_air_h(1,1,0,0.5,100,1d5) =', liq_air_h_si( 1, 1, 0, 0.5d0, 100d0, 1d5), &
' 1.7859992519|6'
write(*,fmt)'liq_air_h(1,0,1,0.5,100,1d5) =', liq_air_h_si( 1, 0, 1, 0.5d0, 100d0, 1d5), &
' 0.811745643965'
write(*,fmt)'liq_air_h(0,2,0,0.5,100,1d5) =', liq_air_h_si( 0, 2, 0, 0.5d0, 100d0, 1d5), &
' 8.91776656830E-02'
write(*,fmt)'liq_air_h(0,1,1,0.5,100,1d5) =', liq_air_h_si( 0, 1, 1, 0.5d0, 100d0, 1d5), &
' 1.55067379031E-04'
write(*,fmt)'liq_air_h(0,0,2,0.5,100,1d5) =', liq_air_h_si( 0, 0, 2, 0.5d0, 100d0, 1d5), &
'-3.83770118470E-06'

write(*,"(A)")' '
write(*,fmt)'liq_air_pottemp(0.5,300,1d4,1d5) =', liq_air_pottemp_si(0.5d0, 300d0, 1d4, 1d5), &
' 348.222379217'
write(*,fmt)'liq_air_potdens(0.5,300,1d4,1d5) =', liq_air_potdensity_si(0.5d0, 300d0, 1d4, 1d5), &
' 1.22550664945'
write(*,fmt)'liq_air_potenth(0.5,300,1d4,1d5) =', liq_air_potenthalpy_si(0.5d0, 300d0, 1d4, 1d5), &
' 655155.797982'
write(*,fmt)'liq_air_h_temp(0.5,100,1d5) =', liq_air_h_temperature_si(0.5d0, 100d0, 1d5), &
' 280.393544899'
write(*,fmt)'liq_air_h_lapserate(0.5,100,1d5) =', liq_air_h_lapserate_si(0.5d0, 100d0, 1d5), &
' 1.55067379031E-04'
write(*,fmt)'liq_air_h_cp(0.5,100,1d5) =', liq_air_h_cp_si(0.5d0, 100d0, 1d5), &
' 3144.21265404'
write(*,fmt)'liq_air_h_kappa_s(0.5,100,1d5) =', liq_air_h_kappa_s_si(0.5d0, 100d0, 1d5), &
' 9.43218607469E-06'
write(*,fmt)'liq_air_h_density(0.5,100,1d5) =', liq_air_h_density_si(0.5d0, 100d0, 1d5), &
' 2.45776980040'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_ICE_AIR_4'
 
write(*,fmt)'liq_ice_air_ifl_si(.99,100) =', liq_ice_air_ifl_si(.99d0, 100d0), & 
' 83234.7314360'
write(*,fmt)'liq_ice_air_iml_si(.99,100) =', liq_ice_air_iml_si(.99d0, 100d0), &
' 81605.5557728'

write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), & 
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 0'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 3496.16306903'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 7356.12943724'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 900361.135280'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 7.74757979404'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 474.974398769'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 38338.9622424'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 10000'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 67931.6010|764'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 112016.075795'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 706.817425301'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.99'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.961024307544'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.994366063923'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.996583352944'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.2'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 0.99'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 0.1'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0, 0d0, 0.5d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 3.30296152581E-03'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 0.3'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 3.30296152581E-03'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 0.2'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 3.39407694837E-03'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 0.4'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.157198087'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.159302793'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.155'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.151724970'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.159992933'

else


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.0 (V1.0)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ENDIF) WILL ONLY BE ACTIVATED IF 
!YOU CHANGE TO VERSION 1.0 (V1.0) IN CONSTANTS_0.F90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4A'

write(*,fmt)'liq_air_massfraction_air(300,1d5) =', liq_air_massfraction_air_si(300d0, 1d5), & 
' 0.977600623582'
write(*,fmt)'liq_air_dewpoint(0.99,1d5) =', liq_air_dewpoint_si(0.99d0, 1d5), & 
' 287.074721927'
write(*,fmt)'liq_air_ict(0.99,300,1d5) =', liq_air_ict_si(0.99d0, 300d0, 1d5), & 
' 284.195916825'
write(*,fmt)'liq_air_icl(0.99,300,1d5) =', liq_air_icl_si(0.99d0, 300d0, 1d5), & 
' 82719.2377442'
write(*,fmt)'liq_air_condenspressure(0.9,300) =', liq_air_condensationpressure_si(0.9d0, 300d0), & 
' 23385.9346276'

write(*,"(A)")' '
write(*,fmt)'liq_air_rh_cct_from_a(0.99,300,1d5) =', liq_air_rh_cct_from_a_si(0.99d0, 300d0, 1d5), & 
' 0.449783278144'
write(*,fmt)'liq_air_rh_wmo_from_a(0.99,300,1d5) =', liq_air_rh_wmo_from_a_si(0.99d0, 300d0, 1d5), & 
' 0.440849494608'
write(*,fmt)'liq_air_a_from_rh_wmo(0.8,300,1d5) =', liq_air_a_from_rh_wmo_si(0.8d0, 300d0, 1d5), & 
' 0.981999860485'
write(*,fmt)'liq_air_a_from_rh_cct(0.8,300,1d5) =', liq_air_a_from_rh_cct_si(0.8d0, 300d0, 1d5), & 
' 0.982129130677'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2458132.60500'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2465665.32464'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2433111.45416'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_enthalpy_evap =', liq_air_enthalpy_evap_si(),& 
' 2434585.86630'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 100'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 145.864613600'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
'-42.0969943912'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_entropy_air =', liq_air_entropy_air_si(),& 
' 296.780112468'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 290.103069429'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 287.074721927'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 300'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_temp =', liq_air_temperature_si(),& 
' 300'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 121541.303927'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 100000'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 223109.793432'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_pressure =', liq_air_pressure_si(),& 
' 100000'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.45116943913'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.20649187941'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 2.57657751339'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_air =', liq_air_density_air_si(),& 
' 1.14587677780'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 1.45116943913d-2'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 1.20649187941d-2'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99,300):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 2.57657751339d-2'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_vap =', liq_air_density_vap_si(),& 
' 2.56669252749d-2'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_eta(0.99d0,100d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_eta(0.99d0,100d0):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 998.795484566'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_p(0.99d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_a_p(0.99d0,1d5):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 999.257178499'

write(*,"(A)")' '
temp = set_liq_air_eq_at_a_t(0.99d0,300d0)
write(*,"(A)")'set_liq_air_equilibrium_at_a_t(0.99d0,300d0):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 996.611605025'

write(*,"(A)")' '
temp = set_liq_air_eq_at_t_p(300d0,1d5)
write(*,"(A)")'set_liq_air_equilibrium_at_t_p(300d0,1d5):'
write(*,fmt)'liq_air_density_liq =', liq_air_density_liq_si(),& 
' 996.556340389'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4B'

write(*,fmt)'liq_air_g(0,0,0,0.5,300,1d5) =', liq_air_g_si( 0, 0, 0, 0.5d0, 300d0, 1d5), &
'-5397.43216342'
write(*,fmt)'liq_air_g(1,0,0,0.5,300,1d5) =', liq_air_g_si( 1, 0, 0, 0.5d0, 300d0, 1d5), &
'-264.763415292'
write(*,fmt)'liq_air_g(0,1,0,0.5,300,1d5) =', liq_air_g_si( 0, 1, 0, 0.5d0, 300d0, 1d5), &
'-343.818233768'
write(*,fmt)'liq_air_g(0,0,1,0.5,300,1d5) =', liq_air_g_si( 0, 0, 1, 0.5d0, 300d0, 1d5), &
' 0.446835172414'
write(*,fmt)'liq_air_g(2,0,0,0.5,300,1d5) =', liq_air_g_si( 2, 0, 0, 0.5d0, 300d0, 1d5), &
' 0'
write(*,fmt)'liq_air_g(1,1,0,0.5,300,1d5) =', liq_air_g_si( 1, 1, 0, 0.5d0, 300d0, 1d5), &
' 98.4884000933'
write(*,fmt)'liq_air_g(1,0,1,0.5,300,1d5) =', liq_air_g_si( 1, 0, 1, 0.5d0, 300d0, 1d5), &
' 0.891663433710'
write(*,fmt)'liq_air_g(0,2,0,0.5,300,1d5) =', liq_air_g_si( 0, 2, 0, 0.5d0, 300d0, 1d5), &
'-14.1012539816'
write(*,fmt)'liq_air_g(0,1,1,0.5,300,1d5) =', liq_air_g_si( 0, 1, 1, 0.5d0, 300d0, 1d5), &
' 2.43241086432d-3'
write(*,fmt)'liq_air_g(0,0,2,0.5,300,1d5) =', liq_air_g_si( 0, 0, 2, 0.5d0, 300d0, 1d5), &
'-4.62469737987d-6'

write(*,fmt)'liq_air_g_cp(0.5,300,1d5) =', liq_air_g_cp_si(0.5d0, 300d0, 1d5), &
' 4230.37619447'
write(*,fmt)'liq_air_g_density(0.5,300,1d5) =', liq_air_g_density_si(0.5d0, 300d0, 1d5), &
' 2.23796169535'
write(*,fmt)'liq_air_g_enth(0.5,300,1d5) =', liq_air_g_enthalpy_si(0.5d0, 300d0, 1d5), &
' 97748.0379670'
write(*,fmt)'liq_air_g_entropy(0.5,300,1d5) =', liq_air_g_entropy_si(0.5d0, 300d0, 1d5), &
' 343.818233768'
write(*,fmt)'liq_air_g_exp(0.5,300, 1d5) =', liq_air_g_expansion_si(0.5d0, 300d0, 1d5), &
' 5.44364234170d-3'
write(*,fmt)'liq_air_g_kap_t(0.5,300,1d5) =', liq_air_g_kappa_t_si(0.5d0, 300d0, 1d5), &
' 1.03498955887d-5'
write(*,fmt)'liq_air_g_lapserate(0.5,300,1d5) =', liq_air_g_lapserate_si(0.5d0, 300d0, 1d5), &
' 1.72496067902d-4'

write(*,"(A)")' '
write(*,fmt)'liq_air_liquidfrac(0.5,300,1d5) =', liq_air_liquidfraction_si(0.5d0, 300d0, 1d5), &
' 0.488543697765'
write(*,fmt)'liq_air_vapourfrac(0.5,300,1d5) =', liq_air_vapourfraction_si(0.5d0, 300d0, 1d5), &
' 1.14563022353d-2'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_AIR_4C'

write(*,fmt)'liq_air_h(0,0,0,0.5,100,1d5) =', liq_air_h_si( 0, 0, 0, 0.5d0, 100d0, 1d5), &
' 26898.2871262'
write(*,fmt)'liq_air_h(1,0,0,0.5,100,1d5) =', liq_air_h_si( 1, 0, 0, 0.5d0, 100d0, 1d5), &
'-1682.2814754|2'
write(*,fmt)'liq_air_h(0,1,0,0.5,100,1d5) =', liq_air_h_si( 0, 1, 0, 0.5d0, 100d0, 1d5), &
' 280.392597112'
write(*,fmt)'liq_air_h(0,0,1,0.5,100,1d5) =', liq_air_h_si( 0, 0, 1, 0.5d0, 100d0, 1d5), &
' 0.406967545475'
write(*,fmt)'liq_air_h(2,0,0,0.5,100,1d5) =', liq_air_h_si( 2, 0, 0, 0.5d0, 100d0, 1d5), &
' 35.6650436194'
write(*,fmt)'liq_air_h(1,1,0,0.5,100,1d5) =', liq_air_h_si( 1, 1, 0, 0.5d0, 100d0, 1d5), &
' 1.78333840654'
write(*,fmt)'liq_air_h(1,0,1,0.5,100,1d5) =', liq_air_h_si( 1, 0, 1, 0.5d0, 100d0, 1d5), &
' 0.811934874970'
write(*,fmt)'liq_air_h(0,2,0,0.5,100,1d5) =', liq_air_h_si( 0, 2, 0, 0.5d0, 100d0, 1d5), &
' 8.91712318143d-2'
write(*,fmt)'liq_air_h(0,1,1,0.5,100,1d5) =', liq_air_h_si( 0, 1, 1, 0.5d0, 100d0, 1d5), &
' 1.55091234445d-4'
write(*,fmt)'liq_air_h(0,0,2,0.5,100,1d5) =', liq_air_h_si( 0, 0, 2, 0.5d0, 100d0, 1d5), &
'-3.83855233558d-6'

write(*,"(A)")' '
write(*,fmt)'liq_air_pottemp(0.5,300,1d4,1d5) =', liq_air_pottemp_si(0.5d0, 300d0, 1d4, 1d5), &
' 348.223146077'
write(*,fmt)'liq_air_potdens(0.5,300,1d4,1d5) =', liq_air_potdensity_si(0.5d0, 300d0, 1d4, 1d5), &
' 1.22518908856'
write(*,fmt)'liq_air_potenth(0.5,300,1d4,1d5) =', liq_air_potenthalpy_si(0.5d0, 300d0, 1d4, 1d5), &
' 655299.468706'
write(*,fmt)'liq_air_h_temp(0.5,100,1d5) =', liq_air_h_temperature_si(0.5d0, 100d0, 1d5), &
' 280.392597112'
write(*,fmt)'liq_air_h_lapserate(0.5,100,1d5) =', liq_air_h_lapserate_si(0.5d0, 100d0, 1d5), &
' 1.55091234445d-4'
write(*,fmt)'liq_air_h_cp(0.5,100,1d5) =', liq_air_h_cp_si(0.5d0, 100d0, 1d5), &
' 3144.42888595'
write(*,fmt)'liq_air_h_kappa_s(0.5,100,1d5) =', liq_air_h_kappa_s_si(0.5d0, 100d0, 1d5), &
' 9.43208464228d-6'
write(*,fmt)'liq_air_h_density(0.5,100,1d5) =', liq_air_h_density_si(0.5d0, 100d0, 1d5), &
' 2.45719839609'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_ICE_AIR_4'
 
write(*,fmt)'liq_ice_air_ifl_si(.99,100) =', liq_ice_air_ifl_si(.99d0, 100d0), & 
' 83240.7163330'
write(*,fmt)'liq_ice_air_iml_si(.99,100) =', liq_ice_air_iml_si(.99d0, 100d0), &
' 81612.0101048'

write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), & 
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 0'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_entropy =', liq_ice_air_entropy_si(), &
' 3496.20368786'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 7358.21843392'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_enthalpy =', liq_ice_air_enthalpy_si(), &
' 900361.142094'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 7.74649959202'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_density =', liq_ice_air_density_si(), &
' 474.974398772'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 38347.8983147'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 10000'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 67931.6010|764'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 112016.572992'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_pressure =', liq_ice_air_pressure_si(), &
' 706.839968365'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.99'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.961015437341'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.99436473982|6'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.996582564510'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_airfraction =', liq_ice_air_airfraction_si(), &
' 0.2'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 0.99'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_dryairfraction =', liq_ice_air_dryairfraction_si(), &
' 0.1'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0, 0d0, 0.5d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 3.30256857013d-3'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0):'
write(*,fmt)'liq_ice_air_solidfraction =', liq_ice_air_solidfraction_si(), &
' 0.3'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 3.30256857013d-3'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_liquidfraction =', liq_ice_air_liquidfraction_si(), &
' 0.2'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 9.99999999d98'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 3.39486285975d-3'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_vapourfraction =', liq_ice_air_vapourfraction_si(), &
' 0.4'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_a(0.99d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_a(0.99d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.157197424'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_p(1d4)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_p(1d4):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.159302793'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_t(273.155d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_t(273.155d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.155'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_eta_wt(0.99d0, 0d0, 0.5d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_eta_wt(0.99d0,0d0,0.5d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.151724933'
 
write(*,"(A)")' '
temp = set_liq_ice_air_eq_at_wa_wl_wi(0.1d0, 0.2d0, 0.3d0)
write(*,"(A)")'set_liq_ice_air_equilibrium_at_wa_wl_wi(0.1d0,0.2d0,0.3d0):'
write(*,fmt)'liq_ice_air_temp =', liq_ice_air_temperature_si(), &
' 273.159992931'

endif

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM LIQ_VAP_4'

write(*,fmt)'liq_vap_vapourpressure(300) =', liq_vap_vapourpressure_si(300d0),&
' 3536.80675227'
write(*,fmt)'liq_vap_boilingtemp(1d4) =', liq_vap_boilingtemperature_si(1d4),&
' 318.956328924'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_pressure_vap =', liq_vap_pressure_vap_si(),&
' 3536.80675227'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_pressure_vap =', liq_vap_pressure_vap_si(),&
' 10000'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_pressure_liq =', liq_vap_pressure_liq_si(),&
' 3536.80675|227'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_pressure_liq =', liq_vap_pressure_liq_si(),&
' 10000.00000|00'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_temp =', liq_vap_temperature_si(),&
' 300'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_temp =', liq_vap_temperature_si(),&
' 318.956328924'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_chempot =', liq_vap_chempot_si(),&
'-5361.84908682'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_chempot =', liq_vap_chempot_si(),&
'-15259.1024273'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_density_liq =', liq_vap_density_liq_si(),&
' 996.513027468'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_density_liq =', liq_vap_density_liq_si(),&
' 989.833275365'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_density_vap =', liq_vap_density_vap_si(),&
' 2.55896736829d-2'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_density_vap =', liq_vap_density_vap_si(),&
' 6.81657223094d-2'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_entropy_liq =', liq_vap_entropy_liq_si(),&
' 393.089029801'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_entropy_liq =', liq_vap_entropy_liq_si(),&
' 649.195605196'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_entropy_vap =', liq_vap_entropy_vap_si(),&
' 8517.38650061'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_entropy_vap =', liq_vap_entropy_vap_si(),&
' 8148.82019425'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_enthalpy_liq =', liq_vap_enthalpy_liq_si(),&
' 112564.859854'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_enthalpy_liq =', liq_vap_enthalpy_liq_si(),&
' 191805.944559'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_enthalpy_vap =', liq_vap_enthalpy_vap_si(),&
' 2549854.10109'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_enthalpy_vap =', liq_vap_enthalpy_vap_si(),&
' 2583858.67179'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_volume_evap =', liq_vap_volume_evap_si(),&
' 39.0772595686'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_volume_evap =', liq_vap_volume_evap_si(),&
' 14.6691196141'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_entropy_evap =', liq_vap_entropy_evap_si(),&
' 8124.29747080'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_entropy_evap =', liq_vap_entropy_evap_si(),&
' 7499.62458905'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_t(300d0)
write(*,"(A)")'set_liq_vap_equilibrium_at_t(300): '
write(*,fmt)'liq_vap_enthalpy_evap =', liq_vap_enthalpy_evap_si(),&
' 2437289.24124'

write(*,"(A)")' '
temp = set_liq_vap_eq_at_p(1d4)
write(*,"(A)")'set_liq_vap_equilibrium_at_p(1d4): '
write(*,fmt)'liq_vap_enthalpy_evap =', liq_vap_enthalpy_evap_si(),&
' 2392052.72723'


if(molar_mass_air_L2000 == molar_mass_air_si)then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.1 (V1.1)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ELSEIF) USES THIS VALUE OF THE MOLAR MASS.
!TO USES THE ORIGINAL VALUE USED BY LEMMON ET AL. (2000) YOU MUST EDIT CONSTANTS_0 TO 
!CORRESPOND TO VERSION 1.0 (V1.0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_AIR_4'

write(*,fmt)'sea_air_massfraction_air(0.035,300,1d5) =', sea_air_massfraction_air_si(0.035d0, 300d0, 1d5), &
' 0.978029483888'
write(*,fmt)'sea_air_vapourpressure(0.035,300,1d5) =', sea_air_vapourpressure_si(0.035d0, 300d0, 1d5), &
' 3485.92986681'
write(*,fmt)'sea_air_chempot_evap(0.9,0.035,300,1d5) =', sea_air_chempot_evap_si(0.9d0,0.035d0,300d0,1d5), &
' 1.45584069071'
write(*,fmt)'sea_air_entropy_air(0.035,300,1d5) =', sea_air_entropy_air_si(0.035d0, 300d0, 1d5), &
' 293.150672957'
write(*,fmt)'sea_air_condens_temp(0.035,0.99,1d5) =', sea_air_condense_temp_si(0.035d0, 0.99d0, 1d5), &
' 287.367456468'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_enthalpy_evap =', sea_air_enthalpy_evap_si(), &
' 2464738.79136'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_enthalpy_evap_si =', sea_air_enthalpy_evap_si(), &
' 2434549.89770'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_density_air_si =', sea_air_density_air_si(), &
' 1.20553988598'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_density_air_si =', sea_air_density_air_si(), &
' 1.14642944448'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_density_vap_si =', sea_air_density_vap_si(), &
' 1.20553988598E-02'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_density_vap_si =', sea_air_density_vap_si(), &
' 2.51876465812E-02'
 

else

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESULTS CORRESPONDING TO VERSION 1.0 (V1.0)
!NOTE THAT LEMMON AND COLLEAGUES HAVE DECIDED THAT IT IS PREFERABLE TO USE THE LATEST 
!MOLAR MASS WITH THEIR MOLAR FORMULATION GIVEN IN LEMMON ET AL. (2000).
!THE PORTION OF THE CODE BELOW (FROM HERE TO ENDIF) WILL ONLY BE ACTIVATED IF 
!YOU CHANGE TO VERSION 1.0 (V1.0) IN CONSTANTS_0.F90.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_AIR_4'

write(*,fmt)'sea_air_massfraction_air(0.035,300,1d5) =', sea_air_massfraction_air_si(0.035d0, 300d0, 1d5), &
' 0.978024405439'
write(*,fmt)'sea_air_vapourpressure(0.035,300,1d5) =', sea_air_vapourpressure_si(0.035d0, 300d0, 1d5), &
' 3485.92802425'
write(*,fmt)'sea_air_chempot_evap(0.9,0.035,300,1d5) =', sea_air_chempot_evap_si(0.9d0,0.035d0,300d0,1d5), &
' 1.45564492467'
write(*,fmt)'sea_air_entropy_air(0.035,300,1d5) =', sea_air_entropy_air_si(0.035d0, 300d0, 1d5), &
' 293.218508134'
write(*,fmt)'sea_air_condens_temp(0.035,0.99,1d5) =', sea_air_condense_temp_si(0.035d0, 0.99d0, 1d5), &
' 287.363869941'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_enthalpy_evap =', sea_air_enthalpy_evap_si(), &
' 2464747.68964'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_enthalpy_evap_si =', sea_air_enthalpy_evap_si(), &
' 2434550.22477'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_density_air_si =', sea_air_density_air_si(), &
' 1.20527399511'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_density_air_si =', sea_air_density_air_si(), &
' 1.14616388271'
 
write(*,"(A)")' '
temp = set_sea_air_eq_at_s_a_p(0.035d0, 0.99d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_a_p(0.035,0.99,1d5):'
write(*,fmt)'sea_air_density_vap_si =', sea_air_density_vap_si(), &
' 1.20527399511E-02'

write(*,"(A)")' '
temp = set_sea_air_eq_at_s_t_p(0.035d0, 300d0, 1d5)
write(*,"(A)")'set_sea_air_equilibrium_at_s_t_p(0.035,300,1d5):'
write(*,fmt)'sea_air_density_vap_si =', sea_air_density_vap_si(), &
' 0.0251876327865'

endif

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_VAP_4'

write(*,fmt)'sea_vap_boilingtemp(0.035,640) =', SEA_VAP_BOILINGTEMPERATURE_SI(0.035d0, 640d0),&
' 274.042416829'
write(*,"(A)")' '
write(*,fmt)'sea_vap_cp_seavap(0.035,274,610) =', sea_vap_cp_seavap_si(0.035d0, 274d0, 610d0), & 
' 756270.43159|3'
write(*,fmt)'sea_vap_dens_seavap(0.035,274,610) =', sea_vap_density_seavap_si(0.035d0, 274d0, 610d0), & 
' 7.27092786882d-3'
write(*,fmt)'sea_vap_enth_seavap(0.035,274,610) =', sea_vap_enthalpy_seavap_si(0.035d0, 274d0, 610d0), & 
' 1661118.41089'
write(*,fmt)'sea_vap_entr_seavap(0.035,274,610) =', sea_vap_entropy_seavap_si(0.035d0, 274d0, 610d0), & 
' 6072.50817709'
write(*,fmt)'sea_vap_exp_seavap(0.035,274,610) =', sea_vap_expansion_seavap_si(0.035d0, 274d0, 610d0), & 
' 0.458863421347'
write(*,fmt)'sea_vap_kap_t_seavap(0.035,274,610) =', sea_vap_kappa_t_seavap_si(0.035d0, 274d0, 610d0), & 
' 1.19990585451d-2'
write(*,fmt)'sea_vap_vapourpres(0.035,274) =', sea_vap_vapourpressure_si(0.035d0, 274d0),&
' 638.044692615'
write(*,fmt)'sea_vap_brinesalinity(272,1d5) =', sea_vap_brinesalinity_si(274d0, 640d0),&
' 2.9439629829|4E-02'
write(*,fmt)'sea_vap_brfr_seavap(.0035,274,640) =', sea_vap_brinefraction_seavap_si(0.0035d0, 274d0, 640d0),&
' 0.118887364425'


write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_density_sea =', sea_vap_density_sea_si(),&
' 1027.87349556'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_density_sea =', sea_vap_density_sea_si(),&
' 1027.87626132'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_density_sea =', sea_vap_density_sea_si(),&
' 1023.42713047'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_density_vap =', sea_vap_density_vap_si(),&
' 5.06324890264d-3'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_density_vap =', sea_vap_density_vap_si(),&
' 5.04855547811d-3'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_density_vap =', sea_vap_density_vap_si(),&
' 5.06403699513d-3'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_enthalpy_evap =', sea_vap_enthalpy_evap_si(),&
' 2498295.32187'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_enthalpy_evap =', sea_vap_enthalpy_evap_si(),&
' 2498395.40101'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_enthalpy_evap =', sea_vap_enthalpy_evap_si(),&
' 2498551.19875'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_enthalpy_sea =', sea_vap_enthalpy_sea_si(),&
' 3465.1189614|4'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_enthalpy_sea =', sea_vap_enthalpy_sea_si(),&
' 3295.9662929|9'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_enthalpy_sea =', sea_vap_enthalpy_sea_si(),&
' 3405.9335373|0'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_enthalpy_vap =', sea_vap_enthalpy_vap_si(),&
' 2502546.89358'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_enthalpy_vap =', sea_vap_enthalpy_vap_si(),&
' 2502469.07187'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_enthalpy_vap =', sea_vap_enthalpy_vap_si(),&
' 2502466.96633'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_entropy_sea =', sea_vap_entropy_sea_si(),&
' 13.0616891215'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_entropy_sea =', sea_vap_entropy_sea_si(),&
' 12.4443983378'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_entropy_sea =', sea_vap_entropy_sea_si(),&
' 14.025681511|2'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_entropy_vap =', sea_vap_entropy_vap_si(),&
' 9140.56256065'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_entropy_vap =', sea_vap_entropy_vap_si(),&
' 9141.68990452'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_entropy_vap =', sea_vap_entropy_vap_si(),&
' 9140.27087793'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_pressure =', sea_vap_pressure_si(),&
' 640'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_pressure =', sea_vap_pressure_si(),&
' 638.044692615'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_pressure =', sea_vap_pressure_si(),&
' 640'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_salinity =', sea_vap_salinity_si(),&
' 0.035'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_salinity =', sea_vap_salinity_si(),&
' 0.035'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_salinity =', sea_vap_salinity_si(),&
' 2.9439629829|4d-2'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_temp =', sea_vap_temperature_si(),&
' 274.042416829'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_temp =', sea_vap_temperature_si(),&
' 274'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_temp =', sea_vap_temperature_si(),&
' 274'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_p(0.035d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_p(0.035,640): '
write(*,fmt)'sea_vap_volume_evap =', sea_vap_volume_evap_si(),&
' 197.500648110'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_s_t(0.035d0, 274d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_s_t(0.035,274): '
write(*,fmt)'sea_vap_volume_evap =', sea_vap_volume_evap_si(),&
' 198.075461154'

write(*,"(A)")' '
temp = set_sea_vap_eq_at_t_p(274d0, 640d0)
write(*,"(A)")'set_sea_vap_equilibrium_at_t_p(274,640): '
write(*,fmt)'sea_vap_volume_evap =', sea_vap_volume_evap_si(),&
' 197.469911653'

write(*,"(A)")' '
write(*,fmt)'sea_vap_g_si(0,0,0,0.035,274,610) =', sea_vap_g_si( 0, 0, 0, 0.035d0, 274d0, 610d0), &
'-2748.82963245'
write(*,fmt)'sea_vap_g_si(1,0,0,0.035,274,610) =', sea_vap_g_si( 1, 0, 0, 0.035d0, 274d0, 610d0), &
' 151028.257424'
write(*,fmt)'sea_vap_g_si(0,1,0,0.035,274,610) =', sea_vap_g_si( 0, 1, 0, 0.035d0, 274d0, 610d0), &
'-6072.50817709'
write(*,fmt)'sea_vap_g_si(0,0,1,0.035,274,610) =', sea_vap_g_si( 0, 0, 1, 0.035d0, 274d0, 610d0), &
' 137.534028399'
write(*,fmt)'sea_vap_g_si(2,0,0,0.035,274,610) =', sea_vap_g_si( 2, 0, 0, 0.035d0, 274d0, 610d0), &
' 0'
write(*,fmt)'sea_vap_g_si(1,1,0,0.035,274,610) =', sea_vap_g_si( 1, 1, 0, 0.035d0, 274d0, 610d0), &
' 14965.0677011'
write(*,fmt)'sea_vap_g_si(1,0,1,0.035,274,610) =', sea_vap_g_si( 1, 0, 1, 0.035d0, 274d0, 610d0), &
'-321.591932572'
write(*,fmt)'sea_vap_g_si(0,2,0,0.035,274,610) =', sea_vap_g_si( 0, 2, 0, 0.035d0, 274d0, 610d0), &
'-2760.11106421'
write(*,fmt)'sea_vap_g_si(0,1,1,0.035,274,610) =', sea_vap_g_si( 0, 1, 1, 0.035d0, 274d0, 610d0), &
' 63.1093348229'
write(*,fmt)'sea_vap_g_si(0,0,2,0.035,274,610) =', sea_vap_g_si( 0, 0, 2, 0.035d0, 274d0, 610d0), &
'-1.65027885871'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM SEA_LIQ_4'

write(*,fmt)'sea_liq_osmoticpressure(0.035,300,1d5) =', sea_liq_osmoticpressure_si(0.035d0, 300d0, 1d5),&
' 2594603.2096|8'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM SEA_ICE_4'

write(*,fmt)'sea_ice_brinesalinity(270,1d5) =', sea_ice_brinesalinity_si(270d0, 1d5),&
' 0.056026415032|2'
write(*,fmt)'sea_ice_brfr_seaice(0.035,270,1d5) =', sea_ice_brinefraction_seaice_si(0.035d0, 270d0, 1d5),&
' 0.62470532836|8'
write(*,fmt)'sea_ice_meltingpr(0.035,270) =', sea_ice_meltingpressure_si(0.035d0, 270d0),&
' 16132047.4385'
write(*,fmt)'sea_ice_freezingtemp(0.035,1d5) =', sea_ice_freezingtemperature_si(0.035d0, 1d5),&
' 271.240373585'
write(*,fmt)'sea_ice_dtfdp(0.035,1d5) =', sea_ice_dtfdp_si(0.035d0, 1d5),&
' 7.48210942879E-08'
write(*,fmt)'sea_ice_dtfds(0.035,1d5) =', sea_ice_dtfds_si(0.035d0, 1d5),&
'-56.8751336296'
write(*,fmt)'sea_ice_cp_seaice(0.035,270,1d5) =', sea_ice_cp_seaice_si(0.035d0, 270d0, 1d5),&
' 62868.901512|6'
write(*,fmt)'sea_ice_dens_seaice(0.035,270,1d5) =', sea_ice_density_seaice_si(0.035d0, 270d0, 1d5),&
' 993.156434117'
write(*,fmt)'sea_ice_enth_seaice(0.035,270,1d5) =', sea_ice_enthalpy_seaice_si(0.035d0, 270d0, 1d5),&
'-135534.28750|3'
write(*,fmt)'sea_ice_entr_seaice(0.035,270,1d5) =', sea_ice_entropy_seaice_si(0.035d0, 270d0, 1d5),&
'-500.44544418|1'
write(*,fmt)'sea_ice_exp_seaice(0.035,270,1d5) =', sea_ice_expansion_seaice_si(0.035d0, 270d0, 1d5),&
'-1.6473132873|8'
write(*,fmt)'sea_ice_kap_t_seaice(0.035,270,1d5) =', sea_ice_kappa_t_seaice_si(0.035d0, 270d0, 1d5),&
' 1.5651344134|8d-9'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035,1d5): '
write(*,fmt)'sea_ice_density_ice =', sea_ice_density_ice_si(),&
' 917.000739687'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035,270): '
write(*,fmt)'sea_ice_density_ice =', sea_ice_density_ice_si(),&
' 918.898527655'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270,1d5): '
write(*,fmt)'sea_ice_density_ice =', sea_ice_density_ice_si(),&
' 917.181167192'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_density_sea =', sea_ice_density_sea_si(),&
' 1028.05199645'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_density_sea =', sea_ice_density_sea_si(),&
' 1035.73670169'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_density_sea =', sea_ice_density_sea_si(),&
' 1045.16805918'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_melt =', sea_ice_enthalpy_melt_si(),&
' 329942.976285'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_enthalpy_melt =', sea_ice_enthalpy_melt_si(),&
' 326829.393605'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_melt =', sea_ice_enthalpy_melt_si(),&
' 328249.119579'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_volume_melt =', sea_ice_volume_melt_si(),&
'-9.10140854473d-5'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_volume_melt =', sea_ice_volume_melt_si(),&
'-9.67135426848d-5'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_volume_melt =', sea_ice_volume_melt_si(),&
'-9.18186917900-5'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_ice =', sea_ice_enthalpy_ice_si(),&
'-337351.999358'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_enthalpy_ice =', sea_ice_enthalpy_ice_si(),&
'-323205.968289'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_ice =', sea_ice_enthalpy_ice_si(),&
'-339929.555499'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_sea =', sea_ice_enthalpy_sea_si(),&
'-7613.19337919'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_enthalpy_sea =', sea_ice_enthalpy_sea_si(),&
' 2832.9491040|7'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_enthalpy_sea =', sea_ice_enthalpy_sea_si(),&
'-12742.8664892'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_entropy_ice =', sea_ice_entropy_ice_si(),&
'-1235.44872812'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_entropy_ice =', sea_ice_entropy_ice_si(),&
'-1247.71314646'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_entropy_ice =', sea_ice_entropy_ice_si(),&
'-1244.97335506'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_entropy_sea =', sea_ice_entropy_sea_si(),&
'-27.9264598103'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_entropy_sea =', sea_ice_entropy_sea_si(),&
'-46.736116956|0'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_entropy_sea =', sea_ice_entropy_sea_si(),&
'-53.166791114|4'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_pressure =', sea_ice_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_pressure =', sea_ice_pressure_si(),&
' 16132047.4385'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_pressure =', sea_ice_pressure_si(),&
' 100000'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_temp =', sea_ice_temperature_si(),&
' 271.240373585'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_temp =', sea_ice_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_temp =', sea_ice_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_p(0.035d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_p(0.035d0,1d5): '
write(*,fmt)'sea_ice_salinity =', sea_ice_salinity_si(),&
' 0.035'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_s_t(0.035d0, 270d0)
write(*,"(A)")'set_sea_ice_equilibrium_at_s_t(0.035d0,270): '
write(*,fmt)'sea_ice_salinity =', sea_ice_salinity_si(),&
' 0.035'

write(*,"(A)")' '
temp = set_sea_ice_eq_at_t_p(270d0, 1d5)
write(*,"(A)")'set_sea_ice_equilibrium_at_t_p(270d0,1d5): '
write(*,fmt)'sea_ice_salinity =', sea_ice_salinity_si(),&
' 0.056026415032|2'

write(*,"(A)")' '
write(*,fmt)'sea_ice_g(0,0,0,0.035,270,1d5) =', sea_ice_g_si( 0, 0, 0, 0.035d0, 270d0, 1d5),&
'-414.0175745|47'
write(*,fmt)'sea_ice_g(1,0,0,0.035,270,1d5) =', sea_ice_g_si( 1, 0, 0, 0.035d0, 270d0, 1d5),&
' 96363.773049|5'
write(*,fmt)'sea_ice_g(0,1,0,0.035,270,1d5) =', sea_ice_g_si( 0, 1, 0, 0.035d0, 270d0, 1d5),&
' 500.44544418|1'
write(*,fmt)'sea_ice_g(0,0,1,0.035,270,1d5) =', sea_ice_g_si( 0, 0, 1, 0.035d0, 270d0, 1d5),&
' 1.00689072300d-3'
write(*,fmt)'sea_ice_g(2,0,0,0.035,270,1d5) =', sea_ice_g_si( 2, 0, 0, 0.035d0, 270d0, 1d5),&
' 0'
write(*,fmt)'sea_ice_g(1,1,0,0.035,270,1d5) =', sea_ice_g_si( 1, 1, 0, 0.035d0, 270d0, 1d5),&
'-1144.02883419'
write(*,fmt)'sea_ice_g(1,0,1,0.035,270,1d5) =', sea_ice_g_si( 1, 0, 1, 0.035d0, 270d0, 1d5),&
'-8.62856321467d-4'
write(*,fmt)'sea_ice_g(0,2,0,0.035,270,1d5) =', sea_ice_g_si( 0, 2, 0, 0.035d0, 270d0, 1d5),&
'-232.84778338|0'
write(*,fmt)'sea_ice_g(0,1,1,0.035,270,1d5) =', sea_ice_g_si( 0, 1, 1, 0.035d0, 270d0, 1d5),&
'-1.65866446694d-5'
write(*,fmt)'sea_ice_g(0,0,2,0.035,270,1d5) =', sea_ice_g_si( 0, 0, 2, 0.035d0, 270d0, 1d5),&
'-1.57591932118d-12'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM SEA_ICE_VAP_4'

temp = set_sea_ice_vap_eq_at_s(0.035d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_s(0.035): '
write(*,fmt)'sea_ice_vap_salinity =', sea_ice_vap_salinity_si(),&
' 0.035'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'sea_ice_vap_salinity =', sea_ice_vap_salinity_si(),&
' 5.61489288506d-2'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_p(500d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_p(500): '
write(*,fmt)'sea_ice_vap_salinity =', sea_ice_vap_salinity_si(),&
' 4.3895587882|8d-2'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_s(0.035d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_s(0.035): '
write(*,fmt)'sea_ice_vap_density_vap =', sea_ice_vap_density_vap_si(),&
' 4.17156419318d-3'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'sea_ice_vap_density_vap =', sea_ice_vap_density_vap_si(),&
' 3.77406140772d-3'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_p(500d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_p(500): '
write(*,fmt)'sea_ice_vap_density_vap =', sea_ice_vap_density_vap_si(),&
' 4.00364833230d-3'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_s(0.035d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_s(0.035): '
write(*,fmt)'sea_ice_vap_temp =', sea_ice_vap_temperature_si(),&
' 271.247815057'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'sea_ice_vap_temp =', sea_ice_vap_temperature_si(),&
' 270'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_p(500d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_p(500): '
write(*,fmt)'sea_ice_vap_temp =', sea_ice_vap_temperature_si(),&
' 270.734430917'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_s(0.035d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_s(0.035): '
write(*,fmt)'sea_ice_vap_pressure =', sea_ice_vap_pressure_si(),&
' 521.950349225'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_t(270d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_t(270): '
write(*,fmt)'sea_ice_vap_pressure =', sea_ice_vap_pressure_si(),&
' 470.059067981'

write(*,"(A)")' '
temp = set_sea_ice_vap_eq_at_p(500d0)
write(*,"(A)")'set_sea_ice_vap_equilibrium_at_p(500): '
write(*,fmt)'sea_ice_vap_pressure =', sea_ice_vap_pressure_si(),&
' 500'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)") 'COMPARISON OF CHECK VALUES FROM ICE_FLU_5'

write(*,fmt)'fit_ice_liq_pressure_si(272.4) =', fit_ice_liq_pressure_si(272.4d0), &
' 0.100213724736d8'

write(*,fmt)'fit_ice_liq_temperature_si(1d7) =', fit_ice_liq_temperature_si(1d7), &
' 0.272401569225d3'

write(*,fmt)'fit_ice_vap_pressure_si(270) =', fit_ice_vap_pressure_si(270d0), &
' 0.470061877574d3'

end subroutine


end module values_mix
