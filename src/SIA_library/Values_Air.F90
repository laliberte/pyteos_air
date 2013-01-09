module values_air

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE AIR ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     AIR_1,       FILE AIR_1.F90
!     AIR_2,       FILE AIR_2.F90
!     AIR_3a,      FILE AIR_3a.F90
!     AIR_3b,      FILE AIR_3b.F90
!     AIR_3c,      FILE AIR_3c.F90
!     AIR_5,       FILE AIR_5.F90

use constants_0
use convert_0
use maths_0
use air_1
use air_2
use air_3a
use air_3b
use air_3c
use air_5

implicit none 
character*16, private :: version = '8 Jun 2010'

character (len = *), parameter :: fmt = "(A,T40,1X,E19.12,3X,A)"

contains

!*******************************************************************
subroutine air_chk_values
!*******************************************************************

implicit none 

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
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF DRY AND MOIST AIR'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_1'
write(*,fmt)'dry_f(0,0,300,0.001) =', dry_f_si( 0, 0, 300d0,0.001d0),&
'-696239.965190'
write(*,fmt)'dry_f(1,0,300,0.001) =', dry_f_si( 1, 0, 300d0,0.001d0),&
'-2124.55145456'
write(*,fmt)'dry_f(0,1,300,0.001) =', dry_f_si( 0, 1, 300d0,0.001d0),&
' 86114714.9596'
write(*,fmt)'dry_f(2,0,300,0.001) =', dry_f_si( 2, 0, 300d0,0.001d0),&
'-2.39242390806'
write(*,fmt)'dry_f(1,1,300,0.001) =', dry_f_si( 1, 1, 300d0,0.001d0),&
' 287049.624545'
write(*,fmt)'dry_f(0,2,300,0.001) =', dry_f_si( 0, 2, 300d0,0.001d0),&
'-86114738036.1'

write(*,"(A)")' '
write(*,fmt)'air_molar_mass(0.5) =', air_molar_mass_si(0.5d0),&
' 0.222142374909E-01'
write(*,fmt)'air_molfraction_air(0.5) =', air_molfraction_air_si(0.5d0),&    !dbg
' 0.383460809718'
write(*,fmt)'air_molfraction_vap(0.5) =', air_molfraction_vap_si(0.5d0),&
' 0.616539190282'
write(*,fmt)'air_massfraction_air(0.5) =', air_massfraction_air_si(0.5d0),&
' 0.616539190282'
write(*,fmt)'air_massfraction_vap(0.5) =', air_massfraction_vap_si(0.5d0),&
' 0.383460809718'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_2'
write(*,fmt)'air_f_mix(0,0,0,0.9,300,1) =', air_f_mix_si(0,0,0,0.9d0,300d0,1d0),&
'-25.8379179182'
write(*,fmt)'air_f_mix(1,0,0,0.9,300,1) =', air_f_mix_si(1,0,0,0.9d0,300d0,1d0),&
' 233.827370434'
write(*,fmt)'air_f_mix(0,1,0,0.9,300,1) =', air_f_mix_si(0,1,0,0.9d0,300d0,1d0),&
' 0.164195952060'
write(*,fmt)'air_f_mix(0,0,1,0.9,300,1) =', air_f_mix_si(0,0,1,0.9d0,300d0,1d0),&
'-26.2357498619'
write(*,fmt)'air_f_mix(2,0,0,0.9,300,1) =', air_f_mix_si(2,0,0,0.9d0,300d0,1d0),&
' 500.273928155'
write(*,fmt)'air_f_mix(1,1,0,0.9,300,1) =', air_f_mix_si(1,1,0,0.9d0,300d0,1d0),&
'-1.53932744055'
write(*,fmt)'air_f_mix(1,0,1,0.9,300,1) =', air_f_mix_si(1,0,1,0.9d0,300d0,1d0),&
' 241.520643317'
write(*,fmt)'air_f_mix(0,2,0,0.9,300,1) =', air_f_mix_si(0,2,0,0.9d0,300d0,1d0),&
'-0.687329742959E-3'
write(*,fmt)'air_f_mix(0,1,1,0.9,300,1) =', air_f_mix_si(0,1,1,0.9d0,300d0,1d0),&
' 0.172192606103'
write(*,fmt)'air_f_mix(0,0,2,0.9,300,1) =', air_f_mix_si(0,0,2,0.9d0,300d0,1d0),&
'-0.795663887493'

write(*,"(A)")' '
write(*,fmt)'air_f(0,0,0,0.9,300,1) =', air_f_si(0,0,0,0.9d0,300d0,1d0),&
'-95019.5943231'
write(*,fmt)'air_f(1,0,0,0.9,300,1) =', air_f_si(1,0,0,0.9d0,300d0,1d0),&
'-205645.554995'
write(*,fmt)'air_f(0,1,0,0.9,300,1) =', air_f_si(0,1,0,0.9d0,300d0,1d0),&
'-940.175394023'
write(*,fmt)'air_f(0,0,1,0.9,300,1) =', air_f_si(0,0,1,0.9d0,300d0,1d0),&
' 91175.3848662'
write(*,fmt)'air_f(2,0,0,0.9,300,1) =', air_f_si(2,0,0,0.9d0,300d0,1d0),&
' 1447768.46379'
write(*,fmt)'air_f(1,1,0,0.9,300,1) =', air_f_si(1,1,0,0.9d0,300d0,1d0),&
' 7443.09771950'
write(*,fmt)'air_f(1,0,1,0.9,300,1) =', air_f_si(1,0,1,0.9d0,300d0,1d0),&
'-48847.9096826'
write(*,fmt)'air_f(0,2,0,0.9,300,1) =', air_f_si(0,2,0,0.9d0,300d0,1d0),&
'-2.96482218054'
write(*,fmt)'air_f(0,1,1,0.9,300,1) =', air_f_si(0,1,1,0.9d0,300d0,1d0),&
' 312.063110700'
write(*,fmt)'air_f(0,0,2,0.9,300,1) =', air_f_si(0,0,2,0.9d0,300d0,1d0),&
'-91421.4440689'

write(*,"(A)")' '
write(*,fmt)'air_f_cp(0.9,300,1) =', air_f_cp_si(0.9d0,300d0,1d0),&
' 1210.74031058'
write(*,fmt)'air_f_cv(0.9,300,1) =', air_f_cv_si(0.9d0,300d0,1d0),&
' 889.446654163'
write(*,fmt)'air_f_enthalpy(0.9,300,1) =', air_f_enthalpy_si(0.9d0,300d0,1d0),&
' 278208.408750'
write(*,fmt)'air_f_entropy(0.9,300,1) =', air_f_entropy_si(0.9d0,300d0,1d0),&
' 940.175394023'
write(*,fmt)'air_f_expansion(0.9,300,1) =', air_f_expansion_si(0.9d0,300d0,1d0),&
' 3.43193033077E-03'
write(*,fmt)'air_f_gibbs_energy(0.9,300,1) =', air_f_gibbs_energy_si(0.9d0,300d0,1d0),&
'-3844.20945693'
write(*,fmt)'air_f_intrnl_energy(0.9,300,1) =', air_f_internal_energy_si(0.9d0,300d0,1d0),&
' 187033.023884'
write(*,fmt)'air_f_kappa_s(0.9,300,1) =', air_f_kappa_s_si(0.9d0,300d0,1d0),&
' 8.07913626816E-06'
write(*,fmt)'air_f_kappa_t(0.9,300,1) =', air_f_kappa_t_si(0.9d0,300d0,1d0),&
' 1.09975521396E-05'
write(*,fmt)'air_f_lapserate(0.9,300,1) =', air_f_lapserate_si(0.9d0,300d0,1d0),&
' 8.50371537341E-04'
write(*,fmt)'air_f_pressure(0.9,300,1) =', air_f_pressure_si(0.9d0,300d0,1d0),&
' 91175.3848662'
write(*,fmt)'air_f_soundspeed(0.9,300,1) =', air_f_soundspeed_si(0.9d0,300d0,1d0),&
' 351.817577078'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3A'

write(*,fmt)'air_g(0,0,0,0.9,300,1d5) =', air_g_si( 0, 0, 0, 0.9d0, 300d0, 1d5),&
' 4577.93065689'
write(*,fmt)'air_g(1,0,0,0.9,300,1d5) =', air_g_si( 1, 0, 0, 0.9d0, 300d0, 1d5),&
'-210141.953243'
write(*,fmt)'air_g(0,1,0,0.9,300,1d5) =', air_g_si( 0, 1, 0, 0.9d0, 300d0, 1d5),&
'-911.170080461'
write(*,fmt)'air_g(0,0,1,0.9,300,1d5) =', air_g_si( 0, 0, 1, 0.9d0, 300d0, 1d5),&
' 0.911504137473'
write(*,fmt)'air_g(2,0,0,0.9,300,1d5) =', air_g_si( 2, 0, 0, 0.9d0, 300d0, 1d5),&
' 1415779.23410'
write(*,fmt)'air_g(1,1,0,0.9,300,1d5) =', air_g_si( 1, 1, 0, 0.9d0, 300d0, 1d5),&
' 7566.34779196'
write(*,fmt)'air_g(1,0,1,0.9,300,1d5) =', air_g_si( 1, 0, 1, 0.9d0, 300d0, 1d5),&
'-0.483353002176'
write(*,fmt)'air_g(0,2,0,0.9,300,1d5) =', air_g_si( 0, 2, 0, 0.9d0, 300d0, 1d5),&
'-4.15449972148'
write(*,fmt)'air_g(0,1,1,0.9,300,1d5) =', air_g_si( 0, 1, 1, 0.9d0, 300d0, 1d5),&
' 3.15111222847E-03'
write(*,fmt)'air_g(0,0,2,0.9,300,1d5) =', air_g_si( 0, 0, 2, 0.9d0, 300d0, 1d5),&
'-9.14462130186E-06'
write(*,fmt)'air_density(0.9,300,1d5) =', air_density_si(0.9d0, 300d0, 1d5),&
' 1.09708772444'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3B'

write(*,fmt)'air_g_chempot_vap(0.9,300,1E5) =', air_g_chempot_vap_si(0.9d0, 300d0, 1d5),&
' 193705.688576'
write(*,fmt)'air_g_contraction(0.9,300,1d5) =', air_g_contraction_si(0.9d0, 300d0, 1d5),&
' 0.530280645260'
write(*,fmt)'air_g_comprfactor(0.9,300,1d5) =', air_g_compressibilityfactor_si(0.9d0, 300d0, 1d5),&
' 0.997825670991'
write(*,fmt)'air_g_cp(0.9,300,1d5) =', air_g_cp_si(0.9d0, 300d0, 1d5),&
' 1246.34991644'
write(*,fmt)'air_g_cv(0.9,300,1d5) =', air_g_cv_si(0.9d0, 300d0, 1d5),&
' 920.600781012'
write(*,fmt)'air_g_enthalpy(0.9,300,1d5) =', air_g_enthalpy_si(0.9d0, 300d0, 1d5),&
' 277928.954795'
write(*,fmt)'air_g_entropy(0.9,300,1d5) =', air_g_entropy_si(0.9d0, 300d0, 1d5),&
' 911.170080461'
write(*,fmt)'air_g_expansion(0.9,300,1d5) =', air_g_expansion_si(0.9d0, 300d0, 1d5),&
' 3.45704654420E-03'
write(*,fmt)'air_g_density(0.9,300,1d5) =', air_g_density_si(0.9d0, 300d0, 1d5),&
' 1.09708772444'
write(*,fmt)'air_g_gibbs_energy(0.9,300,1d5) =', air_g_gibbs_energy_si(0.9d0, 300d0, 1d5),&
' 4577.93065689'
write(*,fmt)'air_g_int_energy(0.9,300,1d5) =', air_g_internal_energy_si(0.9d0, 300d0, 1d5),&
' 186778.541048'
write(*,fmt)'air_g_kappa_s(0.9,300,1d5) =', air_g_kappa_s_si(0.9d0, 300d0, 1d5),&
' 7.41034505449E-06'
write(*,fmt)'air_g_kappa_t(0.9,300,1d5) =', air_g_kappa_t_si(0.9d0, 300d0, 1d5),&
' 1.00324517749E-05'
write(*,fmt)'air_g_lapserate(0.9,300,1d5) =', air_g_lapserate_si(0.9d0, 300d0, 1d5),&
' 7.58481752251E-04'
write(*,fmt)'air_g_soundspeed(0.9,300,1d5) =', air_g_soundspeed_si(0.9d0, 300d0, 1d5),&
' 350.719656182'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3C'

write(*,fmt)'air_h(0,0,0,0.9,900,1d5) =', air_h_si( 0, 0, 0, 0.9d0, 900d0, 1d5),&
' 274592.611783'
write(*,fmt)'air_h(1,0,0,0.9,900,1d5) =', air_h_si( 1, 0, 0, 0.9d0, 900d0, 1d5),&
'-229706.905765'
write(*,fmt)'air_h(0,1,0,0.9,900,1d5) =', air_h_si( 0, 1, 0, 0.9d0, 900d0, 1d5),&
' 297.403043058'
write(*,fmt)'air_h(0,0,1,0.9,900,1d5) =', air_h_si( 0, 0, 1, 0.9d0, 900d0, 1d5),&
' 0.903262695636'
write(*,fmt)'air_h(2,0,0,0.9,900,1d5) =', air_h_si( 2, 0, 0, 0.9d0, 900d0, 1d5),&
' 13963273.0104'
write(*,fmt)'air_h(1,1,0,0.9,900,1d5) =', air_h_si( 1, 1, 0, 0.9d0, 900d0, 1d5),&
' 1676.85098552'
write(*,fmt)'air_h(1,0,1,0.9,900,1d5) =', air_h_si( 1, 0, 1, 0.9d0, 900d0, 1d5),&
'-0.467537679038'
write(*,fmt)'air_h(0,2,0,0.9,900,1d5) =', air_h_si( 0, 2, 0, 0.9d0, 900d0, 1d5),&
' 0.223684689765'
write(*,fmt)'air_h(0,1,1,0.9,900,1d5) =', air_h_si( 0, 1, 1, 0.9d0, 900d0, 1d5),&
' 7.15703143992E-04'
write(*,fmt)'air_h(0,0,2,0.9,900,1d5) =', air_h_si( 0, 0, 2, 0.9d0, 900d0, 1d5),&
'-6.78105152859E-06'

write(*,"(A)")' '
write(*,fmt)'air_pottemp(0.9,300,5d4,1d5) =', air_pottemp_si(0.9d0, 300d0, 5d4, 1d5),&
' 363.653905688'
write(*,fmt)'air_potdens(0.9,300,5d4,1d5) =', air_potdensity_si(0.9d0, 300d0, 5d4, 1d5),&
' 0.903509489711'
write(*,fmt)'air_potenth(0.9,300,5d4,1d5) =', air_potenthalpy_si(0.9d0, 300d0, 5d4, 1d5),&
' 348872.568665'
write(*,fmt)'air_temperature(0.9,900,1d5) =', air_temperature_si(0.9d0, 900d0, 1d5),&
' 297.403043058'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_5'

write(*,fmt)'air_lapserate_moist_c100m(50,20,1000) =',air_lapserate_moist_c100m(50d0, 20d0, 1000d0),&
' 0.971588085046'


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
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF DRY AND MOIST AIR'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_1'
write(*,fmt)'dry_f(0,0,300,0.001) =', dry_f_si( 0, 0, 300d0,0.001d0),&
'-696384.495286'
write(*,fmt)'dry_f(1,0,300,0.001) =', dry_f_si( 1, 0, 300d0,0.001d0),&
'-2124.98673214'
write(*,fmt)'dry_f(0,1,300,0.001) =', dry_f_si( 0, 1, 300d0,0.001d0),&
' 86135114.6608'
write(*,fmt)'dry_f(2,0,300,0.001) =', dry_f_si( 2, 0, 300d0,0.001d0),&
'-2.39299064929'
write(*,fmt)'dry_f(1,1,300,0.001) =', dry_f_si( 1, 1, 300d0,0.001d0),&
' 287117.623822'
write(*,fmt)'dry_f(0,2,300,0.001) =', dry_f_si( 0, 2, 300d0,0.001d0),&
'-86135137748.2'

write(*,"(A)")' '
write(*,fmt)'air_molar_mass(0.5) =', air_molar_mass_si(0.5d0),&
' 2.22122197774E-02'
write(*,fmt)'air_molfraction_air(0.5) =', air_molfraction_air_si(0.5d0),&    !dbg
' 0.383516809814'
write(*,fmt)'air_molfraction_vap(0.5) =', air_molfraction_vap_si(0.5d0),&
' 0.616483190186'
write(*,fmt)'air_massfraction_air(0.5) =', air_massfraction_air_si(0.5d0),&
' 0.616483190186'
write(*,fmt)'air_massfraction_vap(0.5) =', air_massfraction_vap_si(0.5d0),&
' 0.383516809814'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_2'
write(*,fmt)'air_f_mix(0,0,0,0.9,300,1) =', air_f_mix_si(0,0,0,0.9d0,300d0,1d0),&
'-25.8379179182'
write(*,fmt)'air_f_mix(1,0,0,0.9,300,1) =', air_f_mix_si(1,0,0,0.9d0,300d0,1d0),&
' 233.827370434'
write(*,fmt)'air_f_mix(0,1,0,0.9,300,1) =', air_f_mix_si(0,1,0,0.9d0,300d0,1d0),&
' 0.164195952060'
write(*,fmt)'air_f_mix(0,0,1,0.9,300,1) =', air_f_mix_si(0,0,1,0.9d0,300d0,1d0),&
'-26.2357498619'
write(*,fmt)'air_f_mix(2,0,0,0.9,300,1) =', air_f_mix_si(2,0,0,0.9d0,300d0,1d0),&
' 500.273928155'
write(*,fmt)'air_f_mix(1,1,0,0.9,300,1) =', air_f_mix_si(1,1,0,0.9d0,300d0,1d0),&
'-1.53932744055'
write(*,fmt)'air_f_mix(1,0,1,0.9,300,1) =', air_f_mix_si(1,0,1,0.9d0,300d0,1d0),&
' 241.520643317'
write(*,fmt)'air_f_mix(0,2,0,0.9,300,1) =', air_f_mix_si(0,2,0,0.9d0,300d0,1d0),&
'-0.687329742959E-3'
write(*,fmt)'air_f_mix(0,1,1,0.9,300,1) =', air_f_mix_si(0,1,1,0.9d0,300d0,1d0),&
' 0.172192606103'
write(*,fmt)'air_f_mix(0,0,2,0.9,300,1) =', air_f_mix_si(0,0,2,0.9d0,300d0,1d0),&
'-0.795663887493'

write(*,"(A)")' '
write(*,fmt)'air_f(0,0,0,0.9,300,1) =', air_f_si(0,0,0,0.9d0,300d0,1d0),&
'-95024.7900038'
write(*,fmt)'air_f(1,0,0,0.9,300,1) =', air_f_si(1,0,0,0.9d0,300d0,1d0),&
'-205630.937996'
write(*,fmt)'air_f(0,1,0,0.9,300,1) =', air_f_si(0,1,0,0.9d0,300d0,1d0),&
'-940.150652241'
write(*,fmt)'air_f(0,0,1,0.9,300,1) =', air_f_si(0,0,1,0.9d0,300d0,1d0),&
' 91193.7358460'
write(*,fmt)'air_f(2,0,0,0.9,300,1) =', air_f_si(2,0,0,0.9d0,300d0,1d0),&
' 1447791.10862'
write(*,fmt)'air_f(1,1,0,0.9,300,1) =', air_f_si(1,1,0,0.9d0,300d0,1d0),&
' 7443.19342188'
write(*,fmt)'air_f(1,0,1,0.9,300,1) =', air_f_si(1,0,1,0.9d0,300d0,1d0),&
'-48827.5293333'
write(*,fmt)'air_f(0,2,0,0.9,300,1) =', air_f_si(0,2,0,0.9d0,300d0,1d0),&
'-2.96533251321'
write(*,fmt)'air_f(0,1,1,0.9,300,1) =', air_f_si(0,1,1,0.9d0,300d0,1d0),&
' 312.124501058'
write(*,fmt)'air_f(0,0,2,0.9,300,1) =', air_f_si(0,0,2,0.9d0,300d0,1d0),&
'-91439.8037142'

write(*,"(A)")' '
write(*,fmt)'air_f_cp(0.9,300,1) =', air_f_cp_si(0.9d0,300d0,1d0),&
' 1210.95501131'
write(*,fmt)'air_f_cv(0.9,300,1) =', air_f_cv_si(0.9d0,300d0,1d0),&
' 889.599753963'
write(*,fmt)'air_f_enthalpy(0.9,300,1) =', air_f_enthalpy_si(0.9d0,300d0,1d0),&
' 278214.141515'
write(*,fmt)'air_f_entropy(0.9,300,1) =', air_f_entropy_si(0.9d0,300d0,1d0),&
' 940.150652241'
write(*,fmt)'air_f_expansion(0.9,300,1) =', air_f_expansion_si(0.9d0,300d0,1d0),&
' 3.43191318698E-03'
write(*,fmt)'air_f_gibbs_energy(0.9,300,1) =', air_f_gibbs_energy_si(0.9d0,300d0,1d0),&
'-3831.05415778'
write(*,fmt)'air_f_intrnl_energy(0.9,300,1) =', air_f_internal_energy_si(0.9d0,300d0,1d0),&
' 187020.405669'
write(*,fmt)'air_f_kappa_s(0.9,300,1) =', air_f_kappa_s_si(0.9d0,300d0,1d0),&
' 8.07746486644E-06'
write(*,fmt)'air_f_kappa_t(0.9,300,1) =', air_f_kappa_t_si(0.9d0,300d0,1d0),&
' 1.09953341546E-05'
write(*,fmt)'air_f_lapserate(0.9,300,1) =', air_f_lapserate_si(0.9d0,300d0,1d0),&
' 8.50216520413E-04'
write(*,fmt)'air_f_pressure(0.9,300,1) =', air_f_pressure_si(0.9d0,300d0,1d0),&
' 91193.7358460'
write(*,fmt)'air_f_soundspeed(0.9,300,1) =', air_f_soundspeed_si(0.9d0,300d0,1d0),&
' 351.853974518'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3A'

write(*,fmt)'air_g(0,0,0,0.9,300,1d5) =', air_g_si( 0, 0, 0, 0.9d0, 300d0, 1d5),&
' 4574.43347102'
write(*,fmt)'air_g(1,0,0,0.9,300,1d5) =', air_g_si( 1, 0, 0, 0.9d0, 300d0, 1d5),&
'-210115.722120'
write(*,fmt)'air_g(0,1,0,0.9,300,1d5) =', air_g_si( 0, 1, 0, 0.9d0, 300d0, 1d5),&
'-911.203092171'
write(*,fmt)'air_g(0,0,1,0.9,300,1d5) =', air_g_si( 0, 0, 1, 0.9d0, 300d0, 1d5),&
' 0.911688234199'
write(*,fmt)'air_g(2,0,0,0.9,300,1d5) =', air_g_si( 2, 0, 0, 0.9d0, 300d0, 1d5),&
' 1415842.71443'
write(*,fmt)'air_g(1,1,0,0.9,300,1d5) =', air_g_si( 1, 1, 0, 0.9d0, 300d0, 1d5),&
' 7566.48354961'
write(*,fmt)'air_g(1,0,1,0.9,300,1d5) =', air_g_si( 1, 0, 1, 0.9d0, 300d0, 1d5),&
'-0.483164035114'
write(*,fmt)'air_g(0,2,0,0.9,300,1d5) =', air_g_si( 0, 2, 0, 0.9d0, 300d0, 1d5),&
'-4.15491953969'
write(*,fmt)'air_g(0,1,1,0.9,300,1d5) =', air_g_si( 0, 1, 1, 0.9d0, 300d0, 1d5),&
' 3.15167187428E-03'
write(*,fmt)'air_g(0,0,2,0.9,300,1d5) =', air_g_si( 0, 0, 2, 0.9d0, 300d0, 1d5),&
'-9.14645118671E-06'
write(*,fmt)'air_density(0.9,300,1d5) =', air_density_si(0.9d0, 300d0, 1d5),&
' 1.09686619009421'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3B'

write(*,fmt)'air_g_chempot_vap(0.9,300,1E5) =', air_g_chempot_vap_si(0.9d0, 300d0, 1d5),&
' 193678.583379'
write(*,fmt)'air_g_contraction(0.9,300,1d5) =', air_g_contraction_si(0.9d0, 300d0, 1d5),&
' 0.529966294386'
write(*,fmt)'air_g_comprfactor(0.9,300,1d5) =', air_g_compressibilityfactor_si(0.9d0, 300d0, 1d5),&
' 0.997826654404'
write(*,fmt)'air_g_cp(0.9,300,1d5) =', air_g_cp_si(0.9d0, 300d0, 1d5),&
' 1246.47586191'
write(*,fmt)'air_g_cv(0.9,300,1d5) =', air_g_cv_si(0.9d0, 300d0, 1d5),&
' 920.676202551'
write(*,fmt)'air_g_enthalpy(0.9,300,1d5) =', air_g_enthalpy_si(0.9d0, 300d0, 1d5),&
' 277935.361122'
write(*,fmt)'air_g_entropy(0.9,300,1d5) =', air_g_entropy_si(0.9d0, 300d0, 1d5),&
' 911.203092171'
write(*,fmt)'air_g_expansion(0.9,300,1d5) =', air_g_expansion_si(0.9d0, 300d0, 1d5),&
' 3.45696232117E-03'
write(*,fmt)'air_g_density(0.9,300,1d5) =', air_g_density_si(0.9d0, 300d0, 1d5),&
' 1.09686619009'
write(*,fmt)'air_g_gibbs_energy(0.9,300,1d5) =', air_g_gibbs_energy_si(0.9d0, 300d0, 1d5),&
' 4574.43347102'
write(*,fmt)'air_g_int_energy(0.9,300,1d5) =', air_g_internal_energy_si(0.9d0, 300d0, 1d5),&
' 186766.537702'
write(*,fmt)'air_g_kappa_s(0.9,300,1d5) =', air_g_kappa_s_si(0.9d0, 300d0, 1d5),&
' 7.41018952703E-06'
write(*,fmt)'air_g_kappa_t(0.9,300,1d5) =', air_g_kappa_t_si(0.9d0, 300d0, 1d5),&
' 1.00324330660E-05'
write(*,fmt)'air_g_lapserate(0.9,300,1d5) =', air_g_lapserate_si(0.9d0, 300d0, 1d5),&
' 7.58539809057E-04'
write(*,fmt)'air_g_soundspeed(0.9,300,1d5) =', air_g_soundspeed_si(0.9d0, 300d0, 1d5),&
' 350.758752725'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_3C'

write(*,fmt)'air_h(0,0,0,0.9,900,1d5) =', air_h_si( 0, 0, 0, 0.9d0, 900d0, 1d5),&
' 274589.199192'
write(*,fmt)'air_h(1,0,0,0.9,900,1d5) =', air_h_si( 1, 0, 0, 0.9d0, 900d0, 1d5),&
'-229734.926517'
write(*,fmt)'air_h(0,1,0,0.9,900,1d5) =', air_h_si( 0, 1, 0, 0.9d0, 900d0, 1d5),&
' 297.395859294'
write(*,fmt)'air_h(0,0,1,0.9,900,1d5) =', air_h_si( 0, 0, 1, 0.9d0, 900d0, 1d5),&
' 0.903422384672'
write(*,fmt)'air_h(2,0,0,0.9,900,1d5) =', air_h_si( 2, 0, 0, 0.9d0, 900d0, 1d5),&
' 13959323.1016'
write(*,fmt)'air_h(1,1,0,0.9,900,1d5) =', air_h_si( 1, 1, 0, 0.9d0, 900d0, 1d5),&
' 1676.33822835'
write(*,fmt)'air_h(1,0,1,0.9,900,1d5) =', air_h_si( 1, 0, 1, 0.9d0, 900d0, 1d5),&
'-0.467305898316'
write(*,fmt)'air_h(0,2,0,0.9,900,1d5) =', air_h_si( 0, 2, 0, 0.9d0, 900d0, 1d5),&
' 0.223618143111'
write(*,fmt)'air_h(0,1,1,0.9,900,1d5) =', air_h_si( 0, 1, 1, 0.9d0, 900d0, 1d5),&
' 7.15646251375E-04'
write(*,fmt)'air_h(0,0,2,0.9,900,1d5) =', air_h_si( 0, 0, 2, 0.9d0, 900d0, 1d5),&
'-6.78234488400E-06'

write(*,"(A)")' '
write(*,fmt)'air_pottemp(0.9,300,5d4,1d5) =', air_pottemp_si(0.9d0, 300d0, 5d4, 1d5),&
' 363.654373179'
write(*,fmt)'air_potdens(0.9,300,5d4,1d5) =', air_potdensity_si(0.9d0, 300d0, 5d4, 1d5),&
' 0.903326577188'
write(*,fmt)'air_potenth(0.9,300,5d4,1d5) =', air_potenthalpy_si(0.9d0, 300d0, 5d4, 1d5),&
' 348892.581997'
write(*,fmt)'air_temperature(0.9,900,1d5) =', air_temperature_si(0.9d0, 900d0, 1d5),&
' 297.395859294'


write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM AIR_5'

write(*,fmt)'air_lapserate_moist_c100m(50,20,1000) =',air_lapserate_moist_c100m(50d0, 20d0, 1000d0),&
' 0.971359680720'

endif

end subroutine

end module values_air
