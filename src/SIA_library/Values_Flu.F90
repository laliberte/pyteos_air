module values_flu

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE FLU ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     FLU_1,       FILE FLU_1.F90
!     FLU_2,       FILE FLU_2.F90
!     FLU_3A,      FILE FLU_3A.F90
!     FLU_3B,      FILE FLU_3B.F90
!     LIQ_F03_5,   FILE LIQ_F03_5

use constants_0
use convert_0
use flu_1
use flu_2
use flu_3a
use flu_3b
use liq_f03_5

implicit none 
character*16, private :: version = '8 Jun 2010'
character (len = *), parameter :: fmt = "(A,T36,1X,E19.12,3X,A)"

contains

!*******************************************************************
subroutine flu_chk_values
!*******************************************************************

implicit none 

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF PURE WATER'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM FLU_1'
write(*,fmt)'flu_f(0,0,300,1000) =', flu_f_si(0,0,300d0,1000d0), &
'-5351.74115204'
write(*,fmt)'flu_f(1,0,300,1000) =', flu_f_si(1,0,300d0,1000d0), &
'-390.904170767'
write(*,fmt)'flu_f(0,1,300,1000) =', flu_f_si(0,1,300d0,1000d0), &
' 7.83300135597'
write(*,fmt)'flu_f(2,0,300,1000) =', flu_f_si(2,0,300d0,1000d0), &
'-13.6840204925'
write(*,fmt)'flu_f(1,1,300,1000) =', flu_f_si(1,1,300d0,1000d0), &
' 0.639359046588'
write(*,fmt)'flu_f(0,2,300,1000) =', flu_f_si(0,2,300d0,1000d0), &
' 2.24824656167'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM FLU_2'
write(*,fmt)'flu_cp(300,1000) =', flu_cp_si(300d0,1000d0),&
' 4159.37519963'
write(*,fmt)'flu_cv(300,1000) =', flu_cv_si(300d0,1000d0),&
' 4105.20614776'
write(*,fmt)'flu_enthalpy(300,1000) =', flu_enthalpy_si(300d0,1000d0),&
' 119752.511434'
write(*,fmt)'flu_entropy(300,1000) =', flu_entropy_si(300d0,1000d0),&
' 390.904170767'
write(*,fmt)'flu_expansion(300,1000) =', flu_expansion_si(300d0,1000d0),&
' 2.82413312531E-04'
write(*,fmt)'flu_gibbs_energy(300,1000) =', flu_gibbs_energy_si(300d0,1000d0),&
' 2481.26020392'
write(*,fmt)'flu_internal_energy(300,1000) =', flu_internal_energy_si(300d0,1000d0),&
' 111919.510078'
write(*,fmt)'flu_kappa_s(300,1000) =', flu_kappa_s_si(300d0,1000d0),&
' 4.35960581171E-10'
write(*,fmt)'flu_kappa_t(300,1000) =', flu_kappa_t_si(300d0,1000d0),&
' 4.41713172024E-10'
write(*,fmt)'flu_lapserate(300,1000) =', flu_lapserate_si(300d0,1000d0),&
' 2.0369403983|3E-08'
write(*,fmt)'flu_pressure(300,1000) =', flu_pressure_si(300d0,1000d0),&
' 7833001.35597'
write(*,fmt)'flu_soundspeed(300,1000) =', flu_soundspeed_si(300d0,1000d0),&
' 1514.52479780'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM FLU_3A'
write(*,fmt)'liq_g(0,0,300,1E5) =', liq_g_si( 0, 0, 300d0, 1d5),&
'-5265.05045577'
write(*,fmt)'liq_g(1,0,300,1E5) =', liq_g_si( 1, 0, 300d0, 1d5),&
'-393.062433815'
write(*,fmt)'liq_g(0,1,300,1E5) =', liq_g_si( 0, 1, 300d0, 1d5),&
' 1.00345555938E-03'
write(*,fmt)'liq_g(2,0,300,1E5) =', liq_g_si( 2, 0, 300d0, 1d5),&
'-13.9354650734'
write(*,fmt)'liq_g(1,1,300,1E5) =', liq_g_si( 1, 1, 300d0, 1d5),&
' 2.75753316815E-07'
write(*,fmt)'liq_g(0,2,300,1E5) =', liq_g_si( 0, 2, 300d0, 1d5),&
'-4.52072086722E-13'

write(*,"(A)")' '
write(*,fmt)'vap_g(0,0,300,1E3) =', vap_g_si( 0, 0, 300d0, 1d3),&
'-180090.341338'
write(*,fmt)'vap_g(1,0,300,1E3) =', vap_g_si( 1, 0, 300d0, 1d3),&
'-9103.67940087'
write(*,fmt)'vap_g(0,1,300,1E3) =', vap_g_si( 0, 1, 300d0, 1d3),&
' 138.388478069'
write(*,fmt)'vap_g(2,0,300,1E3) =', vap_g_si( 2, 0, 300d0, 1d3),&
'-6.24707163427'
write(*,fmt)'vap_g(1,1,300,1E3) =', vap_g_si( 1, 1, 300d0, 1d3),&
' 0.462704658818'
write(*,fmt)'vap_g(0,2,300,1E3) =', vap_g_si( 0, 2, 300d0, 1d3),&
'-0.138455798864'

write(*,"(A)")' '
write(*,fmt)'liq_density(300,1E5) =', LIQ_DENSITY_SI(300d0, 1d5),&
' 996.556340389'
write(*,fmt)'vap_density(300,1E3) =', VAP_DENSITY_SI(300d0, 1d3),&
' 7.22603510025-03'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM FLU_3B'
write(*,fmt)'liq_cp(300,1e5) =', liq_cp_si(300d0, 1d5),&
' 4180.63952202'
write(*,fmt)'liq_cv(300,1e5) =', liq_cv_si(300d0, 1d5),&
' 4130.17861503'
write(*,fmt)'liq_enthalpy(300,1e5) =', liq_enthalpy_si(300d0, 1d5),&
' 112653.679689'
write(*,fmt)'liq_entropy(300,1e5) =', liq_entropy_si(300d0, 1d5),&
' 393.062433815'
write(*,fmt)'liq_expansion(300,1e5) =', liq_expansion_si(300d0, 1d5),&
' 2.74803716256E-04'
write(*,fmt)'liq_gibbs_energy(300,1e5) =', liq_gibbs_energy_si(300d0, 1d5),&
'-5265.05045577'
write(*,fmt)'liq_internal_energy(300,1e5) =', liq_internal_energy_si(300d0, 1d5),&
' 112553.334133'
write(*,fmt)'liq_kappa_s(300,1e5) =', liq_kappa_s_si(300d0, 1d5),&
' 4.45077521253E-10'
write(*,fmt)'liq_kappa_t(300,1e5) =', liq_kappa_t_si(300d0, 1d5),&
' 4.50515304336E-10'
write(*,fmt)'liq_lapserate(300,1e5) =', liq_lapserate_si(300d0, 1d5),&
' 1.97878804448E-08'
write(*,fmt)'liq_soundspeed(300,1e5) =', liq_soundspeed_si(300d0, 1d5),&
' 1501.52041506'

write(*,"(A)")' '
write(*,fmt)'vap_cp(300,1e3) =', vap_cp_si(300d0, 1d3),&
' 1874.12149028'
write(*,fmt)'vap_cv(300,1e3) =', vap_cv_si(300d0, 1d3),&
' 1410.22845789'
write(*,fmt)'vap_enthalpy(300,1e3) =', vap_enthalpy_si(300d0, 1d3),&
' 2551013.47892'
write(*,fmt)'vap_entropy(300,1e3) =', vap_entropy_si(300d0, 1d3),&
' 9103.67940087'
write(*,fmt)'vap_expansion(300,1e3) =', vap_expansion_si(300d0, 1d3),&
' 3.34352010567E-03'
write(*,fmt)'vap_gibbs_energy(300,1e3) =', vap_gibbs_energy_si(300d0, 1d3),&
'-180090.341338'
write(*,fmt)'vap_internal_energy(300,1e3) =', vap_internal_energy_si(300d0, 1d3),&
' 2412625.00085'
write(*,fmt)'vap_kappa_s(300,1e3) =', vap_kappa_s_si(300d0, 1d3),&
' 7.52840457971E-04'
write(*,fmt)'vap_kappa_t(300,1e3) =', vap_kappa_t_si(300d0, 1d3),&
' 1.00048646242E-03'
write(*,fmt)'vap_lapserate(300,1e3) =', vap_lapserate_si(300d0, 1d3),&
' 7.40674488635E-02'
write(*,fmt)'vap_soundspeed(300,1e3) =', vap_soundspeed_si(300d0, 1d3),&
' 428.744430495'




write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM LIQ_F03_5'
write(*,fmt)'fit_liq_g_f03_si(0,0,300,1E5) =', fit_liq_g_f03_si( 0, 0, 300d0, 1d5),&
'-5265.05056073'
write(*,fmt)'fit_liq_g_f03_si(1,0,300,1E5) =', fit_liq_g_f03_si( 1, 0, 300d0, 1d5),&
'-393.062597709'
write(*,fmt)'fit_liq_g_f03_si(0,1,300,1E5) =', fit_liq_g_f03_si( 0, 1, 300d0, 1d5),&
' 0.100345554745E-02'
write(*,fmt)'fit_liq_g_f03_si(2,0,300,1E5) =', fit_liq_g_f03_si( 2, 0, 300d0, 1d5),&
'-13.9354762020'
write(*,fmt)'fit_liq_g_f03_si(1,1,300,1E5) =', fit_liq_g_f03_si( 1, 1, 300d0, 1d5),&
' 0.275754520492E-06'
write(*,fmt)'fit_liq_g_f03_si(0,2,300,1E5) =', fit_liq_g_f03_si( 0, 2, 300d0, 1d5),&
'-0.452067557155E-12'


write(*,"(A)")' '
write(*,fmt)'fit_liq_cp_f03_si(300,1e5) =', fit_liq_cp_f03_si(300d0, 1d5),&
' 4180.64286060'
write(*,fmt)'fit_liq_density_f03_si(300,1E5) =', fit_liq_density_f03_si(300d0, 1d5),&
' 996.556352243'
write(*,fmt)'fit_liq_expansion_f03_si(300,1e5) =', fit_liq_expansion_f03_si(300d0, 1d5),&
' 0.274804919056E-03'
write(*,fmt)'fit_liq_kappa_t_f03_si(300,1e5) =', fit_liq_kappa_t_f03_si(300d0, 1d5),&
' 0.450510795725E-09'
write(*,fmt)'fit_liq_soundspeed_f03_si(300,1e5) =', fit_liq_soundspeed_f03_si(300d0, 1d5),&
' 1501.52808421'

end subroutine

end module values_flu
