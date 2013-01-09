module values_sea

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE SEAWATER ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     SAL_1,       FILE SAL_1.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_1,       FILE FLU_1.F90
!     FLU_2,       FILE FLU_2.F90
!     FLU_3A,      FILE FLU_3A.F90
!     SEA_3A,      FILE SEA_3A.F90
!     SEA_3B,      FILE SEA_3B.F90
!     SEA_3C,      FILE SEA_3C.F90
!     SEA_3D,      FILE SEA_3D.F90

use constants_0
use convert_0
use sal_1
use sal_2
use flu_1
use flu_2
use flu_3a
use sea_3a
use sea_3b
use sea_3c
use sea_3d
use sea_5a

implicit none 
character*16, private :: version = '8 Jun 2010'
character (len = *), parameter :: fmt = "(A,T47,1X,E19.12,3X,A)"

contains

!==========================================================================
subroutine sea_chk_values
!==========================================================================

implicit none 
character*13 key

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF SEAWATER'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_3A'
write(*,fmt)'sea_g(0,0,0,0.035,300,1d5) =', sea_g_si(0,0,0,0.035d0,300d0,1d5), &
'-5114.9961985|7'
write(*,fmt)'sea_g(1,0,0,0.035,300,1d5) =', sea_g_si(1,0,0,0.035d0,300d0,1d5), &
' 78593.7757371'
write(*,fmt)'sea_g(0,1,0,0.035,300,1d5) =', sea_g_si(0,1,0,0.035d0,300d0,1d5), &
'-374.452240540'
write(*,fmt)'sea_g(0,0,1,0.035,300,1d5) =', sea_g_si(0,0,1,0.035d0,300d0,1d5), &
' 9.77858615182E-04'
write(*,fmt)'sea_g(2,0,0,0.035,300,1d5) =', sea_g_si(2,0,0,0.035d0,300d0,1d5), &
' 2247550.41118'
write(*,fmt)'sea_g(1,1,0,0.035,300,1d5) =', sea_g_si(1,1,0,0.035d0,300d0,1d5), &
' 789.934255688'
write(*,fmt)'sea_g(1,0,1,0.035,300,1d5) =', sea_g_si(1,0,1,0.035d0,300d0,1d5), &
'-7.16682401265E-04'
write(*,fmt)'sea_g(0,2,0,0.035,300,1d5) =', sea_g_si(0,2,0,0.035d0,300d0,1d5), &
'-13.3358324655'
write(*,fmt)'sea_g(0,1,1,0.035,300,1d5) =', sea_g_si(0,1,1,0.035d0,300d0,1d5), &
' 3.04605539768E-07'
write(*,fmt)'sea_g(0,0,2,0.035,300,1d5) =', sea_g_si(0,0,2,0.035d0,300d0,1d5), &
'-4.10945807960E-13'

write(*,"(A)")' '
write(*,fmt)'sea_chempot_h2o(0.035,300,1d5) =', sea_chempot_h2o_si(0.035d0, 300d0, 1d5), &
'-7865.7783493|7'
write(*,fmt)'sea_chempot_rel(0.035,300,1d5) =', sea_chempot_rel_si(0.035d0, 300d0, 1d5), &
' 78593.7757371'
write(*,fmt)'sea_g_contr_t(0.035,300,1d5) =', sea_g_contraction_t_si(0.035d0, 300d0, 1d5), &
' 0.732910044599'
write(*,fmt)'sea_cp(0.035,300,1d5) =', sea_cp_si(0.035d0, 300d0, 1d5), &
' 4000.74973964'
write(*,fmt)'sea_dens(0.035,300,1d5) =', sea_density_si(0.035d0, 300d0, 1d5), &
' 1022.64272613'
write(*,fmt)'sea_enthalpy(0.035,300,1d5) =', sea_enthalpy_si(0.035d0, 300d0, 1d5), &
' 107220.675963'
write(*,fmt)'sea_entropy(0.035,300,1d5) =', sea_entropy_si(0.035d0, 300d0, 1d5), &
' 374.452240540'
write(*,fmt)'sea_g_exp_t(0.035,300,1d5) =', sea_g_expansion_t_si(0.035d0, 300d0, 1d5), &
' 3.11502639583E-04'
write(*,fmt)'sea_gibbs_energy(0.035,300,1d5) =', sea_gibbs_energy_si(0.035d0, 300d0, 1d5), &
'-5114.9961985|7'
write(*,fmt)'sea_int_energy(0.035,300,1d5) =', sea_internal_energy_si(0.035d0, 300d0, 1d5), &
' 107122.890102'
write(*,fmt)'sea_kappa_s(0.035,300,1d5) =', sea_kappa_s_si(0.035d0, 300d0, 1d5), &
' 4.13135667732E-10'
write(*,fmt)'sea_kappa_t(0.035,300,1d5) =', sea_kappa_t_si(0.035d0, 300d0, 1d5), &
' 4.20250741344E-10'
write(*,fmt)'sea_lapserate(0.035,300,1d5) =', sea_lapserate_si(0.035d0, 300d0, 1d5), &
' 2.28411342567E-08'
write(*,fmt)'sea_osm_coeff(0.035,300,1d5) =', sea_osm_coeff_si(0.035d0, 300d0, 1d5), &
' 0.902777495349'
write(*,fmt)'sea_soundspeed(0.035,300,1d5) =', sea_soundspeed_si(0.035d0, 300d0, 1d5), &
' 1538.47940766'
write(*,fmt)'sea_temp_maxdensity(0.01,1d5) =', sea_temp_maxdensity_si(0.01d0, 1d5), &
' 274.950121503'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_3B'
write(*,fmt)'sea_h(0,0,0,0.035,500,1d5) =', sea_h_si(0,0,0,0.035d0,500d0,1d5), &
' 145481.970750'
write(*,fmt)'sea_h(1,0,0,0.035,500,1d5) =', sea_h_si(1,0,0,0.035d0,500d0,1d5), &
' 86860.7980486'
write(*,fmt)'sea_h(0,1,0,0.035,500,1d5) =', sea_h_si(0,1,0,0.035d0,500d0,1d5), &
' 309.557955853'
write(*,fmt)'sea_h(0,0,1,0.035,500,1d5) =', sea_h_si(0,0,1,0.035d0,500d0,1d5), &
' 9.81092930969E-04'
write(*,fmt)'sea_h(2,0,0,0.035,500,1d5) =', sea_h_si(2,0,0,0.035d0,500d0,1d5), &
' 2393730.16716'
write(*,fmt)'sea_h(1,1,0,0.035,500,1d5) =', sea_h_si(1,1,0,0.035d0,500d0,1d5), &
' 72.5298236488'
write(*,fmt)'sea_h(1,0,1,0.035,500,1d5) =', sea_h_si(1,0,1,0.035d0,500d0,1d5), &
'-6.84629317367E-04'
write(*,fmt)'sea_h(0,2,0,0.035,500,1d5) =', sea_h_si(0,2,0,0.035d0,500d0,1d5), &
' 7.72873234085E-02'
write(*,fmt)'sea_h(0,1,1,0.035,500,1d5) =', sea_h_si(0,1,1,0.035d0,500d0,1d5), &
' 2.86305358702E-08'
write(*,fmt)'sea_h(0,0,2,0.035,500,1d5) =', sea_h_si(0,0,2,0.035d0,500d0,1d5), &
'-3.96880481108E-13'

write(*,"(A)")' '
write(*,fmt)'sea_pottemp(0.035,300,1d7,1d5) =', sea_pottemp_si(0.035d0,300d0,1d7,1d5), &
' 299.771869405'
write(*,fmt)'sea_potdens(0.035,300,1d7,1d5) =', sea_potdensity_si(0.035d0,300d0,1d7,1d5), &
' 1022.71520130'
write(*,fmt)'sea_potenthalpy(0.035,300,1d7,1d5) =', sea_potenthalpy_si(0.035d0,300d0,1d7,1d5), &
' 106307.996083'
write(*,fmt)'sea_h_exp_t(0.035,500,1d5) =', sea_h_expansion_t_si(0.035d0, 500d0, 1d5), &
' 3.77581809091E-04'
write(*,fmt)'sea_h_exp_th(0.035,500,1d7,1d5) =', sea_h_expansion_theta_si(0.035d0, 500d0, 1d7, 1d5), &
' 3.84755380181E-04'
write(*,fmt)'sea_h_exp_h(0.035,500,1d7,1d5) =', sea_h_expansion_h_si(0.035d0, 500d0, 1d7, 1d5), &
' 9.60618615640E-08'
write(*,fmt)'sea_h_contr_t(0.035,500,1d5) =', sea_h_contraction_t_si(0.035d0, 500d0, 1d5), &
' 0.725209049049'
write(*,fmt)'sea_h_contr_th(0.035,500,1d7,1d5) =', sea_h_contraction_theta_si(0.035d0, 500d0, 1d7, 1d5), &
' 0.717342103505'
write(*,fmt)'sea_h_contr_h(0.035,500,1d7,1d5) =', sea_h_contraction_h_si(0.035d0, 500d0, 1d7, 1d5), &
' 0.697779873590'
write(*,fmt)'sea_temp(0.035,500,1d5) =', sea_temperature_si(0.035d0, 500d0, 1d5), &
' 309.557955853'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_3C'
write(*,fmt)'sea_eta_entropy(0.035, 1d5, 1d5) =', sea_eta_entropy_si(0.035d0, 1d5, 1d5), &
' 350.310622663'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_temp(0.035,1d5,1d7,1d5,h) =', sea_eta_temperature_si(0.035d0,1d5,1d7,1d5,key), &
' 295.985682129'
key  ='hpot'
write(*,fmt)'sea_eta_temp(0.035,1d5,1d7,1d5,hpot) =', sea_eta_temperature_si(0.035d0,1d5,1d7,1d5,key), &
' 298.413424848'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_pottemp(0.035,1d5,1d7,1d5,h) =', sea_eta_pottemp_si(0.035d0,1d5,1d7,1d5,key), &
' 295.782199115'
key  ='hpot'
write(*,fmt)'sea_eta_pottemp(0.035,1d5,1d7,1d5,hpot) =', sea_eta_pottemp_si(0.035d0,1d5,1d7,1d5,key), &
' 298.194955182'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_dens(0.035,1d5,1d7,1d5,h) =', sea_eta_density_si(0.035d0,1d5,1d7,1d5,key), &
' 1028.10986556'
key  ='hpot'
write(*,fmt)'sea_eta_dens(0.035,1d5,1d7,1d5,hpot) =', sea_eta_density_si(0.035d0,1d5,1d7,1d5,key), &
' 1027.36529798'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_potdens(0.035,1d5,1d7,1d5,h) =', sea_eta_potdensity_si(0.035d0,1d5,1d7,1d5,key), &
' 1023.91737474'
key  ='hpot'
write(*,fmt)'sea_eta_potdens(0.035,1d5,1d7,1d5,hpot) =', sea_eta_potdensity_si(0.035d0,1d5,1d7,1d5,key), &
' 1023.20527737'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_contr_t(0.035,1d5,1d7,1d5,h) =', sea_eta_contraction_t_si(0.035d0,1d5,1d7,1d5,key), &
' 0.728755239644'
key  ='hpot'
write(*,fmt)'sea_eta_contr_t(0.035,1d5,1d7,1d5,hpot) =', sea_eta_contraction_t_si(0.035d0,1d5,1d7,1d5,key), &
' 0.726349317428'
key  ='t'
write(*,fmt)'sea_eta_contr_t(0.035,300,1d7,1d5,t) =', sea_eta_contraction_t_si(0.035d0,300d0,1d7,1d5,key), &
' 0.724913833446'
key  ='tpot'
write(*,fmt)'sea_eta_contr_t(0.035,300,1d7,1d5,tpot) =', sea_eta_contraction_t_si(0.035d0,300d0,1d7,1d5,key), &
' 0.724714253918'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_contr_th(0.035,1d5,1d7,1d5,h) =', sea_eta_contraction_theta_si(0.035d0,1d5,1d7,1d5,key), &
' 0.728499505688'
key  ='hpot'
write(*,fmt)'sea_eta_contr_th(0.035,1d5,1d7,1d5,hpot) =', sea_eta_contraction_theta_si(0.035d0,1d5,1d7,1d5,key), &
' 0.726099732703'
key  ='t'
write(*,fmt)'sea_eta_contr_th(0.035,300,1d7,1d5,t) =', sea_eta_contraction_theta_si(0.035d0,300d0,1d7,1d5,key), &
' 0.724667977117'
key  ='tpot'
write(*,fmt)'sea_eta_contr_th(0.035,300,1d7,1d5,tpot) =', sea_eta_contraction_theta_si(0.035d0,300d0,1d7,1d5,key), &
' 0.724468894946'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_contr_h(0.035,1d5,1d7,1d5,h) =', sea_eta_contraction_h_si(0.035d0,1d5,1d7,1d5,key), &
' 0.718452125957'
key  ='hpot'
write(*,fmt)'sea_eta_contr_h(0.035,1d5,1d7,1d5,hpot) =', sea_eta_contraction_h_si(0.035d0,1d5,1d7,1d5,key), &
' 0.714531922616'
key  ='t'
write(*,fmt)'sea_eta_contr_h(0.035,300,1d7,1d5,t) =', sea_eta_contraction_h_si(0.035d0,300d0,1d7,1d5,key), &
' 0.712069013013'
key  ='tpot'
write(*,fmt)'sea_eta_contr_h(0.035,300,1d7,1d5,tpot) =', sea_eta_contraction_h_si(0.035d0,300d0,1d7,1d5,key), &
' 0.711718411190'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_exp_t(0.035,1d5,1d7,1d5,h) =', sea_eta_expansion_t_si(0.035d0,1d5,1d7,1d5,key), &
' 2.89480851145E-04'
key  ='hpot'
write(*,fmt)'sea_eta_exp_t(0.035,1d5,1d7,1d5,hpot) =', sea_eta_expansion_t_si(0.035d0,1d5,1d7,1d5,key), &
' 3.07242256461E-04'
key  ='t'
write(*,fmt)'sea_eta_exp_t(0.035,300,1d7,1d5,t) =', sea_eta_expansion_t_si(0.035d0,300d0,1d7,1d5,key), &
' 3.18513471410E-04'
key  ='tpot'
write(*,fmt)'sea_eta_exp_t(0.035,300,1d7,1d5,tpot) =', sea_eta_expansion_t_si(0.035d0,300d0,1d7,1d5,key), &
' 3.20122324740E-04'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_exp_th(0.035,1d5,1d7,1d5,h) =', sea_eta_expansion_theta_si(0.035d0,1d5,1d7,1d5,key), &
' 2.91293979902E-04'
key  ='hpot'
write(*,fmt)'sea_eta_exp_th(0.035,1d5,1d7,1d5,hpot) =', sea_eta_expansion_theta_si(0.035d0,1d5,1d7,1d5,key), &
' 3.09134848554E-04'
key  ='t'
write(*,fmt)'sea_eta_exp_th(0.035,300,1d7,1d5,t) =', sea_eta_expansion_theta_si(0.035d0,300d0,1d7,1d5,key), &
' 3.20454167783E-04'
key  ='tpot'
write(*,fmt)'sea_eta_exp_th(0.035,300,1d7,1d5,tpot) =', sea_eta_expansion_theta_si(0.035d0,300d0,1d7,1d5,key), &
' 3.22069710839E-04'

write(*,"(A)")' '
key  ='h'
write(*,fmt)'sea_eta_exp_h(0.035,1d5,1d7,1d5,h) =', sea_eta_expansion_h_si(0.035d0,1d5,1d7,1d5,key), &
' 7.28514646021E-08'
key  ='hpot'
write(*,fmt)'sea_eta_exp_h(0.035,1d5,1d7,1d5,hpot) =', sea_eta_expansion_h_si(0.035d0,1d5,1d7,1d5,key), &
' 7.72876772245E-08'

key  ='t'
write(*,fmt)'sea_eta_exp_h(0.035,300,1d7,1d5,t) =', sea_eta_expansion_h_si(0.035d0,300d0,1d7,1d5,key), &
' 8.01009066333E-08'
key  ='tpot'
write(*,fmt)'sea_eta_exp_h(0.035,300,1d7,1d5,tpot) =', sea_eta_expansion_h_si(0.035d0,300d0,1d7,1d5,key), &
' 8.05023387611E-08'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_3D'
write(*,fmt)'sea_sa(273.15,101325,1028) =', sea_sa_si(273.15d0, 101325d0, 1028d0), &
' 3.50315257709E-02'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SEA_5A'
write(*,fmt)'sea_alpha_ct_si(0.035,300,1E8) =', sea_alpha_ct_si(0.035d0, 300d0, 1d8), &
' 3.91772847589E-04'
write(*,fmt)'sea_alpha_pt0_si(0.035,300,1E8) =', sea_alpha_pt0_si(0.035d0, 300d0, 1d8), &
' 3.92515634064E-04'
write(*,fmt)'sea_alpha_t_si(0.035,300,1E8) =', sea_alpha_t_si(0.035d0, 300d0, 1d8), &
' 3.73608885178E-04'
write(*,fmt)'sea_beta_ct_si(0.035,300,1E8) =', sea_beta_ct_si(0.035d0, 300d0, 1d8), &
' 0.649596383654'
write(*,fmt)'sea_beta_pt0_si(0.035,300,1E8) =', sea_beta_pt0_si(0.035d0, 300d0, 1d8), &
' 0.663973579411'
write(*,fmt)'sea_beta_t_si(0.035,300,1E8) =', sea_beta_t_si(0.035d0, 300d0, 1d8), &
' 0.666238827368'
write(*,fmt)'sea_cabb_ct_si(0.035,300,1d5) =', sea_cabb_ct_si(0.035d0, 300d0, 1d5), &
' 8.612525|67438267E-06'
write(*,fmt)'sea_cabb_pt0_si(0.035,300,1d5) =', sea_cabb_pt0_si(0.035d0, 300d0, 1d5), &
' 8.3387453|7690444E-06'
write(*,fmt)'sea_thrmb_ct_si(0.035,300,1d5) =', sea_thrmb_ct_si(0.035d0, 300d0, 1d5), &
' 1.4810927|1668362E-12'
write(*,fmt)'sea_thrmb_pt0_si(0.035,300,1d5) =', sea_thrmb_pt0_si(0.035d0, 300d0, 1d5), &
' 1.4594101|0702991E-12'
write(*,fmt)'sea_ctmp_from_ptmp0_si(0.035,300) =', sea_ctmp_from_ptmp0_si(0.035d0, 300d0), &
' 300.010069445E+03'
write(*,fmt)'sea_ptmp0_from_ctmp_si(0.035,300.010069445d0) =', sea_ptmp0_from_ctmp_si(0.035d0, 300.010069445d0), &
' 300'

end subroutine

end module values_sea
