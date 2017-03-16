module ice_air_4c

!#########################################################################

!THIS MODULE IMPLEMENTS THE ENTHALPY OF ICE AIR, AS WELL AS ITS PARTIAL 
!DERIVATIVES, DEPENDING ON DRY-AIR FRACTION, ENTROPY AND PRESSURE.

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
!     ICE_LIQ_4,    FILE ICE_LIQ_4.F90
!     ICE_AIR_4B,   FILE ICE_AIR_4B.F90

!NESTED USE ASSOCIATIONS
!     MATHS_0,      FILE MATHS_0.F90
!     AIR_1,        FILE AIR_1.F90
!     FLU_1,        FILE FLU_1.F90
!     ICE_1,        FILE ICE_1.F90
!     AIR_2,        FILE AIR_2.F90
!     FLU_2,        FILE FLU_2.F90
!     ICE_2,        FILE ICE_2.F90
!     AIR_3A,       FILE AIR_3A.F90
!     AIR_3B,       FILE AIR_3B.F90
!     ICE_VAP_4,    FILE ICE_VAP_4.F90
!     ICE_AIR_4A,   FILE ICE_AIR_4A.F90

!#########################################################################

use constants_0
use convert_0
use ice_liq_4
use ice_air_4b

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: ice_air_h_si, ice_air_h_cp_si, ice_air_h_density_si, &
          ice_air_h_kappa_s_si, ice_air_h_lapserate_si, &
          ice_air_h_temperature_si, ice_air_potdensity_si, &
          ice_air_potenthalpy_si, ice_air_pottemp_si, &
          set_it_ctrl_ice_air_pottemp

!CONTROL PARAMETERS OF THE TEMPERATURE ITERATION
integer ctrl_initialized

integer ctrl_mode_pottemp
integer ctrl_loop_maximum

real*8 ctrl_init_pottemp
real*8 ctrl_eps_exit_pottemp

contains

!==========================================================================
function ice_air_h_si(drv_wa, drv_eta, drv_p, wa_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS ENTHALPY OF ICE AIR AS A THERMODYNAMIC POTENTIAL,
!DEPENDING ON AIR MASS FRACTION, ENTROPY AND PRESSURE

!RETURNS ICE_AIR_H_SI AS THE WA-ETA-P DERIVATIVE
!(D/DWA)^DRV_WA (D/DETA)^DRV_ETA (D/DP)^DRV_P H(WA,ETA,P)
!OF THE SPECIFIC ENTHALPY OF ICE AIR, H(WA,ETA,P), IN J/KG

!OUTPUT:
!ICE_AIR_H_SI POTENTIAL RESULTS

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K),
!P_SI      ABSOLUTE PRESSURE IN PA

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A, AND ON THE ITERATION SETTINGS
!      FOR TEMPERATURE BY SET_IT_CTRL_ICE_AIR_POTTEMP OF THIS MODULE

real*8 ice_air_h_si, wa_si, eta_si, p_si
real*8 t

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!ice_air_h_si( 0, 0, 0, 0.5, -600, 1E5) =-164588.277175
!ice_air_h_si( 1, 0, 0, 0.5, -600, 1E5) = 543.382787268
!ice_air_h_si( 0, 1, 0, 0.5, -600, 1E5) = 271.449435548
!ice_air_h_si( 0, 0, 1, 0.5, -600, 1E5) = 0.392073699416
!ice_air_h_si( 2, 0, 0, 0.5, -600, 1E5) = 224938.714580
!ice_air_h_si( 1, 1, 0, 0.5, -600, 1E5) =-177.440249223
!ice_air_h_si( 1, 0, 1, 0.5, -600, 1E5) = 0.781966303025
!ice_air_h_si( 0, 2, 0, 0.5, -600, 1E5) = 0.139971645623
!ice_air_h_si( 0, 1, 1, 0.5, -600, 1E5) = 2.26942461987E-04
!ice_air_h_si( 0, 0, 2, 0.5, -600, 1E5) =-3.57055653968E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!ice_air_h_si( 0, 0, 0, 0.5, -600, 1E5) =-164588.106002
!ice_air_h_si( 1, 0, 0, 0.5, -600, 1E5) = 543.016638396
!ice_air_h_si( 0, 1, 0, 0.5, -600, 1E5) = 271.449994437
!ice_air_h_si( 0, 0, 1, 0.5, -600, 1E5) = 0.391981878510
!ice_air_h_si( 2, 0, 0, 0.5, -600, 1E5) = 224958.525864
!ice_air_h_si( 1, 1, 0, 0.5, -600, 1E5) =-177.457078495
!ice_air_h_si( 1, 0, 1, 0.5, -600, 1E5) = 0.781782661019
!ice_air_h_si( 0, 2, 0, 0.5, -600, 1E5) = 0.139985868894
!ice_air_h_si( 0, 1, 1, 0.5, -600, 1E5) = 2.26912930199E-04
!ice_air_h_si( 0, 0, 2, 0.5, -600, 1E5) =-3.56976697603E-06

integer drv_wa, drv_eta, drv_p

ice_air_h_si = errorreturn

if(drv_wa < 0) return
if(drv_eta < 0) return
if(drv_p < 0) return
if(drv_wa + drv_eta + drv_p > 2) return

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

!COMPUTE TEMPERATURE FROM ENTROPY.
!IF P_SI = IN-SITU PRESSURE, THIS IS IN-SITU TEMPERATURE
!IF P_SI = REFERENCE PRESSURE, THIS IS POTENTIAL TEMPERATURE
t = ice_air_h_temperature_si(wa_si, eta_si, p_si)
if(t == errorreturn) return

ice_air_h_si = ice_air_a_eta_p_derivatives_si(drv_wa, drv_eta, drv_p, wa_si, t, p_si)

end function

!=========================================================================
function ice_air_h_cp_si(wa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!CP(A,ETA,P) = T/(D2H/DETA2)_A_P  HEAT CAPACITY OF ICE AIR IN J/(KG K)

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_cp_si(0.9, -100, 1e5) = 1766.84429536

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_h_cp_si(0.9, -100, 1e5) = 1766.52051488

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 ice_air_h_cp_si, wa_si, eta_si, p_si
real*8 t

ice_air_h_cp_si = errorreturn

t = ice_air_h_temperature_si(wa_si, eta_si, p_si)

if(t == errorreturn) return
if(t <= 0d0) return

ice_air_h_cp_si = ice_air_g_cp_si(wa_si, t, p_si)

end function

!=========================================================================
function ice_air_h_density_si(wa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!D(A,ETA,P) = 1/(DH/DP)_A_ETA DENSITY OF ICE AIR IN KG/M3

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_density_si(0.9, -100, 1e5) = 1.42498657691

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_h_density_si(0.9, -100, 1e5) = 1.4253189599

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 ice_air_h_density_si, wa_si, eta_si, p_si
real*8 t

ice_air_h_density_si = errorreturn

t = ice_air_h_temperature_si(wa_si, eta_si, p_si)
if(t == errorreturn) return
if(t <= 0d0) return

ice_air_h_density_si = ice_air_g_density_si(wa_si, t, p_si)

end function

!=========================================================================
function ice_air_h_kappa_s_si(wa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!KAPPA(A,ETA,P) = - (D2H_DP2)_A_ETA/(DH_P)_A_ETA ADIABATIC COMPRESSIBILITY OF ICE AIR IN 1/PA

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_kappa_s_si(0.9, -100, 1e5) = 8.23024116522E-06

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_h_kappa_s_si(0.9, -100, 1e5) = 8.23031581047E-06

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 ice_air_h_kappa_s_si, wa_si, eta_si, p_si
real*8 d, c

ice_air_h_kappa_s_si = errorreturn

d = ice_air_h_density_si(wa_si, eta_si, p_si)
if(d == errorreturn) return
if(d <= 0d0) return

c = ice_air_h_si(0, 0, 2, wa_si, eta_si, p_si)
if(c == errorreturn) return
if(c >= 0d0) return

ice_air_h_kappa_s_si = -d * c

end function

!=========================================================================
function ice_air_h_lapserate_si(wa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!GAMMA(A,ETA,P) = (D2H/DETA DP)_A  ICE-ADIABATIC LAPSE RATE IN K/PA

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_lapserate_si(0.9, -100, 1e5) = 4.42476802729E-04

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_lapserate_si(0.9, -100, 1e5) = 4.42457786755E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 ice_air_h_lapserate_si, wa_si, eta_si, p_si

ice_air_h_lapserate_si = ice_air_h_si(0, 1, 1, wa_si, eta_si, p_si)

end function

!=========================================================================
function ice_air_h_temperature_si(wa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!T(A,ETA,P) (POTENTIAL) TEMPERATURE OF ICE AIR IN K

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_h_temperature_si(0.9, -100, 1e5) = 270.382806934

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_h_temperature_si(0.9, -100, 1e5) = 270.383680119

!THIS IS THE INVERSE FUNCTION TO
!    ETA_SI = ICE_AIR_G_ENTROPY_SI(WA_SI, T_SI, P_SI)
!IN ICE_AIR_4B

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 ice_air_h_temperature_si, wa_si, eta_si, p_si
real*8 t, eps
integer maxit

ice_air_h_temperature_si = errorreturn

call init_it_ctrl_pottemp

select case(ctrl_mode_pottemp)
  case(0);   t = aux_ice_air_temperature_si(wa_si, eta_si, p_si)
  case default; t = ctrl_init_pottemp
end select

if(t <= 0d0) return
if(t == errorreturn) return

  if(ctrl_loop_maximum == 0) then
     maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
     ice_air_h_temperature_si = t
     return
  elseif(ctrl_loop_maximum > 0) then
     maxit = ctrl_loop_maximum
  else
     return
  endif

eps = ctrl_eps_exit_pottemp
if(eps <= 0d0) return

!RUN THE ITERATION
ice_air_h_temperature_si = pottemp_iteration(wa_si, eta_si, p_si, t, maxit, eps)

end function

!==========================================================================
function ice_air_potdensity_si(a_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL DENSITY OF ICE AIR

!OUTPUT:
!RHO_THETA(A,T,P,PR) POTENTIAL DENSITY OF ICE AIR IN KG/M^3

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_potdensity_si(0.9, 230, 5e4, 1e5) = 1.45013106603

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_potdensity_si(0.9, 230, 5e4, 1e5) = 1.45048110422

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A,
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_ICE_AIR_POTTEMP OF THIS MODULE

real*8 ice_air_potdensity_si, a_si, t_si, p_si, pr_si
real*8 s, v

ice_air_potdensity_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = ice_air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

v = ice_air_h_si(0, 0, 1, a_si, s, pr_si)   !SPECIFIC VOLUME AT REFERENCE PRESSURE
if(v == errorreturn) return
if(v <= 0d0) return

ice_air_potdensity_si = 1d0 / v

end function

!==========================================================================
function ice_air_potenthalpy_si(a_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL ENTHALPY OF ICE AIR

!OUTPUT:
!H_THETA(A,T,P,PR) POTENTIAL ENTHALPY OF ICE AIR IN J/KG

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_potenthalpy_si(0.9, 230, 5e4, 1e5) = -35779.7348519

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_potenthalpy_si(0.9, 230, 5e4, 1e5) = -35781.2564451

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A,
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_ICE_AIR_POTTEMP OF THIS MODULE

real*8 ice_air_potenthalpy_si, a_si, t_si, p_si, pr_si
real*8 s

ice_air_potenthalpy_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = ice_air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

ice_air_potenthalpy_si = ice_air_h_si(0, 0, 0, a_si, s, pr_si)

end function

!==========================================================================
function ice_air_pottemp_si(a_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL TEMPERATURE OF ICE AIR

!OUTPUT:
!THETA(A,T,P,PR) ABSOLUTE POTENTIAL TEMPERATURE OF ICE AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_pottemp_si(0.9, 230, 5E4, 1E5) = 266.106323310

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_pottemp_si(0.9, 230, 5E4, 1E5) = 266.105208871

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A,
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_ICE_AIR_POTTEMP OF THIS MODULE_SI

real*8 ice_air_pottemp_si, a_si, t_si, p_si, pr_si
real*8 s

ice_air_pottemp_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

if(p_si == pr_si) then
  ice_air_pottemp_si = t_si
  return
end if

s = ice_air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

ice_air_pottemp_si = ice_air_h_temperature_si(a_si, s, pr_si)

end function

!=========================================================================
function pottemp_iteration(wa_si, eta_si, p_si, t_si, maxit, eps)
!==========================================================================

!OUTPUT:
!THETA =  ABSOLUTE POTENTIAL TEMPERATURE OF ICE AIR IN K,
!         I.E. THE TEMPERATURE THAT SOLVES THE EQUATION

!INPUTS:
!          ETA_SI = ICE_AIR_G_ENTROPY_SI(WA_SI, THETA, P_SI)
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA
!T_SI      ABSOLUTE (POTENTIAL) TEMPERATURE IN K, INITIAL VALUE
!MAXIT     MAX. NUMBER OF ITERATIONS
!EPS       REQUIRED TOLERANCE IN K

real*8 pottemp_iteration, wa_si, eta_si, p_si, t_si, eps
real*8 s, theta, cp, dt
integer i, maxit

pottemp_iteration = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0 .or. &
  eps <= 0d0 .or. &
  maxit <= 0) return

if(check_limits == 1) then
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) return
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) return
else
  if(t_si <= 0d0 .or. &
     p_si <= 0d0) return
endif

check_limits = check_limits - 1

theta = t_si

do i = 1, maxit

  !GET ENTROPY AND ITS FIRST DERIVATIVE FOR NEWTON ITERATION
  s = ice_air_g_entropy_si(wa_si, theta, p_si)
  cp = ice_air_g_cp_si(wa_si, theta, p_si)
    
  if(s == errorreturn) go to 999
  if(cp == errorreturn) go to 999
  if(cp <= 0d0) go to 999
  
  !NEXT TEMPERATURE IMPROVEMENT STEP
  dt = theta * (eta_si - s) / cp
  theta = theta + dt
  if(theta <= 0d0) go to 999
    
  if(abs(dt) < eps) then
    pottemp_iteration = theta
    go to 999
  end if
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !ICE_LIMITS
  if(theta <= ice_tmin .or. theta > ice_tmax) then
    pottemp_iteration = errorreturn
    return
  endif
  !AIR_LIMITS
  if(theta < dry_air_tmin .or. theta > dry_air_tmax) then
    pottemp_iteration = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_air_a_eta_p_derivatives_si(drv_wa, drv_eta, drv_p, wa_si, t_si, p_si)
!==========================================================================
                                              
!THIS FUNCTION COMPUTES ICE-AIR WA-ETA-P DERIVATIVES OF H FROM WA-T-P DERIVATIVES OF G

real*8 ice_air_a_eta_p_derivatives_si, wa_si, t_si, p_si
real*8 g, g_t, g_p, g_a
real*8 g_aa, g_at, g_ap
real*8 g_tt, g_pp, g_tp

real*8 h
 
integer drv_wa, drv_eta, drv_p

!IN ONE CASE WE DO NOT AT ALL NEED TO COMPUTE THE GIBBS FUNCTION:
if(drv_wa == 0 .and. drv_eta == 1 .and. drv_p == 0) then
  ice_air_a_eta_p_derivatives_si = t_si
  return
end if

ice_air_a_eta_p_derivatives_si = errorreturn
g_tp = errorreturn

select case(drv_wa)

  case(0)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0) 
            g = ice_air_g_si(0, 0, 0, wa_si, t_si, p_si)
            if(g == errorreturn) return
            g_t = ice_air_g_si(0, 1, 0, wa_si, t_si, p_si)
            if(g_t == errorreturn) return
            h = g - t_si * g_t                          !H
          
          case(1)
            g_p = ice_air_g_si(0, 0, 1, wa_si, t_si, p_si)
            if(g_p == errorreturn) return
            h = g_p                                      !DH/DP
                  
          case(2)
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = ice_air_g_si(0, 1, 1, wa_si, t_si, p_si)
            if(g_tp == errorreturn) return
            g_pp = ice_air_g_si(0, 0, 2, wa_si, t_si, p_si)
            if(g_pp == errorreturn) return
            h = (g_tt * g_pp - g_tp**2) / g_tt           !D2H/DP2
                  
          case default; return
        end select
        
      case(1)
        select case(drv_p)
          case(0)
            h = t_si                                     !DH/DETA, HAS ALREADY BEEN HANDLED INITIALLY
          
          case(1)
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = ice_air_g_si(0, 1, 1, wa_si, t_si, p_si)
            if(g_tp == errorreturn) return
            h = -g_tp / g_tt                             !D2H/DETADP
                  
          case default; return
        end select

      case(2)
        select case(drv_p)
          case(0) 
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            h = -1d0 / g_tt                              !D2H/DETA2
                  
          case default; return
        end select

      case default; return
    end select
    
  case(1)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0)
            g_a = ice_air_g_si(1, 0, 0, wa_si, t_si, p_si)
            if(g_a == errorreturn) return
            h = g_a                                      !DH/DA
          
          case(1)
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = ice_air_g_si(1, 1, 0, wa_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_ap = ice_air_g_si(1, 0, 1, wa_si, t_si, p_si)
            if(g_ap == errorreturn) return
            h = (g_tt * g_ap - g_at * g_tp) / g_tt       !D2H/DADP
          
          case default; return
        end select

      case(1)
        select case(drv_p)
        
          case(0)
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = ice_air_g_si(1, 1, 0, wa_si, t_si, p_si)
            if(g_at == errorreturn) return
            h = -g_at / g_tt                             !D2H/DADETA
          
          case default; return
        end select

      case default; return
    end select

  
  case(2)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0) 
            g_tt = ice_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = ice_air_g_si(1, 1, 0, wa_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_aa = ice_air_g_si(2, 0, 0, wa_si, t_si, p_si)
            if(g_aa == errorreturn) return
            h = (g_tt * g_aa - g_at**2) / g_tt           !D2H/DA2
          
          case default; return
        end select

      case default; return
    end select

  case default; return
end select

ice_air_a_eta_p_derivatives_si = h

end function

!=========================================================================
function aux_ice_air_temperature_si(wa_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION ESTIMATES THE TEMPERATURE OF ICE AIR FROM ITS DRY-AIR FRACTION WA_SI IN KG/KG,
!ITS SPECIFIC ENTROPY ETA_SI IN J/(KG K), AND ITS ABSOLUTE PRESSURE, P_SI IN PA

real*8 aux_ice_air_temperature_si, wa_si, eta_si, p_si
real*8 numer, denom, t, tf

real*8 tt, pt, cpa, cpi, sat, sit, ra

tt = tp_temperature_si           !TRIPLE-POINT TEMPERATURE IN K
pt = tp_pressure_iapws95_si      !TRIPLE-POINT PRESSURE IN Pa

cpa = 1003.68997553091d0         !TRIPLE-POINT HEAT CAPACITY OF AIR IN J/(KG K)
cpi = 2096.78431621667d0         !TRIPLE-POINT HEAT CAPACITY OF ICE IN J/(KG K)

sat = 1467.66694249983d0         !TRIPLE-POINT ENTROPY OF AIR IN J/(KG K)
sit = -1220.69433939648d0        !TRIPLE-POINT ENTROPY OF ICE IN J/(KG K)

ra = gas_constant_air_si

aux_ice_air_temperature_si = errorreturn
if(p_si <= 0d0) return
if(wa_si < 0d0 .or. wa_si > 1d0) return

!ASSUME CONSTANT HEAT CAPACITY CP OF ICE AND AIR, USE IDEAL-GAS EQUATION,
!NEGLECT THE VAPOUR FRACTION IN HUMID AIR:

numer = eta_si - sit - wa_si * (sat - sit - ra * log(p_si / pt))

denom = cpi + wa_si * (cpa - cpi)
if(denom == 0d0) return

t = tt * exp(numer / denom)

!IF ESTIMATE IS ABOVE THE FREEZING POINT, RATHER TAKE THE FREEZING POINT
tf = ice_liq_meltingtemperature_si(p_si)
if(tf /= errorreturn .and. tf > 0) then
  if(t > tf) t = tf
end if

aux_ice_air_temperature_si = t

end function

!==========================================================================
subroutine init_it_ctrl_pottemp()

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR TEMPERATURE ITERATION
ctrl_loop_maximum = 100
ctrl_mode_pottemp = 0           !DEFAULT: THETA = AUX_ICE_AIR_TEMPERATURE_SI
ctrl_init_pottemp = 273.15d0
ctrl_eps_exit_pottemp = 0.0001d0  !DEFAULT = 0.1 MK

end subroutine

!==========================================================================
subroutine set_it_ctrl_ice_air_pottemp(key, value)

!THIS SUBROUTINE SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION 
!USED TO COMPUTE POTENTIAL TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0    SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_THETA      0           USE DEFAULT POTENTIAL TEMPERATURE THETA = AUX_TEMPERATURE_MIN_SI TO START
!INIT_THETA      T > 0    USE VALUE T AS POTENTIAL TEMPERATURE TO START

!TOL_THETA       0           USE DEFAULT EXIT ACCURACY FOR POTENTIAL TEMPERATURE (0.1 MK)
!TOL_THETA       EPS > 0  USE EPS AS EXIT ACCURACY FOR POTENTIAL TEMPERATURE

real*8 value
character*13 key

call init_it_ctrl_pottemp

select case (trim(key))

  case ("it_steps");   !ITERATION STEPS
      if(value == 0d0)then
         ctrl_loop_maximum = 100      !DEFAULT = 100
      elseif(value < 0d0) then
         ctrl_loop_maximum = -1
      else
         ctrl_loop_maximum = nint(value)
      endif

  case ("init_theta");   !START THETA
      if(value == 0d0) then
         ctrl_mode_pottemp = 0        !DEFAULT: AUX_TEMPERATURE_MIN_SI
      elseif(value > 0d0) then
         ctrl_mode_pottemp = 1
         ctrl_init_pottemp = value
      endif

  case ("tol_theta");                 !REQUIRED THETA TOLERANCE
      if(value == 0d0) then
         ctrl_eps_exit_pottemp = 0.0001d0 !DEFAULT = 0.1 MK
      elseif(value > 0d0) then
         ctrl_eps_exit_pottemp = value
      endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

end module ice_air_4c
