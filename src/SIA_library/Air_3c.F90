module air_3c

!#########################################################################

!THIS MODULE IMPLEMENTS THE ENTHALPY OF HUMID AIR DEPENDING ON AIR FRACTION, 
!ENTROPY AND PRESSURE, AS WELL AS ITS PARTIAL DERIVATIVES.

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
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     CONVERT_0,    FILE CONVERT_0.F90
!     AIR_2,        FILE AIR_2.F90
!     AIR_3A,       FILE AIR_3A.F90
!     AIR_3B,       FILE AIR_3B.F90

!NESTED USE ASSOCIATIONS
!     MATHS_0,      FILE MATHS_0.F90
!     AIR_1,        FILE AIR_1.F90
!     FLU_1,        FILE FLU_1.F90

!#########################################################################

use constants_0
use convert_0
use air_2
use air_3a
use air_3b

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: air_h_si, air_potdensity_si, air_potenthalpy_si, &
          air_pottemp_si, air_temperature_si, set_it_ctrl_air_pottemp

!CONTROL PARAMETERS OF THE TEMPERATURE ITERATION
integer ctrl_initialized, ctrl_mode_pottemp
integer ctrl_loop_maximum
real*8 ctrl_init_pottemp, ctrl_eps_exit_pottemp

contains

!==========================================================================
function air_h_si(drv_a, drv_eta, drv_p, a_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS ENTHALPY AS A THERMODYNAMIC POTENTIAL OF HUMID AIR, DEPENDING
!ON THE DRY-AIR MASS FRACTION, SPECIFIC ENTROPY AND PRESSURE

!OUTPUT:
!AIR_H_SI AS THE A-ETA-P DERIVATIVE
!(D/DA)^DRV_A (D/DETA)^DRV_ETA (D/DP)^DRV_P H(A,ETA,P)
!OF THE SPECIFIC ENTHALPY OF HUMID AIR, H(A,ETA,P), IN J/KG

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K),
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!air_h_si( 0, 0, 0, 0.9, 900, 1E5) = 274589.199192
!air_h_si( 1, 0, 0, 0.9, 900, 1E5) =-229734.926517
!air_h_si( 0, 1, 0, 0.9, 900, 1E5) = 297.395859294
!air_h_si( 0, 0, 1, 0.9, 900, 1E5) = 0.903422384672
!air_h_si( 2, 0, 0, 0.9, 900, 1E5) = 13959323.1016
!air_h_si( 1, 1, 0, 0.9, 900, 1E5) = 1676.33822835
!air_h_si( 1, 0, 1, 0.9, 900, 1E5) =-0.467305898316
!air_h_si( 0, 2, 0, 0.9, 900, 1E5) = 0.223618143111
!air_h_si( 0, 1, 1, 0.9, 900, 1E5) = 7.15646251375E-04
!air_h_si( 0, 0, 2, 0.9, 900, 1E5) =-6.78234488400E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!air_h_si( 0, 0, 0, 0.9, 900, 1E5) = 274592.611783
!air_h_si( 1, 0, 0, 0.9, 900, 1E5) =-229706.905765
!air_h_si( 0, 1, 0, 0.9, 900, 1E5) = 297.403043058
!air_h_si( 0, 0, 1, 0.9, 900, 1E5) = 0.903262695636
!air_h_si( 2, 0, 0, 0.9, 900, 1E5) = 13963273.0104
!air_h_si( 1, 1, 0, 0.9, 900, 1E5) = 1676.85098552
!air_h_si( 1, 0, 1, 0.9, 900, 1E5) =-0.467537679038
!air_h_si( 0, 2, 0, 0.9, 900, 1E5) = 0.223684689765
!air_h_si( 0, 1, 1, 0.9, 900, 1E5) = 7.15703143992E-04
!air_h_si( 0, 0, 2, 0.9, 900, 1E5) =-6.78105152859E-06


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_AIR_POTTEMP OF THIS MODULE

real*8 air_h_si, a_si, eta_si, p_si
real*8 t
integer drv_a, drv_eta, drv_p

air_h_si = errorreturn

if(drv_a < 0 .or. drv_a > 2) return
if(drv_eta < 0 .or. drv_eta > 2) return
if(drv_p < 0 .or. drv_p > 2) return

if(a_si < 0d0 .or. a_si >= 1d0) return
if(p_si < 0d0) return

!COMPUTE TEMPERATURE FROM ENTROPY.
!IF P_SI = IN-SITU PRESSURE, THIS IS IN-SITU TEMPERATURE
!IF P_SI = REFERENCE PRESSURE, THIS IS POTENTIAL TEMPERATURE
t = air_temperature_si(a_si, eta_si, p_si)
if(t == errorreturn) return
if(t <= 0d0) return

air_h_si = air_a_eta_p_derivatives_si(drv_a, drv_eta, drv_p, a_si, t, p_si)

end function

!==========================================================================
function air_potdensity_si(a_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE POTENTIAL DENSITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE, PRESSURE AND REFERENCE PRESSURE

!OUTPUT:
!RHO_THETA(A,T,P,PR) POTENTIAL DENSITY OF HUMID AIR IN KG/M3

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_POTDENSITY_SI(0.9, 300, 5e4, 1E5) = 0.903326577188

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_POTDENSITY_SI(0.9, 300, 5e4, 1E5) = 0.903509489711

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_AIR_POTTEMP OF THIS MODULE

real*8 air_potdensity_si, a_si, t_si, p_si, pr_si
real*8 s, v

air_potdensity_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

v = air_h_si(0, 0, 1, a_si, s, pr_si) !SPECIFIC VOLUME AT REFERENCE PRESSURE
if(v == errorreturn) return
if(v <= 0d0) return

air_potdensity_si = 1d0 / v

end function

!==========================================================================
function air_potenthalpy_si(a_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE POTENTIAL ENTHALPY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE, PRESSURE AND REFERENCE PRESSURE

!OUTPUT:
!H_THETA(A,T,P,PR) POTENTIAL ENTHALPY OF HUMID AIR IN J/KG

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_POTENTHALPY_SI(0.9, 300, 5e4, 1E5) = 348892.581997

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_POTENTHALPY_SI(0.9, 300, 5e4, 1E5) = 348872.568665

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_AIR_POTTEMP OF THIS MODULE

real*8 air_potenthalpy_si, a_si, t_si, p_si, pr_si
real*8  s

air_potenthalpy_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

air_potenthalpy_si = air_h_si(0, 0, 0, a_si, s, pr_si)

end function

!==========================================================================
function air_pottemp_si(a_si, t_si, p_si, pr_si) 
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE ABSOLUTE POTENTIAL TEMPERATURE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE, PRESSURE AND REFERENCE PRESSURE

!OUTPUT:
!THETA(A,T,P,PR) ABSOLUTE POTENTIAL TEMPERATURE OF HUMID AIR IN K

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_POTTEMP_SI(0.9, 300, 5E4, 1e5) = 363.654373179

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_POTTEMP_SI(0.9, 300, 5E4, 1e5) = 363.653905688

real*8 air_pottemp_si, a_si, t_si, p_si, pr_si
real*8 s

air_pottemp_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

if(p_si == pr_si) then
  air_pottemp_si = t_si
  return
end if

s = air_g_entropy_si(a_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

air_pottemp_si = air_temperature_si(a_si, s, pr_si)

end function

!=========================================================================
function air_temperature_si(a_si, eta_si, p_si)
!==========================================================================

!THIS IS THE INVERSE FUNCTION TO AIR_G_ENTROPY_SI(A_SI, T_SI, P_SI) IN AIR_3B
!COMPUTING THE ABSOLUTE (POTENTIAL) TEMPERATURE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, SPECIFIC ENTROPY AND PRESSURE

!OUTPUT:
!T(A,ETA,P) (POTENTIAL) TEMPERATURE OF HUMID AIR IN K

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!AIR_TEMPERATURE_SI(0.9, 900, 1E5) = 297.395859294

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!AIR_TEMPERATURE_SI(0.9, 900, 1E5) = 297.403043058

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 air_temperature_si, a_si, eta_si, p_si
real*8 t, eps
integer maxit

air_temperature_si = errorreturn

call init_it_ctrl_pottemp

select case(ctrl_mode_pottemp)
  case(0);    t = aux_temperature_si(a_si, eta_si, p_si)  !IDEAL-GAS APPROXIMATION
  case default; t = ctrl_init_pottemp
end select

  if(ctrl_loop_maximum == 0d0) then
     maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
     air_temperature_si = t
     return
  elseif(ctrl_loop_maximum > 0) then
     maxit = ctrl_loop_maximum
  else
     return
  endif

eps = ctrl_eps_exit_pottemp
if(eps <= 0d0) return

!RUN THE ITERATION
air_temperature_si = pottemp_iteration(a_si, eta_si, p_si, t, maxit, eps)

end function

!==========================================================================
function air_a_eta_p_derivatives_si(drv_a, drv_eta, drv_p, a_si, t_si, p_si)
!==========================================================================
                                              
!THIS FUNCTION COMPUTES HUMID AIR A-ETA-P DERIVATIVES OF H FROM A-T-P DERIVATIVES OF G
!THIS FUNCTION COMPUTES A-ETA-P DERIVATIVES OF H FROM A-T-P DERIVATIVES OF G
!AT GIVEN AIR FRACTION, TEMPERATURE AND PRESSURE

!OUTPUT:
!AIR_A_ETA_P_DERIVATEVE_SI = SPECIFIC ENTHAPY OR ITS DERIVATIVES IN J KG-1 (J KG-1 K-1)^-DRV_ETA PA^-DRV_P

!INPUT:
!DRV_A = ORDER OF AIR-FRACTION DERIVATIVE, 0 <= DRV_A <= 2
!DRV_ETA = ORDER OF SPECIFIC ENTROPY DERIVATIVE, 0 <= DRV_ETA + DRV_A <= 2
!DRV_P = ORDER OF PRESSURE DERIVATIVE, 0 <= DRV_P + DRV_ETA + DRV_A <= 2
!A_SI  = MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI  = ABSOLUTE TEMPERATURE IN K
!P_SI  = PRESSURE IN PA

real*8 air_a_eta_p_derivatives_si, a_si, t_si, p_si

real*8 g, g_t, g_p, g_a
real*8  g_aa, g_at, g_ap
real*8  g_tt, g_pp, g_tp

real*8  h

integer drv_a, drv_eta, drv_p

!IN ONE CASE WE DO NOT AT ALL NEED TO COMPUTE THE GIBBS FUNCTION:
if(drv_a == 0 .and. drv_eta == 1 .and. drv_p == 0) then
  air_a_eta_p_derivatives_si = t_si
  return
end if

air_a_eta_p_derivatives_si = errorreturn

g_tp = errorreturn

select case(drv_a)

  case(0)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0) 
            g = air_g_si(0, 0, 0, a_si, t_si, p_si)
            if(g == errorreturn) return
            g_t = air_g_si(0, 1, 0, a_si, t_si, p_si)
            if(g_t == errorreturn) return
            h = g - t_si * g_t                             !H
          
          case(1)
            g_p = air_g_si(0, 0, 1, a_si, t_si, p_si)
            if(g_p == errorreturn) return
            h = g_p                                        !DH/DP
             
          case(2)
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = air_g_si(0, 1, 1, a_si, t_si, p_si)
            if(g_tp == errorreturn) return
            g_pp = air_g_si(0, 0, 2, a_si, t_si, p_si)
            if(g_pp == errorreturn) return
            h = (g_tt * g_pp - g_tp**2) / g_tt             !D2H/DP2
                  
          case default; return
        end select
        
      case(1)
        select case(drv_p)

          case(0) 
            h = t_si                                       !DH/DETA, HAS ALREADY BEEN HANDLED INITIALLY
          
          case(1) 
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = air_g_si(0, 1, 1, a_si, t_si, p_si)
            if(g_tp == errorreturn) return
            h = -g_tp / g_tt                               !D2H/DETADP
             
          case default; return
        end select

      case(2)
        select case(drv_p)
          case(0) 
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            h = -1d0 / g_tt                                !D2H/DETA2
             
          case default; return
        end select

      case default; return
    end select
    
  case(1)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0) 
            g_a = air_g_si(1, 0, 0, a_si, t_si, p_si)
            if(g_a == errorreturn) return
            h = g_a                                        !DH/DA
          
          case(1)
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = air_g_si(1, 1, 0, a_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_ap = air_g_si(1, 0, 1, a_si, t_si, p_si)
            if(g_ap == errorreturn) return
            h = (g_tt * g_ap - g_at * g_tp) / g_tt         !D2H/DADP
         
          case default; return
        end select

      case(1)
        select case(drv_p)
        
          case(0) 
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = air_g_si(1, 1, 0, a_si, t_si, p_si)
            if(g_at == errorreturn) return
            h = -g_at / g_tt                               !D2H/DADETA
         
          case default; return
        end select

      case default; return

    end select

  
  case(2)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0) 
            g_tt = air_g_si(0, 2, 0, a_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = air_g_si(1, 1, 0, a_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_aa = air_g_si(2, 0, 0, a_si, t_si, p_si)
            if(g_aa == errorreturn) return
            h = (g_tt * g_aa - g_at**2) / g_tt             !D2H/DA2
          
          case default; return
        end select

      case default; return
    end select

  case default; return
end select

air_a_eta_p_derivatives_si = h

end function

!==========================================================================
subroutine init_it_ctrl_pottemp()
!==========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_pottemp = 0           !default: theta = ideal-gas approximation
ctrl_init_pottemp = 273.15d0
ctrl_eps_exit_pottemp = 0.0001d0  !default = 0.1 mk

end subroutine

!=========================================================================
function pottemp_iteration(a_si, eta_si, p_si, t_si, maxit, eps)
!==========================================================================

!THIS FUNCTION COMPUTES THE ABSOLUTE (POTENTIAL) TEMPERATURE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, SPECIFIC ENTROPY AND (REFERENCE) PRESSURE

!OUTPUT:
!THETA =   (POTENTIAL) TEMPERATURE OF HUMID AIR IN K,
!          I.E. THE TEMPERATURE T_SI THAT SOLVES ETA_SI = AIR_G_ENTROPY_SI(A_SI, T_SI, P_SI)

!INPUT:
!A_SI      MASS FRACTION OF DRY AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA
!T_SI      ABSOLUTE (POTENTIAL) TEMPERATURE IN K, INITIAL VALUE
!MAXIT     MAX. NUMBER OF ITERATIONS
!EPS       REQUIRED TOLERANCE IN K

real*8 pottemp_iteration, a_si, eta_si, p_si, t_si, eps
real*8  s, theta, cp, dt
integer maxit, i

pottemp_iteration = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return
if(eps <= 0d0) return
if(maxit <= 0) return

check_limits = check_limits - 1

theta = t_si
do i = 1, maxit

  !GET ENTROPY AND ITS FIRST DERIVATIVE FOR NEWTON ITERATION
  s = air_g_entropy_si(a_si, theta, p_si)
  cp = air_g_cp_si(a_si, theta, p_si)
  
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
  !AIR_LIMITS
  if(theta < dry_air_tmin .or. theta > dry_air_tmax) then
    pottemp_iteration = errorreturn
    return
  endif
  if(a_si.ne.1d0 .and. &
  (theta < mix_air_tmin .or. theta > mix_air_tmax)) then
    pottemp_iteration = errorreturn
    return
  endif
endif

end function

!==========================================================================
subroutine set_it_ctrl_air_pottemp(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!POTENTIAL TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0       SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_THETA      0           USE DEFAULT POTENTIAL TEMPERATURE THETA = IDEAL-GAS APPROXIMATION
!INIT_THETA      T > 0       USE VALUE T AS POTENTIAL TEMPERATURE TO START

!TOL_THETA       0           USE DEFAULT EXIT ACCURACY FOR POTENTIAL TEMPERATURE (0.1 MK)
!TOL_THETA       EPS > 0     USE EPS AS EXIT ACCURACY FOR POTENTIAL TEMPERATURE

real*8 value
character*13 key

call init_it_ctrl_pottemp

select case (trim(key))

  case ("it_steps");              !ITERATION STEPS
    if(value == 0d0) then
      ctrl_loop_maximum = 100     !DEFAULT = 100
    elseif(value < 0d0) then
      ctrl_loop_maximum = -1
    else   
      ctrl_loop_maximum = nint(value)
    endif

  case ("init_theta");            !START THETA
    if(value == 0d0) then
      ctrl_mode_pottemp = 0       !DEFAULT: IDEAL-GAS APPROXIMATION
    elseif(value > 0d0) then
      ctrl_mode_pottemp = 1
      ctrl_init_pottemp = value
    endif

  case ("tol_theta");             !REQUIRED THETA TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_pottemp = 0.0001d0 !DEFAULT = 0.1 MK
    elseif(value > 0d0) then
      ctrl_eps_exit_pottemp = value
    endif

  case default;write(*,*)key, " not supported"
end select

end subroutine

!==========================================================================
function aux_temperature_si(a_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES AN IDEAL-GAS APPROX. OF ABSOLUTE (POTENTIAL) TEMPERATURE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, SPECIFIC ENTROPY AND (REFERENCE) PRESSURE
                                    
!OUTPUT:
!T(A,ETA,P) IDEAL-GAS APPROX. OF TEMPERATURE OF HUMID AIR IN K

!INPUT:
!A_SI       MASS FRACTION OF DRY AIR IN KG/KG
!ETA_SI     SPECIFIC ENTROPY IN J/(KG K)
!P_SI       ABSOLUTE (REFERENCE) PRESSURE IN PA

real*8 aux_temperature_si, a_si, eta_si, p_si
real*8 tt, pt, ra, rw, cpa, cpv

real*8 etat, rav, xv
real*8 numer, denom

tt = tp_temperature_si
pt = tp_pressure_iapws95_si

ra = gas_constant_air_si ! R / MA
rw = gas_constant_h2o_si ! R / MW

!PROPERTIES AT THE TRIPLE POINT OF WATER
cpa = 1003.69d0  !HEAT CAPACITY OF DRY AIR, IN J/(KG K)
cpv = 1884.352d0 !HEAT CAPACITY OF VAPOUR, IN J/(KG K)

aux_temperature_si = errorreturn
if(a_si < 0d0 .or. a_si > 1d0) return
if(p_si <= 0d0) return

if(a_si == errorreturn) return
if(eta_si == errorreturn) return
if(p_si == errorreturn) return

etat = air_g_entropy_si(a_si, tt, pt)  !ENTROPY OF HUMID AIR AT THE TRIPLE POINT
if(etat == errorreturn) return
rav = a_si * ra + (1d0 - a_si) * rw    !GAS CONSTANT OF HUMID AIR
xv = air_molfraction_vap_si(a_si)
if(xv == errorreturn) return
if(xv < 0d0 .or. xv > 1d0) return

denom = a_si * (cpa + ra) + (1d0 - a_si) * (cpv + rw)
if(denom == 0d0) return

numer = eta_si - etat + rav * log(p_si / pt)
aux_temperature_si = tt * exp(numer / denom)

end function

end module air_3c
