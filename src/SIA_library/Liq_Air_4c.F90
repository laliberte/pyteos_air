module liq_air_4c

!#########################################################################

!THIS MODULE IMPLEMENTS THE ENTHALPY OF WET AIR, AS WELL AS ITS PARTIAL
!DERIVATIVES, DEPENDING ON AIR FRACTION, ENTROPY AND PRESSURE.

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
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     AIR_3A,      FILE AIR_3A.F90
!     ICE_LIQ_4,   FILE ICE_LIQ_4.F90
!     LIQ_AIR_4A,  FILE LIQ_AIR_4A.F90
!     LIQ_AIR_4B,  FILE LIQ_AIR_4B.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     AIR_1,       FILE AIR_1.F90
!     FLU_1,       FILE FLU_1.F90
!     ICE_1,       FILE ICE_1.F90
!     ICE_2,       FILE ICE_2.F90
!     AIR_2,       FILE AIR_2.F90
!     FLU_2,       FILE FLU_2.F90
!     AIR_3B,      FILE AIR_3A.F90
!     FLU_3A,      FILE FLU_3A.F90
!     LIQ_VAP_4,   FILE LIQ_VAP_4.F90

!#########################################################################

use constants_0
use air_3a
use ice_liq_4
use liq_air_4a
use liq_air_4b

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: liq_air_h_si, liq_air_h_cp_si, liq_air_h_density_si, &
          liq_air_h_kappa_s_si, liq_air_h_lapserate_si, &
          liq_air_h_temperature_si, liq_air_potdensity_si,  &
          liq_air_potenthalpy_si, liq_air_pottemp_si, &
          set_it_ctrl_liq_air_pottemp

!CONTROL PARAMETERS OF THE TEMPERATURE ITERATION
integer ctrl_initialized

integer ctrl_mode_pottemp
integer ctrl_loop_maximum

real*8 ctrl_init_pottemp
real*8 ctrl_eps_exit_pottemp

contains

!==========================================================================
function liq_air_h_si(drv_wa, drv_eta, drv_p, wa_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS ENTHALPY OF WET AIR AS A THERMODYNAMIC POTENTIAL,
!DEPENDING ON AIR MASS FRACTION, ENTROPY AND PRESSURE

!OUTPUT:
!LIQ_AIR_H_SI AS THE WA-ETA-P DERIVATIVE

!(D/DWA)^DRV_WA (D/DETA)^DRV_ETA (D/DP)^DRV_P H(WA,ETA,P)

!OF THE SPECIFIC ENTHALPY OF WET AIR, H(WA,ETA,P), IN J/KG

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF WET AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_SI( 0, 0, 0, 0.5, 100, 1E5) = 26898.2871262
!LIQ_AIR_H_SI( 1, 0, 0, 0.5, 100, 1E5) =-1682.2814754|3
!LIQ_AIR_H_SI( 0, 1, 0, 0.5, 100, 1E5) = 280.392597112
!LIQ_AIR_H_SI( 0, 0, 1, 0.5, 100, 1E5) = 0.406967545475
!LIQ_AIR_H_SI( 2, 0, 0, 0.5, 100, 1E5) = 35.6650436194
!LIQ_AIR_H_SI( 1, 1, 0, 0.5, 100, 1E5) = 1.78333840654
!LIQ_AIR_H_SI( 1, 0, 1, 0.5, 100, 1E5) = 0.811934874970
!LIQ_AIR_H_SI( 0, 2, 0, 0.5, 100, 1E5) = 8.91712318143E-02
!LIQ_AIR_H_SI( 0, 1, 1, 0.5, 100, 1E5) = 1.55091234445E-04
!LIQ_AIR_H_SI( 0, 0, 2, 0.5, 100, 1E5) =-3.83855233558E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_SI( 0, 0, 0, 0.5, 100, 1E5) = 26898.5215492
!LIQ_AIR_H_SI( 1, 0, 0, 0.5, 100, 1E5) =-1681.7936611|3
!LIQ_AIR_H_SI( 0, 1, 0, 0.5, 100, 1E5) = 280.393544899
!LIQ_AIR_H_SI( 0, 0, 1, 0.5, 100, 1E5) = 0.40687293002
!LIQ_AIR_H_SI( 2, 0, 0, 0.5, 100, 1E5) = 35.768970891|5
!LIQ_AIR_H_SI( 1, 1, 0, 0.5, 100, 1E5) = 1.7859992519|6
!LIQ_AIR_H_SI( 1, 0, 1, 0.5, 100, 1E5) = 0.811745643965
!LIQ_AIR_H_SI( 0, 2, 0, 0.5, 100, 1E5) = 8.91776656830E-02
!LIQ_AIR_H_SI( 0, 1, 1, 0.5, 100, 1E5) = 1.55067379031E-04
!LIQ_AIR_H_SI( 0, 0, 2, 0.5, 100, 1E5) =-3.83770118470E-06

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A AND AIR_3A, AND ON THE ITERATION SETTINGS
!      FOR TEMPERATURE BY SET_IT_CTRL_LIQ_AIR_POTTEMP OF THIS MODULE

real*8 liq_air_h_si, wa_si, eta_si, p_si
real*8 t

integer drv_wa, drv_eta, drv_p

liq_air_h_si = errorreturn

if(drv_wa < 0) return
if(drv_eta < 0) return
if(drv_p < 0) return
if(drv_wa + drv_eta + drv_p > 2) return

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

!COMPUTE TEMPERATURE FROM ENTROPY.
!IF P_SI = IN-SITU PRESSURE, THIS IS IN-SITU TEMPERATURE
!IF P_SI = REFERENCE PRESSURE, THIS IS POTENTIAL TEMPERATURE
t = liq_air_h_temperature_si(wa_si, eta_si, p_si)
if(t == errorreturn) return

liq_air_h_si = liq_air_a_eta_p_derivatives_si(drv_wa, drv_eta, drv_p, wa_si, t, p_si)

end function

!=========================================================================
function liq_air_h_cp_si(wa_si, eta_si, p_si)
!=========================================================================
                                 
!OUTPUT:
!CP(A,ETA,P) = T/H_ETA_ETA  HEAT CAPACITY OF WET AIR IN J/(KG K)

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_CP_SI(0.5, 100, 1E5) = 3144.42888595

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_CP_SI(0.5, 100, 1E5) = 3144.21265404

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A, AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 liq_air_h_cp_si, wa_si, eta_si, p_si
real*8 t

liq_air_h_cp_si = errorreturn

t = liq_air_h_temperature_si(wa_si, eta_si, p_si)
if(t == errorreturn) return
if(t <= 0d0) return

liq_air_h_cp_si = liq_air_g_cp_si(wa_si, t, p_si)

end function

!=========================================================================
function liq_air_h_density_si(wa_si, eta_si, p_si)
!=========================================================================
                                 
!OUTPUT:
!D(A,ETA,P) = 1/H_P DENSITY OF WET AIR IN KG/M3

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_DENSITY_SI(0.5, 100, 1E5) = 2.45719839609

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_DENSITY_SI(0.5, 100, 1E5) = 2.45776980040

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A, AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 liq_air_h_density_si, wa_si, eta_si, p_si
real*8 t

liq_air_h_density_si = errorreturn

t = liq_air_h_temperature_si(wa_si, eta_si, p_si)
if(t == errorreturn) return
if(t <= 0d0) return

liq_air_h_density_si = liq_air_g_density_si(wa_si, t, p_si)

end function

!=========================================================================
function liq_air_h_kappa_s_si(wa_si, eta_si, p_si)
!=========================================================================
                                 
!OUTPUT:
!KAPPA(A,ETA,P) = - H_PP/H_P ADIABATIC COMPRESSIBILITY OF WET AIR IN 1/PA

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_KAPPA_S_SI(0.5, 100, 1E5) = 9.43208464228E-06

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_KAPPA_S_SI(0.5, 100, 1E5) = 9.43218607469E-06

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A, AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 liq_air_h_kappa_s_si, wa_si, eta_si, p_si
real*8 d, c

liq_air_h_kappa_s_si = errorreturn

d = liq_air_h_density_si(wa_si, eta_si, p_si)
if(d == errorreturn) return
if(d <= 0d0) return

c = liq_air_h_si(0, 0, 2, wa_si, eta_si, p_si)
if(c == errorreturn) return
if(c >= 0d0) return

liq_air_h_kappa_s_si = -d * c

end function

!=========================================================================
function liq_air_h_lapserate_si(wa_si, eta_si, p_si)
!=========================================================================
                                 
!OUTPUT:
!GAMMA(A,ETA,P) = H_ETA_P  MOIST-ADIABATIC LAPSE RATE IN K/PA

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF WET AIR IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_LAPSERATE_SI(0.5, 100, 1E5) = 1.55091234445E-04

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_LAPSERATE_SI(0.5, 100, 1E5) = 1.55067379031E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A, AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8  liq_air_h_lapserate_si, wa_si, eta_si, p_si

liq_air_h_lapserate_si = liq_air_h_si(0, 1, 1, wa_si, eta_si, p_si)

end function

!=========================================================================
function liq_air_h_temperature_si(wa_si, eta_si, p_si)
!=========================================================================

!THIS IS THE INVERSE FUNCTION FOR TEMPERATURE OF
!    ETA_SI = LIQ_AIR_G_ENTROPY_SI(WA_SI, T_SI, P_SI) IN LIQ_AIR_4B
                                 
!OUTPUT:
!T(A,ETA,P) (POTENTIAL) TEMPERATURE OF WET AIR IN K

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_H_TEMPERATURE_SI(0.5, 100, 1E5) = 280.392597112

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_H_TEMPERATURE_SI(0.5, 100, 1E5) = 280.393544899

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A, AIR_3A AND ON THOSE MADE IN THIS MODULE

real*8 liq_air_h_temperature_si, wa_si, eta_si, p_si
real*8 t1, t2, eps
integer maxit

liq_air_h_temperature_si = errorreturn

call init_it_ctrl_pottemp

select case(ctrl_mode_pottemp)
  case(0);   t1 = aux_temperature_min_si(wa_si, p_si)
  case default; t1 = ctrl_init_pottemp
end select
t2 = aux_temperature_max_si(wa_si, p_si)

  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    liq_air_h_temperature_si = t1
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

eps = ctrl_eps_exit_pottemp
if(eps <= 0d0) return

!RUN THE ITERATION
liq_air_h_temperature_si = temperatureiteration_brent(wa_si, eta_si, p_si, t1, t2, maxit, eps)

end function

!==========================================================================
function liq_air_potdensity_si(a_si, t_si, p_si, pr_si)
!=========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL DENSITY OF WET AIR

!OUTPUT:
!RHO_THETA(A,T,P,PR) POTENTIAL DENSITY OF WET AIR IN KG/M3

!INPUTS:
!A_SI      ABSOLUTE AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_POTDENSITY_SI(0.5, 300, 1E4, 1e5) = 1.22518908856

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_POTDENSITY_SI(0.5, 300, 1E4, 1e5) = 1.22550664945

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_LIQ_AIR_POTTEMP OF THIS MODULE

real*8 liq_air_potdensity_si, a_si, t_si, p_si, pr_si
real*8 s, v

liq_air_potdensity_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = liq_air_g_entropy_si(a_si, t_si, p_si)  !specific entropy in-situ
if(s == errorreturn) return

v = liq_air_h_si(0, 0, 1, a_si, s, pr_si) !specific volume at reference pressure
if(v == errorreturn) return
if(v <= 0d0) return

liq_air_potdensity_si = 1d0 / v

end function

!==========================================================================
function liq_air_potenthalpy_si(a_si, t_si, p_si, pr_si)
!=========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL ENTHALPY OF WET AIR

!OUTPUT:
!H_THETA(A,T,P,PR) POTENTIAL ENTHALPY OF WET AIR IN J/KG

!INPUTS:
!A_SI      ABSOLUTE AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_POTENTHALPY_SI(0.5, 300, 1E4, 1E5) = 655299.468706

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_POTENTHALPY_SI(0.5, 300, 1E4, 1E5) = 655155.797982

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_LIQ_AIR_POTTEMP OF THIS MODULE

real*8 liq_air_potenthalpy_si, a_si, t_si, p_si, pr_si
real*8 s

liq_air_potenthalpy_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = liq_air_g_entropy_si(a_si, t_si, p_si)  !specific entropy in-situ
if(s == errorreturn) return

liq_air_potenthalpy_si = liq_air_h_si(0, 0, 0, a_si, s, pr_si)

end function

!==========================================================================
function liq_air_pottemp_si(a_si, t_si, p_si, pr_si)
!=========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL TEMPERATURE OF WET AIR

!OUTPUT:
!THETA(A,T,P,PR) ABSOLUTE POTENTIAL TEMPERATURE OF WET AIR IN K

!INPUTS:
!A_SI      ABSOLUTE AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_POTTEMP_SI(0.5, 300, 1E4, 1E5) = 348.223146077

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_POTTEMP_SI(0.5, 300, 1E4, 1E5) = 348.222379217

real*8 liq_air_pottemp_si, a_si, t_si, p_si, pr_si
real*8 s

liq_air_pottemp_si = errorreturn

if(a_si < 0d0 .or. a_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

if(p_si == pr_si) then
  liq_air_pottemp_si = t_si
  return
end if

s = liq_air_g_entropy_si(a_si, t_si, p_si)  !specific entropy in-situ
if(s == errorreturn) return

liq_air_pottemp_si = liq_air_h_temperature_si(a_si, s, pr_si)

end function

!==========================================================================
function temperatureiteration_brent(wa_si, eta_si, p_si, t1, t2, &
!=========================================================================
                                    maxit, eps)

!THE FUNCTION COMPUTES THE TEMPERATURE AS A FUNCTION OF AIR FRACTION, ENTROPY AND PRESSURE,
!COMPUTED BY BRENT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/BRENT'S_METHOD

!OUTPUT: 
!        TEMPERATUREITERATION_BRENT: TEMPERATURE IN K
!        THE VALUE ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO LIQ_AIR_G_SI HAS RETURNED AN ERROR
!        - TEMPERATURE HAS TAKEN A ZERO OR NEG_ATIVE VALUE DURING THE ITERATION

!INPUT: 
!       WA_SI: MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!      ETA_SI: SPECIFIC ENTROPY IN J/(KG K)
!        P_SI: ABSOLUTE PRESSURE IN PA
!          T1: INITIAL GUESS FOR TEMPERATURE IN K
!          T2: COUNTERPOINT TEMPERATURE
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF TEMPERATURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE

real*8 temperatureiteration_brent, wa_si, eta_si, p_si, t1, t2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
integer it, maxit

logical mflag

temperatureiteration_brent = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0 .or. &
  t1 <= 0d0 .or. &
  t2 <= 0d0 .or. &
  p_si <= 0d0 .or. &
  eps == 0d0 .or. &
  maxit <= 0) return

a = t1
fa = liq_air_g_si(0, 1, 0, wa_si, a, p_si)
if(fa == ErrorReturn) then !permit excursions to the humid-air region
  fa = air_g_si(0, 1, 0, wa_si, a, p_si)
  If(fa == ErrorReturn) return
end if
fa = -fa - eta_si

b = t2
fb = liq_air_g_si(0, 1, 0, wa_si, b, p_si)
If(fb == ErrorReturn)then !permit excursions to the humid-air region
  fb = air_g_si(0, 1, 0, wa_si, b, p_si)
  if(fb == ErrorReturn) return
end if
fb = -fb - eta_si

if(fa * fb > 0d0) return

if(abs(fa) < abs(fb)) then
  call swap(a, b)
  call swap(fa, fb)
end if

c = a
fc = fa
mflag = .true.

check_limits = check_limits - 1

do it = 1, maxit

  if(fb == 0d0) then
    temperatureiteration_brent = b
    go to 999
  end if

  if(eps > 0d0) then                !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      temperatureiteration_brent = b
      go to 999
    end if
  else                                 !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      temperatureiteration_brent = b
      go to 999
    end if
  end if

  if(fa == fb) go to 999

  if(fa /= fc .and. fb /= fc) then
    s = a * fb * fc / ((fa - fb) * (fa - fc)) +  &
        b * fa * fc / ((fb - fa) * (fb - fc)) +  &
        c * fa * fb / ((fc - fa) * (fc - fb))
  else
    s = b - (b - a) * fb / (fb - fa)
  end if

  if(((3 * a + b) / 4 - s) * (b - s) > 0d0 .or. &
     (mflag .and. abs(s - b) >= 0.5 * abs(b - c)) .or. &
     (.not. mflag .and. abs(s - b) >= 0.5 * (c - d))) then
    s = 0.5 * (a + b)
  else
    mflag = .false.
  end if

  fs = liq_air_g_si(0, 1, 0, wa_si, s, p_si)
  if(fs == ErrorReturn)then !permit excursions to the humid-air region
    fs = air_g_si(0, 1, 0, wa_si, s, p_si)
    if(fs == ErrorReturn) return
  end if
  fs = -fs - eta_si

  d = c
  c = b
  fc = fb

  if(fa * fs < 0d0) then
    b = s
    fb = fs
  else
    a = s
    fa = fs
  end if

  if(abs(fa) < abs(fb)) then
    call swap(a, b)
    call swap(fa, fb)
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(b < flu_tmin .or. b > flu_tmax) then
    temperatureiteration_brent = errorreturn
    return
  endif
  !AIR_LIMITS
  if(b < dry_air_tmin .or. b > dry_air_tmax) &
    temperatureiteration_brent = errorreturn
endif

end function

!==========================================================================
function liq_air_a_eta_p_derivatives_si(drv_wa, drv_eta, drv_p, &
                                        wa_si, t_si, p_si)
!=========================================================================
                                              
!THIS FUNCTION COMPUTES THE WET-AIR PARTIAL WA-ETA-P DERIVATIVES OF H FROM WA-T-P DERIVATIVES OF G

real*8 liq_air_a_eta_p_derivatives_si, wa_si, t_si, p_si
real*8 g, g_t, g_p, g_a
real*8 g_aa, g_at, g_ap
real*8 g_tt, g_pp, g_tp
real*8 h

integer drv_wa, drv_eta, drv_p

!IN ONE CASE WE DO NOT AT ALL NEED TO COMPUTE THE GIBBS FUNCTION:
if(drv_wa == 0 .and. drv_eta == 1 .and. drv_p == 0) then
  liq_air_a_eta_p_derivatives_si = t_si
  return
end if

liq_air_a_eta_p_derivatives_si = errorreturn
g_tp = errorreturn

select case(drv_wa)

  case(0)
    select case(drv_eta)
    
      case(0)
        select case(drv_p)
        
          case(0)
            g = liq_air_g_si(0, 0, 0, wa_si, t_si, p_si)
            if(g == errorreturn) return
            g_t = liq_air_g_si(0, 1, 0, wa_si, t_si, p_si)
            if(g_t == errorreturn) return
            h = g - t_si * g_t                          !H
          
          case(1)
            g_p = liq_air_g_si(0, 0, 1, wa_si, t_si, p_si)
            if(g_p == errorreturn) return
            h = g_p                                      !DH/DP
                  
          case(2)
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = liq_air_g_si(0, 1, 1, wa_si, t_si, p_si)
            if(g_tp == errorreturn) return
            g_pp = liq_air_g_si(0, 0, 2, wa_si, t_si, p_si)
            if(g_pp == errorreturn) return
            h = (g_tt * g_pp - g_tp**2) / g_tt           !D2H/DP2
                  
          case default; return
        end select
        
      case(1)
        select case(drv_p)
          case(0)
            h = t_si                                     !DH/DETA, HAS ALREADY BEEN HANDLED INITIALLY
          
          case(1)
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_tp = liq_air_g_si(0, 1, 1, wa_si, t_si, p_si)
            if(g_tp == errorreturn) return
            h = -g_tp / g_tt                             !D2H/DETADP
                  
          case default; return
        end select

      case(2)
        select case(drv_p)
          case(0)
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
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
            g_a = liq_air_g_si(1, 0, 0, wa_si, t_si, p_si)
            if(g_a == errorreturn) return
            h = g_a                                      !DH/DA
          
          case(1)
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = liq_air_g_si(1, 1, 0, wa_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_ap = liq_air_g_si(1, 0, 1, wa_si, t_si, p_si)
            if(g_ap == errorreturn) return
            h = (g_tt * g_ap - g_at * g_tp) / g_tt       !D2H/DADP
          
          case default; return
        end select

      case(1)
        select case(drv_p)
        
          case(0)
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = liq_air_g_si(1, 1, 0, wa_si, t_si, p_si)
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
            g_tt = liq_air_g_si(0, 2, 0, wa_si, t_si, p_si)
            if(g_tt == errorreturn) return
            if(g_tt == 0d0) return
            g_at = liq_air_g_si(1, 1, 0, wa_si, t_si, p_si)
            if(g_at == errorreturn) return
            g_aa = liq_air_g_si(2, 0, 0, wa_si, t_si, p_si)
            if(g_aa == errorreturn) return
            h = (g_tt * g_aa - g_at**2) / g_tt           !D2H/DA2
          
          case default; return
        end select

      case default; return
    end select

  case default; return
end select

liq_air_a_eta_p_derivatives_si = h

end function

!==========================================================================
function aux_temperature_max_si(wa_si, p_si)
!=========================================================================

!THIS FUNCTION ESTIMATES THE TEMPERATURE OF WET AIR AS A FUNCTION OF
!THE DRY-AIR FRACTION WA_SI AND THE PRESSURE P_SI
!IN THE ASYMPTOTIC LIMIT OF A VANISHING LIQUID FRACTION, WW << 1.

real*8 aux_temperature_max_si, wa_si, p_si, error_flag

aux_temperature_max_si = errorreturn
if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si <= 0d0) return

!MAXIMUM ENTROPY AT THE SATURATION POINT:
error_flag = set_liq_air_eq_at_a_p(wa_si, p_si)
aux_temperature_max_si = liq_air_temperature_si()

end function

!==========================================================================
function aux_temperature_min_si(wa_si, p_si)
!=========================================================================

!THIS FUNCTION ESTIMATES THE TEMPERATURE OF WET AIR AS A FUNCTION OF
!THE DRY-AIR FRACTION WA_SI AND THE PRESSURE P_SI
!IN THE ASYMPTOTIC LIMIT OF A FREEZING LIQUID FRACTION.

real*8 aux_temperature_min_si, wa_si, p_si

aux_temperature_min_si = errorreturn
if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si <= 0d0) return

!MINIMUM ENTROPY AT THE FREEZING POINT:
aux_temperature_min_si = ice_liq_meltingtemperature_si(p_si)

end function

!==========================================================================
subroutine init_it_ctrl_pottemp()
!=========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_pottemp = 0           !default: theta = aux_temperature_min_si to start
ctrl_init_pottemp = 273.15d0
ctrl_eps_exit_pottemp = 0.0001d0  !default = 0.1 mk

end subroutine

!==========================================================================
subroutine set_it_ctrl_liq_air_pottemp(key, value)
!=========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!POTENTIAL TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0    SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_THETA      0           USE DEFAULT POTENTIAL TEMPERATURE THETA = AUX_TEMPERATURE_MIN_SI TO START
!INIT_THETA      T > 0    USE VALUE T AS POTENTIAL TEMPERATURE TO START

!TOL_THETA       0           USE DEFAULT EXIT ACCURACY FOR POTENTIAL TEMPERATURE (0.1 MK)
!TOL_THETA       EPS > 0  USE EPS AS EXIT ACCURACY FOR POTENTIAL TEMPERATURE

character*13 key
real*8 value

call init_it_ctrl_pottemp

select case (trim(key))

  case ("it_steps");                 !ITERATION STEPS
      if(value == 0d0) then
         ctrl_loop_maximum = 100     !DEFAULT = 100
      elseif(value < 0d0) then
         ctrl_loop_maximum = -1
      else
         ctrl_loop_maximum = nint(value)
      endif

  case ("init_theta");               !START THETA
      if(value == 0d0) then
         ctrl_mode_pottemp = 0       !DEFAULT: THETA = AUX_TEMPERATURE_MIN_SI
      elseif(value > 0d0) then
         ctrl_mode_pottemp = 1
         ctrl_init_pottemp = nint(value)
      endif

  case ("tol_theta");                !REQUIRED THETA TOLERANCE
      if(value == 0d0) then
         ctrl_eps_exit_pottemp = 0.0001d0 !DEFAULT = 0.1 MK
      elseif(value > 0d0) then
         ctrl_eps_exit_pottemp = value
      endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine swap(a, b)
!==========================================================================

real*8 a, b, c
c = a
a = b
b = c

end subroutine

end module liq_air_4c
