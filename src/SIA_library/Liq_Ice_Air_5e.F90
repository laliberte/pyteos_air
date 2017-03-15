module liq_ice_air_5e

!#########################################################################

!THIS MODULE IMPLEMENTS THE EQUILIBRIUM PROPERTIES OF WET ICE AIR, 
!I.E., BETWEEN LIQUID WATER, ICE AND HUMID AIR
!
!THIS MODULE ALSO PROVIDES SIMPLE IMPLEMENTATIONS OF DERIVED METEOROLOGICAL
!QUANTITIES LIKE THE EQUIVALENT POTENTIAL TEMPERATURE, ETC.

!#########################################################################

!IMPLEMENTATION IN FORTRAN BY F.B. LALIBERTE AT THE UNIVERSITY OF TORONTO, 2012
!FOR PROBABLE FUTURE PUBLICATION WITH THE DEVELOPERS OF TEOS-10

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
!     AIR_3C,       FILE AIR_3C.F90
!     ICE_LIQ_4,    FILE ICE_LIQ_4.F90
!     ICE_VAP_4,    FILE ICE_VAP_4.F90
!     LIQ_VAP_4,    FILE LIQ_VAP_4.F90

!     LIQ_AIR_4A,   FILE LIQ_AIR_4A.F90
!     LIQ_AIR_4B,   FILE LIQ_AIR_4B.F90
!     LIQ_AIR_4C,   FILE LIQ_AIR_4C.F90

!     ICE_AIR_4A,   FILE ICE_AIR_4A.F90
!     ICE_AIR_4B,   FILE ICE_AIR_4B.F90
!     ICE_AIR_4C,   FILE ICE_AIR_4C.F90

!     LIQ_ICE_AIR_4,FILE LIQ_ICE_AIR_4.F90 

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
use air_3c
use ice_liq_4

use ice_vap_4
use liq_vap_4

use liq_air_4a
use liq_air_4b
use liq_air_4c
use liq_air_4d
use ice_air_4a
use ice_air_4b
use ice_air_4c
use ice_air_4d
use liq_ice_air_4

use liq_ice_air_5a
use liq_ice_air_5b
use liq_ice_air_5c

implicit none
private

character*16, private :: version = 'October 2012'

public :: liq_ice_air_g_temperaturewetbulb_si

!CONTROL PARAMETERS OF THE EQUI TEMPERATURE ITERATION
integer ctrl_initialized

integer ctrl_mode_tempwb
integer ctrl_loop_maximum

real*8 ctrl_init_tempwb
real*8 ctrl_eps_exit_tempwb

contains

function liq_ice_air_g_temperaturewetbulb_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE EQUIVALENT TEMPERATURE OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.
!THIS FUNCTION GIVES AN EXPRESSION THAT CORRESPONDS TO EQUATION (4.5.10) IN EMMANUEL (1994) WHEN
!T_SI IS ABOVE FREEZING AND POTTEMP(A_SI,T_SI,P_SI,PR_SI) IS ABOVE FREEZING.

!OUTPUT:
!TEQ(A,T,P,PR) ABSOLUTE EQUIVALENT TEMPERATURE OF ICE AIR OR LIQUID AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_temperaturewetbulb_si, a_si, t_si, p_si
real*8 t1, t2, eps
integer maxit

liq_ice_air_g_temperaturewetbulb_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

if(a_si==1d0) then
    liq_ice_air_g_temperaturewetbulb_si = t_si
    return
endif

call init_it_ctrl_tempwb

t1=min_liq_ice_air_temperaturewetbulb_si(a_si,t_si,p_si)
t2=max_liq_ice_air_temperaturewetbulb_si(a_si,t_si,p_si)

  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    liq_ice_air_g_temperaturewetbulb_si = t1
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

eps = ctrl_eps_exit_tempwb
if(eps <= 0d0) return

!RUN THE ITERATION
liq_ice_air_g_temperaturewetbulb_si = temperaturewetbulbiteration_brent(a_si, t_si, p_si, t1, t2, maxit, eps)

end function

!=========================================================================
function min_liq_ice_air_temperaturewetbulb_si(wa_si, t_si, p_si)
!==========================================================================
!THIS FUNCTION GIVES A LOWER BOUND ESTIMATE FOR THE WET BULB TEMPERATURE OF WET ICE AIR 
!FROM ITS DRY-AIR FRACTION WA_SI IN KG/KG, ITS ABSOLUTE TEMPERATURE, T_SI IN K AND ITS
!PRESSURE P_SI IN PA

real*8 min_liq_ice_air_temperaturewetbulb_si, wa_si, p_si, t_si

real*8 lv, t_temp, a

min_liq_ice_air_temperaturewetbulb_si = errorreturn

if(t_si <= 0d0) return
if(wa_si < 0d0 .or. wa_si > 1d0) return

if(set_liq_air_eq_at_t_p(tp_temperature_si, sealevel_pressure_si)==errorreturn) return
lv = liq_air_enthalpy_evap_si()

a = liq_ice_air_massfraction_air_si(t_si,p_si)
if(a/=errorreturn.and.a>=wa_si) then
    t_temp = t_si + (1d0-a)*lv/liq_ice_air_g_wb_cp_si(wa_si,t_si,p_si)
else
    t_temp = t_si + (1d0-wa_si)*lv/liq_ice_air_g_wb_cp_si(wa_si,t_si,p_si)
endif

if(t_temp>=mix_air_tmax) then
    min_liq_ice_air_temperaturewetbulb_si = mix_air_tmax
else
    min_liq_ice_air_temperaturewetbulb_si = t_temp
endif

end function

!=========================================================================
function max_liq_ice_air_temperaturewetbulb_si(wa_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION GIVES AN UPPER BOUND ESTIMATE FOR THE WET BULB TEMPERATURE OF WET ICE AIR 
!FROM ITS DRY-AIR FRACTION WA_SI IN KG/KG, ITS ABSOLUTE TEMPERATURE, T_SI IN K AND ITS
!PRESSURE P_SI IN PA

real*8 max_liq_ice_air_temperaturewetbulb_si, wa_si, p_si, t_si

max_liq_ice_air_temperaturewetbulb_si = errorreturn

if(t_si <= 0d0) return
if(wa_si < 0d0 .or. wa_si > 1d0) return

max_liq_ice_air_temperaturewetbulb_si = t_si
end function

!==========================================================================
function temperaturewetbulbiteration_brent(wa_si, t_si, p_si, t1, t2, &
                                    maxit, eps)
!=========================================================================

!THE FUNCTION COMPUTES THE WET BULB TEMPERATURE AS A FUNCTION OF AIR FRACTION, ENTROPY AND PRESSURE,
!COMPUTED BY BRENT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/BRENT'S_METHOD

!OUTPUT: 
!        TEMPERATUREITERATION_BRENT: EQUIVALENT TEMPERATURE IN K
!        THE VALUE ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO LIQ_AIR_G_SI HAS RETURNED AN ERROR
!        - TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT: 
!       WA_SI: MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!        T_SI: TEMPERATURE IN K
!        P_SI: ABSOLUTE PRESSURE IN PA
!          T1: INITIAL GUESS FOR EQUIVALENT TEMPERATURE IN K
!          T2: COUNTERPOINT EQUIVALENT TEMPERATURE
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF TEMPERATURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE

real*8 temperaturewetbulbiteration_brent, wa_si, t_si, p_si, t1, t2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
real*8 eta
integer it, maxit

logical mflag

temperaturewetbulbiteration_brent = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0 .or. &
  t1 <= 0d0 .or. &
  t2 <= 0d0 .or. &
  p_si <= 0d0 .or. &
  eps == 0d0 .or. &
  maxit <= 0) return

a = t1

eta = liq_ice_air_g_wb_entropy_si(wa_si, a, p_si)
if(eta>liq_ice_air_g_entropy_si(wa_si, mix_air_tmax, p_si)) then
    fa = -mix_air_tmax
else
    fa = -liq_ice_air_h_temperature_si(wa_si,eta,p_si)
endif
fa = -fa - t_si

b = t2
eta = liq_ice_air_g_wb_entropy_si(wa_si, b, p_si)
if(eta>liq_ice_air_g_entropy_si(wa_si, mix_air_tmax, p_si)) then
    fb = -mix_air_tmax
else
    fb = -liq_ice_air_h_temperature_si(wa_si,eta,p_si)
endif
fb = -fb - t_si

!write(*,*) a,fa,b,fb
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
    temperaturewetbulbiteration_brent = b
    go to 999
  end if

  if(eps > 0d0) then                !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      temperaturewetbulbiteration_brent = b
      go to 999
    end if
  else                                 !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      temperaturewetbulbiteration_brent = b
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

  eta = liq_ice_air_g_wb_entropy_si(wa_si, s, p_si)
  fs = -liq_ice_air_h_temperature_si(wa_si,eta,p_si)
  fs = -fs - t_si

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
    temperaturewetbulbiteration_brent = errorreturn
    return
  endif
  !AIR_LIMITS
  if(b < dry_air_tmin .or. b > dry_air_tmax) &
    temperaturewetbulbiteration_brent = errorreturn
endif

end function


!==========================================================================
subroutine init_it_ctrl_tempwb()
!=========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_tempwb = 0           !default: minimum temperature is equal to absolute temperature
ctrl_init_tempwb = 273.15d0
ctrl_eps_exit_tempwb = 0.0001d0  !default = 0.1 mk

end subroutine

!==========================================================================
subroutine set_it_ctrl_liq_ice_air_tempwb(key, value)
!=========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!EQUIVALENT TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0    SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_THETA      0           USE DEFAULT EQUIVALENT TEMPERATURE TE=T_SI TO START
!INIT_THETA      T > 0    USE VALUE T AS EQUIVALENT TEMPERATURE TO START

!TOL_THETA       0           USE DEFAULT EXIT ACCURACY FOR EQUIVALENT TEMPERATURE (0.1 MK)
!TOL_THETA       EPS > 0  USE EPS AS EXIT ACCURACY FOR EQUIVALENT TEMPERATURE

character*13 key
real*8 value

call init_it_ctrl_tempwb

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
         ctrl_mode_tempwb = 0       !DEFAULT: TE=T_SI
      elseif(value > 0d0) then
         ctrl_mode_tempwb = 1
         ctrl_init_tempwb = nint(value)
      endif

  case ("tol_theta");                !REQUIRED THETA TOLERANCE
      if(value == 0d0) then
         ctrl_eps_exit_tempwb = 0.0001d0 !DEFAULT = 0.1 MK
      elseif(value > 0d0) then
         ctrl_eps_exit_tempwb = value
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

end module liq_ice_air_5e
