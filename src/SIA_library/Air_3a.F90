module air_3a

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS FUNCTION OF HUMID AIR, DEPENDING ON
!DRY-AIR FRACTION, TEMPERATURE AND PRESSURE, AS WELL AS THEIR PARTIAL 
!DERIVATIVES, COMPUTED NUMERICALLY FROM THE HELMHOLTZ FUNCTION OF HUMID AIR.

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
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     AIR_1,       FILE AIR_1.F90
!     AIR_2,       FILE AIR_2.F90

!NESTED USE ASSOCIATIONS
!     FLU_1,       FILE FLU_1.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use air_1
use air_2

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: air_density_si, air_g_si, get_it_ctrl_airdensity, &
          set_it_ctrl_airdensity

!CONTROL PARAMETERS OF THE DENSITY ITERATION
integer ctrl_initialized
integer ctrl_mode_air
integer ctrl_loop_maximum
integer ctrl_method_air

real*8 ctrl_density_air
real*8 ctrl_eps_exit_air
real*8 ctrl_density2_air

contains 

!==========================================================================
function air_g_si(drv_a, drv_t, drv_p, a_si, t_si, p_si)
!==========================================================================

!THIS ROUTINE IMPLEMENTS THE GIBBS FUNCTION OF HUMID AIR COMPUTED FROM ITS 
!HELMHOLTZ FUNCTION AND ITS FIRST AND SECOND DERIVATIVES WITH RESPECT TO THE 
!AIR FRACTION A_SI, THE ABSOLUTE TEMPERATURE T_SI AND THE ABSOLUTE PRESSURE P_SI

!NOTE: THIS GIBBS FUNCTION IS NOT DESIGNED FOR LIQUID OR CRITICAL AIR STATES
!NOTE: THE ACCURACY OF THIS FUNCTIONS DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!OUTPUT: 
!AIR_G_SI = SPECIFIC GIBBS ENERGY OR ITS DERIVATIVES IN J KG-1 K^-DRV_T PA^-DRV_P

!INPUT:  
!DRV_A = ORDER OF AIR-FRACTION DERIVATIVE, 0 <= DRV_A <= 2
!DRV_T = ORDER OF TEMPERATURE DERIVATIVE, 0 <= DRV_T + DRV_A <= 2
!DRV_P = ORDER OF PRESSURE DERIVATIVE, 0 <= DRV_P + DRV_T + DRV_A <= 2
!A_SI  = MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI  = ABSOLUTE TEMPERATURE IN K
!P_SI  = ABSOLUTE PRESSURE IN PA


!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!air_g_si( 0, 0, 0, 0.9, 300, 1E5) = 4574.43347102
!air_g_si( 1, 0, 0, 0.9, 300, 1E5) =-210115.722120
!air_g_si( 0, 1, 0, 0.9, 300, 1E5) =-911.203092171
!air_g_si( 0, 0, 1, 0.9, 300, 1E5) = 0.911688234199
!air_g_si( 2, 0, 0, 0.9, 300, 1E5) = 1415842.71443
!air_g_si( 1, 1, 0, 0.9, 300, 1E5) = 7566.48354961
!air_g_si( 1, 0, 1, 0.9, 300, 1E5) =-0.483164035114
!air_g_si( 0, 2, 0, 0.9, 300, 1E5) =-4.15491953969
!air_g_si( 0, 1, 1, 0.9, 300, 1E5) = 3.15167187428E-03
!air_g_si( 0, 0, 2, 0.9, 300, 1E5) =-9.14645118671E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!air_g_si( 0, 0, 0, 0.9, 300, 1E5) = 4577.93065689
!air_g_si( 1, 0, 0, 0.9, 300, 1E5) =-210141.953243
!air_g_si( 0, 1, 0, 0.9, 300, 1E5) =-911.170080461
!air_g_si( 0, 0, 1, 0.9, 300, 1E5) = 0.911504137473
!air_g_si( 2, 0, 0, 0.9, 300, 1E5) = 1415779.23410
!air_g_si( 1, 1, 0, 0.9, 300, 1E5) = 7566.34779196
!air_g_si( 1, 0, 1, 0.9, 300, 1E5) =-0.483353002176
!air_g_si( 0, 2, 0, 0.9, 300, 1E5) =-4.15449972148
!air_g_si( 0, 1, 1, 0.9, 300, 1E5) = 3.15111222847E-03
!air_g_si( 0, 0, 2, 0.9, 300, 1E5) =-9.14462130186E-06

real*8 air_g_si, a_si, t_si, p_si
integer drv_a, drv_t, drv_p
real*8 g, d

air_g_si = errorreturn

if (a_si < 0d0 .or. a_si > 1d0) return
if (p_si <= 0d0) return

if(check_limits == 1) then
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) return
  if(a_si /= 1d0 .and. (t_si < mix_air_tmin .or. t_si > mix_air_tmax)) return
else
  if (t_si <= 0d0) return
endif

d = air_density_si(a_si, t_si, p_si)  !NUMERICAL INVERSE FUNCTION TO p = -df/dv

if(d == errorreturn) return
if(d <= 0d0) return

if(check_limits == 1) then
  if(d <= dry_air_dmin .or. d > dry_air_dmax) return
endif

g = air_a_t_p_derivative_si(drv_a, drv_t, drv_p, a_si, t_si, d)
if(g == errorreturn) return

air_g_si = g

end function

!==========================================================================
function air_density_si(a_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS HUMID-AIR DENSITY AS A FUNCTION OF AIR FRACTION, TEMPERATURE AND PRESSURE
!FROM NUMERICAL ITERATION OF THE HELMHOLTZ FUNCTION DERIVATIVE, P = - (DF/DV)_A,T

!OUTPUT:
!AIR_DENSITY_SI: DENSITY OF HUMID AIR IN KG/M3

!INPUT:
!A_SI:  MASS FRACTION OF DRY AIR IN KG/KG
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!air_density_si(0.9, 300, 1E5) = 1.09686619009

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!air_density_si(0.9, 300, 1E5) = 1.09708772444

real*8 air_density_si, a_si, t_si, p_si
real*8 d, d2, eps
integer maxit

air_density_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_it_ctrl_airdensity

!CONSIDER THE CONTROL SETTINGS FOR THE ITERATION
select case(ctrl_mode_air)
  case (1); d = ctrl_density_air
  case (0); d = aux_density_ideal_si(a_si, t_si, p_si)
  case (-1); d = aux_density_ideal_si(a_si, t_si, p_si)
  case default; return
end select

if(d == errorreturn) return
if(d <= 0d0) return

if(ctrl_loop_maximum == 0) then
   maxit = 100
elseif(ctrl_loop_maximum == (-1)) then
   air_density_si = d
   return
elseif(ctrl_loop_maximum > 0) then
   maxit = ctrl_loop_maximum
else
   return
endif

eps = ctrl_eps_exit_air
if(eps == 0d0) return

if(ctrl_method_air > 1) then
  !SPECIFY THE SECOND POINT FOR SECANT OR BRENT METHOD
  if(ctrl_density2_air == 0d0) then
    d2 = densityiteration_newton(a_si, t_si, p_si, d, 1, 1d99)
    if(d2 == errorreturn) return
    if(ctrl_method_air == 2) then !BRENT
      d2 = d + 2d0 * (d2 - d)
    end if
  elseif(ctrl_density2_air == (-1d0)) then
    d2 = densityiteration_newton(a_si, t_si, p_si, d, 1, 1d99)
    if(d2 == errorreturn) return
  elseif(ctrl_density2_air == (-2d0)) then
    d2 = densityiteration_newton(a_si, t_si, p_si, d, 1, 1d99)
    if(d2 == errorreturn) return
    d2 = d + 0.5d0 * (d2 - d)
  elseif(ctrl_density2_air == (-3d0)) then
    d2 = densityiteration_newton(a_si, t_si, p_si, d, 1, 1d99)
    if(d2 == errorreturn) return
    d2 = d + 2d0 * (d2 - d)
  elseif(ctrl_density2_air > 0d0) then
    d2 = ctrl_density2_air
  else 
    return
  endif
end if

!RUN THE ITERATION
select case(ctrl_method_air)
  case(0); d = densityiteration_newton(a_si, t_si, p_si, d, maxit, eps)
  case(1); d = densityiteration_newton(a_si, t_si, p_si, d, maxit, eps)
  case(2); d = densityiteration_brent(a_si, t_si, p_si, d, d2, maxit, eps)
  case(3); d = densityiteration_secant(a_si, t_si, p_si, d, d2, maxit, eps)
  case default; return
end select

if(d == errorreturn) return
if(d <= 0) return

air_density_si = d

end function

!==========================================================================
function densityiteration_newton(a, t, p, d, maxit, eps)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF HUMID AIR AS A FUNCTION OF AIR FRACTION, TEMPERATURE AND PRESSURE
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY NEWTON ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/NEWTON%27S_METHOD

!OUTPUT: 
!DENSITYITERATION_NEWTON: DENSITY OF HUMID AIR IN KG/M3
!THE VALUE ERRORRETURN IS RETURNED IF
!    - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!    - THE FUNCTION CALL TO AIR_F_SI HAS RETURNED AN ERROR
!    - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  
!A: MASS FRACTION OF DRY AIR IN KG/KG
!T: ABSOLUTE TEMPERATURE IN K
!P: ABSOLUTE PRESSURE IN PA
!D: INITIAL GUESS FOR DENSITY IN KG/M3
!MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!EPS: REQUIRED ACCURACY OF DENSITY
!     EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!     EPS < 0: RELATIVE DENSITY TOLERANCE

real*8 densityiteration_newton, a, t, p, d, eps
real*8 dd, f_d, f_dd
integer it, maxit

densityiteration_newton = errorreturn

if(a < 0d0 .or. a > 1d0) return
if (p <= 0d0) return

if(check_limits == 1) then
  if(t < dry_air_tmin .or. t > dry_air_tmax) return
  if(a /= 1d0 .and. (t < mix_air_tmin .or. t > mix_air_tmax)) return
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

do it = 1, maxit

  f_d = air_f_si(0, 0, 1, a, t, d)
  if(f_d == errorreturn) goto 999
  
  f_dd = air_f_si(0, 0, 2, a, t, d)
  if(f_dd == errorreturn) goto 999
  
  dd = d * (2d0 * f_d + d * f_dd)
  if(dd <= 0d0) goto 999

  dd = (p - d**2 * f_d) / dd
  d = d + dd

  if(d <= 0d0) goto 999

  if(eps > 0d0) then            !ABSOLUTE LIMIT
    if(abs(dd) < eps) then
      densityiteration_newton = d
      goto 999
    end if
  else                             !RELATIVE LIMIT
    if(abs(dd) < -eps * d) then
      densityiteration_newton = d
      goto 999
    end if
  end if
  
enddo

999 continue

check_limits = check_limits + 1

end function

!==========================================================================
function densityiteration_brent(a_si, t, p, d1, d2, maxit, eps)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF HUMID AIR AS A FUNCTION OF AIR FRACTION, TEMPERATURE AND PRESSURE
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY BRENT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/BRENT'S_METHOD

!OUTPUT: 
!DENSITYITERATION_BRENT: DENSITY OF HUMID AIR IN KG/M3
!THE VALUE ERRORRETURN IS RETURNED IF
!    - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!    - THE FUNCTION CALL TO AIR_F_SI HAS RETURNED AN ERROR
!    - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  
!A_SI: MASS FRACTION OF DRY AIR IN KG/KG
!T: ABSOLUTE TEMPERATURE IN K
!P: ABSOLUTE PRESSURE IN PA
!D1: INITIAL GUESS FOR DENSITY IN KG/M3
!D2: COUNTERPOINT DENSITY
!MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!EPS: REQUIRED ACCURACY OF DENSITY
!     EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!     EPS < 0: RELATIVE DENSITY TOLERANCE

real*8 densityiteration_brent, a_si, t, p, d1, d2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
integer maxit, it
logical mflag

densityiteration_brent = errorreturn

d = errorreturn

If (d1 < d2) then
  a = d1
  b = d2
else
  b = d1
  a = d2
endif

if(a_si < 0d0 .or. a_si > 1d0) return
if(p <= 0d0) return

if(check_limits == 1) then
  if(t < dry_air_tmin .or. t > dry_air_tmax) return
  if(a /=1d0 .and. (t < mix_air_tmin .or. t > mix_air_tmax)) return
  if(a == 0d0) a= 1d-6
  if(b > dry_air_dmax) b = dry_air_dmax
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

fa = a**2 * air_f_si(0, 0, 1, a_si, t, a) - p
fb = b**2 * air_f_si(0, 0, 1, a_si, t, b) - p
if(fa * fb > 0) go to 999

if(abs(fa) < abs(fb)) then
  call swap(a, b)
  call swap(fa, fb)
end if

c = a
fc = fa
mflag = .true.

do it = 1, maxit

  if(fb == 0d0) then
    densityiteration_brent = b
    go to 999
  end if
  if(eps > 0d0) then          !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      densityiteration_brent = b
      go to 999
    end if
  else                           !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      densityiteration_brent = b
      go to 999
    end if
  end if

  if(fa == fb) go to 999

  if(fa /= fc .and. fb /= fc) then
    s = a * fb * fc / ((fa - fb) * (fa - fc)) + &
        b * fa * fc / ((fb - fa) * (fb - fc)) + &
        c * fa * fb / ((fc - fa) * (fc - fb))
  else
    s = b - (b - a) * fb / (fb - fa)
  end if

  if(((3d0 * a + b) / 4d0 - s) * (b - s) > 0d0 .or. &
    (mflag .and. abs(s - b) >= 0.5d0 * abs(b - c)) .or. &
    (.not. mflag .and. abs(s - b) >= 0.5d0 * (c - d))) then
    s = 0.5d0 * (a + b)
  else
    mflag = .false.
  end if

  fs = s**2 * air_f_si(0, 0, 1, a_si, t, s) - p
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

end function

!==========================================================================
function densityiteration_secant(a, t, p, d, d2, maxit, eps)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF HUMID AIR AS A FUNCTION AIR FRACTION, TEMPERATURE AND PRESSURE
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY SECANT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/SECANT_METHOD

!OUTPUT: 
!DENSITYITERATION_SECANT: DENSITY OF HUMID AIR IN KG/M3
!THE VALUE ERRORRETURN IS RETURNED IF
!    - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!    - THE FUNCTION CALL TO AIR_F_SI HAS RETURNED AN ERROR
!    - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  
!A: MASS FRACTION OF DRY AIR IN KG/KG
!T: ABSOLUTE TEMPERATURE IN K
!P: ABSOLUTE PRESSURE IN PA
!D: INITIAL GUESS FOR DENSITY IN KG/M3
!D2: COUNTERPOINT DENSITY
!MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!EPS: REQUIRED ACCURACY OF DENSITY
!     EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!     EPS < 0: RELATIVE DENSITY TOLERANCE

real*8 densityiteration_secant, a, t, p, d, d2, eps
real*8 dd, d1, p1, p2
integer it, maxit

densityiteration_secant = errorreturn

if(a < 0d0 .or. a > 1d0) return
if(p <= 0d0) return

if(check_limits == 1) then
  if(t < dry_air_tmin .or. t > dry_air_tmax) return
  if(a /=1d0 .and. (t < mix_air_tmin .or. t > mix_air_tmax)) return
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

p2 = d2**2 * air_f_si(0, 0, 1, a, t, d2)
if(p2 == errorreturn) go to 999

do it = 1, maxit
  d1 = d2
  p1 = p2
  d2 = d
  p2 = d**2 * air_f_si(0, 0, 1, a, t, d)

  if(p2 == errorreturn) go to 999
  if(p2 == p1) go to 999
  
  dd = -(d2 - d1) * (p2 - p) / (p2 - p1)
  d = d + dd

  if(d <= 0) go to 999

  if(eps > 0d0) then           !ABSOLUTE LIMIT
    if(abs(dd) < eps) then
      densityiteration_secant = d
      go to 999
    end if
  else                            !RELATIVE LIMIT
    if(abs(dd) < -eps * d) then
      densityiteration_secant = d
      go to 999
    end if
  end if
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  if(d <= dry_air_dmin .or. d > dry_air_dmax) then
    densityiteration_secant = errorreturn
    return
  endif
endif

end function

!==========================================================================
function air_a_t_p_derivative_si(drv_a, drv_t, drv_p, a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION CONPUTES A-T-P DERIVATIVES OF G FROM A-T-D DERIVATIVES OF F
!AT GIVEN AIR FRACTION, TEMPERATURE AND DENSITY

!OUTPUT:
!AIR_A_T_P_DERIVATEVE_SI = SPECIFIC GIBBS ENERGY OR ITS DERIVATIVES IN J KG-1 K^-DRV_T PA^-DRV_P

!INPUT:
!DRV_A = ORDER OF AIR-FRACTION DERIVATIVE, 0 <= DRV_A <= 2
!DRV_T = ORDER OF TEMPERATURE DERIVATIVE, 0 <= DRV_T + DRV_A <= 2
!DRV_P = ORDER OF PRESSURE DERIVATIVE, 0 <= DRV_P + DRV_T + DRV_A <= 2
!A_SI  = MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI  = ABSOLUTE TEMPERATURE IN K
!D_SI  = DENSITY OF HUMID AIR IN KG/M3

real*8 air_a_t_p_derivative_si, a_si, t_si, d_si
real*8 g
real*8 n
real*8 f, f_t, f_d, f_tt, f_td, f_dd, f_at, f_ad, f_aa
integer drv_a, drv_t, drv_p

air_a_t_p_derivative_si = errorreturn

n = errorreturn

if(drv_a < 0) return
if(drv_t < 0) return
if(drv_p < 0) return
if(drv_a + drv_t + drv_p > 2) return

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(d_si <= 0d0) return

if(drv_a + drv_t + drv_p == 2) then !COMPUTE COMMON DENOMINATOR OF EACH 2ND DERIVATIVE
  f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
  if(f_d == errorreturn) return
  f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
  if(f_dd == errorreturn) return
  n = 2d0 * f_d + d_si * f_dd
  if(n == 0d0) return
  n = 1d0 / n
end if

g = errorreturn

select case(drv_a)
  case(0)
    select case(drv_t)
      case(0)
        select case(drv_p)
          case(0)
            f = air_f_si(0, 0, 0, a_si, t_si, d_si)
            if(f == errorreturn) return
            f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
            if(f_d == errorreturn) return
            g = f + d_si * f_d                              !G
          case(1); g = 1d0 / d_si                           !G_P
          case(2); g = -n / d_si**3                        !G_PP
        end select
      case(1)
        select case(drv_p)
          case(0) 
            f_t = air_f_si(0, 1, 0, a_si, t_si, d_si)       !G_T
            if(f_t == errorreturn) return
            g = f_t
          case(1) 
            f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
            if(f_td == errorreturn) return
            g = f_td * n / d_si                             !G_TP
        end select
      case(2)
        select case(drv_p)
          case (0)
            f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
            if(f_tt == errorreturn) return
            f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
            if(f_td == errorreturn) return
            g = f_tt - d_si * n * f_td**2                    !G_TT
        end select
    end select

  case(1)
    select case(drv_t)
      case(0)
        select case(drv_p)
          case(0)
            g = air_f_si(1, 0, 0, a_si, t_si, d_si)          !G_A
          case(1) 
            f_ad = air_f_si(1, 0, 1, a_si, t_si, d_si)
            if(f_ad == errorreturn) return
            g = f_ad * n / d_si                              !G_AP
        end select
      case(1)
        select case(drv_p)
          case(0) 
            f_at = air_f_si(1, 1, 0, a_si, t_si, d_si)
            if(f_at == errorreturn) return
            f_ad = air_f_si(1, 0, 1, a_si, t_si, d_si)
            if(f_ad == errorreturn) return
            f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
            if(f_td == errorreturn) return
            g = f_at - d_si * f_ad * f_td * n                !G_AT
        end select
    end select

  case(2)
    select case(drv_t)
      case(0)
        select case(drv_p)
          case(0)
            f_aa = air_f_si(2, 0, 0, a_si, t_si, d_si)
            if(f_aa == errorreturn) return
            f_ad = air_f_si(1, 0, 1, a_si, t_si, d_si)
            if(f_ad == errorreturn) return
            g = f_aa - d_si * f_ad**2 * n                    !G_AA
        end select
    end select
end select

air_a_t_p_derivative_si = g

end function

!==========================================================================
subroutine set_it_ctrl_airdensity(key, value)
!==========================================================================

!THIS SUBROUTINE SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!HUMID-AIR DENSITY FROM AIR-FRACTION, TEMPERATURE AND PRESSURE

!KEY              VALUE
!IT_STEPS         0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS         N > 0       SET ITERATION NUMBER TO N
!IT_STEPS        -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_AIR_DENS    0           USE DEFAULT DENSITY TO START (IDEAL GAS)
!INIT_AIR_DENS   -1           USE IDEAL-GAS DENSITY TO START
!INIT_AIR_DENS    D > 0       USE VALUE D AS DENSITY TO START

!TOL_AIR_DENS     0           USE DEFAULT EXIT ACCURACY FOR DENSITY (0.1 PPM)
!TOL_AIR_DENS     EPS         USE EPS AS EXIT ACCURACY FOR DENSITY (EPS < 0 MEANS RELATIVE ERROR)

!METHOD_AIR       0           USE DEFAULT ITERATION METHOD (NOW: NEWTON METHOD) FOR DENSITY
!METHOD_AIR       1           USE NEWTON METHOD FOR DENSITY
!METHOD_AIR       2           USE BRENT METHOD FOR DENSITY
!METHOD_AIR       3           USE SECANT METHOD FOR DENSITY

!DENS2_AIR        0           USE DEFAULT COUNTERPOINT FOR BRENT/SECANT METHOD FOR DENSITY
!                             BRENT: 2 * NEWTON STEP, SECANT: 1 * NEWTON STEP
!DENS2_AIR       -1           USE NEWTON STEP AS THE FIRST COUNTERPOINT FOR DENSITY
!DENS2_AIR       -2           USE 0.5 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR DENSITY
!DENS2_AIR       -3           USE 2 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR DENSITY
!DENS2_AIR        D > 0       USE D AS THE FIRST COUNTERPOINT FOR DENSITY

character*13 key
real*8 value

call init_it_ctrl_airdensity

select case (trim(key))

  case ("it_steps");            !ITERATION STEPS
    if(value == 0d0) then
       ctrl_loop_maximum = 100  !DEFAULT = 100
    elseif(value < 0d0) then
       ctrl_loop_maximum = -1
    else   
       ctrl_loop_maximum = nint(value)
    endif

  case ("init_air_dens");       !START DENSITY
    if(value == 0d0) then
       ctrl_mode_air = 0        !DEFAULT = IDEAL GAS
    elseif(value < -1d0) then
       !return
    elseif(value < 0d0) then
       ctrl_mode_air = nint(value)
    else   
       ctrl_mode_air = 1
       ctrl_density_air = value
    endif

  case ("tol_air_dens");        !REQUIRED DENSITY TOLERANCE
    if(value == 0d0)then
       ctrl_eps_exit_air = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
       ctrl_eps_exit_air = value
    endif

  case ("method_air");
    if(value == 0d0 .or. value == 1d0 .or. value == 2d0 .or. value == 3d0) &
       ctrl_method_air = nint(value)

  case ("dens2_air");
    if(value >= (-3d0)) ctrl_density2_air = value

  case default; write(*,*)key, " not supported"
end select

end subroutine

!==========================================================================
function get_it_ctrl_airdensity(key)
!==========================================================================

!THIS FUNCTION RETURNS CONTROL PARAMETERS AS SET FOR THE NEWTON OR OTHER ITERATION
!METHOD USED TO COMPUTE HUMID-AIR DENSITY FROM  AIR FRACTION, TEMPERATURE AND PRESSURE 

real*8 get_it_ctrl_airdensity
character*13 key

get_it_ctrl_airdensity = errorreturn

call init_it_ctrl_airdensity

select case (trim(key))

  case ("it_steps");            !MAX. ITERATION STEPS
    get_it_ctrl_airdensity = ctrl_loop_maximum

  case ("init_air_dens");       !INITIAL HUMID AIR DENSITY
    if(ctrl_mode_air == 1) then
      get_it_ctrl_airdensity = ctrl_density_air
    else
      get_it_ctrl_airdensity = ctrl_mode_air
    end if

  case ("tol_air_dens");        !REQUIRED HUMID AIR DENSITY TOLERANCE
    get_it_ctrl_airdensity = ctrl_eps_exit_air

  case ("method_air");          !SELECTED ITERATION METHOD FOR HUMID AIR
    get_it_ctrl_airdensity = ctrl_method_air

  case ("dens2_air");           !COUNTERPOINT VALUE FOR HUMID AIR (IRRELEVANT FOR NEWTON)
    get_it_ctrl_airdensity = ctrl_density2_air

  case default;
    get_it_ctrl_airdensity = errorreturn
    
end select

end function

!==========================================================================
function aux_density_ideal_si(a_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF IDEAL-GAS HUMID AIR AS A FUNCTION OF
!AIR FRACTION, TEMPERATURE AND PRESSURE

!OUTPUT:  
!AUX_DENSITY_IDEAL_SI: DENSITY IN KG/M3

!INPUT:   
!A_SI: MASS FRACTION OF DRY AIR IN KG/KG
!T_SI: ABSOLUTE TEMPERATURE IN K
!p_si: ABSOLUTE PRESSURE IN PA

real*8 aux_density_ideal_si, a_si, t_si, p_si
real*8 ma, r

r = gas_constant_molar_si        !MOLAR GAS CONSTANT R IN J/(MOL K)

aux_density_ideal_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

ma = air_molar_mass_si(a_si)

aux_density_ideal_si = ma * p_si / (r * t_si)

end function

!==========================================================================
subroutine init_it_ctrl_airdensity()
!==========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_mode_air = 0
ctrl_loop_maximum = 100
ctrl_density_air = 1d0
ctrl_eps_exit_air = -0.0000001d0 !RELATIVE, 0.1 PPM

!SET DEFAULT VALUES FOR ALTERNATIVE ITERATION METHODS
ctrl_method_air = 0 !DEFAULT = NEWTON

!SET DEFAULT COUNTERPOINT VALUES FOR ALTERNATIVE ITERATION METHODS
ctrl_density2_air = 0d0  !DEFAULT = .5 * NEWTON STEP

end subroutine

!==========================================================================
subroutine swap(a, b)
!==========================================================================

real*8 a, b, c
c = a
a = b
b = c
end subroutine

end module air_3a
