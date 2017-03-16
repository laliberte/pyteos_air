module flu_3a

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS FUNCTIONS OF LIQUID WATER AND VAPOUR
!DEPENDING ON TEMPERATURE AND PRESSURE, AS WELL AS THEIR PARTIAL DERIVATIVES,
!COMPUTED NUMERICALLY FROM THE HELMHOLTZ FUNCTION OF FLUID WATER, IAPWS-95.

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., !WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     FLU_1,       FILE FLU_1.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use flu_1

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: get_it_ctrl_density, liq_density_si, liq_g_si, &
          set_it_ctrl_density, vap_density_si, vap_g_si

integer, private :: ctrl_mode_liquid, ctrl_mode_vapour, ctrl_loop_maximum, &
         ctrl_method_vapour, ctrl_method_liquid, ctrl_initialized 
real*8, private :: ctrl_density_liquid, ctrl_density_vapour, ctrl_eps_exit_liquid,  &
        ctrl_eps_exit_vapour, ctrl_density2_liquid, ctrl_density2_vapour

!Coefficients of IF97
!Coefficients of region 1
integer, private, save :: i1i(34), j1i(34)
real*8,  private, save :: n1i(34)

!Coefficients of region 2, ideal part
integer, private, save :: j0i(9)
real*8,  private, save :: n0i(9)

!Coefficients of region 2, residual part
integer, private, save :: iri(43), jri(43)
real*8,  private, save :: nri(43)

!Coefficients of the F03 Gibbs function
integer, private, save :: maxt = 7, maxp = 6
real*8, private, save :: gc03(0:7, 0:6)

contains

!==========================================================================
function liq_g_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS IMPLEMENTS THE GIBBS FUNCTION OF LIQUID WATER COMPUTED FROM IAPWS-95
!AND ITS FIRST AND SECOND DERIVATIVES WITH RESPECT TO THE ABSOLUTE TEMPERATURE T_SI
!AND THE ABSOLUTE PRESSURE P_SI

!OUTPUT: LIQ_G_SI = SPECIFIC GIBBS ENERGY IN J/KG OR ITS DERIVATIVE

!INPUT:  DRV_T = ORDER OF TEMPERATURE DERIVATIVE, 0 <= DRV_T <= 2
!        DRV_P = ORDER OF PRESSURE DERIVATIVE, 0 <= DRV_P + DRV_T <= 2
!        T_SI = ABSOLUTE TEMPERATURE, IN K
!        P_SI = ABSOLUTE PRESSURE, IN PA
!
!CHECK VALUES WITH DEFAULT SETTINGS:
!liq_g_si( 0, 0, 300, 1E5) =-5265.05045577
!liq_g_si( 1, 0, 300, 1E5) =-393.062433815
!liq_g_si( 0, 1, 300, 1E5) = 1.00345555938E-03
!liq_g_si( 2, 0, 300, 1E5) =-13.9354650734
!liq_g_si( 1, 1, 300, 1E5) = 2.75753316815E-07
!liq_g_si( 0, 2, 300, 1E5) =-4.52072086722E-13

!NOTE: THE ACCURACY OF THIS FUNCTIONS DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
integer drv_t, drv_p
real*8 liq_g_si, t_si, p_si
real*8 g, d

liq_g_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)  !NUMERICAL INVERSE FUNCTION OF IAPWS-95

if(d == errorreturn) return
if(d <= 0d0) return

g = flu_t_p_derivative_si(drv_t, drv_p, t_si, d)
if(g == errorreturn) return

liq_g_si = g

end function

!==========================================================================
function vap_g_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS IMPLEMENTS THE GIBBS FUNCTION OF WATER VAPOUR COMPUTED FROM IAPWS-95
!AND ITS FIRST AND SECOND DERIVATIVES WITH RESPECT TO THE ABSOLUTE TEMPERATURE T_SI
!AND THE ABSOLUTE PRESSURE P_SI

!OUTPUT: LIQ_G_SI = SPECIFIC GIBBS ENERGY IN J/KG OR ITS DERIVATIVE

!INPUT:  DRV_T = ORDER OF TEMPERATURE DERIVATIVE, 0 <= DRV_T <= 2
!        DRV_P = ORDER OF PRESSURE DERIVATIVE, 0 <= DRV_P + DRV_T <= 2
!        T_SI = ABSOLUTE TEMPERATURE, IN K
!        P_SI = ABSOLUTE PRESSURE, IN PA

!CHECK VALUES WITH DEFAULT SETTINGS:
!vap_g_si( 0, 0, 300, 1E3) =-180090.341338
!vap_g_si( 1, 0, 300, 1E3) =-9103.67940087
!vap_g_si( 0, 1, 300, 1E3) = 138.388478069
!vap_g_si( 2, 0, 300, 1E3) =-6.24707163427
!vap_g_si( 1, 1, 300, 1E3) = 0.462704658818
!vap_g_si( 0, 2, 300, 1E3) =-0.138455798864

!NOTE: THE ACCURACY OF THIS FUNCTIONS DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
integer drv_t, drv_p
real*8 vap_g_si, t_si, p_si
real*8 g, d

vap_g_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)  !NUMERICAL INVERSE FUNCTION OF IAPWS-95

if(d == errorreturn) return
if(d <= 0) return

g = flu_t_p_derivative_si(drv_t, drv_p, t_si, d)
if(g == errorreturn) return

vap_g_si = g

end function

!==========================================================================
function liq_density_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER AS A FUNCTION OF TEMPERATURE AND PRESSURE

!OUTPUT: 
!LIQ_DENSITY_SI = DENSITY IN KG/M3

!INPUT:  T_SI = ABSOLUTE TEMPERATURE, IN K
!        P_SI = ABSOLUTE PRESSURE, IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_density_si(300, 1E5) = 996.556340389

implicit none
integer maxit
real*8 liq_density_si, t_si, p_si
real*8 d, d2, eps

liq_density_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

!AVOID NUMERICAL PROBLEMS AT THE VERY CRITICAL POINT
if (t_si == cp_temperature_si .and. p_si == cp_pressure_si) then
  liq_density_si = cp_density_si
  return
end if

call init_it_ctrl_density

!CONSIDER THE CONTROL SETTINGS FOR THE ITERATION
select case (ctrl_mode_liquid)
  case (1);  d = ctrl_density_liquid
  case (0)
     if(t_si <= 623.15d0 .and. p_si <= 100000000d0) then  !dbg 'or' changed to 'and'
        !IF-97 AT SUBCRITICAL TEMPERATURE, LIQUID
        d = aux_liq_density_if97_si(t_si, p_si)
     elseif (t_si <= cp_temperature_si .and. p_si <= 16529000d0) then
        !IF-97 AT SUBCRITICAL TEMPERATURE, SUPERHEATED LIQUID
        d = aux_liq_density_if97_si(t_si, p_si)
     elseif (t_si <= 1073.15d0 .and. p_si <= 16529000d0) then
        !IF-97 AT SUBCRITICAL PRESSURE, FLUID
        d = aux_vap_density_if97_si(t_si, p_si)
     elseif (t_si <= 650d0 .and. p_si <= 35000000d0) then
        !CUBIC EOF IN THE CRITICAL REGION
        d = aux_liq_density_critical_si(t_si, p_si)
     elseif (t_si <= 650d0) then
        !DENSE FLUID
        d = 1000d0
     else
        !IDEAL GAS ANYWHERE ELSE
        d = aux_density_ideal_si(t_si, p_si)
     end if
  case (-1); d = aux_liq_density_if97_si(t_si, p_si)
  case (-2); d = aux_density_eos80_si(t_si, p_si)
  case (-3); d = aux_liq_density_f03_si(t_si, p_si)
  case default; return
end select

if(d == errorreturn) return
if(d <= 0d0) return

  if(ctrl_loop_maximum == 0)then
     maxit = 100
  elseif(ctrl_loop_maximum == (-1))then
     liq_density_si = d
     return
  elseif(ctrl_loop_maximum > 0)then
     maxit = ctrl_loop_maximum
  else
     return
  endif

eps = ctrl_eps_exit_liquid
if(eps == 0d0) return

if(ctrl_method_liquid > 1) then
  !SPECIFY THE SECOND POINT FOR SECANT OR BRENT METHOD
  if(ctrl_density2_liquid == 0)then
     d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
     if(d2 == errorreturn) return
     if(ctrl_method_liquid == 2) then !BRENT
        d2 = d + 2d0 * (d2 - d)
     end if
   elseif(ctrl_density2_liquid == (-1))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
   elseif(ctrl_density2_liquid == (-2))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
      d2 = d + 0.5d0 * (d2 - d)
   elseif(ctrl_density2_liquid == (-3))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
      d2 = d + 2d0 * (d2 - d)
   elseif(ctrl_density2_liquid > 0)then
      d2 = ctrl_density2_liquid
   else
      return
   endif
end if

select case (ctrl_method_liquid)
  case (0); d = densityiteration_newton(t_si, p_si, d, maxit, eps)
  case (1); d = densityiteration_newton(t_si, p_si, d, maxit, eps)
  case (2); d = densityiteration_brent(t_si, p_si, d, d2, maxit, eps)
  case (3); d = densityiteration_secant(t_si, p_si, d, d2, maxit, eps)
  case default; return
end select

if(d == errorreturn) return
if(d <= 0d0) return

!AVOID ACCIDENTAL VAPOUR DENSITY
if (t_si < cp_temperature_si .and. p_si < cp_pressure_si) then
  if (d < cp_density_si) return
end if

liq_density_si = d

end function

!==========================================================================
function vap_density_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF WATER VAPOUR AS A FUNCTION OF TEMPERATURE AND PRESSURE

!OUTPUT: VAP_DENSITY_SI = DENSITY IN KG/M3

!INPUT:  T_SI = ABSOLUTE TEMPERATURE, IN K
!        P_SI = ABSOLUTE PRESSURE, IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_density_si(300, 1e3) = 7.22603510025E-03

implicit none
integer maxit
real*8 vap_density_si, t_si, p_si
real*8 d, d2, eps

vap_density_si = errorreturn

d = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

!AVOID NUMERICAL PROBLEMS RIGHT AT THE CRITICAL POINT
if (t_si == cp_temperature_si .and. p_si == cp_pressure_si) then
  vap_density_si = cp_density_si
  return
end if

call init_it_ctrl_density

!CONSIDER THE CONTROL SETTINGS FOR THE ITERATION

select case (ctrl_mode_liquid)
  case (1); d = ctrl_density_vapour
  case (0);
     if (t_si <= 623.15d0 .and. p_si <= cp_pressure_si) then
       !IF-97 AT SUBCRITICAL PRESSURE, VAPOUR
       d = aux_vap_density_if97_si(t_si, p_si)
     elseif (t_si <= 1073.15d0 .and. p_si <= 16529000d0) then
       !IF-97 AT SUBCRITICAL PRESSURE, FLUID
       d = aux_vap_density_if97_si(t_si, p_si)
     elseif (t_si <= 623.15d0 .and. p_si <= 100000000d0) then
       !IF-97 AT SUBCRITICAL TEMPERATURE, LIQUID
       d = aux_liq_density_if97_si(t_si, p_si)
     elseif (t_si <= 650d0 .and. p_si <= cp_pressure_si) then
       !CUBIC EOF IN THE CRITICAL REGION, VAPOUR
       d = aux_vap_density_critical_si(t_si, p_si)
     elseif (t_si <= 650d0 .and. p_si <= 35000000d0) then
       !CUBIC EOF IN THE CRITICAL REGION, LIQUID SINCE P > PC
       d = aux_liq_density_critical_si(t_si, p_si)
     elseif (t_si <= 650d0) then
       !DENSE FLUID
       d = 1000d0
     else
       !IDEAL GAS ANYWHERE ELSE
       d = aux_density_ideal_si(t_si, p_si)
     end if
  case (-1); d = aux_vap_density_if97_si(t_si, p_si)
  case (-2); d = aux_density_ideal_si(t_si, p_si)
  case default; return
end select

if(d == errorreturn) return
if(d <= 0d0) return

if(ctrl_loop_maximum == 0)then
   maxit = 100
elseif(ctrl_loop_maximum == (-1))then
   vap_density_si = d
   return
elseif(ctrl_loop_maximum > 0)then
   maxit = ctrl_loop_maximum
else
   return
endif

eps = ctrl_eps_exit_vapour
if(eps == 0d0) return

if(ctrl_method_vapour > 1) then
  !SPECIFY THE SECOND POINT FOR SECANT OR BRENT METHOD
  if(ctrl_density2_vapour == 0)then
     d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
     if(d2 == errorreturn) return
     if(ctrl_method_vapour == 2) then !BRENT
        d2 = d + 2d0 * (d2 - d)
     end if
   elseif(ctrl_density2_vapour == (-1))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
   elseif(ctrl_density2_vapour == (-2))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
      d2 = d + 0.5d0 * (d2 - d)
   elseif(ctrl_density2_vapour == (-3))then
      d2 = densityiteration_newton(t_si, p_si, d, 1, 1d+99)
      if(d2 == errorreturn) return
      d2 = d + 2d0 * (d2 - d)
   elseif(ctrl_density2_vapour > 0)then
      d2 = ctrl_density2_vapour
   else
      return
   endif
end if

!RUN THE ITERATION
select case (ctrl_method_vapour)
  case (0); d = densityiteration_newton(t_si, p_si, d, maxit, eps)
  case (1); d = densityiteration_newton(t_si, p_si, d, maxit, eps)
  case (2); d = densityiteration_brent(t_si, p_si, d, d2, maxit, eps)
  case (3); d = densityiteration_secant(t_si, p_si, d, d2, maxit, eps)
  case default; return
end select


if(d == errorreturn) return
if(d <= 0d0) return

!AVOID ACCIDENTAL LIQUID DENSITY
if (t_si < cp_temperature_si .and. p_si < cp_pressure_si) then
 if (d > cp_density_si) return
end if

vap_density_si = d

end function

!==========================================================================
function densityiteration_newton(t, p, d, maxit, eps)
!==========================================================================

!THE FUNCTION RETURNS THE DENSITY AS A FUNCTION OF TEMPERATURE AND PRESSURE,
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY NEWTON ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/NEWTON%27S_METHOD

!OUTPUT: 
!        DENSITYITERATION_NEWTON: DENSITY IN KG/M3
!        THE VALUE ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  
!        T: ABSOLUTE TEMPERATURE IN K
!        P: ABSOLUTE PRESSURE IN PA
!        D: INITIAL GUESS FOR DENSITY IN KG/M3
!    MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!      EPS: REQUIRED ACCURACY OF DENSITY
!           EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!           EPS < 0: RELATIVE DENSITY TOLERANCE


implicit none
integer it, maxit
real*8 densityiteration_newton, t, p, d, eps
real*8 dd, f_d, f_dd

densityiteration_newton = errorreturn

if(p <= 0d0) return

if(check_limits == 1) then
  if(t < flu_tmin .or. t > flu_tmax) return
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

do it = 1, maxit

  f_d = flu_f_si(0, 1, t, d)
  if(f_d == errorreturn) goto 999
  
  f_dd = flu_f_si(0, 2, t, d)
  if(f_dd == errorreturn) goto 999
  
  dd = d * (2d0 * f_d + d * f_dd)
  if(dd == 0d0) goto 999
  
  dd = (p - d ** 2 * f_d) / dd
  d = d + dd

  if(d <= 0d0) goto 999

  if(eps > 0d0) then           !ABSOLUTE LIMIT
    if(abs(dd) < eps) then
      densityiteration_newton = d
      goto 999
    end if
  else                            !RELATIVE LIMIT
    if(abs(dd) < -eps * d) then
      densityiteration_newton = d
      goto 999
    end if
  end if
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d <= flu_dmin .or. d > flu_dmax) then
    densityiteration_newton = errorreturn
    return
  endif
endif

end function

!==========================================================================
function densityiteration_brent(t, p, d1, d2, maxit, eps)
!==========================================================================

!THE FUNCTION RETURNS THE DENSITY AS A FUNCTION OF TEMPERATURE AND PRESSURE,
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY BRENT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/BRENT!S_METHOD

!OUTPUT: 
!        DENSITYITERATION_BRENT: DENSITY IN KG/M3
!        THE VALUE ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  T: ABSOLUTE TEMPERATURE IN K
!        P: ABSOLUTE PRESSURE IN PA
!       D1: INITIAL GUESS FOR DENSITY IN KG/M3
!       D2: COUNTERPOINT DENSITY
!    MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!      EPS: REQUIRED ACCURACY OF DENSITY
!           EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!           EPS < 0: RELATIVE DENSITY TOLERANCE


implicit none
integer it, maxit
real*8 densityiteration_brent, t, p, d1, d2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
logical mflag

densityiteration_brent = errorreturn

d = errorreturn

if(p <= 0d0) return

if(check_limits == 1) then
  if(t < flu_tmin .or. t > flu_tmax) return
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

a = d1
fa = a ** 2 * flu_f_si(0, 1, t, a) - p
b = d2
fb = b ** 2 * flu_f_si(0, 1, t, b) - p
if(fa * fb > 0d0) go to 999

if(abs(fa) < abs(fb)) then
  call swap (a, b)
  call swap (fa, fb)
end if

c = a
fc = fa
mflag = .true.

do it = 1, maxit

  if(fb == 0d0) then
    densityiteration_brent = b
    goto 999
  end if
  if(eps > 0d0) then           !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      densityiteration_brent = b
      goto 999
    end if
  else                            !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      densityiteration_brent = b
      goto 999
    end if
  end if

  if(fa == fb) goto 999

  if(fa /= fc .and. fb /= fc) then
    s = a * fb * fc / ((fa - fb) * (fa - fc)) +  &
        b * fa * fc / ((fb - fa) * (fb - fc)) +  &
        c * fa * fb / ((fc - fa) * (fc - fb))
  else
    s = b - (b - a) * fb / (fb - fa)
  end if

  if(((3d0 * a + b) / 4d0 - s) * (b - s) > 0d0 .or. &
     (mflag .and. abs(s - b) >= 0.5d0 * abs(b - c)) .or.  &
     (.not.mflag .and. abs(s - b) >= 0.5d0 * (c - d))) then
    s = 0.5d0 * (a + b)
  else
    mflag = .false.
  end if

  fs = s ** 2 * flu_f_si(0, 1, t, s) - p
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
    call swap (a, b)
    call swap (fa, fb)
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d <= flu_dmin .or. d > flu_dmax) then
    densityiteration_brent = errorreturn
    return
  endif
endif

end function

!==========================================================================
function densityiteration_secant(t, p, d, d2, maxit, eps)
!==========================================================================

!THE FUNCTION RETURNS THE DENSITY AS A FUNCTION OF TEMPERATURE AND PRESSURE,
!COMPUTED FROM THE HELMHOLTZ FUNCTION BY SECANT ITERATION
!HTTP://EN.WIKIPEDIA.ORG/WIKI/SECANT_METHOD

!OUTPUT: 
!        DENSITYITERATION_SECANT: DENSITY IN KG/M3
!        THE VALUE ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!INPUT:  
!        T: ABSOLUTE TEMPERATURE IN K
!        P: ABSOLUTE PRESSURE IN PA
!        D: INITIAL GUESS FOR DENSITY IN KG/M3
!       D2: COUNTERPOINT DENSITY
!    MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!      EPS: REQUIRED ACCURACY OF DENSITY
!           EPS > 0: ABSOLUTE DENSITY TOLERANCE IN KG/M3
!           EPS < 0: RELATIVE DENSITY TOLERANCE


implicit none
integer it, maxit
real*8 densityiteration_secant, t, p, d, d2, eps
real*8 dd, d1, p1, p2

densityiteration_secant = errorreturn

if(p <= 0d0) return

if(check_limits == 1) then
  if(t < flu_tmin .or. t > flu_tmax) return
else
  if(t <= 0d0) return
endif

check_limits = check_limits - 1

p2 = d2 ** 2 * flu_f_si(0, 1, t, d2)
if(p2 == errorreturn) go to 999

do it = 1, maxit
  d1 = d2
  p1 = p2
  d2 = d
  p2 = d ** 2 * flu_f_si(0, 1, t, d)

  if(p2 == errorreturn) goto 999
  if(p2 == p1) goto 999
  
  dd = -(d2 - d1) * (p2 - p) / (p2 - p1)
  d = d + dd

  if(d <= 0d0) goto 999

  if(eps > 0d0) then           !ABSOLUTE LIMIT
    if(abs(dd) < eps) then
      densityiteration_secant = d
      goto 999
    end if
  else                            !RELATIVE LIMIT
    if(abs(dd) < -eps * d) then
      densityiteration_secant = d
      goto 999
    end if
  end if
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  if(d <= flu_dmin .or. d > flu_dmax) then
    densityiteration_secant = errorreturn
    return
  endif
endif

end function

!==========================================================================
function flu_t_p_derivative_si(drv_t, drv_p, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES T-P DERIVATIVES OF G FROM T-D DERIVATIVES OF F

implicit none
integer drv_t, drv_p
real*8 flu_t_p_derivative_si, t_si, d_si

real*8 g

real*8 f, f_t, f_d
real*8 f_tt, f_td, f_dd

flu_t_p_derivative_si = errorreturn

if(t_si <= 0d0) return
if(d_si <= 0d0) return

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0) 
        f = flu_f_si(0, 0, t_si, d_si)          !G
        if(f == errorreturn) return
        f_d = flu_f_si(0, 1, t_si, d_si)
        if(f_d == errorreturn) return
        g = f + d_si * f_d

      case (1) 
        g = 1d0 / d_si                          !G_P

      case (2) 
        f_d = flu_f_si(0, 1, t_si, d_si)        !G_PP
        if(f_d == errorreturn) return
        f_dd = flu_f_si(0, 2, t_si, d_si)
        if(f_dd == errorreturn) return
        g = d_si ** 3 * (2d0 * f_d + d_si * f_dd)
        if(g == 0d0) return
        g = -1d0 / g

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0)
        f_t = flu_f_si(1, 0, t_si, d_si)        !G_T
        if(f_t == errorreturn) return
        g = f_t

      case (1)
        f_d = flu_f_si(0, 1, t_si, d_si)        !G_TP
        if(f_d == errorreturn) return
        f_td = flu_f_si(1, 1, t_si, d_si)
        if(f_td == errorreturn) return
        f_dd = flu_f_si(0, 2, t_si, d_si)
        if(f_dd == errorreturn) return
        g = d_si * (2d0 * f_d + d_si * f_dd)
        if(g == 0d0) return
        g = f_td / g

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0)
        f_d = flu_f_si(0, 1, t_si, d_si)        !G_TT
        if(f_d == errorreturn) return
        f_tt = flu_f_si(2, 0, t_si, d_si)
        if(f_tt == errorreturn) return
        f_td = flu_f_si(1, 1, t_si, d_si)
        if(f_td == errorreturn) return
        f_dd = flu_f_si(0, 2, t_si, d_si)
        if(f_dd == errorreturn) return
        g = 2d0 * f_d + d_si * f_dd
        if(g == 0d0) return
        g = f_tt - d_si * f_td ** 2 / g

      case default; return
    end select

  case default; return
end select

flu_t_p_derivative_si = g

end function

!==========================================================================
function polyf03_gyz(drv_y, drv_z, y, z)
!==========================================================================

!RETURNS THE VALUE OF THE POLYNOMIAL DERIVATIVE
!(D/DY)^DRV_Y (D/DZ)^DRV_Z SUM(J,K) GC03(J,K)*Y^J*Z^K

implicit none

integer drv_y, drv_z
integer j, jmax
integer k, kmax
integer l
real*8 polyf03_gyz, g
real*8 y, z
real*8 yj, zk
real*8 c

g = 0d0
if(y == 0d0) then 
   jmax = drv_y 
else 
   jmax = maxt
endif

if(z == 0d0) then 
   kmax = drv_z 
else 
   kmax = maxp
endif

yj = 1d0
do j = drv_y, jmax   !LOOP OVER POWERS OF Y
  zk = 1d0
  do k = drv_z, kmax    !LOOP OVER POWERS OF Z
    if(gc03(j, k) /= 0d0) then
      c = gc03(j, k) * yj * zk
      do l = 1, drv_y            !FACTORS FROM Y-DERIVATIVES
        c = c * dfloat(j - l + 1)
      enddo
      do l = 1, drv_z            !FACTORS FROM Z-DERIVATIVES
        c = c * dfloat(k - l + 1)
      enddo
      g = g + c
    end if
    if(k < kmax) zk = zk * z
  enddo
  if(j < jmax) yj = yj * y
enddo

polyf03_gyz = g

end function
 
!=========================================================================
function gamma_0(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_0 AS GIVEN IN TABLE 13 OF IF-97

implicit none

integer i, k, drv_t, drv_p
real*8 tau, psi
real*8 gamma_0
real*8 g, pwrt

gamma_0 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

g = 0d0
if(drv_t == 0) then
  if(drv_p == 0) then
    g = log(psi)
  else
    g = 1d0 / psi
    do k = 2, drv_p
      g = dfloat(1 - k) * g / psi
    enddo
  end if
end if

if(drv_p == 0) then
  do i = 1, 9
    pwrt = tau ** (j0i(i) - drv_t)
    do k = 0, drv_t - 1
      pwrt = pwrt * dfloat(j0i(i) - k)
    enddo
    g = g + n0i(i) * pwrt
  enddo
end if

gamma_0 = g

end function

!=========================================================================
function gamma_1(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN TABLE 4 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 tau, psi
real*8 gamma_1
real*8 g, pp, tt
real*8 pwrt, pwrp

gamma_1 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

pp = 7.1d0 - psi
tt = tau - 1.222d0

g = 0d0
do i = 1, 34

  if(tt == 0d0) then
    if(j1i(i) == drv_t) then
      pwrt = 1d0
    elseif(j1i(i) > drv_t) then
      pwrt = 0d0
    else
      return
    endif
  else
    pwrt = tt ** (j1i(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(j1i(i) - k)
  enddo

  if(pp == 0d0) then
    if(i1i(i) == drv_p) then
      pwrp = 1d0
    elseif(i1i(i) > drv_p) then
      pwrp = 0d0
    else
      return
    endif
  else
    pwrp = pp ** (i1i(i) - drv_p)
  end if
  do k = 0, drv_p - 1
    pwrp = -pwrp * dfloat(i1i(i) - k)
  enddo

  g = g + n1i(i) * pwrp * pwrt

enddo

gamma_1 = g

end function

!=========================================================================
function gamma_2(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN EQ. 15 OF IF-97

implicit none

integer drv_t, drv_p
real*8 tau, psi
real*8 gamma_2
real*8 g0, gr

call init_if97_2

gamma_2 = errorreturn

g0 = gamma_0(drv_t, drv_p, tau, psi)
if(g0 == errorreturn) return

gr = gamma_r(drv_t, drv_p, tau, psi)
if(gr == errorreturn) return

gamma_2 = g0 + gr

end function

!=========================================================================
function gamma_r(drv_t, drv_p, tau, psi)
!=========================================================================
                         
!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_R AS GIVEN IN TABLE 14 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 gamma_r
real*8 tau, psi
real*8 g, tt
real*8 pwrt, pwrp

gamma_r = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

tt = tau - 0.5d0

g = 0d0
do i = 1, 43

  if(tt == 0d0) then
    if(jri(i) == drv_t) then
      pwrt = 1d0
    elseif(jri(i) > drv_t) then
      pwrt = 0d0
    else
      return
    endif
  else
    pwrt = tt ** (jri(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(jri(i) - k)
  enddo

  pwrp = psi ** (iri(i) - drv_p)
  do k = 0, drv_p - 1
    pwrp = pwrp * dfloat(iri(i) - k)
  enddo

  g = g + nri(i) * pwrp * pwrt

enddo

gamma_r = g

end function

!==========================================================================
function aux_density_eos80_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER COMPUTED FROM THE
!INTERNATIONAL EQUATION OF STATE OF SEAWATER 1980, EOS-80, AS A FUNCTION OF TEMPERATURE
!AND PRESSURE

!OUTPUT:  
!AUX_DENSITY_EOS80_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none
real*8 aux_density_eos80_si, t_si, p_si
real*8 aw, bw, kw, rw
real*8 t, p

aux_density_eos80_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

!R.C.MILLARD JR.: INTERNATIONAL OCEANOGRAPHIC TABLES VOL.4
!UNESCO TECHNICAL PAPERS IN MARINE SCIENCE 40
!UNESCO 1987

t = 1.00024d0*(t_si - 273.15d0) !CONVERT TO DEGC(T90) THEN USE SAUNDER'S FORMULA
p = (p_si - 101325d0)/1d5

bw = 0.0000850935d0 + t * (-0.00000612293d0 + t * 0.000000052787d0)
aw = 3.239908d0 + t * (0.00143713d0 + t * (0.000116092d0 - t * 0.000000577905d0))
kw = 0.01360477d0 - t * 0.00005155288d0
kw = 19652.21d0 + t * (148.4206d0 + t * (-2.327105d0 + t * kw))
kw = kw + p * (aw + p * bw)
if(kw == 0d0) return

kw = 1d0 - p / kw
if(kw == 0d0) return

rw = 0.0001001685d0 + t * (-0.000001120083d0 + t * 0.000000006536332d0)
rw = 999.842594d0 + t * (0.06793952d0 + t * (-0.00909529d0 + t * rw))

aux_density_eos80_si = rw / kw

end function

!==========================================================================
function aux_density_ideal_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF IDEAL-GAS VAPOUR AS A FUNCTION OF TEMPERATURE
!AND PRESSURE

!OUTPUT:  
!AUX_DENSITY_IDEAL_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA


implicit none
real*8 aux_density_ideal_si, t_si, p_si, r

r = gas_constant_h2o_si         !SPECIFIC GAS CONSTANT OF WATER IN J/(KG K)

aux_density_ideal_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

aux_density_ideal_si = p_si / (r * t_si)

end function

!==========================================================================
subroutine aux_density_critical_si(t_si, p_si, d_liq, d_vap)
!==========================================================================
                                      
!RETURNS APPROXIMATE LIQUID AND/OR VAPOUR DENSITY FROM A CUBIC EQUATION OF STATE AT GIVEN T, P
!USABLE IN THE RANGE 620 K < T < 650 K, 10 MPA < P < 35 MPA
!SUPERCRITICAL RANGE:  RETURN VALID D_LIQ = D_VAP
!2-PHASE REGION:       RETURN VALID D_LIQ /= D_VAP          (INCLUDING THE METASTABLE REGION)
!LIQUID REGION:        RETURN VALID D_LIQ,  INVALID D_VAP
!VAPOUR REGION:        RETURN VALID D_VAP,  INVALID D_LIQ

implicit none
real*8 pr, tr
real*8 r, s, t
real*8 d1, d2, d3
real*8 p1, p2
real*8 dc, pc, tc
real*8 a0, a1, a2, a3
real*8 a10, a20, a11, a21, a12, a22, a03, a13, a23
real*8 d_liq, d_vap, t_si, p_si

!THE CRITICAL PROPERTIES OF IAPWS-95
dc = cp_density_si
pc = cp_pressure_si
tc = cp_temperature_si

a10 = -7.60041479494879d0
a20 = 118.661872386874d0
a11 = -17.463827264079d0
a21 = 186.040087842884d0
a12 = 0.69701967809328d0
a22 = 25.5059905941023d0
a03 = -0.602044738250314d0
a13 = 30.8633119943879d0
a23 = 14.4873846518829d0

d_liq = errorreturn
d_vap = errorreturn
if(p_si <= 0d0) return
if(t_si <= 0d0) return

!REDUCED P AND T
pr = pc / p_si - 1d0
tr = t_si / tc - 1d0

!CUBIC EQ. PR = A0 + A1*DR + A2*DR^2 + A3*DR^3
!WITH REDUCED DENSITY DR = D/DC - 1
a0 = (a10 + a20 * tr) * tr
a1 = (a11 + a21 * tr) * tr
a2 = (a12 + a22 * tr) * tr
a3 = a03 + (a13 + a23 * tr) * tr  !A3 < 0 HOLDS FOR T < 659.6 K

if(a3 >= 0d0) return

if(tr < 0d0) then

  !GET THE PRESSURE RANGE OF 2-PHASE SOLUTIONS AT THE GIVEN TEMPERATURE.
  !SOLUTIONS OF DP/DV = 0, I.E. THE SPINODAL CURVE:
  r = a2 * a2 - 3d0 * a1 * a3
  if(r < 0d0) return
  r = sqrt(r)
  d1 = -(a2 + r) / (3d0 * a3)
  d2 = -(a2 - r) / (3d0 * a3)
  !THE PRESSURE RANGE IS
  p1 = a0 + d1 * (a1 + d1 * (a2 + d1 * a3)) !HIGHEST PRESSURE FOR SUBCOOLED VAPOUR (LOWEST REDUCED PRESSURE)
  p2 = a0 + d2 * (a1 + d2 * (a2 + d2 * a3)) !LOWEST PRESSURE FOR SUPERHEATED LIQUID (HIGHEST REDUCED PRESSURE)

else

  !ONE FLUID STATE
  p1 = pr
  p2 = pr
  
end if

!COEFFICIENTS OF THE CUBIC EQ.   DR^3 + R * DR^2 + S * DR + T = 0
r = a2 / a3
s = a1 / a3
t = (a0 - pr) / a3

select case (get_cubicroots(r, s, t, d1, d2, d3))

  case (1)
     if(d1 <= -1d0) return
     if(pr >= p1) then
       !GET THE VAPOUR DENSITY
       d_vap = (d1 + 1d0) * dc
     end if
     if(pr <= p2) then
       !GET THE LIQUID DENSITY
       d_liq = (d1 + 1d0) * dc
     end if

  case (2) 

  case (3) 
     call sort3up (d1, d2, d3)     !3 SOLUTIONS, MIN IS VAPOUR, MAX IS LIQUID
     if(d1 <= -1d0) return
     if(d3 <= -1d0) return
     d_vap = (d1 + 1d0) * dc
     d_liq = (d3 + 1d0) * dc
          
end select

end subroutine

!==========================================================================
function aux_liq_density_critical_si(t_si, p_si)
!==========================================================================

!RETURNS THE APPROXIMATE LIQUID DENSITY FROM THE CUBIC EQUATION

implicit none
real*8 aux_liq_density_critical_si, t_si, p_si
real*8 d_liq, d_vap

call aux_density_critical_si (t_si, p_si, d_liq, d_vap)
aux_liq_density_critical_si = d_liq

end function

!==========================================================================
function aux_vap_density_critical_si(t_si, p_si)
!==========================================================================

!RETURNS THE APPROXIMATE VAPOUR DENSITY FROM THE CUBIC EQUATION

implicit none
real*8 aux_vap_density_critical_si, t_si, p_si
real*8 d_liq, d_vap

call aux_density_critical_si (t_si, p_si, d_liq, d_vap)
aux_vap_density_critical_si = d_vap

end function

!==========================================================================
function aux_liq_density_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!AUX_LIQ_DENSITY_F03_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 aux_liq_density_f03_si
real*8 t_si, p_si
real*8 d

aux_liq_density_f03_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = aux_liq_g_f03_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

aux_liq_density_f03_si = 1d0 / d

end function

!==========================================================================
function aux_liq_g_f03_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE GIBBS FUNCTION OF PURE WATER AS DEFINED IN
!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58 (2003) 43-115

implicit none

integer drv_t, drv_p
real*8 aux_liq_g_f03_si
real*8 t_si, p_si
real*8 t0, tu, p0, pu
real*8 y, z, g

t0 = celsius_temperature_si     !IN K
tu = 40d0                       !IN K
p0 = sealevel_pressure_si       !IN Pa
pu = 100000000d0                !IN Pa

aux_liq_g_f03_si = errorreturn

if(drv_t < 0) return
if(drv_p < 0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_f03

y = (t_si - t0) / tu
z = (p_si - p0) / pu

g = polyf03_gyz(drv_t, drv_p, y, z)
if(g == errorreturn) return

aux_liq_g_f03_si = g / (tu ** drv_t * pu ** drv_p)

end function

!=========================================================================
function aux_liq_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 1 (LIQUID) IN IAPWS-IF97

!OUTPUT:  
!AUX_LIQ_DENSITY_IF97_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 aux_liq_density_if97_si
real*8 t_si, p_si
real*8 d

aux_liq_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = aux_liq_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

aux_liq_density_if97_si = 1d0 / d

end function
 
!=========================================================================
function aux_vap_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF WATER VAPOUR AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 2 (VAPOUR) IN IAPWS-IF97

!OUTPUT:  
!AUX_VAP_DENSITY_IF97_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 aux_vap_density_if97_si
real*8 t_si, p_si
real*8 d

aux_vap_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = aux_vap_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

aux_vap_density_if97_si = 1d0 / d

end function

!=========================================================================
function aux_liq_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 1 (LIQUID) IN IAPWS-IF97
!
!OUTPUT: 
!AUX_LIQ_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!DRV_P: ORDER OF THE PRESSURE DERIVATIVE (0-2)
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 aux_liq_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 1386d0
pu = 16530000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_SI = 461.51805

aux_liq_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_1

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0); g = gamma_1(0, 0, tau, psi)
        if(g == errorreturn) return
        aux_liq_g_if97_si = rt * g

      case (1); g_p = gamma_1(0, 1, tau, psi)
        if(g_p == errorreturn) return
        aux_liq_g_if97_si = rt * g_p / pu

      case (2); g_pp = gamma_1(0, 2, tau, psi)
        if(g_pp == errorreturn) return
        aux_liq_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0)
        g = gamma_1(0, 0, tau, psi)
        if(g == errorreturn) return
        g_t = gamma_1(1, 0, tau, psi)
        if(g_t == errorreturn) return
        aux_liq_g_if97_si = r * (g - tau * g_t)

      case (1) 
        g_p = gamma_1(0, 1, tau, psi)
        if(g_p == errorreturn) return
        g_tp = gamma_1(1, 1, tau, psi)
        if(g_tp == errorreturn) return
        aux_liq_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
        g_tt = gamma_1(2, 0, tau, psi)
        if(g_tt == errorreturn) return
        aux_liq_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function

!=========================================================================
function aux_vap_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 2 (VAPOUR) IN IAPWS-IF97
!
!OUTPUT: 
!AUX_VAP_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!DRV_P: ORDER OF THE PRESSURE DERIVATIVE (0-2)
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 aux_vap_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 540d0
pu = 1000000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_SI = 461.51805

aux_vap_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_2

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0) 
        g = gamma_2(0, 0, tau, psi)
        if(g == errorreturn) return
        aux_vap_g_if97_si = rt * g

      case (1) 
        g_p = gamma_2(0, 1, tau, psi)
        if(g_p == errorreturn) return
        aux_vap_g_if97_si = rt * g_p / pu

      case (2) 
        g_pp = gamma_2(0, 2, tau, psi)
        if(g_pp == errorreturn) return
        aux_vap_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0) 
        g = gamma_2(0, 0, tau, psi)
        if(g == errorreturn) return
        g_t = gamma_2(1, 0, tau, psi)
        if(g_t == errorreturn) return
        aux_vap_g_if97_si = r * (g - tau * g_t)

      case (1) 
        g_p = gamma_2(0, 1, tau, psi)
        if(g_p == errorreturn) return
        g_tp = gamma_2(1, 1, tau, psi)
        if(g_tp == errorreturn) return
        aux_vap_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
        g_tt = gamma_2(2, 0, tau, psi)
        if(g_tt == errorreturn) return
        aux_vap_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function

!==========================================================================
subroutine init_it_ctrl_density
!==========================================================================

implicit none

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_mode_liquid = 0
ctrl_mode_vapour = 0
ctrl_loop_maximum = 100
ctrl_density_liquid = 1000
ctrl_density_vapour = 1
ctrl_eps_exit_liquid = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_vapour = -0.0000001d0 !RELATIVE, 0.1 PPM

!SET DEFAULT VALUES FOR ALTERNATIVE ITERATION METHODS
ctrl_method_liquid = 0 !DEFAULT = NEWTON
ctrl_method_vapour = 0 !DEFAULT = NEWTON

!SET DEFAULT COUNTERPOINT VALUES FOR ALTERNATIVE ITERATION METHODS
ctrl_density2_liquid = 0d0  !DEFAULT = .5 * NEWTON STEP
ctrl_density2_vapour = 0d0  !DEFAULT = .5 * NEWTON STEP

end subroutine

!==========================================================================
subroutine set_it_ctrl_density(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!IAPWS-95 DENSITY FROM PRESSURE

!KEY                    VALUE
!IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY (IF97) TO START IN THE SUBCRITICAL REGION (T <= 623.16 K OR P <=16.529 MPA)
!                                  USE DEFAULT CUBIC EOS TO START IN THE CRITICAL REGION (623.16 K < T < 650 K AND 16.529 MPA < P < 35 MPA)
!                                  USE DEFAULT 1000 KG/M^3 TO START IN THE REMAINING SUPERCRITICAL REGION
!INIT_LIQ_DENS'       -1           USE IF97 LIQUID DENSITY TO START
!INIT_LIQ_DENS'       -2           USE EOS80 LIQUID DENSITY TO START
!INIT_LIQ_DENS'       -3           USE F03 LIQUID DENSITY TO START
!INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!INIT_VAP_DENS'        0           USE DEFAULT VAPOUR DENSITY (IF97) TO START IN THE SUBCRITICAL REGION (T <= 623.16 K OR P <=16.529 MPA)
!                                  USE DEFAULT CUBIC EOS TO START IN THE CRITICAL REGION (623.16 K < T < 650 K AND 16.529 MPA < P < 35 MPA)
!                                  USE DEFAULT 1000 KG/M^3 TO START IN THE REMAINING SUPERCRITICAL REGION
!INIT_VAP_DENS'       -1           USE IF97 VAPOUR DENSITY TO START
!INIT_VAP_DENS'       -2           USE IDEAL-GAS VAPOUR DENSITY TO START
!INIT_VAP_DENS'        D > 0       USE VALUE D AS VAPOUR DENSITY TO START

!TOL_LIQ_DENS '        0           USE DEFAULT EXIT ACCURACY FOR LIQUID DENSITY (0.1 PPM)
!TOL_LIQ_DENS '        EPS         USE EPS AS EXIT ACCURACY FOR LIQUID DENSITY (EPS < 0 MEANS RELATIVE ERROR)

!TOL_VAP_DENS '        0           USE DEFAULT EXIT ACCURACY FOR VAPOUR DENSITY (0.1 PPM)
!TOL_VAP_DENS '        EPS         USE EPS AS EXIT ACCURACY FOR VAPOUR DENSITY (EPS < 0 MEANS RELATIVE ERROR)

!METHOD_LIQ   '        0           USE DEFAULT ITERATION METHOD (NOW: NEWTON METHOD) FOR LIQUID
!METHOD_LIQ   '        1           USE NEWTON METHOD FOR LIQUID
!METHOD_LIQ   '        2           USE BRENT METHOD FOR LIQUID
!METHOD_LIQ   '        3           USE SECANT METHOD FOR LIQUID

!DENS2_LIQ    '        0           USE DEFAULT COUNTERPOINT FOR BRENT/SECANT METHOD FOR LIQUID
!                                  BRENT: 2 * NEWTON STEP, SECANT: 1 * NEWTON STEP
!DENS2_LIQ    '       -1           USE NEWTON STEP AS THE FIRST COUNTERPOINT FOR LIQUID
!DENS2_LIQ    '       -2           USE 0.5 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR LIQUID
!DENS2_LIQ    '       -3           USE 2 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR LIQUID
!DENS2_LIQ    '        D > 0       USE D AS THE FIRST COUNTERPOINT DENSITY FOR LIQUID

!METHOD_VAP   '        0           USE DEFAULT ITERATION METHOD (NOW: NEWTON METHOD) FOR VAPOUR
!METHOD_VAP   '        1           USE NEWTON METHOD FOR VAPOUR
!METHOD_VAP   '        2           USE BRENT METHOD FOR VAPOUR
!METHOD_VAP   '        3           USE SECANT METHOD FOR VAPOUR

!DENS2_VAP    '        0           USE DEFAULT COUNTERPOINT FOR BRENT/SECANT METHOD FOR VAPOUR
!                                  BRENT: 2 * NEWTON STEP, SECANT: 1 * NEWTON STEP
!DENS2_VAP    '       -1           USE NEWTON STEP AS THE FIRST COUNTERPOINT FOR VAPOUR
!DENS2_VAP    '       -2           USE 0.5 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR VAPOUR
!DENS2_VAP    '       -3           USE 2 * NEWTON STEP AS THE FIRST COUNTERPOINT FOR VAPOUR
!DENS2_VAP    '        D > 0       USE D AS THE FIRST COUNTERPOINT DENSITY FOR VAPOUR

implicit none
character*13 key
real*8 value

call init_it_ctrl_density

select case (trim(key))

  case ('it_steps')             !ITERATION STEPS
    if(value == 0d0)then
       ctrl_loop_maximum = 100  !DEFAULT = 100
    elseif(value < 0d0)then
       ctrl_loop_maximum = -1
    else
       ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')        !START LIQUID DENSITY
    if(value == 0d0)then
       ctrl_mode_liquid = 0     !DEFAULT = IF97
    elseif(value < (-3d0))then
       !exit sub
    elseif(value < 0d0)then
       ctrl_mode_liquid = nint(value)
    else
       ctrl_mode_liquid = 1
       ctrl_density_liquid = value
    endif

  case ('init_vap_dens')        !START VAPOUR DENSITY
    if(value == 0d0)then
       ctrl_mode_vapour = 0     !DEFAULT = IF97
    elseif(value < (-2d0))then
       !exit sub
    elseif(value < 0d0)then
       ctrl_mode_vapour = nint(value)
    else
       ctrl_mode_vapour = 1
       ctrl_density_vapour = value
    endif

  case ('tol_liq_dens')         !REQUIRED LIQUID DENSITY TOLERANCE
    if(value == 0d0)then
       ctrl_eps_exit_liquid = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
       ctrl_eps_exit_liquid = value
    endif

  case ('tol_vap_dens')         !REQUIRED VAPOUR DENSITY TOLERANCE
    if(value == 0d0)then
       ctrl_eps_exit_vapour = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
       ctrl_eps_exit_vapour = value
    endif

  case ('method_liq')
      if(value == 0d0 .or. value == 1d0 .or. value == 2d0 .or. &
      value == 3d0) ctrl_method_liquid = nint(value)

  case ('dens2_liq')
    if(value >= -3d0) ctrl_density2_liquid = value

  case ('method_vap')
      if(value == 0d0 .or. value == 1d0 .or. value == 2d0 .or. &
      value == 3d0) ctrl_method_vapour = nint(value)

  case ('dens2_vap')
      if(value >= -3d0) ctrl_density2_vapour = value

  case default; write(*,*)key,' not supported'
  
end select

end subroutine

!==========================================================================
function get_it_ctrl_density(key)
!==========================================================================

!THIS FUNCTION RETURNS CONTROL PARAMETERS AS SET FOR THE NEWTON ITERATION
!USED TO COMPUTE IAPWS-95 DENSITY FROM PRESSURE

implicit none
character*13 key
real*8 get_it_ctrl_density

get_it_ctrl_density = errorreturn

call init_it_ctrl_density

select case (trim(key))

  case ('it_steps')   !ITERATION STEPS
    get_it_ctrl_density = ctrl_loop_maximum

  case ('init_liq_dens')   !START LIQUID DENSITY
    if(ctrl_mode_liquid == 1) then
      get_it_ctrl_density = ctrl_density_liquid
    else
      get_it_ctrl_density = ctrl_mode_liquid
    end if

  case ('init_vap_dens')   !START VAPOUR DENSITY
    if(ctrl_mode_vapour == 1) then
      get_it_ctrl_density = ctrl_density_vapour
    else
      get_it_ctrl_density = ctrl_mode_vapour
    end if

  case ('tol_liq_dens')      !REQUIRED LIQUID DENSITY TOLERANCE
    get_it_ctrl_density = ctrl_eps_exit_liquid

  case ('tol_vap_dens')      !REQUIRED VAPOUR DENSITY TOLERANCE
    get_it_ctrl_density = ctrl_eps_exit_vapour

  case ('method_liq')        !SELECTED ITERATION METHOD FOR LIQUID
    get_it_ctrl_density = ctrl_method_liquid

  case ('dens2_liq')         !COUNTERPOINT VALUE FOR LIQUID (IRRELEVANT FOR NEWTON)
    get_it_ctrl_density = ctrl_density2_liquid

  case ('method_vap')        !SELECTED ITERATION METHOD FOR VAPOUR
    get_it_ctrl_density = ctrl_method_vapour

  case ('dens2_vap')         !COUNTERPOINT VALUE FOR VAPOUR (IRRELEVANT FOR NEWTON)
    get_it_ctrl_density = ctrl_density2_vapour

  case default;
    get_it_ctrl_density = errorreturn
    
end select

end function

!==========================================================================
subroutine init_f03
!==========================================================================

! G(T,P) = SUM GC03(J,K) * (T/40 C)^J * (P/100 MPA)^K
!          J,K

implicit none

real*8 g, g_t, g_p
real*8 tt, t0, tu, pt, p0, pu, yt, zt

if(gc03(2, 0) == -12357.785933039d0) return

tt = tp_temperature_si
t0 = celsius_temperature_si
tu = 40d0

pt = tp_pressure_iapws95_si !IAPWS-95 TRIPLE POINT PRESSURE
p0 = sealevel_pressure_si
pu = 100000000d0

yt = (tt - t0) / tu
zt = (pt - p0) / pu

gc03 = 0d0

gc03(2, 0) = -12357.785933039d0
gc03(3, 0) = 736.741204151612d0
gc03(4, 0) = -148.185936433658d0
gc03(5, 0) = 58.0259125842571d0
gc03(6, 0) = -18.9843846514172d0
gc03(7, 0) = 3.05081646487967d0
gc03(0, 1) = 100015.695367145d0
gc03(1, 1) = -270.983805184062d0
gc03(2, 1) = 1455.0364540468d0
gc03(3, 1) = -672.50778314507d0
gc03(4, 1) = 397.968445406972d0
gc03(5, 1) = -194.618310617595d0
gc03(6, 1) = 63.5113936641785d0
gc03(7, 1) = -9.63108119393062d0
gc03(0, 2) = -2544.5765420363d0
gc03(1, 2) = 776.153611613101d0
gc03(2, 2) = -756.558385769359d0
gc03(3, 2) = 499.360390819152d0
gc03(4, 2) = -301.815380621876d0
gc03(5, 2) = 120.520654902025d0
gc03(6, 2) = -22.2897317140459d0
gc03(0, 3) = 284.517778446287d0
gc03(1, 3) = -196.51255088122d0
gc03(2, 3) = 273.479662323528d0
gc03(3, 3) = -239.545330654412d0
gc03(4, 3) = 152.196371733841d0
gc03(5, 3) = -55.2723052340152d0
gc03(6, 3) = 8.17060541818112d0
gc03(0, 4) = -33.3146754253611d0
gc03(1, 4) = 28.9796526294175d0
gc03(1, 4) = 28.9796526294175d0
gc03(2, 4) = -55.5604063817218d0
gc03(2, 4) = -55.5604063817218d0
gc03(3, 4) = 48.8012518593872d0
gc03(4, 4) = -26.3748377232802d0
gc03(5, 4) = 6.48190668077221d0
gc03(0, 5) = 4.20263108803084d0
gc03(1, 5) = -2.13290083518327d0
gc03(2, 5) = 4.34420671917197d0
gc03(3, 5) = -1.66307106208905d0
gc03(0, 6) = -0.546428511471039d0

!IAPWS-95 REFERENCE STATE CONDITION
!ENERGY = 0 AND ENTROPY = 0 AT THE TRIPLE POINT:
!THESE ARE ONLY USED AS A FIRST GUESS HERE SO HIGH ACCURACY IS NOT REQUIRED.
gc03(0, 0) = 101.342743139674d0
gc03(1, 0) =  5.90578347909402d0

!NEVERTHELESS, THE ACCURACY CAN BE IMPROVED FOLLOWING KIYOSHI'S SUGGESTION.
!ADJUSTMENT OF COEFFICIENTS

!WRITE(*,*)"GC00, GC10: FIRST ESTIMATES = ", GC03(0,0), GC03(1,0)

g = polyf03_gyz(0, 0, yt, zt)
g_t = polyf03_gyz(1, 0, yt, zt) / tu
g_p = polyf03_gyz(0, 1, yt, zt) / pu

!WRITE(*,*)'G,G_T,G_P=',G,G_T,G_P
gc03(0, 0) = gc03(0, 0) - g + pt * g_p
gc03(1, 0) = gc03(1, 0) - g_t * tu

!WRITE(*,*)"GC00, GC10: SECOND ESTIMATES = ", GC03(0,0), GC03(1,0)

end subroutine
 
!=========================================================================
subroutine init_if97_1
!=========================================================================

implicit none
integer i

if(n1i(1) == 0.14632971213167d0) return

!TABLE 2. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE DIMENSIONLESS GIBBS FREE ENERGY
!FOR REGION 1, EQ. (7)
!I     II            JI            NI
i = 1; i1i(i)  = 0;  j1i(i) = -2;  n1i(i) = 0.14632971213167d0
i = 2; i1i(i)  = 0;  j1i(i) = -1;  n1i(i) = -0.84548187169114d0
i = 3; i1i(i)  = 0;  j1i(i) = 0;   n1i(i) = -3.756360367204d0
i = 4; i1i(i)  = 0;  j1i(i) = 1;   n1i(i) = 3.3855169168385d0
i = 5; i1i(i)  = 0;  j1i(i) = 2;   n1i(i) = -0.95791963387872d0
i = 6; i1i(i)  = 0;  j1i(i) = 3;   n1i(i) = 0.15772038513228d0
i = 7; i1i(i)  = 0;  j1i(i) = 4;   n1i(i) = -0.016616417199501d0
i = 8; i1i(i)  = 0;  j1i(i) = 5;   n1i(i) = 8.1214629983568d-04
i = 9; i1i(i)  = 1;  j1i(i) = -9;  n1i(i) = 2.8319080123804d-04
i = 10; i1i(i) = 1;  j1i(i) = -7;  n1i(i) = -6.0706301565874d-04
i = 11; i1i(i) = 1;  j1i(i) = -1;  n1i(i) = -0.018990068218419d0
i = 12; i1i(i) = 1;  j1i(i) = 0;   n1i(i) = -0.032529748770505d0
i = 13; i1i(i) = 1;  j1i(i) = 1;   n1i(i) = -0.021841717175414d0
i = 14; i1i(i) = 1;  j1i(i) = 3;   n1i(i) = -5.283835796993d-05  
i = 15; i1i(i) = 2;  j1i(i) = -3;  n1i(i) = -4.7184321073267d-04
i = 16; i1i(i) = 2;  j1i(i) = 0;   n1i(i) = -3.0001780793026d-04
i = 17; i1i(i) = 2;  j1i(i) = 1;   n1i(i) = 4.7661393906987d-05
i = 18; i1i(i) = 2;  j1i(i) = 3;   n1i(i) = -4.4141845330846d-06
i = 19; i1i(i) = 2;  j1i(i) = 17;  n1i(i) = -7.2694996297594d-16
i = 20; i1i(i) = 3;  j1i(i) = -4;  n1i(i) = -3.1679644845054d-05
i = 21; i1i(i) = 3;  j1i(i) = 0;   n1i(i) = -2.8270797985312d-06
i = 22; i1i(i) = 3;  j1i(i) = 6;   n1i(i) = -8.5205128120103d-10
i = 23; i1i(i) = 4;  j1i(i) = -5;  n1i(i) = -2.2425281908d-06
i = 24; i1i(i) = 4;  j1i(i) = -2;  n1i(i) = -6.5171222895601d-07
i = 25; i1i(i) = 4;  j1i(i) = 10;  n1i(i) = -1.4341729937924d-13
i = 26; i1i(i) = 5;  j1i(i) = -8;  n1i(i) = -4.0516996860117d-07
i = 27; i1i(i) = 8;  j1i(i) = -11; n1i(i) = -1.2734301741641d-09
i = 28; i1i(i) = 8;  j1i(i) = -6;  n1i(i) = -1.7424871230634d-10
i = 29; i1i(i) = 21; j1i(i) = -29; n1i(i) = -6.8762131295531d-19
i = 30; i1i(i) = 23; j1i(i) = -31; n1i(i) = 1.4478307828521d-20
i = 31; i1i(i) = 29; j1i(i) = -38; n1i(i) = 2.6335781662795d-23
i = 32; i1i(i) = 30; j1i(i) = -39; n1i(i) = -1.1947622640071d-23
i = 33; i1i(i) = 31; j1i(i) = -40; n1i(i) = 1.8228094581404d-24
i = 34; i1i(i) = 32; j1i(i) = -41; n1i(i) = -9.3537087292458d-26

end subroutine
 
!=========================================================================
subroutine init_if97_2
!=========================================================================

implicit none
integer i

if(n0i(1) == -9.6927686500217d0) return

!TABLE 10. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE IDEAL-GAS PART GAMMA_0 OF THE
!DIMENSIONLESS GIBBS FREE ENERGY FOR REGION 2, EQ. (16)
!I     JI           NI
i = 1; j0i(i) = 0;  n0i(i) = -9.6927686500217d0
i = 2; j0i(i) = 1;  n0i(i) = 10.086655968018d0
i = 3; j0i(i) = -5; n0i(i) = -0.005608791128302d0
i = 4; j0i(i) = -4; n0i(i) = 0.071452738081455d0
i = 5; j0i(i) = -3; n0i(i) = -0.40710498223928d0
i = 6; j0i(i) = -2; n0i(i) = 1.4240819171444d0
i = 7; j0i(i) = -1; n0i(i) = -4.383951131945d0
i = 8; j0i(i) = 2;  n0i(i) = -0.28408632460772d0
i = 9; j0i(i) = 3;  n0i(i) = 0.021268463753307d0
!IF(EQ. (16) IS INCORPORATED INTO EQ. (18), INSTEAD OF THE-VALUES FOR N0I(1) AND N0I(1)
!GIVEN ABOVE, THE FOLLOWING VALUES
!N0I(1) = -0.96937268393049D1
!N0I(2) = 0.10087275970006D2
!SHOULD BE USED


!TABLE 11. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE
!RESIDUAL PART GAMMA_R OF THE DIMENSIONLESS GIBBS FREE ENERGY FOR
!REGION 2, EQ. (17)
!I      II           JI           NI
i = 1;  iri(i) = 1;  jri(i) = 0;  nri(i) = -1.7731742473213d-03
i = 2;  iri(i) = 1;  jri(i) = 1;  nri(i) = -0.017834862292358d0
i = 3;  iri(i) = 1;  jri(i) = 2;  nri(i) = -0.045996013696365d0
i = 4;  iri(i) = 1;  jri(i) = 3;  nri(i) = -0.057581259083432d0
i = 5;  iri(i) = 1;  jri(i) = 6;  nri(i) = -0.05032527872793d0
i = 6;  iri(i) = 2;  jri(i) = 1;  nri(i) = -3.3032641670203d-05
i = 7;  iri(i) = 2;  jri(i) = 2;  nri(i) = -1.8948987516315d-04
i = 8;  iri(i) = 2;  jri(i) = 4;  nri(i) = -3.9392777243355d-03
i = 9;  iri(i) = 2;  jri(i) = 7;  nri(i) = -0.043797295650573d0
i = 10; iri(i) = 2;  jri(i) = 36; nri(i) = -2.6674547914087d-05
i = 11; iri(i) = 3;  jri(i) = 0;  nri(i) = 2.0481737692309d-08
i = 12; iri(i) = 3;  jri(i) = 1;  nri(i) = 4.3870667284435d-07
i = 13; iri(i) = 3;  jri(i) = 3;  nri(i) = -3.227767723857d-05
i = 14; iri(i) = 3;  jri(i) = 6;  nri(i) = -1.5033924542148d-03
i = 15; iri(i) = 3;  jri(i) = 35; nri(i) = -0.040668253562649d0
i = 16; iri(i) = 4;  jri(i) = 1;  nri(i) = -7.8847309559367d-10
i = 17; iri(i) = 4;  jri(i) = 2;  nri(i) = 1.2790717852285d-08
i = 18; iri(i) = 4;  jri(i) = 3;  nri(i) = 4.8225372718507d-07
i = 19; iri(i) = 5;  jri(i) = 7;  nri(i) = 2.2922076337661d-06
i = 20; iri(i) = 6;  jri(i) = 3;  nri(i) = -1.6714766451061d-11
i = 21; iri(i) = 6;  jri(i) = 16; nri(i) = -2.1171472321355d-03
i = 22; iri(i) = 6;  jri(i) = 35; nri(i) = -23.895741934104d0
i = 23; iri(i) = 7;  jri(i) = 0;  nri(i) = -5.905956432427d-18
i = 24; iri(i) = 7;  jri(i) = 11; nri(i) = -1.2621808899101d-06
i = 25; iri(i) = 7;  jri(i) = 25; nri(i) = -0.038946842435739d0
i = 26; iri(i) = 8;  jri(i) = 8;  nri(i) = 1.1256211360459d-11
i = 27; iri(i) = 8;  jri(i) = 36; nri(i) = -8.2311340897998d0
i = 28; iri(i) = 9;  jri(i) = 13; nri(i) = 1.9809712802088d-08
i = 29; iri(i) = 10; jri(i) = 4;  nri(i) = 1.0406965210174d-19
i = 30; iri(i) = 10; jri(i) = 10; nri(i) = -1.0234747095929d-13
i = 31; iri(i) = 10; jri(i) = 14; nri(i) = -1.0018179379511d-09
i = 32; iri(i) = 16; jri(i) = 29; nri(i) = -8.0882908646985d-11
i = 33; iri(i) = 16; jri(i) = 50; nri(i) = 0.10693031879409d0
i = 34; iri(i) = 18; jri(i) = 57; nri(i) = -0.33662250574171d0
i = 35; iri(i) = 20; jri(i) = 20; nri(i) = 8.9185845355421d-25
i = 36; iri(i) = 20; jri(i) = 35; nri(i) = 3.0629316876232d-13
i = 37; iri(i) = 20; jri(i) = 48; nri(i) = -4.2002467698208d-06
i = 38; iri(i) = 21; jri(i) = 21; nri(i) = -5.9056029685639d-26
i = 39; iri(i) = 22; jri(i) = 53; nri(i) = 3.7826947613457d-06
i = 40; iri(i) = 23; jri(i) = 39; nri(i) = -1.2768608934681d-15
i = 41; iri(i) = 24; jri(i) = 26; nri(i) = 7.3087610595061d-29
i = 42; iri(i) = 24; jri(i) = 40; nri(i) = 5.5414715350778d-17
i = 43; iri(i) = 24; jri(i) = 58; nri(i) = -9.436970724121d-07

end subroutine

!==========================================================================
subroutine swap(a, b)
!==========================================================================

implicit none
real*8 a, b, c
c = a
a = b
b = c

end subroutine

!==========================================================================
subroutine sort3up(d1, d2, d3)
!==========================================================================

!SORT D1, D2, D3 INCREASING

implicit none
real*8 d1, d2, d3

if(d2 > d3) call swap(d2, d3)
if(d1 > d2) call swap(d1, d2)
if(d2 > d3) call swap(d2, d3)

end subroutine


end module flu_3a

