module sea_liq_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIUM PROPERTIES OF LIQUID WATER 
!WITH SEAWATER COMPUTED FROM IAPWS-95 AND IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

!RELEASE ON THE IAPWS FORMULATION 2008 FOR THE THERMODYNAMIC PROPERTIES OF SEAWATER
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     FLU_1,       FILE FLU_1.F90
!     SAL_1,       FILE SAL_1.F90
!     FLU_2,       FILE FLU_2.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_3A,      FILE FLU_3A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90

!=========================================================================

use constants_0
use flu_1
use sal_1
use flu_2
use sal_2
use flu_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_liq_osmoticpressure_si, set_it_ctrl_sea_liq, &
          set_sea_liq_eq_at_s_t_p

!CONTROL PARAMETERS OF THE SEAWATER-LIQUID ITERATION
real*8, private :: ctrl_init_d_liq, ctrl_eps_exit_p

integer, private :: ctrl_initialized, ctrl_loop_maximum, ctrl_mode_liquid

!PROPERTIES OF THE CURRENT SEAWATER-LIQUID EQUILIBRIUM STATE

integer, private :: equi_sea_liq_done
real*8, private :: equi_sea_liq_d_liq, equi_sea_liq_d_sea, equi_sea_liq_s, &
        equi_sea_liq_t, equi_sea_liq_p_liq, equi_sea_liq_p_sea

contains

!==========================================================================
function sea_liq_osmoticpressure_si(sa_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE OSMOTIC PRESSURE IN PA AT THE ABSOLUTE TEMPERATURE T_SI IN K,
!THE ABSOLUTE PRESSURE OF LIQUID WATER P_SI IN PA AND THE SEAWATER SALINITY SA_SI IN KG/KG

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_LIQ_OSMOTICPRESSURE_SI(0.035, 300, 1E5) = 2594603.20968

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE
 
implicit none
real*8 sea_liq_osmoticpressure_si, sa_si, t_si, p_si

sea_liq_osmoticpressure_si = errorreturn

if(set_sea_liq_eq_at_s_t_p(sa_si, t_si, p_si) == errorreturn) return

sea_liq_osmoticpressure_si = equi_sea_liq_p_sea - equi_sea_liq_p_liq

end function

!==========================================================================
function set_sea_liq_eq_at_s_t_p(sa_si, t_si, p_liq_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-WATER EQUILIBRIUM AT GIVEN TEMPERATURE AND SALINITY
!BY NEWTON ITERATION

implicit none

integer maxit
real*8 set_sea_liq_eq_at_s_t_p, sa_si, t_si, p_liq_si
real*8 eps
real*8 d_liq
real*8 d_sea, p_sea

if(equi_sea_liq_done == isok .and. sa_si == equi_sea_liq_s &
  .and. t_si == equi_sea_liq_t .and. p_liq_si == equi_sea_liq_p_liq) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_liq_eq_at_s_t_p = isok
  return
end if

call clear_sea_liq_state !WASTE ANY PREVIOUS STATE

set_sea_liq_eq_at_s_t_p = errorreturn

if(sa_si <= 0d0) return
if(sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_liq_si <= 0d0) return

d_liq = liq_density_si(t_si, p_liq_si)
if(d_liq == errorreturn) return

call init_it_ctrl_sea_liq

!SET INITIAL DENSITY GUESS FOR WATER UNDER SEAWATER PRESSURE
select case (ctrl_mode_liquid)
  case (0);  d_sea = d_liq
  case (1);  d_sea = ctrl_init_d_liq
  case default; return
end select
if(d_sea <= 0d0) return
if(d_sea == errorreturn) return

!SET MAX. ITERATION NUMBER
  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1))then
    p_sea = d_sea ** 2 * flu_f_si(0, 1, t_si, d_sea)
    call set_sea_liq_state (d_liq, d_sea, sa_si, t_si, p_liq_si, p_sea)
    set_sea_liq_eq_at_s_t_p = isok
    return
  elseif(ctrl_loop_maximum > 0)then
    maxit = ctrl_loop_maximum
  else
    return
  endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_liq_iteration_at_s_t_d(sa_si, t_si, d_liq, maxit, eps, d_sea, p_sea) == errorreturn) then
  return
end if

call set_sea_liq_state (d_liq, d_sea, sa_si, t_si, p_liq_si, p_sea)

set_sea_liq_eq_at_s_t_p = isok

end function

!==========================================================================
function sea_liq_iteration_at_s_t_d(sa_si, t_si, d_liq_si, maxit, eps, d_sea_si, p_sea_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-WATER PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN SALINITY, SA_SI, AND TEMPERATURE, T_SI,
!AND LIQUID DENSITY, D_LIQ_SI, FROM AN INITIAL GUESS FOR THE DENSITY, D_SEA_SI.
!D_LIQ_SI IS THE DENSITY IN KG/M3 OF LIQUID WATER UNDER THE PRESSURE OF THE PURE WATER COMPARTMENT.
!D_SEA_SI IS THE DENSITY IN KG/M3 OF LIQUID WATER UNDER THE PRESSURE OF THE SEAWATER COMPARTMENT.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT: 
!        SEA_LIQ_ITERATION_AT_S_T_D = ISOK IF SUCCESSFULLY DONE
!        SEA_LIQ_ITERATION_AT_S_T_D = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
! D_SEA_SI: DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SEA_SI) 
! P_SEA_SI: PRESSURE IN PA ON SEAWATER = PRESSURE ON WATER + OSMOTIC PRESSURE

!INPUT: 
!       SA_SI: ABSOLUTE SALINITY IN KG/KG
!        T_SI: ABSOLUTE TEMPERATURE IN K
!    D_LIQ_SI: DENSITY OF LIQUID WATER IN KG/M3
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!              EPS < 0: RELATIVE PRESSURE TOLERANCE
!    D_SEA_SI: INITIAL GUESS FOR DENSITY OF WATER IN KG/M3 UNDER OSMOTIC PRESSURE

implicit none

integer it, maxit

real*8 sea_liq_iteration_at_s_t_d, sa_si, t_si, d_liq_si, eps, d_sea_si, p_sea_si
real*8 p, d, dd
real*8 mus, p_old
real*8 gl, gw
real*8 gs_p, gs_sp
real*8 f, f_d, f_dd

sea_liq_iteration_at_s_t_d = errorreturn

if(d_sea_si <= 0d0) then
  d_sea_si = errorreturn
  p_sea_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    d_sea_si = errorreturn
    p_sea_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax) then
    d_sea_si = errorreturn
    p_sea_si = errorreturn
    return
  endif
else
  if((d_liq_si <= 0d0) .or. &
    (t_si <= 0d0) .or. &
    (sa_si <= 0d0) .or. (sa_si >= 1d0)) then
    d_sea_si = errorreturn
    p_sea_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

gw = flu_gibbs_energy_si(t_si, d_liq_si)
if(gw == errorreturn) go to 999
d = d_sea_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  f = flu_f_si(0, 0, t_si, d)
  if(f == errorreturn) go to 999
  f_d = flu_f_si(0, 1, t_si, d)
  if(f_d == errorreturn) go to 999
  f_dd = flu_f_si(0, 2, t_si, d)
  if(f_dd == errorreturn) go to 999
  p_old = p
  p = d ** 2 * f_d
  gl = f + d * f_d

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(sa_si, t_si, p) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  gs_p = sal_g_si(0, 0, 1, sa_si, t_si, p)
  if(gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, sa_si, t_si, p)
  if(gs_sp == errorreturn) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0 .and. abs(p - p_old) < eps) .or. &
       (eps < 0 .and. abs(p - p_old) < -eps * abs(p))) then
      if(p <= 0) go to 999
      sea_liq_iteration_at_s_t_d = isok
      p_sea_si = p
      d_sea_si = d
      go to 999
    end if
  end if

  if(it == maxit) go to 999

  !DENSITY ITERATION STEP
  dd = (2d0 * f_d + d * f_dd) * (1d0 + d * (gs_p - sa_si * gs_sp))
  if(dd == 0d0) go to 999
  dd = (gw - gl - mus) / dd

  !UPDATE LIQUID DENSITY
  d = d + dd
  if(d <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_sea_si <= flu_dmin .or. d_sea_si > flu_dmax) then
    d_sea_si = errorreturn
    p_sea_si = errorreturn
    sea_liq_iteration_at_s_t_d = errorreturn
    return
  endif
  !SAL_LIMITS
  if(p_sea_si < sal_pmin .or. p_sea_si > sal_pmax) then
    d_sea_si = errorreturn
    p_sea_si = errorreturn
    sea_liq_iteration_at_s_t_d = errorreturn
    return
  endif
endif

end function

!==========================================================================
subroutine set_it_ctrl_sea_liq(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!SEAWATER-WATER EQUILIBRIUM

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VALUES TO COMPUTE EQULIBRIUM

!'INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY TO START =
!                                   LIQUID DENSITY IAPWS95, D = LIQ_DENSITY_SI(T,P)
!'INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!'TOL_LIQ_PRESS'        0           USE DEFAULT EXIT ACCURACY FOR SEAWATER PRESSURE (0.1 PPM)
!'TOL_LIQ_PRESS'        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

implicit none

character*13 key
real*8 value

call init_it_ctrl_sea_liq

call clear_sea_liq_state

select case (trim(key))

  case ('it_steps')   !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100  !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')   !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_liquid = 0    !DEFAULT = PURE WATER DENSITY
    elseif(value < 0d0)then
      !ignore it
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ('tol_liq_press')      !REQUIRED LIQUID PRESSURE TOLERANCE
    if(value == 0d0)then
      ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine clear_sea_liq_state()
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

equi_sea_liq_done = 0

end subroutine

!==========================================================================
subroutine init_it_ctrl_sea_liq()
!==========================================================================

implicit none

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_mode_liquid = 0
ctrl_loop_maximum = 100
ctrl_init_d_liq = 1000
ctrl_eps_exit_p = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

!==========================================================================
subroutine set_sea_liq_state(d_liq, d, s, t, p_liq, p)
!==========================================================================

implicit none
real*8 d_liq, d, s, t, p_liq, p

equi_sea_liq_d_liq = d_liq    !EQUILIBRIUM DENSITY OF LIQUID WATER
equi_sea_liq_d_sea = d        !EQUILIBRIUM DENSITY OF LIQUID WATER UNDER SEAWATER PRESSURE
equi_sea_liq_s = s            !EQUILIBRIUM SALINITY
equi_sea_liq_t = t            !EQUILIBRIUM TEMPERATURE
equi_sea_liq_p_liq = p_liq    !EQUILIBRIUM PRESSURE OF LIQUID WATER
equi_sea_liq_p_sea = p        !EQUILIBRIUM PRESSURE OF SEAWATER

end subroutine

end module sea_liq_4
