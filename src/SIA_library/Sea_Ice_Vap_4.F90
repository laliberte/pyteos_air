module sea_ice_vap_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIUM PROPERTIES OF SEA ICE WITH VAPOUR
!COMPUTED FROM IAPWS-95, IAPWS-06 AND IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

!RELEASE ON AN EQUATION OF STATE FOR H2O ICE IH
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!WITNEY, UK, SEPTEMBER 2006

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
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     FLU_1,         FILE FLU_1.F90
!     ICE_1,         FILE ICE_1.F90
!     SAL_1,         FILE SAL_1.F90
!     SAL_2,         FILE SAL_2.F90

!#########################################################################

use constants_0
use maths_0
use flu_1
use ice_1
use sal_1
use sal_2

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_ice_vap_density_vap_si, sea_ice_vap_pressure_si, &
          sea_ice_vap_salinity_si, sea_ice_vap_temperature_si,  &
          set_it_ctrl_sea_ice_vap, set_sea_ice_vap_eq_at_s, &
          set_sea_ice_vap_eq_at_t, set_sea_ice_vap_eq_at_p

!INTEGER, PRIVATE :: CTRL_INITIALIZED, CTRL_LOOP_MAXIMUM, CTRL_MODE_LIQUID
integer, private :: ctrl_loop_maximum

!CONTROL PARAMETERS OF THE SEA-ICE-VAP ITERATION
integer :: &
ctrl_initialized, ctrl_mode_salinity, &
ctrl_mode_temperature, ctrl_mode_pressure, &
ctrl_mode_liquid, ctrl_mode_vapour

real*8 :: ctrl_init_s, ctrl_init_t, ctrl_init_p, ctrl_init_d_liq, &  
  ctrl_init_d_vap, ctrl_eps_exit_t, ctrl_eps_exit_p 

!PROPERTIES OF THE CURRENT SEA-ICE-VAPOUR EQUILIBRIUM STATE
integer :: equi_sea_ice_vap_done  !FLAG = ISOK IF EQUILIBRIUM WAS SUCCESSFULLY COMPUTED

real*8 :: equi_sea_ice_vap_t, equi_sea_ice_vap_p, equi_sea_ice_vap_s, &  
  equi_sea_ice_vap_d_liq, equi_sea_ice_vap_d_vap 

contains

!==========================================================================
function set_sea_ice_vap_eq_at_s(sa_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEA-ICE-VAPOUR EQUILIBRIUM AT GIVEN SALINITY SA_SI IN KG/KG
!BY NEWTON ITERATION

real*8 :: set_sea_ice_vap_eq_at_s, sa_si
real*8 :: tt, pt, dt_liq, dt_vap
real*8 :: t, p, dl, dv, eps
integer maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dt_vap = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR

if (equi_sea_ice_vap_done == isok .and. &
   sa_si == equi_sea_ice_vap_s) then
  !the requested state has already been computed earlier
  set_sea_ice_vap_eq_at_s = isok
  return
end if

call clear_sea_ice_vap_state !waste any previous state

set_sea_ice_vap_eq_at_s = errorreturn

if (sa_si < 0d0) return
if (sa_si >= 1d0) return

if (sa_si == 0) then  !PURE WATER: SET IAPWS-95 TRIPLE POINT
  call set_sea_ice_vap_state(sa_si, tt, pt, dt_vap, dt_liq)
  set_sea_ice_vap_eq_at_s = isok
  return
end if

call init_it_ctrl_sea_ice_vap

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case (0);  t = aux_temperature_si(sa_si)
  case (-1); t = tt
  case (1);  t = ctrl_init_t
  case default; return
end select
if (t <= 0) return
if (t == errorreturn) return

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case (0);  p = aux_pressure_si(sa_si)
  case (-1); p = pt
  case (1);  p = ctrl_init_p
  case default; return
end select
if (p <= 0d0) return
if (p == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case (0);  dl = dt_liq
  case (-1); dl = dt_liq
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select
if (dl <= 0d0) return
if (dl == errorreturn) return

!SET INITIAL VAPOUR DENSITY GUESS
select case(ctrl_mode_vapour)
  case (0);  dv = aux_density_ideal_si(t, p)
  case (-1); dv = dt_vap
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select
if (dv <= 0d0) return
if (dv == errorreturn) return

!SET MAX. ITERATION NUMBER
select case(ctrl_loop_maximum)
  case (0); maxit = 100
  case (-1); call set_sea_ice_vap_state(sa_si, t, p, dv, dl)
           set_sea_ice_vap_eq_at_s = isok
           return
  case (1:); maxit = ctrl_loop_maximum
  case default; return
end select

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_t
if (eps == 0d0) return

!RUN ITERATION LOOP
if (sea_ice_vap_iteration_at_s(sa_si, maxit, eps, t, p, dv, dl) &
   == errorreturn) return

call set_sea_ice_vap_state(sa_si, t, p, dv, dl)

set_sea_ice_vap_eq_at_s = isok

end function

!==========================================================================
function set_sea_ice_vap_eq_at_t(t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEA-ICE-VAPOUR EQUILIBRIUM AT GIVEN ABSOLUTE TEMPERATURE T_SI IN K
!BY NEWTON ITERATION

real*8 :: set_sea_ice_vap_eq_at_t, t_si
real*8 :: tt, pt, dt_liq, dt_vap
real*8 :: s, p, dl, dv, eps
integer :: maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dt_vap = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR
!

if (equi_sea_ice_vap_done == isok .and. &
  t_si == equi_sea_ice_vap_t) then
  !the requested state has already been computed earlier
  set_sea_ice_vap_eq_at_t = isok
  return
end if

call clear_sea_ice_vap_state !waste any previous state

set_sea_ice_vap_eq_at_t = errorreturn

if (t_si <= 0d0) return
if (t_si > tt) return

if (t_si == tt) then  !pure water: set iapws-95 triple point
  call set_sea_ice_vap_state(0d0, tt, pt, dt_vap, dt_liq)
  set_sea_ice_vap_eq_at_t = isok
  return
end if

call init_it_ctrl_sea_ice_vap

!SET INITIAL SALINITY GUESS
select case(ctrl_mode_salinity)
  case (0);  s = aux_salinity_of_t_si(t_si)
  case (-1); s = so_salinity_si
  case (1);  s = ctrl_init_s
  case default; return
end select
if (s <= 0d0 .or. s >= 1d0) return
if (s == errorreturn) return

!SET INITIAL PRESSURE GUESS
select case(ctrl_mode_pressure)
  case (0);  p = aux_pressure_si(s)
  case (-1); p = pt
  case (1);  p = ctrl_init_p
  case default; return
end select
if (p <= 0d0) return
if (p == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case (0);  dl = dt_liq
  case (-1); dl = dt_liq
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select
if (dl <= 0d0) return
if (dl == errorreturn) return

!SET INITIAL VAPOUR DENSITY GUESS
select case(ctrl_mode_vapour)
  case (0);  dv = aux_density_ideal_si(t_si, p)
  case (-1); dv = dt_vap
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select
if (dv <= 0d0) return
if (dv == errorreturn) return

!SET MAX. ITERATION NUMBER
select case(ctrl_loop_maximum)
  case (0); maxit = 100
  case (-1); call set_sea_ice_vap_state(s, t_si, p, dv, dl)
             set_sea_ice_vap_eq_at_t = isok
             return
  case (1:); maxit = ctrl_loop_maximum
  case default; return
end select

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if (eps == 0d0) return

!RUN ITERATION LOOP
if (sea_ice_vap_iteration_at_t(t_si, maxit, eps, s, p, dv, dl) &
   == errorreturn) return

call set_sea_ice_vap_state(s, t_si, p, dv, dl)

set_sea_ice_vap_eq_at_t = isok

end function

!==========================================================================
function set_sea_ice_vap_eq_at_p(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEA-ICE-VAPOUR EQUILIBRIUM AT GIVEN ABSOLUTE PRESSURE P_SI IN PA
!BY NEWTON ITERATION

real*8 :: set_sea_ice_vap_eq_at_p, p_si
real*8 :: tt, pt, dt_liq, dt_vap, s, t, dl, dv, eps
integer :: maxit

tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dt_liq = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dt_vap = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR

if (equi_sea_ice_vap_done == isok .and. &
   p_si == equi_sea_ice_vap_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_ice_vap_eq_at_p = isok
  return
end if

call clear_sea_ice_vap_state !waste any previous state

set_sea_ice_vap_eq_at_p = errorreturn

if (p_si <= 0d0) return
if (p_si > pt) return

if (p_si == pt) then  !pure water: set iapws-95 triple point
  call set_sea_ice_vap_state(0d0, tt, pt, dt_vap, dt_liq)
  set_sea_ice_vap_eq_at_p = isok
  return
end if

call init_it_ctrl_sea_ice_vap

!SET INITIAL SALINITY GUESS
select case(ctrl_mode_salinity)
  case (0);  s = aux_salinity_of_p_si(p_si)
  case (-1); s = so_salinity_si
  case (1);  s = ctrl_init_s
  case default; return
end select
if (s <= 0d0 .or. s >= 1d0) return
if (s == errorreturn) return

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case (0);  t = aux_temperature_si(s)
  case (-1); t = tt
  case (1);  t = ctrl_init_t
  case default; return
end select
if (t <= 0d0) return
if (t == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case(ctrl_mode_liquid)
  case (0);  dl = dt_liq
  case (-1); dl = dt_liq
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select
if (dl <= 0d0) return
if (dl == errorreturn) return

!SET INITIAL VAPOUR DENSITY GUESS
select case(ctrl_mode_vapour)
  case (0);  dv = aux_density_ideal_si(t, p_si)
  case (-1); dv = dt_vap
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select
if (dv <= 0d0) return
if (dv == errorreturn) return

!SET MAX. ITERATION NUMBER
select case(ctrl_loop_maximum)
  case (0); maxit = 100
  case (-1); call set_sea_ice_vap_state(s, t, p_si, dv, dl)
           set_sea_ice_vap_eq_at_p = isok
           return
  case (1:); maxit = ctrl_loop_maximum
  case default; return
end select

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_t
if (eps == 0d0) return

!RUN ITERATION LOOP
if (sea_ice_vap_iteration_at_p(p_si, maxit, eps, s, t, dv, dl) &
   == errorreturn) return

call set_sea_ice_vap_state(s, t, p_si, dv, dl)

set_sea_ice_vap_eq_at_p = isok

end function

!==========================================================================
function sea_ice_vap_salinity_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE BRINE SALINITY IN KG/KG OF SEA ICE AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_VAP_EQUILIBRIUM_TO_S,
!SET_SEA_ICE_VAP_EQUILIBRIUM_TO_T OR SET_SEA_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_VAP_EQ_AT_S 0.035
!SEA_ICE_VAP_SALINITY_SI = 0.035

!SET_SEA_ICE_VAP_EQ_AT_T 270
!SEA_ICE_VAP_SALINITY_SI = 5.61489288506E-02

!SET_SEA_ICE_VAP_EQ_AT_P 500
!SEA_ICE_VAP_SALINITY_SI = 4.38955878828E-02

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 :: sea_ice_vap_salinity_si

sea_ice_vap_salinity_si = errorreturn

if (equi_sea_ice_vap_done /= isok) return

sea_ice_vap_salinity_si = equi_sea_ice_vap_s

end function

!==========================================================================
function sea_ice_vap_density_vap_si()
!==========================================================================

!THIS FUNCTION RETURNS THE VAPOUR DENSITY IN KG M-3 OF SEA ICE AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_VAP_EQUILIBRIUM_TO_S,
!SET_SEA_ICE_VAP_EQUILIBRIUM_TO_T OR SET_SEA_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_VAP_EQ_AT_S 0.035
!SEA_ICE_VAP_DENSITY_VAP_SI = 4.17156419318E-03

!SET_SEA_ICE_VAP_EQ_AT_T 270
!SEA_ICE_VAP_DENSITY_VAP_SI = 3.77406140772E-03

!SET_SEA_ICE_VAP_EQ_AT_P 500
!SEA_ICE_VAP_DENSITY_VAP_SI = 4.00364833230E-03

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 :: sea_ice_vap_density_vap_si

sea_ice_vap_density_vap_si = errorreturn

if (equi_sea_ice_vap_done /= isok) return

sea_ice_vap_density_vap_si = equi_sea_ice_vap_d_vap

end function

!==========================================================================
function sea_ice_vap_temperature_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE TEMPERATURE IN K OF SEA ICE AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_VAP_EQUILIBRIUM_TO_S,
!SET_SEA_ICE_VAP_EQUILIBRIUM_TO_T OR SET_SEA_ICE_VAP_EQUILIBRIUM_TO_P

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_VAP_EQ_AT_S 0.035
!SEA_ICE_VAP_TEMPERATURE_SI = 271.247815057

!SET_SEA_ICE_VAP_EQ_AT_T 270
!SEA_ICE_VAP_TEMPERATURE_SI = 270

!SET_SEA_ICE_VAP_EQ_AT_P 500
!SEA_ICE_VAP_TEMPERATURE_SI = 270.734430917

real*8 :: sea_ice_vap_temperature_si

sea_ice_vap_temperature_si = errorreturn

if (equi_sea_ice_vap_done /= isok) return

sea_ice_vap_temperature_si = equi_sea_ice_vap_t

end function

!==========================================================================
function sea_ice_vap_pressure_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE PRESSURE IN PA OF SEA ICE AT EQUILIBRIUM WITH VAPOUR,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_VAP_EQUILIBRIUM_TO_S,
!SET_SEA_ICE_VAP_EQUILIBRIUM_TO_T OR SET_SEA_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_VAP_EQ_AT_S 0.035
!SEA_ICE_VAP_PRESSURE_SI = 521.950349225

!SET_SEA_ICE_VAP_EQ_AT_T 270
!SEA_ICE_VAP_PRESSURE_SI = 470.059067981

!SET_SEA_ICE_VAP_EQ_AT_P 500
!SEA_ICE_VAP_PRESSURE_SI = 500

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

real*8 :: sea_ice_vap_pressure_si

sea_ice_vap_pressure_si = errorreturn

if (equi_sea_ice_vap_done /= isok) return

sea_ice_vap_pressure_si = equi_sea_ice_vap_p

end function

!==========================================================================
function sea_ice_vap_iteration_at_s(sa_si, maxit, eps, t_si, &
                                    p_si, d_vap_si, d_liq_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN SALINITY, SA_SI,
!FROM INITIAL GUESSES FOR THE TEMPERATURE, T_SI, PRESSURE, P_SI,
!THE VAPOUR DENSITY, DVAP_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN TEMPERATURE.

!OUTPUT: 
!        SEA_ICE_VAP_ITERATION_AT_S = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_VAP_ITERATION_AT_S = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, TEMPERATURE OR PRESSURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - S-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
!  T_SI: TEMPERATURE IN K OF THE SEAWATER-ICE EQUILIBRIUM
!  P_SI: PRESSURE IN PA OF THE SEAWATER-ICE EQUILIBRIUM
!  DVAP_SI: DENSITY OF VAPOUR IN KG/M3 AT (T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID IN KG/M3 AT (T_SI, P_SI)

!INPUT: 
!       SA_SI: ABSOLUTE SALINITY IN KG/KG
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF TEMPERATURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE
!        T_SI: INITIAL GUESS FOR ABSOLUTE (FREEZING) TEMPERATURE IN K
!        P_SI: INITIAL GUESS FOR ABSOLUTE PRESSURE IN PA
!        DVAP_SI: INITIAL GUESS FOR DENSITY OF VAPOUR IN KG/M3 AT (T_SI, P_SI)
!        DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID IN KG/M3 AT (T_SI, P_SI)

real*8 :: sea_ice_vap_iteration_at_s, sa_si, eps, t_si, p_si, d_vap_si, d_liq_si
real*8 :: dv, dl, t, p, ddv, ddl, dt, dp, gi, gi_t, gi_p, &
       gs, gs_s, gs_t, gs_st, gs_p, gs_sp, fv, fv_t, fv_td, &
       fv_d, fv_dd, fl, fl_t, fl_td, fl_d, fl_dd
real*8 :: a(4, 4), b(4), x(4)
integer :: maxit, it

sea_ice_vap_iteration_at_s = errorreturn

if(d_vap_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  t_si <= 0d0 .or. &
  p_si <= 0d0) then
  t_si = errorreturn
  p_si = errorreturn
  d_liq_si = errorreturn 
  d_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn 
    d_vap_si = errorreturn
    return
  endif
else
  if ((sa_si <= 0d0) .or. (sa_si >= 1d0)) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn 
    d_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
p = p_si
dl = d_liq_si
dv = d_vap_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF WATER VAPOUR FOR 4D NEWTON ITERATION
  fv = flu_f_si(0, 0, t, dv)
  if (fv == errorreturn) go to 999
  fv_t = flu_f_si(1, 0, t, dv)
  if (fv_t == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t, dv)
  if (fv_d == errorreturn) go to 999
  fv_td = flu_f_si(1, 1, t, dv)
  if (fv_td == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t, dv)
  if (fv_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if (fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if (fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if (fl_d == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if (fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if (fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p)
  if (gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p)
  if (gi_t == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t, p)
  if (gi_p == errorreturn) go to 999

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  gs = sal_g_si(0, 0, 0, sa_si, t, p)
  if (gs == errorreturn) go to 999
  gs_s = sal_g_si(1, 0, 0, sa_si, t, p)
  if (gs_s == errorreturn) go to 999
  gs_t = sal_g_si(0, 1, 0, sa_si, t, p)
  if (gs_t == errorreturn) go to 999
  gs_st = sal_g_si(1, 1, 0, sa_si, t, p)
  if (gs_st == errorreturn) go to 999
  gs_p = sal_g_si(0, 0, 1, sa_si, t, p)
  if (gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, sa_si, t, p)
  if (gs_sp == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = -fl_t - dl * fl_td + fv_t + dv * fv_td - gs_t + sa_si * gs_st
  a(1, 2) = -gs_p + sa_si * gs_sp
  a(1, 3) = 2d0 * fv_d + dv * fv_dd
  a(1, 4) = -2d0 * fl_d - dl * fl_dd
  a(2, 1) = fv_t + dv * fv_td - gi_t
  a(2, 2) = -gi_p
  a(2, 3) = 2d0 * fv_d + dv * fv_dd
  a(2, 4) = 0d0
  a(3, 1) = -dv * fv_td
  a(3, 2) = 1d0 / dv
  a(3, 3) = -2d0 * fv_d - dv * fv_dd
  a(3, 4) = 0d0
  a(4, 1) = -dl * fl_td
  a(4, 2) = 1d0 / dl
  a(4, 3) = 0d0
  a(4, 4) = -2d0 * fl_d - dl * fl_dd

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + dl * fl_d - fv - dv * fv_d + gs - sa_si * gs_s
  b(2) = gi - fv - dv * fv_d
  b(3) = dv * fv_d - p / dv
  b(4) = dl * fl_d - p / dl

  !SOLVE EQUATIONS
  if (matrix_solve(a, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  dt = x(1)
  dp = x(2)
  ddv = x(3)
  ddl = x(4)
  
  !UPDATE UNKNOWNS
  dv = dv + ddv
  if (dv <= 0d0) go to 999
  dl = dl + ddl
  if (dl <= 0d0) go to 999
  t = t + dt
  if (t <= 0d0) go to 999
  p = p + dp
  if (p <= 0d0) go to 999
  
  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if ((eps > 0d0 .and. abs(dt) < eps) .or. &
     (eps < 0d0 .and. abs(dt) < -eps * abs(t))) then
    sea_ice_vap_iteration_at_s = isok
    t_si = t
    p_si = p
    d_vap_si = dv
    d_liq_si = dl
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_s = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax) then
    t_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_s = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_ice_vap_iteration_at_t(t_si, maxit, eps, sa_si, &
                                    p_si, d_vap_si, d_liq_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN TEMPERATURE, T_SI,
!FROM INITIAL GUESSES FOR THE SALINITY, SA_SI, PRESSURE, P_SI,
!THE VAPOUR DENSITY, DVAP_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT: 
!        SEA_ICE_VAP_ITERATION_AT_T = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_VAP_ITERATION_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, SALNITY OR PRESSURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - S-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
!  SA_SI: ABSOLUTE SALINITY IN KG/KG AT THE SEAWATER-ICE EQUILIBRIUM
!  P_SI:  ABSOLUTE PRESSURE IN PA AT THE SEAWATER-ICE EQUILIBRIUM
!  DVAP_SI: DENSITY OF VAPOUR IN KG/M3 AT (T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID IN KG/M3 AT (T_SI, P_SI)

!INPUT:  
!        T_SI: ABSOLUTE TEMPERATURE IN K
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN K
!              EPS < 0: RELATIVE PRESSURE TOLERANCE
!       SA_SI: INITIAL GUESS FOR ABSOLUTE SALINITY IN KG/KG
!        P_SI: INITIAL GUESS FOR ABSOLUTE PRESSURE IN PA
!        DVAP_SI: INITIAL GUESS FOR DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!        DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)

real*8 :: sea_ice_vap_iteration_at_t, t_si, eps, sa_si, p_si, d_vap_si, d_liq_si
real*8 :: dv, dl, s, p, ddv, ddl, ds, dp, &
          gi, gi_p, gs, gs_s, gs_ss, gs_p, gs_sp, &
          fv, fv_d, fv_dd, fl, fl_d, fl_dd
real*8 :: a(4, 4), b(4), x(4)
integer :: maxit, it

sea_ice_vap_iteration_at_t = errorreturn

if(d_vap_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  p_si <= 0d0 .or. &
  sa_si <= 0d0 .or. sa_si >= 1d0) then
  sa_si = errorreturn
  p_si = errorreturn
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    sa_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
else
  if (t_si <= 0d0) then
    sa_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

s = sa_si
p = p_si
dl = d_liq_si
dv = d_vap_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF WATER VAPOUR FOR 4D NEWTON ITERATION
  fv = flu_f_si(0, 0, t_si, dv)
  if (fv == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t_si, dv)
  if (fv_d == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t_si, dv)
  if (fv_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t_si, dl)
  if (fl == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t_si, dl)
  if (fl_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t_si, dl)
  if (fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t_si, p)
  if (gi == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t_si, p)
  if (gi_p == errorreturn) go to 999

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  gs = sal_g_si(0, 0, 0, s, t_si, p)
  if (gs == errorreturn) go to 999
  gs_s = sal_g_si(1, 0, 0, s, t_si, p)
  if (gs_s == errorreturn) go to 999
  gs_ss = sal_g_si(2, 0, 0, s, t_si, p)
  if (gs_ss == errorreturn) go to 999
  gs_p = sal_g_si(0, 0, 1, s, t_si, p)
  if (gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, s, t_si, p)
  if (gs_sp == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = s * gs_ss
  a(1, 2) = -gs_p + s * gs_sp
  a(1, 3) = 2d0 * fv_d + dv * fv_dd
  a(1, 4) = -2d0 * fl_d - dl * fl_dd
  a(2, 1) = 0d0
  a(2, 2) = -gi_p
  a(2, 3) = 2d0 * fv_d + dv * fv_dd
  a(2, 4) = 0d0
  a(3, 1) = 0d0
  a(3, 2) = 1d0 / dv
  a(3, 3) = -2d0 * fv_d - dv * fv_dd
  a(3, 4) = 0d0
  a(4, 1) = 0d0
  a(4, 2) = 1d0 / dl
  a(4, 3) = 0d0
  a(4, 4) = -2d0 * fl_d - dl * fl_dd

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + dl * fl_d - fv - dv * fv_d + gs - s * gs_s
  b(2) = gi - fv - dv * fv_d
  b(3) = dv * fv_d - p / dv
  b(4) = dl * fl_d - p / dl

  !SOLVE EQUATIONS
  if (matrix_solve(a, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  ds = x(1)
  dp = x(2)
  ddv = x(3)
  ddl = x(4)
  
  !UPDATE UNKNOWNS
  dv = dv + ddv
  if (dv <= 0d0) go to 999
  dl = dl + ddl
  if (dl <= 0d0) go to 999
  s = s + ds
  if (s <= 0d0) go to 999
  if (s >= 1d0) go to 999
  p = p + dp
  if (p <= 0d0) go to 999
  
  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if ((eps > 0d0 .and. abs(dp) < eps) .or. &
     (eps < 0d0 .and. abs(dp) < -eps * abs(p))) then
    sea_ice_vap_iteration_at_t = isok
    sa_si = s
    p_si = p
    d_vap_si = dv
    d_liq_si = dl
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    sa_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_t = errorreturn
    return
  endif
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) then
    sa_si = errorreturn
    p_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_ice_vap_iteration_at_p(p_si, maxit, eps, sa_si, t_si, &
                                            d_vap_si, d_liq_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE-VAPOUR PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER AT GIVEN PRESSURE, P_SI,
!FROM INITIAL GUESSES FOR THE SALINITY, SA_SI, TEMPERATURE, T_SI,
!THE VAPOUR DENSITY, DVAP_SI AND THE LIQUID DENSITY, DLIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN TEMPERATURE.

!OUTPUT: 
!        SEA_ICE_VAP_ITERATION_AT_P = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_VAP_ITERATION_AT_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY, SALNITY OR PRESSURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - S-T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
!    SA_SI: ABSOLUTE SALINITY IN KG/KG AT THE SEAWATER-ICE EQUILIBRIUM
!     T_SI: ABSOLUTE TEMPERATURE IN K AT THE SEAWATER-ICE EQUILIBRIUM
!  DVAP_SI: DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!  DLIQ_SI: DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)

!INPUT:  
!        P_SI: ABSOLUTE PRESSURE IN PA
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF TEMPERATURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE
!       SA_SI: INITIAL GUESS FOR ABSOLUTE SALINITY IN KG/KG
!        T_SI: INITIAL GUESS FOR ABSOLUTE TEMPERATURE IN K
!     DVAP_SI: INITIAL GUESS FOR DENSITY OF WATER VAPOUR IN KG/M3 AT (T_SI, P_SI)
!     DLIQ_SI: INITIAL GUESS FOR DENSITY OF LIQUID WATER IN KG/M3 AT (T_SI, P_SI)

real*8 :: sea_ice_vap_iteration_at_p, p_si, eps, sa_si, t_si, d_vap_si, d_liq_si
real*8 :: dv, dl, s, t, ddv, ddl, ds, dt, &
       gi, gi_t, gs, gs_s, gs_ss, gs_t, gs_st, &
       fv, fv_t, fv_td, fv_d, fv_dd, fl, fl_t, fl_td, fl_d, fl_dd
real*8 :: a(4, 4), b(4), x(4)
integer :: maxit, it

sea_ice_vap_iteration_at_p = errorreturn

if(d_vap_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  t_si <= 0d0 .or. &
  sa_si <= 0d0 .or. sa_si >= 1d0) then
  sa_si = errorreturn
  t_si = errorreturn
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !SAL_LIMITS
  if(p_si < sal_pmin .or. p_si > sal_pmax) then
    sa_si = errorreturn
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
else
  if (p_si <= 0d0) then
    sa_si = errorreturn
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

s = sa_si
t = t_si
dl = d_liq_si
dv = d_vap_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF WATER VAPOUR FOR 4D NEWTON ITERATION
  fv = flu_f_si(0, 0, t, dv)
  if (fv == errorreturn) go to 999
  fv_t = flu_f_si(1, 0, t, dv)
  if (fv_t == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t, dv)
  if (fv_d == errorreturn) go to 999
  fv_td = flu_f_si(1, 1, t, dv)
  if (fv_td == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t, dv)
  if (fv_dd == errorreturn) go to 999

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF LIQUID WATER FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if (fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if (fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if (fl_d == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if (fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if (fl_dd == errorreturn) go to 999

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p_si)
  if (gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p_si)
  if (gi_t == errorreturn) go to 999

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  gs = sal_g_si(0, 0, 0, s, t, p_si)
  if (gs == errorreturn) go to 999
  gs_s = sal_g_si(1, 0, 0, s, t, p_si)
  if (gs_s == errorreturn) go to 999
  gs_ss = sal_g_si(2, 0, 0, s, t, p_si)
  if (gs_ss == errorreturn) go to 999
  gs_t = sal_g_si(0, 1, 0, s, t, p_si)
  if (gs_t == errorreturn) go to 999
  gs_st = sal_g_si(1, 1, 0, s, t, p_si)
  if (gs_st == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = s * gs_ss
  a(1, 2) = -fl_t - dl * fl_td + fv_t + dv * fv_td - gs_t + sa_si * gs_st
  a(1, 3) = 2d0 * fv_d + dv * fv_dd
  a(1, 4) = -2d0 * fl_d - dl * fl_dd
  a(2, 1) = 0d0
  a(2, 2) = fv_t + dv * fv_td - gi_t
  a(2, 3) = 2d0 * fv_d + dv * fv_dd
  a(2, 4) = 0d0
  a(3, 1) = 0d0
  a(3, 2) = -dv * fv_td
  a(3, 3) = -2d0 * fv_d - dv * fv_dd
  a(3, 4) = 0d0
  a(4, 1) = 0d0
  a(4, 2) = -dl * fl_td
  a(4, 3) = 0d0
  a(4, 4) = -2d0 * fl_d - dl * fl_dd

  !RIGHT-HAND SIDES, MUST VANISH AT EQUILIBRIUM
  b(1) = fl + dl * fl_d - fv - dv * fv_d + gs - s * gs_s
  b(2) = gi - fv - dv * fv_d
  b(3) = dv * fv_d - p_si / dv
  b(4) = dl * fl_d - p_si / dl

  !SOLVE EQUATIONS
  if (matrix_solve(a, b, x, 4) /= 0) go to 999 !MATRIX SINGULAR
  ds = x(1)
  dt = x(2)
  ddv = x(3)
  ddl = x(4)
  
  !UPDATE UNKNOWNS
  dv = dv + ddv
  if (dv <= 0d0) go to 999
  dl = dl + ddl
  if (dl <= 0d0) go to 999
  s = s + ds
  if (s <= 0d0) go to 999
  if (s >= 1d0) go to 999
  t = t + dt
  if (t <= 0d0) go to 999
  
  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if ((eps > 0d0 .and. abs(dt) < eps) .or. &
     (eps < 0d0 .and. abs(dt) < -eps * abs(t))) then
    sea_ice_vap_iteration_at_p = isok
    sa_si = s
    t_si = t
    d_vap_si = dv
    d_liq_si = dl
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    sa_si = errorreturn
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_p = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax) then
    sa_si = errorreturn
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    sea_ice_vap_iteration_at_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_temperature_si(sa_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED FREEZING TEMPERATURE IN K AT GIVEN SALINITY IN KG/KG

real*8 :: aux_temperature_si, sa_si, tt

tt = tp_temperature_si

aux_temperature_si = tt - 60d0 * sa_si

end function

!==========================================================================
function aux_salinity_of_t_si(t_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED SALINITY IN KG/KG AT GIVEN TEMPERATURE IN K

real*8 :: aux_salinity_of_t_si, t_si, tt

tt = tp_temperature_si

aux_salinity_of_t_si = (tt - t_si) / 60d0

end function

!==========================================================================
function aux_pressure_si(sa_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED FREEZING PRESSURE IN PA AT GIVEN SALINITY IN KG/KG

real*8 :: aux_pressure_si, sa_si, pt

pt = tp_pressure_iapws95_si

aux_pressure_si = pt - 3000d0 * sa_si

end function

!==========================================================================
function aux_salinity_of_p_si(p_si)
!==========================================================================

!FUNCTION RETURNS THE ESTIMATED FREEZING SALINITY IN KG/KG AT GIVEN PRESSURE IN PA

real*8 :: aux_salinity_of_p_si, p_si, pt

pt = tp_pressure_iapws95_si

aux_salinity_of_p_si = (pt - p_si) / 3000d0

end function

!==========================================================================
function aux_density_ideal_si(t_si, p_si)
!==========================================================================

!FUNCTION RETURNS THE IDEAL-GAS DENSITY IN KG M-3 AT GIVEN
!ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA

real*8 :: aux_density_ideal_si, t_si, p_si, pt, rw

pt = tp_pressure_iapws95_si
rw = gas_constant_h2o_si

aux_density_ideal_si = errorreturn
if (t_si <= 0d0) return
if (p_si <= 0d0) return

aux_density_ideal_si = p_si / (rw * t_si)

end function

!==========================================================================
subroutine init_it_ctrl_sea_ice_vap()
!==========================================================================

real*8 :: s_norm, tt, pt, dlt, dvt

s_norm = so_salinity_si  !KG/KG OF KCL NORMALIZED SEAWATER
tt = tp_temperature_si
pt = tp_pressure_iapws95_si
dlt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER
dvt = tp_density_vap_iapws95_si  !TRIPLE POINT DENSITY OF WATER VAPOUR

if (ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR THE ITERATION
ctrl_mode_liquid = 0
ctrl_mode_vapour = 0
ctrl_mode_salinity = 0
ctrl_mode_temperature = 0
ctrl_mode_pressure = 0
ctrl_loop_maximum = 100
ctrl_init_s = s_norm
ctrl_init_t = tt
ctrl_init_p = pt
ctrl_init_d_liq = dlt
ctrl_init_d_vap = dvt
!
ctrl_eps_exit_p = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_t = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

!==========================================================================
subroutine set_it_ctrl_sea_ice_vap(key, value)

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!SEA-ICE-VAPOUR EQUILIBRIUM

!KEY              VALUE
!IT_STEPS         0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS         N > 0       SET ITERATION NUMBER TO N
!IT_STEPS        -1           DO NOT ITERATE, USE INITIAL VALUES TO COMPUTE EQULIBRIUM

!INIT_LIQ_DENS    0           USE DEFAULT LIQUID DENSITY TO START, D = TRIPLE POINT LIQUID DENSITY
!INIT_LIQ_DENS   -1           USE TRIPLE POINT LIQUID DENSITY TO START
!INIT_LIQ_DENS    D > 0       USE VALUE D AS LIQUID DENSITY TO START

!INIT_VAP_DENS    0           USE DEFAULT VAPOUR DENSITY TO START, D = AUX_DENSITY_IDEAL_SI(T, P)
!INIT_VAP_DENS   -1           USE TRIPLE POINT VAPOUR DENSITY TO START
!INIT_VAP_DENS    D > 0       USE VALUE D AS VAPOUR DENSITY TO START

!INIT_BRINE_SA    0           USE DEFAULT BRINE SALINITY TO START (FITTED GUESS FUNCTION AUX_SALINITY)
!INIT_BRINE_SA   -1           USE NORMAL SALINITY TO START
!INIT_BRINE_SA    S > 0       USE VALUE S AS TEMPERATURE TO START

!INIT_TEMP        0           USE DEFAULT TEMPERATURE TO START (FITTED GUESS FUNCTION AUX_TEMPERATURE)
!INIT_TEMP       -1           USE TIPLE POINT TEMPERATURE TO START
!INIT_TEMP        T > 0       USE VALUE T AS TEMPERATURE TO START

!INIT_PRESS       0           USE DEFAULT PRESSURE TO START (FITTED GUESS FUNCTION AUX_PRESSURE)
!INIT_PRESS      -1           USE TRIPLE POINT PRESSURE TO START
!INIT_PRESS       P > 0       USE VALUE P AS PRESSURE TO START

!TOL_TEMP         0           USE DEFAULT EXIT ACCURACY FOR FREEZING TEMPERATURE (0.1 PPM)
!TOL_TEMP         EPS         USE EPS AS EXIT ACCURACY FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_PRESS        0           USE DEFAULT EXIT ACCURACY FOR MELTING PRESSURE (0.1 PPM)
!TOL_PRESS        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

real*8 :: value
character*13 key

call init_it_ctrl_sea_ice_vap

call clear_sea_ice_vap_state

select case (trim(key))

  case ("it_steps");             !ITERATION STEPS
    if (value == 0d0) then
      ctrl_loop_maximum = 100   !DEFAULT = 100
    elseif (value < 0d0) then
      ctrl_loop_maximum = -1
    else   
      ctrl_loop_maximum = nint(value)
    endif

  case ("init_vap_dens");        !STARTING VAPOUR DENSITY
     if (value == 0d0) then
       ctrl_mode_vapour = 0     !DEFAULT = IDEAL GAS
     elseif (value < -1d0) then
       !ignore it
     elseif (value < 0d0) then
       ctrl_mode_vapour = nint(value)
     else
       ctrl_mode_vapour = 1
       ctrl_init_d_vap = value
     endif

  case ("init_liq_dens");        !STARTING LIQUID DENSITY
    if (value == 0d0) then
      ctrl_mode_liquid = 0      !DEFAULT = TRIPLE POINT
    elseif (value < -1d0) then
      !ignore it
    elseif (value < 0d0) then
      ctrl_mode_liquid = nint(value)
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ("init_brine_sa");        !STARTING BRINE SALINITY
    if (value == 0d0) then
      ctrl_mode_salinity = 0    !DEFAULT = AUX FUNCTION
    elseif (value < -1d0) then
      !ignore it
    elseif (value < 0d0) then
      ctrl_mode_salinity = nint(value)
    else
      ctrl_mode_salinity = 1
      ctrl_init_s = value
    endif

  case ("init_temp");            !STARTING TEMPERATURE
    if (value == 0d0) then
      ctrl_mode_temperature = 0 !DEFAULT = AUX FUNCTION
    elseif (value < -1d0) then
      !ignore it
    elseif (value < 0d0) then
      ctrl_mode_temperature = nint(value)
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif

  case ("init_press");           !STARTING PRESSURE
    if (value == 0d0) then
      ctrl_mode_pressure = 0    !DEFAULT = AUX FUNCTION
    elseif (value < -1d0) then
      !ignore it
    elseif (value < 0d0) then
      ctrl_mode_pressure = nint(value)
    else
      ctrl_mode_pressure = 1
      ctrl_init_p = value
    endif

  case ("tol_temp");             !REQUIRED TEMPERATURE TOLERANCE
    if (value == 0d0) then
      ctrl_eps_exit_t = -0.0000001d0  !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_t = value
    endif

  case ("tol_press");            !REQUIRED PRESSURE TOLERANCE
    if (value == 0d0) then
      ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine clear_sea_ice_vap_state()
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_sea_ice_vap_done = 0

end subroutine

!==========================================================================
subroutine set_sea_ice_vap_state(sa_si, t_si, p_si, d_vap_si, d_liq_si)
!==========================================================================
                                  
!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

real*8 :: sa_si, t_si, p_si, d_vap_si, d_liq_si

equi_sea_ice_vap_done = isok

equi_sea_ice_vap_s = sa_si         !TRIPLE-POINT SALINITY
equi_sea_ice_vap_t = t_si          !TRIPLE-POINT TEMPERATURE
equi_sea_ice_vap_p = p_si          !TRIPLE-POINT PRESSURE
equi_sea_ice_vap_d_liq = d_liq_si   !DENSITY OF LIQUID PURE WATER
equi_sea_ice_vap_d_vap = d_vap_si   !DENSITY OF WATER VAPOUR

end subroutine

end module sea_ice_vap_4
