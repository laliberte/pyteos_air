module liq_vap_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIUM PROPERTIES OF FLUID WATER
!COMPUTED FROM IAPWS-95:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL .AND. SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER .AND. STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

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
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     MATHS_0,      FILE MATHS_0.F90
!     FLU_1,        FILE FLU_1.F90
!     FLU_2,        FILE FLU_2.F90
!     FLU_3A,       FILE FLU_3A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,    FILE CONVERT_0.F90

!#########################################################################

use constants_0
use maths_0
use flu_1
use flu_2
use flu_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: chk_iapws95_table8, &
          liq_vap_boilingtemperature_si, liq_vap_chempot_si, &
          liq_vap_density_liq_si, liq_vap_density_vap_si, &
          liq_vap_enthalpy_evap_si, liq_vap_enthalpy_liq_si, &
          liq_vap_enthalpy_vap_si, liq_vap_entropy_evap_si, &
          liq_vap_entropy_liq_si, liq_vap_entropy_vap_si, &
          liq_vap_pressure_liq_si, liq_vap_pressure_vap_si, &
          liq_vap_temperature_si, liq_vap_vapourpressure_si, &
          liq_vap_volume_evap_si, set_liq_vap_eq_at_p, &
          set_liq_vap_eq_at_t, set_it_ctrl_liq_vap 

integer, private :: equi_liq_vap_done
integer, private :: ctrl_mode_liquid, ctrl_mode_vapour, &
        ctrl_loop_maximum, ctrl_initialized, ctrl_mode_temperature
real*8, private :: ctrl_init_d_liq, ctrl_init_d_vap, ctrl_init_t, ctrl_eps_exit_p_vap
real*8, private :: equi_liq_vap_t, equi_liq_vap_p_liq, &
        equi_liq_vap_p_vap, equi_liq_vap_g_liq, equi_liq_vap_g_vap, &
        equi_liq_vap_d_liq, equi_liq_vap_d_vap
real*8, private :: d_liq_si, d_vap_si, g_liq_si, g_vap_si, p_liq_si, p_vap_si
real*8, private :: t_si

contains

!==========================================================================
function liq_vap_vapourpressure_si(t_si)
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE OF WATER VAPOUR IN PA AT THE  VAPOUR-LIQUID EQUILIBRIUM
!AS A FUNCTION OF ABSOLUTE TEMPERATURE T_SI IN K

!CHECK VALUE WITH DEFAULT SETTINGS: liq_vap_vapourpressure_si(300) = 3536.80675227

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_vapourpressure_si
real*8 t_si

liq_vap_vapourpressure_si = errorreturn

if(set_liq_vap_eq_at_t(t_si) == errorreturn) then
  return
endif

liq_vap_vapourpressure_si = equi_liq_vap_p_vap

end function

!==========================================================================
function liq_vap_boilingtemperature_si(p_si)
!==========================================================================

!RETURNS THE ABSOLUTE TEMPERATURE IN K AT THE VAPOUR-LIQUID EQUILIBRIUM
!AS A FUNCTION OF ABSOLUTE PRESSURE P_SI IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_vap_boilingtemperature_si(1E4) = 318.956328924

implicit none

real*8 liq_vap_boilingtemperature_si
real*8 p_si

liq_vap_boilingtemperature_si = errorreturn

if(set_liq_vap_eq_at_p(p_si) == errorreturn) return

liq_vap_boilingtemperature_si = equi_liq_vap_t

end function

!==========================================================================
function liq_vap_pressure_vap_si()
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE OF WATER VAPOUR IN PA AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_liq_vap_eq_at_t 300
!liq_vap_pressure_vap_si = 3536.80675227

!set_liq_vap_eq_at_p 1e4
!liq_vap_pressure_vap_si = 10000

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_pressure_vap_si

liq_vap_pressure_vap_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_pressure_vap_si = equi_liq_vap_p_vap

end function

!==========================================================================
function liq_vap_pressure_liq_si()
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE OF WATER VAPOUR IN PA AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_liq_vap_eq_at_t 300
!liq_vap_pressure_liq_si = 3536.80675|227

!set_liq_vap_eq_at_p 1e4
!liq_vap_pressure_liq_si = 10000.00000|00

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_pressure_liq_si

liq_vap_pressure_liq_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_pressure_liq_si = equi_liq_vap_p_liq

end function

!==========================================================================
function liq_vap_temperature_si()
!==========================================================================

!RETURNS THE ABSOLUTE TEMPERATUE IN K AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_TEMPERATURE_SI = 300

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_TEMPERATURE_SI = 318.956328924

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_temperature_si

liq_vap_temperature_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_temperature_si = equi_liq_vap_t

end function

!==========================================================================
function liq_vap_chempot_si()
!==========================================================================

!RETURNS THE CHEMICAL POTENTIAL OF WATER VAPOUR IN J/KG AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_CHEMPOT_SI = -5361.84908682

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_CHEMPOT_SI = -15259.1024273

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_chempot_si

liq_vap_chempot_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_chempot_si = equi_liq_vap_g_vap

end function

!==========================================================================
function liq_vap_density_liq_si()
!==========================================================================

!RETURNS THE DENSITY OF LIQUID WATER IN KG/M3 AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_DENSITY_LIQ_SI = 996.513027468

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_DENSITY_LIQ_SI = 989.833275365

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_density_liq_si

liq_vap_density_liq_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_density_liq_si = equi_liq_vap_d_liq

end function

!==========================================================================
function liq_vap_density_vap_si()
!==========================================================================

!RETURNS THE DENSITY OF WATER VAPOUR IN KG/M3 AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_DENSITY_VAP_SI = 2.55896736829e-02

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_DENSITY_VAP_SI = 6.81657223094e-02

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_density_vap_si

liq_vap_density_vap_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_density_vap_si = equi_liq_vap_d_vap

end function

!==========================================================================
function liq_vap_entropy_liq_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY OF LIQUID WATER IN J/(KG K) AT THE  VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTROPY_LIQ_SI = 393.089029801

!set_liq_vap_eq_at_p 1e4
!liq_vap_entropy_liq_si = 649.195605196

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_entropy_liq_si

liq_vap_entropy_liq_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_entropy_liq_si = flu_entropy_si(equi_liq_vap_t, equi_liq_vap_d_liq)

end function

!==========================================================================
function liq_vap_entropy_vap_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY OF WATER VAPOUR IN J/(KG K) AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTROPY_VAP_SI = 8517.38650061

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_ENTROPY_VAP_SI = 8148.82019425

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_entropy_vap_si

liq_vap_entropy_vap_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_entropy_vap_si = flu_entropy_si(equi_liq_vap_t, equi_liq_vap_d_vap)

end function

!==========================================================================
function liq_vap_enthalpy_liq_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY OF LIQUID WATER IN J/KG AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTHALPY_LIQ_SI = 112564.859854

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_ENTHALPY_LIQ_SI = 191805.944559

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_enthalpy_liq_si

liq_vap_enthalpy_liq_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_enthalpy_liq_si = flu_enthalpy_si(equi_liq_vap_t, equi_liq_vap_d_liq)

end function

!==========================================================================
function liq_vap_enthalpy_vap_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY OF WATER VAPOUR IN J/KG AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTHALPY_VAP_SI = 2549854.10109

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_ENTHALPY_VAP_SI = 2583858.67179

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_enthalpy_vap_si

liq_vap_enthalpy_vap_si = errorreturn

if(equi_liq_vap_done /= isok) return

liq_vap_enthalpy_vap_si = flu_enthalpy_si(equi_liq_vap_t, equi_liq_vap_d_vap)

end function

!==========================================================================
function liq_vap_volume_evap_si()
!==========================================================================

!RETURNS THE SPECIFIC EVAPORATION VOLUME IN M3/KG AT THE VAPOUR-LIQUID EQUILIBRIUM
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_VOLUME_EVAP_SI = 39.0772595686

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_VOLUME_EVAP_SI = 14.6691196141

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_volume_evap_si

liq_vap_volume_evap_si = errorreturn

if(equi_liq_vap_done /= isok) return

if(equi_liq_vap_d_liq <= 0) return
if(equi_liq_vap_d_vap <= 0) return

liq_vap_volume_evap_si = 1d0 / equi_liq_vap_d_vap - 1d0 / equi_liq_vap_d_liq

end function

!==========================================================================
function liq_vap_entropy_evap_si()
!==========================================================================

!RETURNS THE SPECIFIC EVAPORATION ENTROPY IN J/(KG K) OF WATER
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTROPY_EVAP_SI = 8124.29747080

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_ENTROPY_EVAP_SI = 7499.62458905

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_entropy_evap_si
real*8 sl, sv

liq_vap_entropy_evap_si = errorreturn

if(equi_liq_vap_done /= isok) return

sl = flu_entropy_si(equi_liq_vap_t, equi_liq_vap_d_liq)
if(sl == errorreturn) return
sv = flu_entropy_si(equi_liq_vap_t, equi_liq_vap_d_vap)
if(sv == errorreturn) return

liq_vap_entropy_evap_si = sv - sl

end function

!==========================================================================
function liq_vap_enthalpy_evap_si()
!==========================================================================

!RETURNS THE SPECIFIC EVAPORATION ENTHALPY IN J/KG OF WATER
!AFTER EXECUTION OF SET_LIQ_VAP_EQUILIBRIUM_TO_T .OR. SET_LIQ_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_LIQ_VAP_EQ_AT_T 300
!LIQ_VAP_ENTHALPY_EVAP_SI = 2437289.24124

!SET_LIQ_VAP_EQ_AT_P 1E4
!LIQ_VAP_ENTHALPY_EVAP_SI = 2392052.72723

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 liq_vap_enthalpy_evap_si
real*8 hl, hv

liq_vap_enthalpy_evap_si = errorreturn

if(equi_liq_vap_done /= isok) return

hl = flu_enthalpy_si(equi_liq_vap_t, equi_liq_vap_d_liq)
if(hl == errorreturn) return
hv = flu_enthalpy_si(equi_liq_vap_t, equi_liq_vap_d_vap)
if(hv == errorreturn) return

liq_vap_enthalpy_evap_si = hv - hl

end function

!==========================================================================
function set_liq_vap_eq_at_p(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE LIQUID-VAPOUR EQUILIBRIUM AT GIVEN PRESSURE IN PA
!BY 3-DIMENSIONAL NEWTON ITERATION

implicit none

integer maxit
real*8 set_liq_vap_eq_at_p
real*8 t_si, p_si
real*8 pl, pv
real*8 gl, gv
real*8 dl, dv
real*8 eps
real*8 dc, pc, tc, gc, dl_triple, dv_triple, t_triple

!CRITICAL POINT
dc = cp_density_si
pc = cp_pressure_si
tc = cp_temperature_si
gc = cp_chempot_si

!TRIPLE POINT
dl_triple = tp_density_liq_iapws95_si
dv_triple = tp_density_vap_iapws95_si
t_triple = tp_temperature_si

if(equi_liq_vap_done == isok .and. &
  p_si == equi_liq_vap_p_vap) then
  !the requested state has already been computed earlier
  set_liq_vap_eq_at_p = isok
  return
end if

call clear_liq_vap_state !waste any previous state

set_liq_vap_eq_at_p = errorreturn

if(p_si <= 0d0) return
if(p_si > pc) return

if(p_si == pc) then
  call set_liq_vap_state(tc, pc, pc, gc, gc, dc, dc)
  set_liq_vap_eq_at_p = isok
  return
end if

call init_it_ctrl_liq_vap

!SET INITIAL TEMPERATURE GUESS
select case (ctrl_mode_temperature)
  case (0); t_si = aux_boil_temperature(p_si)
  case (-1); t_si = t_triple
  case (1); t_si = ctrl_init_t
  case default; return
end select

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0); !default liquid density:
    if(t_si < 350d0) then
      dl = aux_liq_density_triplepoint(t_si)
    else
      dl = aux_liq_density_critical(t_si)
    end if
  case (-1); dl = dl_triple
  case (1); dl = ctrl_init_d_liq
  case default; return
end select

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0)
    if(t_si < 550d0) then
      dv = aux_vap_density_triplepoint(t_si)
    else
      dv = aux_vap_density_critical(t_si)
    end if
  case (-1); dv = dv_triple
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  set_liq_vap_eq_at_p = isok
  pl = flu_pressure_si(t_si, dl)
  pv = flu_pressure_si(t_si, dv)
  gl = flu_gibbs_energy_si(t_si, dl)
  gv = flu_gibbs_energy_si(t_si, dv)
  call set_liq_vap_state(t_si, pl, pv, gl, gv, dl, dv)
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else 
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p_vap
if(eps == 0) return

!RUN ITERATION LOOP
if(liq_vap_iteration_at_p(p_si, maxit, eps, t_si, dl, dv, gl, gv, pl, pv) == errorreturn) then
  return
end if

call set_liq_vap_state(t_si, pl, pv, gl, gv, dl, dv)

set_liq_vap_eq_at_p = isok

end function

!==========================================================================
function set_liq_vap_eq_at_t(t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE LIQUID-VAPOUR EQUILIBRIUM AT GIVEN TEMPERATURE IN K
!BY 2-DIMENSIONAL NEWTON ITERATION

implicit none

integer maxit
real*8 set_liq_vap_eq_at_t
real*8 t_si
real*8 pl, pv
real*8 gl, gv
real*8 dl, dv
real*8 eps
real*8 dc, pc, tc, gc, dl_triple, dv_triple

!CRITICAL POINT
dc = cp_density_si
pc = cp_pressure_si
tc = cp_temperature_si
gc = cp_chempot_si

!TRIPLE POINT
dl_triple = tp_density_liq_iapws95_si
dv_triple = tp_density_vap_iapws95_si

if(equi_liq_vap_done == isok .and.  &
   t_si == equi_liq_vap_t) then
  !the requested state has already been computed earlier
  set_liq_vap_eq_at_t = isok
  return
end if

call clear_liq_vap_state !waste any previous state

set_liq_vap_eq_at_t = errorreturn

if(t_si <= 0) return
if(t_si > tc) return

if(t_si == tc) then
  call set_liq_vap_state(t_si, pc, pc, gc, gc, dc, dc)
  set_liq_vap_eq_at_t = isok
  return
end if

call init_it_ctrl_liq_vap

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  !default liquid density:
    if(t_si < 350d0) then
      dl = aux_liq_density_triplepoint(t_si)
    else
      dl = aux_liq_density_critical(t_si)
    end if
  case (-1); dl = dl_triple
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0)
    if(t_si < 550d0) then
      dv = aux_vap_density_triplepoint(t_si)
    else
      dv = aux_vap_density_critical(t_si)
    end if
  case (-1); dv = dv_triple
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  set_liq_vap_eq_at_t = isok
  pl = flu_pressure_si(t_si, dl)
  pv = flu_pressure_si(t_si, dv)
  gl = flu_gibbs_energy_si(t_si, dl)
  gv = flu_gibbs_energy_si(t_si, dv)
  call set_liq_vap_state(t_si, pl, pv, gl, gv, dl, dv)
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else 
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(liq_vap_iteration_at_t(t_si, maxit, eps, dl, dv, gl, gv, pl, pv) == errorreturn) &
  return

call set_liq_vap_state(t_si, pl, pv, gl, gv, dl, dv)

set_liq_vap_eq_at_t = isok

end function

!==========================================================================
function liq_vap_iteration_at_p(p_si, maxit, eps, t_si, d_liq_si, &
         d_vap_si, g_liq_si, g_vap_si, p_liq_si, p_vap_si)
!==========================================================================

!THIS FUNCTION RETURNS THE LIQUID-VAPOUR PHASE EQUILIBRIUM FROM EQUAL PRESSURES,
!TEMPERATURES .AND. CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN PRESSURE, P_SI,
!FROM INITIAL GUESSES FOR THE LIQUID DENSITY, D_LIQ_SI, THE VAPOUR DENSITY, D_VAP_SI,
!AND THE TEMPERATURE, T_SI
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PRESSURE

!OUTPUT: 
!        LIQ_VAP_ITERATION_AT_T = ISOK IF SUCCESSFUL
!        LIQ_VAP_ITERATION_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITIES .OR. TEMPERATURE HAVE TAKEN A ZERO .OR. NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
! D_LIQ_SI: LIQUID DENSITY IN KG/M3 AT THE LIQUID-VAPOUR EQUILIBRIUM
! D_VAP_SI: VAPOUR DENSITY IN KG/M3 AT THE LIQUID-VAPOUR EQUILIBRIUM
! G_LIQ_SI: LIQUID CHEMICAL POTENTIAL IN J/KG AT THE LIQUID-VAPOUR EQUILIBRIUM
! G_VAP_SI: VAPOUR CHEMICAL POTENTIAL IN J/KG AT THE LIQUID-VAPOUR EQUILIBRIUM
! P_LIQ_SI: LIQUID PRESSURE IN PA AT THE LIQUID-VAPOUR EQUILIBRIUM
! P_VAP_SI: VAPOUR PRESSURE IN PA AT THE LIQUID-VAPOUR EQUILIBRIUM

!INPUT: 
!       P_SI: ABSOLUTE PRESSURE IN PA
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!   D_LIQ_SI: INITIAL GUESS OF LIQUID DENSITY IN KG/M3
!   D_VAP_SI: INITIAL GUESS OF VAPOUR DENSITY IN KG/M3

implicit none
real*8 liq_vap_iteration_at_p

integer it, maxit
real*8 eps, t_si, d_liq_si, d_vap_si,  &
       g_liq_si, g_vap_si, p_liq_si, p_vap_si
real*8 p_si
real*8 dl, dv
real*8 pl, pv
real*8 gl, gv
real*8 t
real*8 fl, fl_d, fl_dd
real*8 fv, fv_d, fv_dd
real*8 fl_t, fl_td
real*8 fv_t, fv_td
real*8 a(3, 3), b(3), x(3)
real*8 ddl, ddv, dt, pv_old
real*8 dc, pc, tc, gc

!THE CRITICAL PROPERTIES OF IAPWS-95
dc = cp_density_si
pc = cp_pressure_si
tc = cp_temperature_si
gc = cp_chempot_si

liq_vap_iteration_at_p = errorreturn

pv = errorreturn

if(p_si == pc) then !equlibrium = critical point
  t_si = tc
  d_liq_si = dc
  d_vap_si = dc
  g_liq_si = gc
  g_vap_si = gc
  p_liq_si = pc
  p_vap_si = pc
  liq_vap_iteration_at_p = isok
  return
end if

if((d_liq_si < dc) .or. &
  (d_vap_si <= 0) .or. (d_vap_si > dc) .or. &
  (p_si <=  0d0 .or. p_si > pc)) then
  t_si = errorreturn
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  g_liq_si = errorreturn
  g_vap_si = errorreturn
  p_liq_si = errorreturn
  p_vap_si = errorreturn
  return
endif

check_limits = check_limits - 1

dl = d_liq_si
dv = d_vap_si
t = t_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 3D-NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fv = flu_f_si(0, 0, t, dv)
  if(fv == errorreturn) go to 999

  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t, dv)
  if(fv_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t, dv)
  if(fv_dd == errorreturn) go to 999
  
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fv_t = flu_f_si(1, 0, t, dv)
  if(fv_t == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999
  fv_td = flu_f_si(1, 1, t, dv)
  if(fv_td == errorreturn) go to 999

  pv_old = pv  !KEEP PREVIOUS PV TO CHECK ACCURACY
  !LIQUID AND VAPOUR PRESSURES AND CHEMICAL POTENTIALS
  pl = dl ** 2 * fl_d    !PRESSURE OF LIQUID  !THIS OFTEN TAKES NEGATIVE VALUES, TEMPORARILY
  pv = dv ** 2 * fv_d    !PRESSURE OF VAPOUR
  gl = fl + dl * fl_d    !CHEM. POT. OF LIQUID
  gv = fv + dv * fv_d    !CHEM. POT. OF VAPOUR

  if(it > 0) then
    !CHECK ABSOLUTE .OR. RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or. &
       (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      liq_vap_iteration_at_p = isok
      p_liq_si = pl
      p_vap_si = pv
      g_liq_si = gl
      g_vap_si = gv
      d_liq_si = dl
      d_vap_si = dv
      t_si = t
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = dl * (2d0 * fl_d + dl * fl_dd)
  a(1, 2) = 0d0
  a(1, 3) = dl ** 2 * fl_td
  a(2, 1) = 0d0
  a(2, 2) = dv * (2d0 * fv_d + dv * fv_dd)
  a(2, 3) = dv ** 2 * fv_td
  a(3, 1) = (2d0 * fl_d + dl * fl_dd)
  a(3, 2) = -(2d0 * fv_d + dv * fv_dd)
  a(3, 3) = fl_t + dl * fl_td - fv_t - dv * fv_td

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = p_si - pl
  b(2) = p_si - pv
  b(3) = gv - gl

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 3) /= 0) go to 999 !MATRIX SINGULAR
  ddl = x(1)
  ddv = x(2)
  dt = x(3)
  
  !UPDATE DENSITIES & TEMPERATURE
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  dv = dv + ddv
  if(dv <= 0d0) go to 999
  t = t + dt
  if(t <= 0d0) go to 999
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
    d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    g_liq_si = errorreturn
    g_vap_si = errorreturn
    p_liq_si = errorreturn
    p_vap_si = errorreturn
    liq_vap_iteration_at_p = errorreturn
    return
  endif
!dbg  !AIR_LIMITS
!dbg  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) then
!dbg    t_si = errorreturn
!dbg    d_liq_si = errorreturn
!dbg    d_vap_si = errorreturn
!dbg    g_liq_si = errorreturn
!dbg    g_vap_si = errorreturn
!dbg    p_liq_si = errorreturn
!dbg    p_vap_si = errorreturn
!dbg    liq_vap_iteration_at_p = errorreturn
!dbg    return
!dbg  endif
endif

end function

!==========================================================================
function liq_vap_iteration_at_t(t_si, maxit, eps, d_liq_si, &
                 d_vap_si, g_liq_si, g_vap_si, p_liq_si, p_vap_si)
!==========================================================================

!THIS FUNCTION RETURNS THE LIQUID-VAPOUR PHASE EQUILIBRIUM FROM EQUAL PRESSURES
!AND CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN TEMPERATURE, T_SI,
!FROM INITIAL GUESSES FOR THE LIQUID DENSITY, D_LIQ_SI, .AND. THE VAPOUR DENSITY, D_VAP_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PRESSURE

!OUTPUT: 
!        LIQ_VAP_ITERATION_AT_T = ISOK IF SUCCESSFUL
!        LIQ_VAP_ITERATION_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITIES HAVE TAKEN A ZERO .OR. NEGATIVE VALUE DURING THE ITERATION
! D_LIQ_SI: LIQUID DENSITY IN KG/M3 AT THE LIQUID-VAPOUR EQUILIBRIUM
! D_VAP_SI: VAPOUR DENSITY IN KG/M3 AT THE LIQUID-VAPOUR EQUILIBRIUM
! G_LIQ_SI: LIQUID CHEMICAL POTENTIAL IN J/KG AT THE LIQUID-VAPOUR EQUILIBRIUM
! G_VAP_SI: VAPOUR CHEMICAL POTENTIAL IN J/KG AT THE LIQUID-VAPOUR EQUILIBRIUM
! P_LIQ_SI: LIQUID PRESSURE IN PA AT THE LIQUID-VAPOUR EQUILIBRIUM
! P_VAP_SI: VAPOUR PRESSURE IN PA AT THE LIQUID-VAPOUR EQUILIBRIUM

!INPUT: 
!       T_SI: ABSOLUTE TEMPERATURE IN K
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!   D_LIQ_SI: INITIAL GUESS OF LIQUID DENSITY IN KG/M3
!   D_VAP_SI: INITIAL GUESS OF VAPOUR DENSITY IN KG/M3

implicit none
real*8 liq_vap_iteration_at_t

integer it, maxit
real*8 t_si, eps, d_liq_si, d_vap_si, g_vap_si, g_liq_si, &
       p_liq_si, p_vap_si
real*8 dl, dv
real*8 pl, pv
real*8 gl, gv
real*8 fv, fv_d, fv_dd
real*8 fl, fl_d, fl_dd
real*8 a(2, 2), b(2), x(2)
real*8 ddl, ddv, pv_old
real*8 dc, pc, tc, gc

!THE CRITICAL PROPERTIES OF IAPWS-95
dc = cp_density_si
pc = cp_pressure_si
tc = cp_temperature_si
gc = cp_chempot_si

liq_vap_iteration_at_t = errorreturn

pv = errorreturn

if(t_si == tc) then    !EQULIBRIUM = CRITICAL POINT
  d_liq_si = dc
  d_vap_si = dc
  g_liq_si = gc
  g_vap_si = gc
  p_liq_si = pc
  p_vap_si = pc
  liq_vap_iteration_at_t = isok
  return
end if

if(t_si > tc .or. &
  d_liq_si < dc .or. &
  d_vap_si <= 0d0 .or. d_vap_si > dc) then  !no equilibrium in the supercritical range
  d_liq_si = errorreturn
  d_vap_si = errorreturn
  g_liq_si = errorreturn
  g_vap_si = errorreturn
  p_liq_si = errorreturn
  p_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    g_liq_si = errorreturn
    g_vap_si = errorreturn
    p_liq_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
else
  if(t_si <= 0d0)then
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    g_liq_si = errorreturn
    g_vap_si = errorreturn
    p_liq_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
endif

If(d_liq_si < dc) return
If(d_vap_si <= 0) return
If(d_vap_si > dc) return

check_limits = check_limits - 1

dl = d_liq_si
dv = d_vap_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 2D-NEWTON ITERATION
  fl = flu_f_si(0, 0, t_si, dl)
  if(fl == errorreturn) go to 999
  fv = flu_f_si(0, 0, t_si, dv)
  if(fv == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t_si, dl)
  if(fl_d == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t_si, dv)
  if(fv_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t_si, dl)
  if(fl_dd == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t_si, dv)
  if(fv_dd == errorreturn) go to 999
  
  pv_old = pv            !KEEP PREVIOUS PV TO CHECK ACCURACY
  !LIQUID AND VAPOUR PRESSURES AND CHEMICAL POTENTIALS
  pl = dl ** 2 * fl_d    !PRESSURE OF LIQUID  !THIS OFTEN TAKES NEGATIVE VALUES TEMPORARILY
  pv = dv ** 2 * fv_d    !PRESSURE OF VAPOUR
  gl = fl + dl * fl_d    !CHEM. POT. OF LIQUID
  gv = fv + dv * fv_d    !CHEM. POT. OF VAPOUR

  if(it > 0) then
    !CHECK ABSOLUTE .OR. RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or. &
       (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      liq_vap_iteration_at_t = isok
      p_liq_si = pl
      p_vap_si = pv
      g_liq_si = gl
      g_vap_si = gv
      d_liq_si = dl
      d_vap_si = dv
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(2, 1) = 2d0 * fl_d + dl * fl_dd
  a(1, 1) = dl * a(2, 1)
  a(2, 2) = -2d0 * fv_d - dv * fv_dd
  a(1, 2) = dv * a(2, 2)
  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = pv - pl
  b(2) = gv - gl

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR
  ddl = x(1)
  ddv = x(2)

  !UPDATE DENSITIES
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  dv = dv + ddv
  if(dv <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    d_liq_si = errorreturn
    d_vap_si = errorreturn
    g_liq_si = errorreturn
    g_vap_si = errorreturn
    p_liq_si = errorreturn
    p_vap_si = errorreturn
    liq_vap_iteration_at_t = errorreturn
    return
  endif
!dbg  !AIR_LIMITS
!dbg  if(d_vap_si <= dry_air_dmin .or. d_vap_si > dry_air_dmax) then
!dbg    d_liq_si = errorreturn
!dbg    d_vap_si = errorreturn
!dbg    g_liq_si = errorreturn
!dbg    g_vap_si = errorreturn
!dbg    p_liq_si = errorreturn
!dbg    p_vap_si = errorreturn
!dbg    liq_vap_iteration_at_t = errorreturn
!dbg    return
!dbg  endif
endif

end function

!==========================================================================
function aux_liq_density_critical(t_si)
!==========================================================================

!COMPUTES THE LIQUID DENSITY ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_T (300 K < T < TC WITH RMS = 1 K)

! T/TC - 1 = A1 * (D/DC-1)^3 + A2 * (D/DC-1)^6 + A3 * (D/DC-1)^9

implicit none

real*8 aux_liq_density_critical, t_si
real*8 tc, dc
real*8 a1, a2, a3
real*8 r, s, t
real*8 x1, x2, x3
real*8 d

tc = cp_temperature_si
dc = cp_density_si

a1 = -7.34017329598858d-02
a2 = 5.70516487711065d-03
a3 = -4.31313846955949d-04

aux_liq_density_critical = errorreturn

!SOLVE       X^3 + R * X^2 + S * X + T = 0
r = a2 / a3
s = a1 / a3
t = (1d0 - t_si / tc) / a3

if(get_cubicroots(r, s, t, x1, x2, x3) /= 1) return

if(x1 /= 0d0) then
  x1 = sign(1d0, x1) * abs(x1) ** (1d0 / 3d0)
end if
d = (x1 + 1d0) * dc
if(d < dc) return

aux_liq_density_critical = d

end function

!==========================================================================
function aux_vap_density_critical(t_si)
!==========================================================================

!COMPUTES THE VAPOUR DENSITY ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_T (550 K < T < TC WITH RMS = 0.4 K)

! T/TC - 1 = A1 * (1-D/DC)^4 + A2 * (1-D/DC)^8 + A3 * (1-D/DC)^12

implicit none

real*8 aux_vap_density_critical, t_si
real*8 tc, dc
real*8 a1, a2, a3
real*8 r, s, t
real*8 x1, x2, x3
real*8 d

tc = cp_temperature_si
dc = cp_density_si
a1 = -0.237216002118091d0
a2 = 0.186593118426901d0
a3 = -0.258472040504799d0

aux_vap_density_critical = errorreturn

!SOLVE       X^3 + R * X^2 + S * X + T = 0
r = a2 / a3
s = a1 / a3
t = (1d0 - t_si / tc) / a3

if(get_cubicroots(r, s, t, x1, x2, x3) /= 1) return

if(x1 < 0d0) return
if(x1 > 0d0) then
  x1 = abs(x1) ** 0.25d0
end if
d = (1d0 - x1) * dc
if(d <= 0d0) return

aux_vap_density_critical = d

end function

!==========================================================================
function aux_liq_density_triplepoint(t_si)
!==========================================================================

!COMPUTES THE LIQUID DENSITY ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_T (TT < T < 350 K WITH RMS = 0.002 KG/M3)

! (D/DT-1) = A1 * (T/TT-1) + A2 * (T/TT-1)^2 + A3 * (T/TT-1)^3 + A4 * (T/TT-1)^4 + A5 * (T/TT-1)^5

implicit none

real*8 aux_liq_density_triplepoint, t_si
real*8 tt, dt
real*8 a1, a2, a3, a4, a5
real*8 tr
real*8 d

tt = tp_temperature_si
dt = tp_density_liq_iapws95_si

a1 = 1.80066818428501d-02
a2 = -0.648994409718973d0
a3 = 1.56594764908347d0
a4 = -3.18116999660964d0
a5 = 2.98590977093295d0

tr = t_si / tt - 1d0
d = dt*(a1*tr + a2*tr**2 + a3*tr**3 + a4*tr**4 + a5*tr**5 + 1)

aux_liq_density_triplepoint = d

end function

!==========================================================================
function aux_vap_density_triplepoint(t_si)
!==========================================================================

!COMPUTES THE VAPOUR DENSITY ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_T (TT < T < 550 K WITH RMS = 0.01 IN LN(D/DT))

! LN(D/DT) = A1 * (TT/T-1) + A2 * (TT/T-1)^2 + A3 * (TT/T-1)^3

implicit none
real*8 aux_vap_density_triplepoint, t_si
real*8 tt, dt
real*8 a1, a2, a3
real*8 tr
real*8 d

tt = tp_temperature_si
dt = tp_density_vap_iapws95_si

a1 = -19.2235086866063d0
a2 = -6.15770193302955d0
a3 = -4.965736126494d0

tr = tt / t_si - 1d0
d = dt * exp(a1 * tr + a2 * tr ** 2 + a3 * tr ** 3)

aux_vap_density_triplepoint = d

end function

!==========================================================================
function aux_boil_temperature(p_si)
!==========================================================================

!COMPUTES THE BOILING TEMPERATURE ON THE SATURATION CURVE FROM A CORRELATION POLYNOMIAL
!AS A FIRST GUESS FOR LIQ_VAP_ITERATION_AT_P (TT < T < 640 K WITH RMS = 0.01 IN LN(P/PT))

! CLAUSIUS-CLAPEYRON TYPE EQUATION:
! LN(P/PT) = A1 * (TT/T-1) + A2 * (TT/T-1)^2

implicit none
real*8 aux_boil_temperature, p_si
real*8 tt, pt
real*8 a1, a2
real*8 p, q
real*8 tr

tt = tp_temperature_si
pt = tp_pressure_iapws95_si

a1 = -19.8731005709116d0
a2 = -3.08975437352998d0

p = a1 / a2
q = -log(p_si / pt) / a2

tr = -0.5d0*p + sqrt(0.25d0*p**2 - q)

aux_boil_temperature = tt / (tr + 1d0)

end function

!==========================================================================
subroutine set_liq_vap_state(t, pl, pv, gl, gv, dl, dv)
!==========================================================================

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none
real*8 t, pl, pv, gl, gv, dl, dv

equi_liq_vap_done = isok

equi_liq_vap_t = t
equi_liq_vap_p_liq = pl
equi_liq_vap_p_vap = pv
equi_liq_vap_g_liq = gl
equi_liq_vap_g_vap = gv
equi_liq_vap_d_liq = dl
equi_liq_vap_d_vap = dv

return
end subroutine

!==========================================================================
subroutine clear_liq_vap_state()
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

equi_liq_vap_done = 0

return
end subroutine

!==========================================================================
subroutine init_it_ctrl_liq_vap()
!==========================================================================

implicit none

real*8 dl_triple, dv_triple, t_triple

!TRIPLE POINT
dl_triple = tp_density_liq_iapws95_si
dv_triple = tp_density_vap_iapws95_si
t_triple = tp_temperature_si

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES .AND. MODES FOR DENSITY ITERATION
ctrl_mode_liquid = 0
ctrl_mode_vapour = 0
ctrl_mode_temperature = 0
ctrl_loop_maximum = 100
ctrl_init_d_liq = dl_triple
ctrl_init_d_vap = dv_triple
ctrl_init_t = t_triple
ctrl_eps_exit_p_vap = -0.0000001d0 !RELATIVE, 0.1 PPM

return
end subroutine

!==========================================================================
subroutine set_it_ctrl_liq_vap(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!IAPWS-95 VAPOUR PRESSURE FROM TEMPERATURE

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VAPOUR DENSITY TO COMPUTE VAPOUR PRESSURE

!'INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY TO START:
!                                   T < 350 K: LIQUID DENSITY CORRELATION FUNCTION, D = AUX_LIQ_DENSITY_TRIPLEPOINT(T)
!                                   T > 350 K: LIQUID DENSITY CORRELATION FUNCTION, D = AUX_LIQ_DENSITY_CRITICAL(T)
!'INIT_LIQ_DENS'       -1           USE TRIPLE POINT LIQUID DENSITY TO START
!'INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!'INIT_VAP_DENS'        0           USE DEFAULT VAPOUR DENSITY TO START
!                                   T < 550 K: VAPOUR DENSITY CORRELATION FUNCTION, D = AUX_VAP_DENSITY_TRIPLEPOINT(T)
!                                   T > 550 K: VAPOUR DENSITY CORRELATION FUNCTION, D = AUX_VAP_DENSITY_CRITICAL(T)
!'INIT_VAP_DENS'       -1           USE TRIPLE POINT VAPOUR DENSITY TO START
!'INIT_VAP_DENS'        D > 0       USE VALUE D AS VAPOUR DENSITY TO START

!'INIT_TEMP    '        0           USE DEFAULT TEMPERATURE TO START ( = AUX_BOIL_TEMPERATURE(P))
!'INIT_TEMP    '       -1           USE TRIPLE POINT TEMPERATURE TO START
!'INIT_TEMP    '        T > 0       USE VALUE T AS TEMPERATURE TO START

!'TOL_VAP_PRESS'        0           USE DEFAULT EXIT ACCURACY FOR VAPOUR PRESSURE (0.1 PPM)
!'TOL_VAP_PRESS'        EPS         USE EPS AS EXIT ACCURACY FOR VAPOUR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

implicit none

character*13 key
real*8 value

call init_it_ctrl_liq_vap

call clear_liq_vap_state

select case (trim(key))

  case ('it_steps')             !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100   !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')        !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_liquid = 0      !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then !IGNORE IT
    elseif(value < 0d0)then
      ctrl_mode_liquid = nint(value)
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ('init_vap_dens')        !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_vapour = 0      !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then !IGNORE IT
      elseif(value < 0d0)then
        ctrl_mode_vapour = nint(value)
      else
        ctrl_mode_vapour = 1
        ctrl_init_d_vap = value
      endif

  case ('init_temp')            !START TEMPERATURE
    if(value == 0d0)then
      ctrl_mode_temperature = 0 !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then !IGNORE IT
    elseif(value < 0d0)then
      ctrl_mode_temperature = nint(value)
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif

  case ('tol_vap_press')        !REQUIRED VAPOUR PRESSURE TOLERANCE
    if(value == 0d0)then
      ctrl_eps_exit_p_vap = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p_vap = value
    endif

  case default; write(*,*)key,' not supported'

end select

return
end subroutine

!==========================================================================
subroutine chk_iapws95_table8
!==========================================================================

implicit none
integer it
real*8 pl, pv, rhol, rhov
real*8 t, zt(3)
real*8 g, gt, gp
real*8 temp

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of IAPWS-95 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "TABLE A3 OF OS2008 AND TABLE 8 OF IAPWS-95 WITH P CALCULATED INTERNALLY"
write(*,*) "PROPERTIES OF H2O IN THE TWO-PHASE REGION AT SELECTED VALUES OF TEMPERATURE" 

zt(1) = 275d0
zt(2) = 450d0
zt(3) = 625d0

do it = 1,3

  t = zt(it)

#ifdef interactive
  write(*,*)' '
  write(*,*)'ENTER 1 TO SHOW RESULTS FOR T = ',t,'; 0 TO EXIT'
  read(*,*)icheck
  if(icheck /= 1) go to 999
#endif

  !if(it == 1)rhol = 9.99887406119850d2; rhov = 5.50664918504123d-3  
  !if(it == 2)rhol = 8.90341249761673d2; rhov = 4.81200360125671d0
  !if(it == 3)rhol = 5.67090385146353d2; rhov = 1.18290280451157d2
  temp = liq_vap_vapourpressure_si(t)
  rhol = liq_density_si(t,temp)
  rhov = vap_density_si(t,temp)

  pl=flu_pressure_si(t, rhol)
  g =  flu_f_si(0, 0, t, rhol) + pl/rhol
  gt =  flu_f_si(1, 0, t, rhol)
  gp =  1d0/rhol

  write(*,*) " "
  write(*,*) " "
  write(*,*) " "
  write(*,*) "PROPERTIES AT T =", zt(it), " K, P =", pl," Pa"
  write(*,*) " "

  write(*,*) "QUANTITY       TABLE VALUE             COMPUTED VALUE         &
  &UNITS"
  write(*,*) " "

  if(it == 1)then
    write(*,'(" pl             6.9845117e8       ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           999.887406120     ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             7.759722016e3     ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         2.83094669595e1   ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 2)then
    write(*,'(" pl             9.32203563628e5   ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           890.34124976167   ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             7.4916158501217e5 ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         2.10865844688447e3",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 3)then
    write(*,'(" pl             1.6908269318578e7 ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           5.670903851464e2  ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             1.6862697594697e6 ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         3.8019468301114e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  endif

  write(*,*)' '

  pv=flu_pressure_si(t, rhov)
  g = flu_f_si(0, 0, t, rhov) + pv/rhov
  gt =  flu_f_si(1, 0, t, rhov)

  if(it == 1)then
    write(*,'(" pv             6.98451166701e8   ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           5.506649185041e-3 ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.5042899500405e6 ",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         9.1066012052322e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 2)then
    write(*,'(" pv             9.322035636282e5  ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           4.812003601257    ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.77441077988962e6",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         6.6092122132788e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 3)then
    write(*,'(" pv             1.6908269318578e7 ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           1.182902804512e2  ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.5507162456235e6 ",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         5.1850612079574e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  endif

go to 999
999 write(*,*)' '
 
enddo

end subroutine

end module liq_vap_4
