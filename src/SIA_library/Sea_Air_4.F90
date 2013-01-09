module sea_air_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILIBRIUM PROPERTIES OF SEAWATER WITH HUMID AIR

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
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     AIR_1,         FILE AIR_1.F90
!     FLU_1,         FILE FLU_1.F90
!     SAL_1,         FILE SAL_1.F90
!     AIR_2,         FILE AIR_2.F90
!     FLU_2,         FILE FLU_2.F90
!     SAL_2,         FILE SAL_2.F90
!     AIR_3A,        FILE AIR_3A.F90
!     AIR_3B,        FILE AIR_3B.F90
!     FLU_3A,        FILE FLU_3A.F90
!     LIQ_VAP_4,     FILE LIQ_VAP_4.F90
!     LIQ_AIR_4A,    FILE LIQ_AIR_4A.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use air_1
use flu_1
use sal_1
use air_2
use flu_2
use sal_2
use air_3a
use air_3b
use flu_3a
use sea_3a
use liq_vap_4
use liq_air_4a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_air_chempot_evap_si, sea_air_condense_temp_si, &
          sea_air_density_air_si, sea_air_density_vap_si, &
          sea_air_enthalpy_evap_si, sea_air_entropy_air_si, &
          sea_air_massfraction_air_si, sea_air_vapourpressure_si, &
          set_it_ctrl_sea_air, set_sea_air_eq_at_s_a_p, &
          set_sea_air_eq_at_s_t_p

!CONTROL PARAMETERS OF THE VAPOUR PRESSURE ITERATION
integer ctrl_initialized

!SELECT INITIAL GUESSING METHOD
integer ctrl_mode_humidity
integer ctrl_mode_temperature
integer ctrl_mode_pressure

!STORE CUSTOMIZED INITIAL VALUES HERE
integer ctrl_loop_maximum
real*8 ctrl_init_t
real*8 ctrl_init_a
real*8 ctrl_init_p

!SET TOLERANCE OF THE ITERATION LOOP FOR
real*8 ctrl_eps_exit_p_vap   !vapour partial pressure
real*8 ctrl_eps_exit_temp    !temperature
real*8 ctrl_eps_exit_press    !pressure

!PROPERTIES OF THE CURRENT LIQUID-VAPOUR EQUILIBRIUM STATE
integer equi_sea_air_done

!STORAGE OF THE ITERATION RESULTS = THE EQUILIBRIUM STATE (SATURATED AIR)
real*8 equi_sea_air_s        !salinity
real*8 equi_sea_air_a        !dry air mass fraction of humid air
real*8 equi_sea_air_t        !temperature
real*8 equi_sea_air_p        !pressure
real*8 equi_sea_air_d_liq    !liquid water density
real*8 equi_sea_air_d_air    !humid air density

contains

!==========================================================================
function sea_air_massfraction_air_si(sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE AIR MASS FRACTION IN KG/KG OF SEA AIR

!OUTPUT: 
!A_SI  MASS FRACTION OF DRY AIR IN SEA AIR, IN KG/KG

!INPUTS: 
!SA_SI  ABSOLUTE SALINITY IN KG/KG
! T_SI  ABSOLUTE TEMPERATURE IN K
! P_SI  ABSOLUTE PRESSURE OF HUMID AIR IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!SEA_AIR_MASSFRACTION_AIR_SI(0.035, 300, 1E5) = 0.978024405439

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!SEA_AIR_MASSFRACTION_AIR_SI(0.035, 300, 1E5) = 0.978029483888

real*8 sea_air_massfraction_air_si, sa_si, t_si, p_si

sea_air_massfraction_air_si = errorreturn

if(set_sea_air_eq_at_s_t_p(sa_si, t_si, p_si) == errorreturn) return

sea_air_massfraction_air_si = equi_sea_air_a

end function

!==========================================================================
function sea_air_vapourpressure_si(sa_si, t_si, p_si)
!==========================================================================
 
!RETURNS THE PARTIAL PRESSURE OF VAPOUR IN HUMID AIR AT EQUILIBRIUM WITH SEAWATER,
!FROM ABSOLUTE SALINITY, TEMPERATURE AND TOTAL PRESSURE
 
!INPUT:  SA_SI ABSOLUTE SALINITY IN KG/KG
!        T_SI  ABSOLUTE TEMPERATURE IN K
!        P_SI  ABSOLUTE (TOTAL) PRESSURE IN PA
!OUTPUT: VAPOUR PRESSURE IN PA
 
!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!SEA_AIR_VAPOURPRESSURE_SI(0.035,300,1E5) = 3485.92802425
 
!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!SEA_AIR_VAPOURPRESSURE_SI(0.035,300,1E5) = 3485.92986681
 
real*8 sa_si, t_si, p_si, sea_air_vapourpressure_si
real*8 a_sat, xv
 
sea_air_vapourpressure_si = errorreturn
 
if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return
 
!GET THE AIR FRACTION AT EQUILIBRIUM WITH SEAWATER
a_sat = sea_air_massfraction_air_si(sa_si, t_si, p_si)
if(a_sat == errorreturn) return
if(a_sat <= 0d0 .or. a_sat > 1d0) return
 
!CONVERT AIR MASS FRACTION TO VAPOUR MOLE FRACTION
xv = air_molfraction_vap_si(a_sat)
if(xv == errorreturn) return
if(xv <= 0d0 .or. xv > 1d0) return
 
!COMPUTE THE PARTIAL PRESSURE
sea_air_vapourpressure_si = xv * p_si
 
end function

!==========================================================================
function sea_air_chempot_evap_si(a_si, sa_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE CHEMICAL POTENTIAL DIFFERENCE X_SA = LN(PHI(A_SI)/PHI(SA_SI))
!COMPUTED FROM THE RELATIVE FUGACITY OF HUMID AIR, PHI(A), AND OF SEA AIR, PHI(SA)
!AS A MEASURE FOR THE DEVIATION FROM AIR-SEA EQUILIBRIUM (ONSAGER FORCE)

!INPUT:  A_SI:  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
!        SA_SI: ABSOLUTE SALINITY OF SEAWATER
!        T_SI:  ABSOLUTE TEMPERATURE IN K
!        P_SI:  ABSOLUTE PRESSURE IN PA
!OUTPUT: CHEMICAL POTENTIAL DIFFERENCE /(R_W T), UNITLESS

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!SEA_AIR_CHEMPOT_EVAP_SI(0.9,0.035,300,1E5) = 1.45564492467

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!SEA_AIR_CHEMPOT_EVAP_SI(0.9,0.035,300,1E5) = 1.45584069071

real*8 sea_air_chempot_evap_si, a_si, sa_si, t_si, p_si
real*8 mu_v, mu_w

sea_air_chempot_evap_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0)return
if(sa_si < 0d0 .or. sa_si >= 1d0)return
if(t_si <= 0d0)return
if(p_si <= 0d0)return

!CHEMICAL POTENTIAL OF VAPOUR IN HUMID AIR
mu_v = air_g_chempot_vap_si(a_si, t_si, p_si)
if(mu_v == errorreturn)return

!CHEMICAL POTENTIAL OF WATER IN SEAWATER
mu_w = sea_chempot_h2o_si(sa_si, t_si, p_si)
if(mu_w == errorreturn)return

sea_air_chempot_evap_si = (mu_v - mu_w) / (gas_constant_h2O_si * t_si)

end function

!==========================================================================
function sea_air_condense_temp_si(sa_si, a_si, p_si)
!==========================================================================

!RETURNS THE SEAWATER CONDENSATION TEMPERATURE OF UNSATURATED HUMID AIR

!OUTPUT: 
!T_SI  ABSOLUTE CONDENSATION TEMPERATURE IN K

!INPUT: 
!SA_SI  ABSOLUTE SALINITY IN KG/KG
! A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR, IN KG/KG
! P_SI  ABSOLUTE PRESSURE OF HUMID AIR IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!SEA_AIR_CONDENSE_TEMP_SI(0.035, 0.99, 1E5) = 287.363869941

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!SEA_AIR_CONDENSE_TEMP_SI(0.035, 0.99, 1E5) = 287.367456468

real*8 sea_air_condense_temp_si, sa_si, a_si, p_si

sea_air_condense_temp_si = errorreturn

if(sa_si == 0d0) then
  sea_air_condense_temp_si = liq_air_dewpoint_si(a_si, p_si)
  return
end if

if(set_sea_air_eq_at_s_a_p(sa_si, a_si, p_si) == errorreturn) return

sea_air_condense_temp_si = equi_sea_air_t

end function

!==========================================================================
function sea_air_density_air_si()
!==========================================================================

!RETURNS THE HUMID AIR DENSITY IN KG/M3 OF SEA AIR
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_SEA_AIR_EQUILIBRIUM_TO_S_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_DENSITY_AIR_SI = 1.20527399511

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_DENSITY_AIR_SI = 1.14616388271

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_DENSITY_AIR_SI = 1.20553988598

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_DENSITY_AIR_SI = 1.14642944448

real*8 sea_air_density_air_si

sea_air_density_air_si = errorreturn

if(equi_sea_air_done /= isok) return

sea_air_density_air_si = equi_sea_air_d_air

end function

!==========================================================================
function sea_air_density_vap_si()
!==========================================================================

!RETURNS THE VAPOUR DENSITY IN KG/M3 OF SEA AIR
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_SEA_AIR_EQUILIBRIUM_TO_S_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_DENSITY_VAP_SI = 1.20527399511E-02

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_DENSITY_VAP_SI = 0.0251876327865

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_DENSITY_VAP_SI = 1.20553988598E-02

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_DENSITY_VAP_SI = 2.51876465812E-02

real*8 sea_air_density_vap_si

sea_air_density_vap_si = errorreturn

if(equi_sea_air_done /= isok) return

sea_air_density_vap_si = equi_sea_air_d_air * (1d0 - equi_sea_air_a)

end function

!==========================================================================
function sea_air_enthalpy_evap_si()
!==========================================================================

!RETURNS THE LATENT HEAT IN J/KG OF SEA AIR
!AFTER SETTING THE EQUILIBRIUM STATE BY CALLING SET_SEA_AIR_EQUILIBRIUM_TO_S_T_P ETC

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_ENTHALPY_EVAP_SI = 2464747.68964

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_ENTHALPY_EVAP_SI = 2434550.22477

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!SET_SEA_AIR_EQ_AT_S_A_P(0.035, 0.99, 1E5)
!SEA_AIR_ENTHALPY_EVAP_SI = 2464738.79136

!SET_SEA_AIR_EQ_AT_S_T_P(0.035, 300, 1E5)
!SEA_AIR_ENTHALPY_EVAP_SI = 2434549.89770

real*8 sea_air_enthalpy_evap_si
real*8 a, t, dh, dl
real*8 hh, hh_a, hl
real*8 fh_a, fh_d
real*8 fh_at, fh_ad
real*8 fh_td, fh_dd
real*8 p_a, p_d
real*8 h_a, h_d
real*8 s, p, hs, hs_s
real*8 gs, gs_s, gs_t, gs_st

sea_air_enthalpy_evap_si = errorreturn

if(equi_sea_air_done /= isok) return

!GET THE EQUILIBRIUM PROPERTIES
s = equi_sea_air_s
if(s == errorreturn) return
if(s < 0d0 .or. s >= 1d0) return
a = equi_sea_air_a
if(a == errorreturn) return
if(a < 0d0 .or. a > 1d0) return
t = equi_sea_air_t
if(t == errorreturn) return
if(t <= 0d0) return
p = equi_sea_air_p
if(p == errorreturn) return
if(p <= 0d0) return
dh = equi_sea_air_d_air
if(dh == errorreturn) return
if(dh <= 0d0) return
dl = equi_sea_air_d_liq
if(dl == errorreturn) return
if(dl <= 0d0) return

!ENTHALPY OF HUMID AIR
hh = air_f_enthalpy_si(a, t, dh)
if(hh == errorreturn) return

!DERIVATIVE HH_A = (DH/DA)_T_P OF ENTHALPY OF HUMID AIR
fh_a = air_f_si(1, 0, 0, a, t, dh)
if(fh_a == errorreturn) return
fh_d = air_f_si(0, 0, 1, a, t, dh)
if(fh_d == errorreturn) return
fh_at = air_f_si(1, 1, 0, a, t, dh)
if(fh_at == errorreturn) return
fh_ad = air_f_si(1, 0, 1, a, t, dh)
if(fh_ad == errorreturn) return
fh_td = air_f_si(0, 1, 1, a, t, dh)
if(fh_td == errorreturn) return
fh_dd = air_f_si(0, 0, 2, a, t, dh)
if(fh_dd == errorreturn) return

p_a = dh * fh_ad                            !DP/DA AT CONST T,D
p_d = dh * (2d0 * fh_d + dh * fh_dd)        !DP/DD AT CONST A,T
if(p_d == 0d0) return
h_a = fh_a + dh * fh_ad - t * fh_at         !DH/DA AT CONST T,D
h_d = 2d0 * fh_d + dh * fh_dd - t * fh_td   !DH/DD AT CONST A,T

hh_a = h_a - h_d * p_a / p_d                !DH/DA AT CONST T,P

!ENTHALPY OF LIQUID WATER
hl = flu_enthalpy_si(t, dl)
if(hl == errorreturn) return

if(s == 0d0) then !treat the pure-water case separately
  !compute the latent heat of water
  sea_air_enthalpy_evap_si = hh - a * hh_a - hl
  return
end if

!SALINE ENTHALPY
gs = sal_g_si(0, 0, 0, s, t, p)
if(gs == errorreturn) return
gs_t = sal_g_si(0, 1, 0, s, t, p)
if(gs_t == errorreturn) return
hs = gs - t * gs_t
gs_s = sal_g_si(1, 0, 0, s, t, p)
if(gs_s == errorreturn) return
gs_st = sal_g_si(1, 1, 0, s, t, p)
if(gs_st == errorreturn) return
hs_s = gs_s - t * gs_st

!COMPUTE THE LATENT HEAT OF SEAWATER
sea_air_enthalpy_evap_si = hh - a * hh_a - hl - hs + s * hs_s

end function

!==========================================================================
function sea_air_entropy_air_si(sa_si, t_si, p_si)
!==========================================================================

!RETURNS THE HUMID-AIR ENTROPY IN J/(KG K) OF SEA AIR

!OUTPUT: 
!ETA_SI  ENTROPY OF HUMID AIR, IN J/(KG K)

!INPUTS: 
!SA_SI  ABSOLUTE SALINITY IN KG/KG OF THE SEAWATER FRACTION
! T_SI  ABSOLUTE TEMPERATURE IN K
! P_SI  ABSOLUTE PRESSURE OF HUMID AIR IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!SEA_AIR_ENTROPY_AIR_SI(0.035, 300, 1E5) = 293.218508134

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!SEA_AIR_ENTROPY_AIR_SI(0.035, 300, 1E5) = 293.150672957

real*8 sea_air_entropy_air_si
real*8 sa_si, t_si, p_si

sea_air_entropy_air_si = errorreturn

if (set_sea_air_eq_at_s_t_p(sa_si, t_si, p_si) == errorreturn) return

sea_air_entropy_air_si = air_g_entropy_si(equi_sea_air_a, t_si, p_si)

end function

!==========================================================================
function set_sea_air_eq_at_s_a_p(sa_si, a_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE EQUILIBRIUM STATE BETWEEN HUMID AIR AND SEAWATER,
!IN PARTICULAR, THE SEA-AIR TEMPERATURE AT WHICH HUMID AIR CONDENSATES
!AT THE SEA SURFACE (ALREADY BEFORE THE DEWPOINT IS REACHED)

!INPUTS: 
!SA_SI  ABSOLUTE SALINITY IN KG/KG
! A_SI  MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
! P_SI  ABSOLUTE PRESSURE OF HUMID AIR IN PA

real*8 set_sea_air_eq_at_s_a_p, sa_si, a_si, p_si
real*8 t_si, dl, dh, eps

integer maxit

if(equi_sea_air_done == isok .and. &
   sa_si == equi_sea_air_s .and. &
   a_si == equi_sea_air_a .and. &
   p_si == equi_sea_air_p) then
  !the requested state has already been computed earlier
  set_sea_air_eq_at_s_a_p = isok
  return
end if

call clear_sea_air_state !waste any previous state

set_sea_air_eq_at_s_a_p = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(a_si < 0d0 .or. a_si >= 1d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_air

!SET INITIAL TEMPERATURE GUESS
select case(ctrl_mode_temperature)
  case(0); t_si = aux_temperature_si(a_si, p_si)
  case(1); t_si = ctrl_init_t
  case default; return
end select
if(t_si == errorreturn) return
if(t_si <= 0d0) return

!SET INITIAL LIQUID DENSITY
dl = liq_density_si(t_si, p_si)
if(dl == errorreturn) return
if(dl <= 0d0) return

!SET INITIAL HUMID AIR DENSITY
dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0) then
  maxit = 100
elseif(ctrl_loop_maximum == (-1)) then
  set_sea_air_eq_at_s_a_p = isok
  call set_sea_air_state(sa_si, a_si, t_si, p_si, dl, dh)
  return
elseif(ctrl_loop_maximum > 0) then 
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_air_iteration_at_s_a_p(sa_si, a_si, p_si, maxit, eps, t_si, dl, dh) == errorreturn) then
  return
end if

call set_sea_air_state(sa_si, a_si, t_si, p_si, dl, dh)

set_sea_air_eq_at_s_a_p = isok

end function

!==========================================================================
function set_sea_air_eq_at_s_t_p(sa_si, t_si, p_si)
!==========================================================================
                                             
!THIS FUNCTION COMPUTES THE EQUILIBRIUM STATE BETWEEN HUMID AIR AND
!SEAWATER, IN PARTICULAR, THE SPECIFIC HUMIDITY OF SEA AIR.

!INPUTS: 
!SA_SI  ABSOLUTE SALINITY IN KG/KG
! T_SI  ABSOLUTE TEMPERATURE IN K
! P_SI  ABSOLUTE PRESSURE OF HUMID AIR IN PA

real*8 set_sea_air_eq_at_s_t_p, sa_si, t_si, p_si
real*8 a_si, dl, dh, eps

integer maxit

if(equi_sea_air_done == isok .and. &
  sa_si == equi_sea_air_s .and. &
  t_si == equi_sea_air_t .and. &
  p_si == equi_sea_air_p) then
  !the requested state has already been computed earlier
  set_sea_air_eq_at_s_t_p = isok
  return
end if

call clear_sea_air_state !waste any previous state

set_sea_air_eq_at_s_t_p = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_air

!SET INITIAL AIR FRACTION GUESS: USE PURE WATER
select case(ctrl_mode_humidity)
  case(0)
    a_si = aux_airfraction_si(t_si, p_si)
    if(a_si == errorreturn) return
    if(a_si < 0d0 .or. a_si > 1d0) return
  case(1)
    a_si = ctrl_init_a
  case default; return
end select
if(a_si == errorreturn) return
if(a_si < 0d0 .or. a_si > 1d0) return

!LIQUID DENSITY CAN BE COMPUTED FROM T, P
dl = liq_density_si(t_si, p_si)
if(dl == errorreturn) return
if(dl <= 0d0) return

dh = air_density_si(a_si, t_si, p_si)
if(dh == errorreturn) return
if(dh <= 0d0) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0) then
  maxit = 100
elseif(ctrl_loop_maximum == (-1)) then
  set_sea_air_eq_at_s_t_p = isok
  call set_sea_air_state(sa_si, a_si, t_si, p_si, dl, dh)
  return
elseif(ctrl_loop_maximum > 0) then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT FOR THE VAPOUR PARTIAL PRESSURE
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_air_iteration_at_s_t_p(sa_si, t_si, p_si, dl, maxit, eps, a_si, dh) == errorreturn) then
  return
end if

call set_sea_air_state(sa_si, a_si, t_si, p_si, dl, dh)

set_sea_air_eq_at_s_t_p = isok

end function

!==========================================================================
function sea_air_iteration_at_s_t_p(sa_si, t_si, p_si, d_liq_si, maxit, eps, a_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE SEAWATER-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURES,
!TEMPERATURES AND CHEMICAL POTENTIAL OF WATER OF THE TWO PHASES AT GIVEN SALINITY, SA_SI,
!PRESSURE, P_SI, AND  TEMPERATURE, T_SI, FROM INITIAL GUESSES FOR THE HUMID-AIR DENSITY,
!D_AIR_SI, AND THE AIR FRACTION, A_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PARTIAL PRESSURE

!OUTPUT: 
!        SEA_AIR_ITERATION_AT_S_T_P = ISOK IF SUCCESSFUL
!        SEA_AIR_ITERATION_AT_S_T_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO SAL_G_SI, FLU_F_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - DENSITIES OR AIR FRACTION HAVE TAKEN AN INVALID VALUE DURING THE ITERATION
!     A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
! D_AIR_SI: HUMID-AIR DENSITY AT THE EQUILIBRIUM

!INPUTS: 
!      SA_SI: ABSOLUTE SALINITY IN KG/KG
!       P_SI: ABSOLUTE PRESSURE IN PA
!       T_SI: ABSOLUTE TEMPERATURE IN K
!   D_LIQ_SI: LIQUID PURE-WATER DENSITY AT T_SI, P_SI
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       A_SI: INITIAL GUESS OF THE MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY

real*8 sea_air_iteration_at_s_t_p, sa_si, t_si, p_si, d_liq_si, eps, a_si, d_air_si
real*8 ah, fh_a, fh_aa, fh_ad
real*8 gl, gss
real*8 fh, fh_d, fh_dd
real*8 dh
real*8 pv, pv_old
real*8 a(2, 2), b(2), x(2)

integer it, maxit

sea_air_iteration_at_s_t_p = errorreturn
pv = errorreturn

!CHECK FOR INVALID INPUT VALUES
if((a_si < 0d0 .or. a_si > 1d0) .or. &
  (d_air_si <= 0d0)) then
  a_si = errorreturn
  d_air_si = errorreturn
  return
endif

if(check_limits == 1) then
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
else
  if(d_liq_si <= 0d0 .or. &
    sa_si < 0d0 .or. sa_si >= 1d0 .or. &
    t_si <= 0d0 .or. &
    p_si <= 0d0) then
    a_si = errorreturn
    d_air_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

ah = a_si
dh = d_air_si

!CHEMICAL POTENTIAL OF WATER
gl = flu_gibbs_energy_si(t_si, d_liq_si)
if(gl == errorreturn) go to 999
gss = sal_chempot_h2o_si(sa_si, t_si, p_si)
if(gss == errorreturn) go to 999

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 2D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, ah, t_si, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, ah, t_si, dh)
  if(fh_a == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, ah, t_si, dh)
  if(fh_d == errorreturn) go to 999
  fh_aa = air_f_si(2, 0, 0, ah, t_si, dh)
  if(fh_aa == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, ah, t_si, dh)
  if(fh_ad == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, ah, t_si, dh)
  if(fh_dd == errorreturn) go to 999
  
  !VAPOUR PARTIAL PRESSURE FOR ERROR CHECKING
  pv_old = pv
  pv = air_molfraction_vap_si(ah)
  if(pv == errorreturn) go to 999
  pv = p_si * pv
  
  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or. &
      (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      sea_air_iteration_at_s_t_p = isok
      d_air_si = dh
      a_si = ah
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = -ah * fh_aa
  a(1, 2) = fh_d - ah * fh_ad - p_si / dh**2
  a(2, 1) = dh * fh_ad
  a(2, 2) = 2d0 * fh_d + dh * fh_dd

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gl - fh - p_si / dh + ah * fh_a + gss
  b(2) = p_si / dh - dh * fh_d
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !matrix singular

  !UPDATE AIR FRACTION & DENSITY
  ah = ah + x(1)
  if(ah < 0d0 .or. ah > 1d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !AIR_LIMITS
  if(a_si < 0d0 .or. a_si > 1d0 .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    a_si = errorreturn
    d_air_si = errorreturn
    sea_air_iteration_at_s_t_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_air_iteration_at_s_a_p(sa_si, a_si, p_si, maxit, eps, t_si, &
                                    d_liq_si, d_air_si)
!==========================================================================

!THIS FUNCTION RETURNS THE SEAWATER-HUMID-AIR PHASE EQUILIBRIUM FROM EQUAL PRESSURES,
!TEMPERATURES AND CHEMICAL POTENTIAL OF WATER OF THE TWO PHASES AT GIVEN PRESSURE, P_SI,
!SALINITY, SA_SI, AND AIR FRACTION, A_SI, FROM INITIAL GUESSES FOR THE LIQUID DENSITY,
!D_LIQ_SI, THE HUMID-AIR DENSITY, D_AIR_SI, AND THE TEMPERATURE, T_SI

!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN THE TEMPERATURE

!OUTPUT: 
!        SEA_AIR_ITERATION_AT_S_A_P = ISOK IF SUCCESSFUL
!        SEA_AIR_ITERATION_AT_S_A_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - FUNCTION CALLS TO SAL_G_SI, FLU_F_SI OR AIR_F_SI HAVE RETURNED AN ERROR
!        - DENSITIES OR TEMPERATURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
! D_LIQ_SI: LIQUID DENSITY AT T_SI, P_SI
! D_AIR_SI: HUMID-AIR DENSITY AT THE SATURATION EQUILIBRIUM

!INPUTS: 
!      SA_SI: ABSOLUTE SALINITY IN KG/KG
!       P_SI: ABSOLUTE PRESSURE IN PA
!       A_SI: MASS FRACTION OF DRY AIR IN HUMID AIR, KG/KG
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF TEMPERATURE
!             EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN PA
!             EPS < 0: RELATIVE TEMPERATURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!   D_LIQ_SI: INITIAL GUESS OF LIQUID DENSITY AT P_SI, T_SI
!   D_AIR_SI: INITIAL GUESS OF HUMID-AIR DENSITY


real*8 sea_air_iteration_at_s_a_p, sa_si, a_si, p_si, eps, &
       t_si, d_liq_si, d_air_si
real*8 fh, fh_a, fh_t, fh_d
real*8 fh_ad, fh_at, fh_td, fh_dd
real*8 fl, fl_d, fl_t, fl_td, fl_dd
real*8 dl, dh, t
real*8 gss, gs_t, gs_st
real*8 a(3, 3), b(3), x(3)

integer it, maxit

sea_air_iteration_at_s_a_p = errorreturn

!CHECK FOR INVALID INPUT VALUES
if(d_air_si <= 0d0 .or. &
  d_liq_si <= 0d0 .or. &
  a_si < 0d0 .or. a_si >= 1d0 .or. &
  t_si <= 0d0) then
  t_si = errorreturn
  d_liq_si = errorreturn
  d_air_si = errorreturn
  return
endif

if(check_limits == 1) then
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax .or. &
    p_si <  sal_pmin .or. p_si > sal_pmax) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    return
  endif
else
  if((sa_si < 0d0 .or. sa_si >= 1d0) .or. &
    (p_si <= 0d0)) then
    t_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
dh = d_air_si
dl = d_liq_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR 3D-NEWTON ITERATION
  !HUMID AIR:
  fh = air_f_si(0, 0, 0, a_si, t, dh)
  if(fh == errorreturn) go to 999
  fh_a = air_f_si(1, 0, 0, a_si, t, dh)
  if(fh_a == errorreturn) go to 999
  fh_t = air_f_si(0, 1, 0, a_si, t, dh)
  if(fh_t == errorreturn) go to 999
  fh_d = air_f_si(0, 0, 1, a_si, t, dh)
  if(fh_d == errorreturn) go to 999
  fh_ad = air_f_si(1, 0, 1, a_si, t, dh)
  if(fh_ad == errorreturn) go to 999
  fh_at = air_f_si(1, 1, 0, a_si, t, dh)
  if(fh_at == errorreturn) go to 999
  fh_td = air_f_si(0, 1, 1, a_si, t, dh)
  if(fh_td == errorreturn) go to 999
  fh_dd = air_f_si(0, 0, 2, a_si, t, dh)
  if(fh_dd == errorreturn) go to 999

  !LIQUID WATER:
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fl_td = flu_f_si(0, 2, t, dl)
  if(fl_td == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999
  
  !SEAWATER
  gss = sal_chempot_h2o_si(sa_si, t, p_si)
  if(gss == errorreturn) go to 999
  gs_t = sal_g_si(0, 1, 0, sa_si, t, p_si)
  if(gs_t == errorreturn) go to 999
  gs_st = sal_g_si(1, 1, 0, sa_si, t, p_si)
  if(gs_st == errorreturn) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = fh_t - a_si * fh_at - fl_t - gs_t + sa_si * gs_st
  a(1, 2) = fh_d - a_si * fh_ad - p_si / dh**2
  a(1, 3) = -fl_d + p_si / dl**2
  a(2, 1) = dh * fh_td
  a(2, 2) = 2d0 * fh_d + dh * fh_dd
  a(2, 3) = 0d0
  a(3, 1) = dl * fl_td
  a(3, 2) = 0d0
  a(3, 3) = 2d0 * fl_d + dl * fl_dd

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = fl - fh + a_si * fh_a + p_si * (1d0 / dl - 1d0 / dh) + gss
  b(2) = p_si / dh - dh * fh_d
  b(3) = p_si / dl - dl * fl_d
  
  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 3) /= 0) go to 999 !matrix singular

  !UPDATE TEMPERATURE & DENSITIES
  t = t + x(1)
  if(t <= 0d0) go to 999
  dh = dh + x(2)
  if(dh <= 0d0) go to 999
  dl = dl + x(3)
  if(dl <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(x(1)) < eps) .or. &
    (eps < 0d0 .and. abs(x(1)) < -eps * t)) then
    sea_air_iteration_at_s_a_p = isok
    d_air_si = dh
    d_liq_si = dl
    t_si = t
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !AIR_LIMITS
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax .or. &
  d_air_si <= dry_air_dmin .or. d_air_si > dry_air_dmax) then
    a_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    sea_air_iteration_at_s_a_p = errorreturn
    return
  endif
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    a_si = errorreturn
    d_liq_si = errorreturn
    d_air_si = errorreturn
    sea_air_iteration_at_s_a_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_airfraction_si(t_si, p_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE AIR MASS FRACTION OF SATURATED
!HUMID AIR AT GIVEN ABSOLUTE TEMPERATURE, T_SI, IN K AND ABSOLUTE
!PRESSURE, P_SI, IN PA

real*8 aux_airfraction_si, t_si, p_si
real*8 dh, dv, pv

aux_airfraction_si = errorreturn
if(t_si <= 0d0) return
if(p_si <= 0d0) return

!AS AN ESTIMATE, ADD THE DENSITIES OF DRY AIR AND VAPOUR

!VAPOUR PRESSURE
pv = liq_vap_vapourpressure_si(t_si)
if(pv == errorreturn) return
if(pv > p_si) return

if(pv == p_si) then
  aux_airfraction_si = 0d0
  return
end if

!VAPOUR DENSITY
dv = liq_vap_density_vap_si()
if(dv == errorreturn) return
if(dv <= 0d0) return

!DRY AIR DENSITY AT REMAINING PARTIAL PRESSURE
dh = air_density_si(1d0, t_si, p_si - pv)
if(dh == errorreturn) return
if(dh <= 0d0) return

!AIR FRACTION OF THE IDEAL MIXTURE
aux_airfraction_si = dh / (dh + dv)
           
end function

!==========================================================================
function aux_temperature_si(a_si, p_si)
!==========================================================================

!FUNCTION RETURNS AN ESTIMATE OF THE TEMPERATURE IN K OF SATURATED
!HUMID AIR AT GIVEN AIR MASS FRACTION, A_SI, IN KG/KG AND ABSOLUTE
!PRESSURE, P_SI, IN PA

real*8 aux_temperature_si, a_si, p_si
real*8 pv, t

aux_temperature_si = errorreturn
if(a_si < 0d0 .or. a_si >= 1d0) return
if(p_si <= 0d0) return

!AS AN ESTIMATE, IGNORE THE AIR PRESENT

!PARTIAL VAPOUR PRESSURE
pv = p_si * air_molfraction_vap_si(a_si)
!RELATED SATURATION TEMPERATURE
t = liq_vap_boilingtemperature_si(pv)
if(t == errorreturn) return
if(t <= 0d0) return

aux_temperature_si = t

end function

!==========================================================================
subroutine set_sea_air_state(s, a, t, p, dl, dh)
!==========================================================================

real*8 s, a, t, p, dl, dh

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR
equi_sea_air_done = isok

equi_sea_air_s = s
equi_sea_air_a = a
equi_sea_air_t = t
equi_sea_air_p = p
equi_sea_air_d_liq = dl
equi_sea_air_d_air = dh

end subroutine

!==========================================================================
subroutine clear_sea_air_state()
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_sea_air_done = 0

end subroutine

!==========================================================================
subroutine init_it_ctrl_sea_air()
!==========================================================================

real*8 t_triple, p_triple

!TRIPLE POINT VALUES
t_triple = tp_temperature_si
p_triple = tp_pressure_iapws95_si

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR ITERATION
ctrl_mode_pressure = 0
ctrl_mode_temperature = 0
ctrl_mode_humidity = 0

ctrl_loop_maximum = 100
ctrl_init_t = t_triple
ctrl_init_p = p_triple
ctrl_init_a = 1d0 - 0.001d0

ctrl_eps_exit_p_vap = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_temp =  -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_press = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

!==========================================================================
subroutine set_it_ctrl_sea_air(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!SATURATED HUMID AIR PROPERTIES

!KEY              VALUE
!IT_STEPS         0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS         N > 0    SET ITERATION NUMBER TO N
!IT_STEPS        -1           DO NOT ITERATE, USE INITIAL VALUES AS THE SATURATION STATE

!INIT_AIR         0           USE DEFAULT AIR FRACTION TO START ( = AUX_AIRFRACTION_SI(T, P))
!INIT_AIR         A > 0    USE VALUE A AS AIR FRACTION TO START

!INIT_TEMP        0           USE DEFAULT TEMPERATURE TO START ( = AUX_TEMPERATURE_SI(A, P))
!INIT_TEMP        T > 0    USE VALUE T AS TEMPERATURE TO START

!INIT_PRESS       0           USE DEFAULT PRESSURE TO START ( = AUX_PRESSURE_SI(A, T))
!INIT_PRESS       P > 0    USE VALUE P AS PRESSURE TO START

!TOL_VAP_PRESS    0           USE DEFAULT EXIT ACCURACY FOR VAPOUR PRESSURE (0.1 PPM)
!TOL_VAP_PRESS    EPS         USE EPS AS EXIT ACCURACY FOR VAPOUR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_TEMP         0           USE DEFAULT EXIT ACCURACY FOR TEMPERATURE (0.1 PPM)
!TOL_TEMP         EPS         USE EPS AS EXIT ACCURACY FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!TOL_PRESS        0           USE DEFAULT EXIT ACCURACY FOR PRESSURE (0.1 PPM)
!TOL_PRESS        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

real*8 value
character*13 key

call init_it_ctrl_sea_air

call clear_sea_air_state

select case (trim(key))

  case ("it_steps");                !ITERATION STEPS
    if(value == 0d0) then
      ctrl_loop_maximum = 100    !DEFAULT = 100
    elseif(value < 0d0) then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ("init_air");                !START AIR FRACTION
    if(value == 0d0) then
      ctrl_mode_humidity = 0     !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_humidity = 1
      ctrl_init_a = value
    endif
    
  case ("init_temp");               !START TEMPERATURE
    if(value == 0d0) then
      ctrl_mode_temperature = 0  !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif
    
  case ("init_press");              !START PRESSURE
    if(value == 0d0) then
      ctrl_mode_pressure = 0     !DEFAULT = AUX POLYNOMIAL
    else
      ctrl_mode_pressure = 1
      ctrl_init_p = value
    endif

  case ("tol_vap_press");           !REQUIRED VAPOUR PRESSURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_p_vap = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p_vap = value
    endif

  case ("tol_temp");                !REQUIRED TEMPERATURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_temp = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_temp = value
    endif

  case ("tol_press");               !REQUIRED TOTAL PRESSURE TOLERANCE
    if(value == 0d0) then
      ctrl_eps_exit_press = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_press = value

    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

end module sea_air_4





