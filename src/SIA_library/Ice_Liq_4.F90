Module ice_liq_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIA PROPERTIES OF ICE WITH LIQUID WATER
!COMPUTED FROM IAPWS-95 AND IAPWS-06:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

!RELEASE ON AN EQUATION OF STATE FOR H2O ICE IH
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!WITNEY, UK, SEPTEMBER 2006

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
!     CONSTANTS_0,     FILE CONSTANTS_0.F90
!     MATHS_0,         FILE MATHS_0.F90
!     FLU_1,           FILE FLU_1.F90
!     ICE_1,           FILE ICE_1.F90
!     FLU_2,           FILE FLU_2.F90
!     ICE_2,           FILE ICE_2.F90

!#########################################################################

use constants_0
use maths_0
use flu_1
use ice_1
use flu_2
use ice_2

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: ice_liq_chempot_si, ice_liq_density_ice_si, &
          ice_liq_density_liq_si, ice_liq_enthalpy_ice_si, &
          ice_liq_enthalpy_liq_si, ice_liq_enthalpy_melt_si, &
          ice_liq_entropy_ice_si, ice_liq_entropy_liq_si, &
          ice_liq_entropy_melt_si, ice_liq_meltingpressure_si, &
          ice_liq_meltingtemperature_si, ice_liq_pressure_liq_si, &
          ice_liq_temperature_si, ice_liq_volume_melt_si, &
          set_ice_liq_eq_at_p, set_ice_liq_eq_at_t, &
          set_it_ctrl_ice_liq 

real*8, private :: equi_ice_liq_done, equi_ice_liq_t, equi_ice_liq_p, &
        equi_ice_liq_g_ice, equi_ice_liq_g_liq, equi_ice_liq_d_ice, &
        equi_ice_liq_d_liq

integer, private :: ctrl_mode_temperature, ctrl_mode_liquid, ctrl_mode_ice, &
         ctrl_loop_maximum, ctrl_initialized
real*8, private :: ctrl_init_t, ctrl_init_d_liq, ctrl_eps_exit_p

contains

!==========================================================================
function ice_liq_chempot_si()
!==========================================================================

!RETURNS THE CHEMICAL POTENTIAL IN J/KG LIQUID WATER IN LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_chempot_si = 9972.8817138|1

!set_ice_liq_eq_at_t 270
!ice_liq_chempot_si = 38870.060519|2

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_chempot_si

ice_liq_chempot_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_chempot_si = equi_ice_liq_g_liq

return

end function

!==========================================================================
function ice_liq_density_ice_si()
!==========================================================================

!RETURNS THE DENSITY  IN KG/M3 OF ICE AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_density_ice_si = 917.896690831

!set_ice_liq_eq_at_t 270
!ice_liq_density_ice_si = 921.359428514

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_density_ice_si

ice_liq_density_ice_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_density_ice_si = equi_ice_liq_d_ice

return

end function

!==========================================================================
function ice_liq_density_liq_si()
!==========================================================================

!RETURNS THE DENSITY IN KG/M3 OF WATER AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_density_liq_si = 1004.79353660

!set_ice_liq_eq_at_t 270
!ice_liq_density_liq_si = 1019.05568894

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_density_liq_si

ice_liq_density_liq_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_density_liq_si = equi_ice_liq_d_liq

return

end function

!==========================================================================
function ice_liq_enthalpy_ice_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY J/KG OF ICE AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_enthalpy_ice_si = -324602.983812

!set_ice_liq_eq_at_t 270
!ice_liq_enthalpy_ice_si = -299055.938629

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_enthalpy_ice_si

ice_liq_enthalpy_ice_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_enthalpy_ice_si = ice_enthalpy_si(equi_ice_liq_t, equi_ice_liq_p)

return

end function

!==========================================================================
function ice_liq_enthalpy_liq_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY IN J/KG OF WATER AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_enthalpy_liq_si = 6945.9270048|3

!set_ice_liq_eq_at_t 270
!ice_liq_enthalpy_liq_si = 26110.7481094

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_enthalpy_liq_si

ice_liq_enthalpy_liq_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_enthalpy_liq_si = flu_enthalpy_si(equi_ice_liq_t, equi_ice_liq_d_liq)

return

end function

!==========================================================================
function ice_liq_enthalpy_melt_si()
!==========================================================================

!RETURNS THE SPECIFIC MELTING ENTHALPY OF WATER IN J/KG
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_enthalpy_melt_si = 331548.910817

!set_ice_liq_eq_at_t 270
!ice_liq_enthalpy_melt_si = 325166.686739

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_enthalpy_melt_si
real*8 hi, hl

ice_liq_enthalpy_melt_si = errorreturn

if(equi_ice_liq_done /= isok) return

hi = ice_enthalpy_si(equi_ice_liq_t, equi_ice_liq_p)
if(hi == errorreturn) return
hl = flu_enthalpy_si(equi_ice_liq_t, equi_ice_liq_d_liq)
if(hl == errorreturn) return

ice_liq_enthalpy_melt_si = hl - hi

return

end function

!==========================================================================
function ice_liq_entropy_ice_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY IN J/(KG K) OF ICE AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_entropy_ice_si = -1228.24464139

!set_ice_liq_eq_at_t 270
!ice_liq_entropy_ice_si = -1251.57777462

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_entropy_ice_si

ice_liq_entropy_ice_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_entropy_ice_si = ice_entropy_si(equi_ice_liq_t, equi_ice_liq_p)

return

end function

!==========================================================================
function ice_liq_entropy_liq_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY J/(KG K) OF WATER AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_entropy_liq_si = -11.112101272|3

!set_ice_liq_eq_at_t 270
!ice_liq_entropy_liq_si = -47.2567126291

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_entropy_liq_si

ice_liq_entropy_liq_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_entropy_liq_si = flu_entropy_si(equi_ice_liq_t, equi_ice_liq_d_liq)

return

end function

!==========================================================================
function ice_liq_entropy_melt_si()
!==========================================================================

!RETURNS THE SPECIFIC MELTING ENTROPY IN J/(KG K) OF WATER
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_entropy_melt_si = 1217.13254011

!set_ice_liq_eq_at_t 270
!ice_liq_entropy_melt_si = 1204.32106199

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_entropy_melt_si
real*8 si, sl

ice_liq_entropy_melt_si = errorreturn

if(equi_ice_liq_done /= isok) return

si = ice_entropy_si(equi_ice_liq_t, equi_ice_liq_p)
if(si == errorreturn) return
sl = flu_entropy_si(equi_ice_liq_t, equi_ice_liq_d_liq)
if(sl == errorreturn) return

ice_liq_entropy_melt_si = sl - si

return

end function

!==========================================================================
function ice_liq_meltingpressure_si(t_si)
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE OF WATER IN LIQUID-ICE EQUILIBRIUM IN PA
!AS A FUNCTION OF ABSOLUTE TEMPERATURE T_SI IN K

!CHECK VALUE WITH DEFAULT SETTINGS: ice_liq_meltingpressure_si(270) = 39313338.8825

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_meltingpressure_si
real*8 t_si

ice_liq_meltingpressure_si = errorreturn

if(set_ice_liq_eq_at_t(t_si) == errorreturn) return

ice_liq_meltingpressure_si = equi_ice_liq_p

return

end function

!==========================================================================
function ice_liq_meltingtemperature_si(p_si)
!==========================================================================

!RETURNS THE ABSOLUTE TEMPERATURE OF LIQUID-ICE EQUILIBRIUM IN K
!AS A FUNCTION OF ABSOLUTE PRESSURE P_SI IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: ice_liq_meltingtemperature_si(1E7) = 272.401648869

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_meltingtemperature_si
real*8 p_si

ice_liq_meltingtemperature_si = errorreturn

if(set_ice_liq_eq_at_p(p_si) == errorreturn) return

ice_liq_meltingtemperature_si = equi_ice_liq_t

return

end function

!==========================================================================
function ice_liq_pressure_liq_si()
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE IN PA OF WATER AT THE  LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_pressure_liq_si = 10000000.0070

!set_ice_liq_eq_at_t 270
!ice_liq_pressure_liq_si = 39313338.8825

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_pressure_liq_si

ice_liq_pressure_liq_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_pressure_liq_si = equi_ice_liq_p

return

end function

!==========================================================================
function ice_liq_temperature_si()
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE OF WATER IN LIQUID-ICE EQUILIBRIUM IN PA
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_temperature_si = 272.401648869

!set_ice_liq_eq_at_t 270
!ice_liq_temperature_si = 270

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_temperature_si

ice_liq_temperature_si = errorreturn

if(equi_ice_liq_done /= isok) return

ice_liq_temperature_si = equi_ice_liq_t

return

end function

!==========================================================================
function ice_liq_volume_melt_si()
!==========================================================================

!RETURNS THE SPECIFIC MELTING VOLUME IN M3/KG AT THE LIQUID-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_LIQ_EQUILIBRIUM_TO_T OR SET_ICE_LIQ_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_liq_eq_at_p 1e7
!ice_liq_volume_melt_si = -9.42178903288E-05

!set_ice_liq_eq_at_t 270
!ice_liq_volume_melt_si = -1.04052121182E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_liq_volume_melt_si

ice_liq_volume_melt_si = errorreturn

if(equi_ice_liq_done /= isok) return

if(equi_ice_liq_d_ice <= 0d0) return
if(equi_ice_liq_d_liq <= 0d0) return

ice_liq_volume_melt_si = 1d0 / equi_ice_liq_d_liq - 1d0 / equi_ice_liq_d_ice

return

end function

!==========================================================================
function set_ice_liq_eq_at_p(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ICE-LIQUID EQUILIBRIUM AT GIVEN PRESSURE
!BY 2-DIMENSIONAL NEWTON ITERATION

implicit none

integer maxit
real*8 p_si
real*8 set_ice_liq_eq_at_p
real*8  t, pl, gi, gl, di, dl, eps
real *8 dl_triple, t_triple

!TRIPLE POINT
dl_triple = tp_density_liq_iapws95_si
t_triple = tp_temperature_si

if(equi_ice_liq_done == isok .and. &
   p_si == equi_ice_liq_p) then
  !the requested state has already been computed earlier
  set_ice_liq_eq_at_p = isok
  return
end if

call clear_ice_liq_state !waste any previous state

set_ice_liq_eq_at_p = errorreturn

if(p_si <= 0d0) return

call init_it_ctrl_ice_liq

!SET INITIAL TEMPERATURE GUESS
select case (ctrl_mode_temperature)
  case (0);  t = aux_temperature_correlation(p_si)
             if(t == errorreturn) return
  case (-1); t = t_triple
  case (1);  t = ctrl_init_t
  case default; return
end select

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  dl = aux_liq_density_correlation(t)
  case (-1); dl = dl_triple
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
   maxit = 100
elseif (ctrl_loop_maximum == (-1))then
   pl = flu_pressure_si(t, dl)
   gi = ice_g_si(0, 0, t, pl)
   gl = flu_gibbs_energy_si(t, dl)
   di = ice_g_si(0, 1, t, pl)
   if(di == errorreturn .or. di <= 0d0) return
   di = 1d0 / di
   call set_ice_liq_state (t, pl, gi, gl, di, dl)
   set_ice_liq_eq_at_p = isok
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
if(ice_liq_iteration_at_p(p_si, maxit, eps, t, di, dl, gi, gl, pl) == errorreturn) then
  return
end if

call set_ice_liq_state (t, pl, gi, gl, di, dl)

set_ice_liq_eq_at_p = isok

return

end function

!==========================================================================
function set_ice_liq_eq_at_t(t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ICE-LIQUID EQUILIBRIUM AT GIVEN TEMPERATURE
!BY NEWTON ITERATION

implicit none

integer maxit
real*8 t_si
real*8 set_ice_liq_eq_at_t

real*8 pl, gi, gl, di, dl, eps
real*8 t_triple, dl_triple

!TRIPLE POINT
t_triple = tp_temperature_si
dl_triple = tp_density_liq_iapws95_si

if(equi_ice_liq_done == isok .and. &
  t_si == equi_ice_liq_t) then
  !the requested state has already been computed earlier
  set_ice_liq_eq_at_t = isok
  return
end if

call clear_ice_liq_state !waste any previous state

set_ice_liq_eq_at_t = errorreturn

if(t_si <= 0d0) return
if(t_si > t_triple) return

call init_it_ctrl_ice_liq

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  dl = aux_liq_density_correlation(t_si)
  case (-1); dl = dl_triple
  case (1);  dl = ctrl_init_d_liq
  case default; return
end select

!SET MAX. ITERATION NUMBER
if (ctrl_loop_maximum == 0)then
   maxit = 100
elseif (ctrl_loop_maximum == (-1))then
   pl = flu_pressure_si(t_si, dl)
   gi = ice_g_si(0, 0, t_si, pl)
   gl = flu_gibbs_energy_si(t_si, dl)
   di = ice_g_si(0, 1, t_si, pl)
   if(di == errorreturn .or. di <= 0d0) return
   di = 1d0 / di
   call set_ice_liq_state (t_si, pl, gi, gl, di, dl)
   set_ice_liq_eq_at_t = isok
   return
elseif (ctrl_loop_maximum > 0)then
   maxit = ctrl_loop_maximum
else
   return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_liq_iteration_at_t(t_si, maxit, eps, di, dl, gi, gl, pl) == errorreturn) then
  return
end if

call set_ice_liq_state (t_si, pl, gi, gl, di, dl)

set_ice_liq_eq_at_t = isok

return

end function

!==========================================================================
function ice_liq_iteration_at_p(p_si, maxit, eps, t_si, d_ice_si, d_liq_si, &
                                g_ice_si, g_liq_si, p_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-LIQUID PHASE EQUILIBRIUM FROM EQUAL PRESSURES,
!TEMPERATURES AND CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN PRESSURE, P_SI,
!FROM INITIAL GUESSES FOR THE LIQUID DENSITY, D_LIQ_SI, AND THE TEMPERATURE, T_SI
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN LIQUID PRESSURE

!OUTPUT: ICE_LIQ_ITERATION_AT_T = ISOK IF SUCCESSFUL
!        ICE_LIQ_ITERATION_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI OR TO ICE_G_SI HAS RETURNED AN ERROR
!        - DENSITY OR TEMPERATURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
! D_ICE_SI: ICE DENSITY IN KG/M3 AT THE ICE-LIQUID EQUILIBRIUM
! D_LIQ_SI: LIQUID DENSITY IN KG/M3 AT THE ICE-LIQUID EQUILIBRIUM
! G_ICE_SI: ICE CHEMICAL POTENTIAL IN J/KG AT THE ICE-LIQUID EQUILIBRIUM
! G_LIQ_SI: LIQUID CHEMICAL POTENTIALIN J/KG AT THE ICE-LIQUID EQUILIBRIUM
! P_LIQ_SI: LIQUID PRESSURE IN PA AT THE ICE-LIQUID EQUILIBRIUM

!INPUT: P_SI: ABSOLUTE PRESSURE IN PA
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF LIQUID PRESSURE
!             EPS > 0: ABSOLUTE LIQUID PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE LIQUID PRESSURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!   D_LIQ_SI: INITIAL GUESS OF LIQUID DENSITY IN KG/M3

implicit none

integer it, maxit
real*8 ice_liq_iteration_at_p
real*8 p_si, eps, t_si, d_ice_si, d_liq_si,  &
       g_ice_si, g_liq_si, p_liq_si
real*8 vi, dl, pl, gl, gi, gi_t, t, fl, fl_d, fl_dd, &
       fl_t, fl_td, a(2, 2), b(2), x(2), ddl, dt, pl_old 
real*8 tt

tt = tp_temperature_si !TRIPLE POINT TEMPERATURE IN K

ice_liq_iteration_at_p = errorreturn

if(t_si <= 0d0 .or. &
  d_liq_si <= 0d0) then
  t_si = errorreturn
  d_ice_si = errorreturn
  d_liq_si = errorreturn
  g_ice_si = errorreturn
  g_liq_si = errorreturn
  p_liq_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

dl = d_liq_si
t = t_si

do it = 0, maxit

  !DERIVATIVES OF THE GIBBS & HELMHOLTZ FUNCTIONS FOR 2D-NEWTON ITERATION
  fl = flu_f_si(0, 0, t, dl)
  if(fl == errorreturn) go to 999
  
  fl_d = flu_f_si(0, 1, t, dl)
  if(fl_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t, dl)
  if(fl_dd == errorreturn) go to 999
  
  fl_t = flu_f_si(1, 0, t, dl)
  if(fl_t == errorreturn) go to 999
  fl_td = flu_f_si(1, 1, t, dl)
  if(fl_td == errorreturn) go to 999

  pl_old = pl             !KEEP PREVIOUS PL TO CHECK ACCURACY
  !LIQUID AND ICE PRESSURE + LIQUID CHEMICAL POTENTIALS
  pl = dl ** 2 * fl_d     !PRESSURE OF LIQUID
  gl = fl + dl * fl_d     !CHEM. POT. OF LIQUID
  
  gi = ice_g_si(0, 0, t, pl)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, pl)
  if(gi_t == errorreturn) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pl - pl_old) < eps) .or. &
       (eps < 0 .and. abs(pl - pl_old) < -eps * pl)) then
      p_liq_si = pl
      g_ice_si = gi
      g_liq_si = gl
      vi = ice_g_si(0, 1, t, pl)
      if(vi == errorreturn .or. vi <= 0d0) go to 999
      d_ice_si = 1d0 / vi
      d_liq_si = dl
      t_si = t
      ice_liq_iteration_at_p = isok
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = dl * (2d0 * fl_d + dl * fl_dd)
  a(1, 2) = dl ** 2 * fl_td
  a(2, 1) = (2d0 * fl_d + dl * fl_dd)
  a(2, 2) = -gi_t + fl_t + dl * fl_td

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = p_si - pl
  b(2) = gi - gl

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR
  ddl = x(1)
  dt = x(2)
  
  !UPDATE DENSITY & TEMPERATURE
  dl = dl + ddl
  if(dl <= 0d0) go to 999
  t = t + dt
  if(t <= 0d0) go to 999
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  p_liq_si <= 0d0 .or. &
  d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    ice_liq_iteration_at_p = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  d_ice_si <= 0d0) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    ice_liq_iteration_at_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_liq_iteration_at_t(t_si, maxit, eps, d_ice_si, d_liq_si,  &
                                g_ice_si, g_liq_si, p_liq_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-LIQUID PHASE EQUILIBRIUM FROM EQUAL PRESSURES
!AND CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN TEMPERATURE, T_SI,
!FROM AN INITIAL GUESS FOR THE LIQUID DENSITY, D_LIQ_SI.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN LIQUID PRESSURE

!OUTPUT: ICE_LIQ_ITERATION = ISOK IF SUCCESSFULLY DONE
!        ICE_LIQ_ITERATION = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI HAS RETURNED AN ERROR
!        - DENSITIES OR PRESSURES HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
! D_LIQ_SI: LIQUID DENSITY IN KG/M3 AT THE ICE-LIQUID EQUILIBRIUM
! D_ICE_SI: ICE DENSITY IN KG/M3 AT THE ICE-LIQUID EQUILIBRIUM
! G_LIQ_SI: LIQUID CHEMICAL POTENTIAL IN J/KG AT THE ICE-LIQUID EQUILIBRIUM
! G_ICE_SI: ICE CHEMICAL POTENTIAL IN J/KG AT THE ICE-LIQUID EQUILIBRIUM
! P_LIQ_SI: LIQUID PRESSURE IN PA AT THE ICE-LIQUID EQUILIBRIUM

!INPUT: T_SI: ABSOLUTE TEMPERATURE IN K
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF PRESSURE
!             EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN KG/M3
!             EPS < 0: RELATIVE PRESSURE TOLERANCE
!   D_LIQ_SI: INITIAL GUESS OF LIQUID DENSITY IN KG/M3

implicit none

integer it, maxit
real*8 t_si, eps
real*8 d_ice_si, d_liq_si,g_ice_si, g_liq_si, p_liq_si
real*8 ice_liq_iteration_at_t
real*8 vi, dl, pl, gi, gl, fl, fl_d, fl_dd, ddl, pl_old

ice_liq_iteration_at_t = errorreturn

If(t_si <= 0) return

if(d_liq_si <= 0d0) then
  d_ice_si = errorreturn
  d_liq_si = errorreturn
  g_ice_si = errorreturn
  g_liq_si = errorreturn
  p_liq_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    return
  endif
else
  if(t_si <= 0d0) then
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

dl = d_liq_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR NEWTON ITERATION
  fl = flu_f_si(0, 0, t_si, dl)
  if(fl == errorreturn) go to 999
  fl_d = flu_f_si(0, 1, t_si, dl)
  if(fl_d == errorreturn) go to 999
  fl_dd = flu_f_si(0, 2, t_si, dl)
  if(fl_dd == errorreturn) go to 999
  
  pl_old = pl            !KEEP PREVIOUS PL TO CHECK ACCURACY
  !LIQUID AND ICE PRESSURE + LIQUID CHEMICAL POTENTIALS
  pl = dl ** 2 * fl_d    !PRESSURE OF LIQUID
  gl = fl + dl * fl_d    !CHEM. POT. OF LIQUID
  
  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t_si, pl)
  if(gi == errorreturn) go to 999
  vi = ice_g_si(0, 1, t_si, pl)
  if(vi == errorreturn) go to 999
  if(vi <= 0d0) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pl - pl_old) < eps) .or. &
       (eps < 0 .and. abs(pl - pl_old) < -eps * pl)) then
      ice_liq_iteration_at_t = isok
      p_liq_si = pl
      g_ice_si = gi
      g_liq_si = gl
      d_ice_si = 1d0 / vi
      d_liq_si = dl
      go to 999
    end if
  end if

  if(it == maxit) go to 999
  
  !DO ITERATION STEP
  ddl = (2d0 * fl_d + dl * fl_dd) * (vi * dl - 1d0)
  if(ddl == 0d0) go to 999
  ddl = (gl - gi) / ddl
  
  !UPDATE LIQUID DENSITY
  dl = dl + ddl
  if(dl <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_liq_si <= flu_dmin .or. d_liq_si > flu_dmax .or. &
  p_liq_si <= 0d0) then
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    ice_liq_iteration_at_t = errorreturn
    return
  endif
  !ICE_LIMITS
  if(d_ice_si <= 0d0) then
    d_ice_si = errorreturn
    d_liq_si = errorreturn
    g_ice_si = errorreturn
    g_liq_si = errorreturn
    p_liq_si = errorreturn
    ice_liq_iteration_at_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_liq_density_correlation(t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE DENSITY IN KG/M3 OF THE FREEZING LIQUID
!AS A FUNCTION OF ABSOLUTE TEMPERATURE T_SI IN K FROM AN APPROXIMATE
!CORRELATION FIT BETWEEN 252 AND 273 K WITH RMS = 1.2 E-4 IN D/DT

implicit none

real*8 aux_liq_density_correlation
real*8 t_si, tt, dt, a1, a2, a3, tr, dr

tt = tp_temperature_si          !TRIPLE POINT TEMPERATURE IN K
dt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER

a1 = -1.78582981492113d0
a2 = -12.2325084306734d0
a3 = -52.8236936433529d0

aux_liq_density_correlation = errorreturn

if(t_si <= 0d0) return

tr = t_si / tt - 1d0
dr = a1 * tr + a2 * tr ** 2 + a3 * tr ** 3

aux_liq_density_correlation = dt * (dr + 1d0)

return

end function

!==========================================================================
function aux_temperature_correlation(p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE MELTING TEMPERATURE IN K OF ICE AS AN APPROXIMATE
!FUNCTION OF THE ABSOLUTE PRESSURE IN PA FROM A CORRELATION FIT
!BETWEEN 273 K AND 252 K WITH RMS = 1.4E-5 IN T/TT

implicit none

real*8 aux_temperature_correlation
real*8 p_si
real*8 pt, tt, a1, a2, tr, pr

pt = tp_pressure_iapws95_si  !IAPWS-95 TRIPLE POINT PRESSURE IN PA
tt = tp_temperature_si       !TRIPLE POINT TEMPERATURE OF WATER

a1 = -1.67329759176351d-07
a2 = -2.02262929999658d-13

aux_temperature_correlation = errorreturn

if(p_si <= 0d0) return

pr = p_si / pt - 1d0
tr = a1 * pr + a2 * pr ** 2

aux_temperature_correlation = tt * (tr + 1d0)

return

end function

!==========================================================================
subroutine init_it_ctrl_ice_liq()
!==========================================================================

implicit none

real*8 dl_triple, t_triple 

!TRIPLE POINT
dl_triple = tp_density_liq_iapws95_si
t_triple = tp_temperature_si

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR ITERATION
ctrl_mode_ice = 0
ctrl_mode_liquid = 0
ctrl_mode_temperature = 0
ctrl_loop_maximum = 100
ctrl_init_d_liq = dl_triple
ctrl_init_t = t_triple
ctrl_eps_exit_p = -0.0000001d0 !relative, 0.1 ppm

end subroutine

!==========================================================================
subroutine set_it_ctrl_ice_liq(key, value)
!==========================================================================

implicit none

character*13 key
real*8 value

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!ICE-LIQUID EQUILIBRIUM

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL LIQUID DENSITY TO COMPUTE LIQUID PRESSURE

!'INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY TO START =
!                                   LIQUID DENSITY CORRELATION FUNCTION, D = AUX_LIQ_DENSITY_CORRELATION(T)
!'INIT_LIQ_DENS'       -1           USE TRIPLE POINT LIQUID DENSITY TO START
!'INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!'INIT_TEMP    '        0           USE DEFAULT TEMPERATURE TO START ( = AUX_MELT_TEMPERATURE(P))
!'INIT_TEMP    '       -1           USE TRIPLE POINT TEMPERATURE TO START
!'INIT_TEMP    '        T > 0       USE VALUE T AS TEMPERATURE TO START

!'TOL_LIQ_PRESS'        0           USE DEFAULT EXIT ACCURACY FOR LIQUID PRESSURE (0.1 PPM)
!'TOL_LIQ_PRESS'        EPS         USE EPS AS EXIT ACCURACY FOR LIQUID PRESSURE (EPS < 0 MEANS RELATIVE ERROR)


call init_it_ctrl_ice_liq

call clear_ice_liq_state

select case (trim(key))

  case ('it_steps')                                  !ITERATION STEPS
    if(value == 0d0)then; ctrl_loop_maximum = 100    !DEFAULT = 100
    elseif(value < 0d0)then; ctrl_loop_maximum = -1
    else;   ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')                             !START LIQUID DENSITY
    if(value == 0d0)then; ctrl_mode_liquid = 0       !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then; !ignore it
    elseif(value < 0d0)then;  ctrl_mode_liquid = nint(value)
    else; ctrl_mode_liquid = 1; ctrl_init_d_liq = value
    endif

  case ('init_temp')                                 !START TEMPERATURE
    if(value == 0d0)then; ctrl_mode_temperature = 0  !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then;                     !IGNORE IT
    elseif(value < 0d0)then; ctrl_mode_temperature = nint(value)
    else; ctrl_mode_temperature = 1; ctrl_init_t = value
    endif

  case ('tol_liq_press')                             !REQUIRED LIQUID PRESSURE TOLERANCE
    if(value == 0d0)then; ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else; ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine set_ice_liq_state(t, pl, gi, gl, di, dl)
!==========================================================================

implicit none

real*8 t, pl, gi, gl, di, dl

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_ice_liq_done = isok

equi_ice_liq_t = t
equi_ice_liq_p = pl
equi_ice_liq_g_ice = gi
equi_ice_liq_g_liq = gl
equi_ice_liq_d_ice = di
equi_ice_liq_d_liq = dl

end subroutine

!==========================================================================
subroutine clear_ice_liq_state()
!==========================================================================

implicit none

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

equi_ice_liq_done = 0

end subroutine

end module ice_liq_4
