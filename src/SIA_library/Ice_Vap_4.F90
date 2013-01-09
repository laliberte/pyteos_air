module ice_vap_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIUM PROPERTIES OF ICE WITH WATER VAPOUR
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
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     MATHS_0,      FILE MATHS_0.F90
!     FLU_1,        FILE FLU_1.F90
!     ICE_1,        FILE ICE_1.F90
!     FLU_2,        FILE FLU_2.F90
!     ICE_2,        FILE ICE_2.F90

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

public :: ice_vap_chempot_si, ice_vap_density_ice_si, &
          ice_vap_density_vap_si, ice_vap_enthalpy_ice_si, &
          ice_vap_enthalpy_subl_si, ice_vap_enthalpy_vap_si, &
          ice_vap_entropy_ice_si, ice_vap_entropy_subl_si, &
          ice_vap_entropy_vap_si, ice_vap_pressure_vap_si, &
          ice_vap_sublimationpressure_si, ice_vap_sublimationtemp_si, &
          ice_vap_temperature_si, ice_vap_volume_subl_si, &
          set_ice_vap_eq_at_p, set_ice_vap_eq_at_t, &
          set_it_ctrl_ice_vap

integer, private :: equi_ice_vap_done
real*8,  private :: equi_ice_vap_p_vap, equi_ice_vap_t, &
         equi_ice_vap_g_ice,  equi_ice_vap_g_vap, equi_ice_vap_d_ice,  &
         equi_ice_vap_d_vap
integer, private :: ctrl_mode_ice, ctrl_initialized, ctrl_mode_vapour, &
         ctrl_mode_temperature, ctrl_loop_maximum
real*8,  private :: ctrl_init_t, ctrl_eps_exit_p_vap, ctrl_init_d_vap

contains

!==========================================================================
function ice_vap_chempot_si()
!==========================================================================

!RETURNS THE CHEMICAL POTENTIAL IN J/KG OF WATER VAPOUR AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_chempot_si = -3895.26747392

!set_ice_vap_eq_at_p 100
!ice_vap_chempot_si = -26421.2820403

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_chempot_si

ice_vap_chempot_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_chempot_si = equi_ice_vap_g_vap

end function

!==========================================================================
function ice_vap_density_ice_si()
!==========================================================================

!RETURNS THE DENSITY IN KG/M3 OF ICE AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_density_ice_si = 917.170465733

!SET_ICE_VAP_EQUILIBRIUM_TO_P 100
!ICE_VAP_DENSITY_ICE_SI = 919.600269745

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_density_ice_si

ice_vap_density_ice_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_density_ice_si = equi_ice_vap_d_ice

end function

!==========================================================================
function ice_vap_density_vap_si()
!==========================================================================

!RETURNS THE DENSITY OF WATER VAPOUR IN KG/M3 AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_density_vap_si = 3.77406140772E-03

!set_ice_vap_eq_at_p 100
!ice_vap_density_vap_si = 8.57185487853E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_density_vap_si

ice_vap_density_vap_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_density_vap_si = equi_ice_vap_d_vap

end function

!==========================================================================
function ice_vap_enthalpy_ice_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY OF ICE IN J/KG AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ICE_VAP_ENTHALPY_ICE_SI = -340033.434649

!SET_ICE_VAP_EQUILIBRIUM_TO_P 100
!ICE_VAP_ENTHALPY_ICE_SI = -374576.247867

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_enthalpy_ice_si

ice_vap_enthalpy_ice_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_enthalpy_ice_si = ice_enthalpy_si(equi_ice_vap_t, equi_ice_vap_p_vap)

end function

!==========================================================================
function ice_vap_enthalpy_vap_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTHALPY OF WATER VAPOUR IN J/KG AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_enthalpy_vap_si = 2495132.21977

!set_ice_vap_eq_at_p 100
!ice_vap_enthalpy_vap_si = 2463525.19629

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_enthalpy_vap_si

ice_vap_enthalpy_vap_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_enthalpy_vap_si = flu_enthalpy_si(equi_ice_vap_t, equi_ice_vap_d_vap)

end function

!==========================================================================
function ice_vap_entropy_ice_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY OF ICE IN J/(KG K) AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_entropy_ice_si = -1244.95617472

!set_ice_vap_eq_at_p 100
!ice_vap_entropy_ice_si = -1377.09771247

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_entropy_ice_si

ice_vap_entropy_ice_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_entropy_ice_si = ice_entropy_si(equi_ice_vap_t, equi_ice_vap_p_vap)

end function

!==========================================================================
function ice_vap_entropy_vap_si()
!==========================================================================

!RETURNS THE SPECIFIC ENTROPY OF WATER VAPOUR IN J/(KG K) AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_entropy_vap_si = 9255.65736018

!set_ice_vap_eq_at_p 100
!ice_vap_entropy_vap_si = 9848.77406912

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_entropy_vap_si

ice_vap_entropy_vap_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_entropy_vap_si = flu_entropy_si(equi_ice_vap_t, equi_ice_vap_d_vap)

end function

!==========================================================================
function ice_vap_pressure_vap_si()
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE IN PA OF WATER VAPOUR AT THE VAPOUR-ICE EQUILIBRIUM 
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_pressure_vap_si = 470.059067981

!set_ice_vap_eq_at_p 100
!ice_vap_pressure_vap_si = 100

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_pressure_vap_si

ice_vap_pressure_vap_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_pressure_vap_si = equi_ice_vap_p_vap

end function

!==========================================================================
function ice_vap_sublimationpressure_si(t_si)
!==========================================================================

!RETURNS THE ABSOLUTE PRESSURE IN PA OF WATER VAPOUR AT THE VAPOUR-ICE EQUILIBRIUM
!AS A FUNCTION OF ABSOLUTE TEMPERATURE T_SI IN K

!CHECK VALUE WITH DEFAULT SETTINGS: ice_vap_sublimationpressure_si(270) = 470.059067981

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_sublimationpressure_si
real*8 t_si

ice_vap_sublimationpressure_si = errorreturn

if(set_ice_vap_eq_at_t(t_si) == errorreturn) return

ice_vap_sublimationpressure_si = equi_ice_vap_p_vap

end function

!==========================================================================
function ice_vap_sublimationtemp_si(p_si)
!==========================================================================

!RETURNS THE ABSOLUTE TEMPERATURE IN K AT THE VAPOUR-ICE EQUILIBRIUM
!AS A FUNCTION OF ABSOLUTE PRESSURE P_SI IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: ice_vap_sublimationtemp_si(100) = 252.817910215

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none
 
real*8 ice_vap_sublimationtemp_si
real*8 p_si

ice_vap_sublimationtemp_si = errorreturn

if(set_ice_vap_eq_at_p(p_si) == errorreturn) return

ice_vap_sublimationtemp_si = equi_ice_vap_t

end function

!==========================================================================
function ice_vap_temperature_si()
!==========================================================================

!RETURNS THE ABSOLUTE TEMPERATURE OF WATER VAPOUR IN K AT THE VAPOUR-ICE EQUILIBRIUM
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_temperature_si = 270

!set_ice_vap_eq_at_p 100
!ice_vap_temperature_si = 252.817910215

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_temperature_si

ice_vap_temperature_si = errorreturn

if(equi_ice_vap_done /= isok) return

ice_vap_temperature_si = equi_ice_vap_t

end function

!==========================================================================
function ice_vap_volume_subl_si()
!==========================================================================

!RETURNS THE SPECIFIC SUBLIMATION VOLUME IN M3/KG OF WATER
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_volume_subl_si = 264.965451558

!set_ice_vap_eq_at_p 100
!ice_vap_volume_subl_si = 1166.60755699

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_volume_subl_si

ice_vap_volume_subl_si = errorreturn

if(equi_ice_vap_done /= isok) return

if(equi_ice_vap_d_ice <= 0d0) return
if(equi_ice_vap_d_vap <= 0d0) return

ice_vap_volume_subl_si = 1d0 / equi_ice_vap_d_vap - 1d0 / equi_ice_vap_d_ice

end function

!==========================================================================
function ice_vap_entropy_subl_si()
!==========================================================================

!RETURNS THE SPECIFIC SUBLIMATION ENTROPY IN J/(KG K) OF WATER
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_entropy_subl_si = 10500.6135349

!set_ice_vap_eq_at_p 100
!ice_vap_entropy_subl_si = 11225.8717816

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_entropy_subl_si
real*8 si, sv

ice_vap_entropy_subl_si = errorreturn

if(equi_ice_vap_done /= isok) return

si = ice_entropy_si(equi_ice_vap_t, equi_ice_vap_p_vap)
if(si == errorreturn) return
sv = flu_entropy_si(equi_ice_vap_t, equi_ice_vap_d_vap)
if(sv == errorreturn) return

ice_vap_entropy_subl_si = sv - si

end function

!==========================================================================
function ice_vap_enthalpy_subl_si()
!==========================================================================

!RETURNS THE SPECIFIC SUBLIMATION ENTHALPY IN J/KG OF WATER
!AFTER EXECUTION OF SET_ICE_VAP_EQUILIBRIUM_TO_T OR SET_ICE_VAP_EQUILIBRIUM_TO_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!set_ice_vap_eq_at_t 270
!ice_vap_enthalpy_subl_si = 2835165.65442

!set_ice_vap_eq_at_p 100
!ice_vap_enthalpy_subl_si = 2838101.44416

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS IN THIS MODULE

implicit none

real*8 ice_vap_enthalpy_subl_si
real*8 hi, hv

ice_vap_enthalpy_subl_si = errorreturn

if(equi_ice_vap_done /= isok) return

hi = ice_enthalpy_si(equi_ice_vap_t, equi_ice_vap_p_vap)
if(hi == errorreturn) return
hv = flu_enthalpy_si(equi_ice_vap_t, equi_ice_vap_d_vap)
if(hv == errorreturn) return

ice_vap_enthalpy_subl_si = hv - hi

end function

!==========================================================================
function set_ice_vap_eq_at_p(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ICE-VAPOUR EQUILIBRIUM AT GIVEN PRESSURE
!BY 2-DIMENSIONAL NEWTON ITERATION

implicit none

integer maxit
real*8 set_ice_vap_eq_at_p
       
real*8 dv_triple, t_triple
real*8 t_si, p_si
real*8 pv
real*8 gi, gv
real*8 di, dv
real*8 eps

!TRIPLE POINT
dv_triple = tp_density_vap_iapws95_si
t_triple = tp_temperature_si

if(equi_ice_vap_done == isok .and.  &
   p_si == equi_ice_vap_p_vap) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_ice_vap_eq_at_p = isok
  return
end if

call clear_ice_vap_state !WASTE ANY PREVIOUS STATE

set_ice_vap_eq_at_p = errorreturn

if(p_si <= 0) return

call init_it_ctrl_ice_vap

!SET INITIAL TEMPERATURE GUESS
select case (ctrl_mode_temperature)
  case (0); t_si = aux_temperature_correlation(p_si)
            if(t_si == errorreturn) return
  case (-1); t_si = t_triple
  case (1);  t_si = ctrl_init_t
  case default; return
end select

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0);  dv = aux_vap_density_correlation(t_si)
  case (-1); dv = dv_triple
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then; maxit = 100
elseif(ctrl_loop_maximum == (-1))then 
   pv = flu_pressure_si(t_si, dv)
   gi = ice_g_si(0, 0, t_si, pv)
   gv = flu_gibbs_energy_si(t_si, dv)
   di = ice_g_si(0, 1, t_si, pv)
   if(di == errorreturn .or. di <= 0d0) return
   di = 1d0 / di
   call set_ice_vap_state(t_si, pv, gi, gv, di, dv)
   set_ice_vap_eq_at_p = isok
   return
elseif(ctrl_loop_maximum > 0)then; maxit = ctrl_loop_maximum
else; return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_vap_iteration_at_p(p_si, maxit, eps, t_si, di, dv, gi, gv, pv) == errorreturn) then
  return
end if

call set_ice_vap_state(t_si, pv, gi, gv, di, dv)

set_ice_vap_eq_at_p = isok

end function

!==========================================================================
function set_ice_vap_eq_at_t(t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ICE-VAPOUR EQUILIBRIUM AT GIVEN TEMPERATURE
!BY NEWTON ITERATION

implicit none

integer maxit
real*8 set_ice_vap_eq_at_t
real*8 t_triple, dv_triple
real*8 t_si, pv
real*8 gi, gv
real*8 di, dv
real*8 eps

!TRIPLE POINT
t_triple = tp_temperature_si
dv_triple = tp_density_vap_iapws95_si

if(equi_ice_vap_done == isok .and.  &
   t_si == equi_ice_vap_t) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_ice_vap_eq_at_t = isok
  return
end if

call clear_ice_vap_state !WASTE ANY PREVIOUS STATE

set_ice_vap_eq_at_t = errorreturn

if(t_si <= 0d0) return
if(t_si > t_triple) return

call init_it_ctrl_ice_vap

!SET INITIAL VAPOUR DENSITY GUESS
select case (ctrl_mode_vapour)
  case (0); dv = aux_vap_density_correlation(t_si)
  case (-1); dv = dv_triple
  case (1);  dv = ctrl_init_d_vap
  case default; return
end select

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then; maxit = 100
elseif(ctrl_loop_maximum == (-1))then 
   pv = flu_pressure_si(t_si, dv)
   gi = ice_g_si(0, 0, t_si, pv)
   gv = flu_gibbs_energy_si(t_si, dv)
   di = ice_g_si(0, 1, t_si, pv)
   if(di == errorreturn .or. di <= 0d0) return
   di = 1d0 / di
   call set_ice_vap_state(t_si, pv, gi, gv, di, dv)
   set_ice_vap_eq_at_t = isok
   return
elseif(ctrl_loop_maximum > 0)then; maxit = ctrl_loop_maximum
else; return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p_vap
if(eps == 0d0) return

!RUN ITERATION LOOP
if(ice_vap_iteration_at_t(t_si, maxit, eps, di, dv, gi, gv, pv) == errorreturn) then
  return
end if

call set_ice_vap_state(t_si, pv, gi, gv, di, dv)

set_ice_vap_eq_at_t = isok

end function

!==========================================================================
function ice_vap_iteration_at_p(p_si, maxit, eps, t_si, d_ice_si,  &
                                d_vap_si, g_ice_si, g_vap_si, p_vap_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-VAPOUR PHASE EQUILIBRIUM FROM EQUAL PRESSURES,
!TEMPERATURES AND CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN PRESSURE, P_SI,
!FROM INITIAL GUESSES FOR THE VAPOUR DENSITY, D_VAP_SI, AND THE TEMPERATURE, T_SI
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PRESSURE

!OUTPUT: 
!        ICE_VAP_ITERATION_AT_T = ISOK IF SUCCESSFUL
!        ICE_VAP_ITERATION_AT_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO FLU_F_SI OR TO ICE_G_SI HAS RETURNED AN ERROR
!        - DENSITY OR TEMPERATURE HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!     T_SI: ABSOLUTE TEMPERATURE IN K
! D_ICE_SI: ICE DENSITY IN KG/M3 AT THE ICE-VAPOUR EQUILIBRIUM
! D_VAP_SI: VAPOUR DENSITY IN KG/M3 AT THE ICE-VAPOUR EQUILIBRIUM
! G_ICE_SI: ICE CHEMICAL POTENTIAL IN J/KG AT THE ICE-VAPOUR EQUILIBRIUM
! G_VAP_SI: VAPOUR CHEMICAL POTENTIAL J/KG AT IN THE ICE-VAPOUR EQUILIBRIUM
! P_VAP_SI: VAPOUR PRESSURE IN PA AT THE ICE-VAPOUR EQUILIBRIUM

!INPUT: 
!       P_SI: ABSOLUTE PRESSURE IN PA
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE VAPOUR PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE VAPOUR PRESSURE TOLERANCE
!       T_SI: INITIAL GUESS OF ABSOLUTE TEMPERATURE IN K
!   D_VAP_SI: INITIAL GUESS OF VAPOUR DENSITY IN KG/M3

implicit none

integer it, maxit
real*8 p_si, t_si, eps, d_ice_si, d_vap_si, g_ice_si, g_vap_si, p_vap_si
real*8 ice_vap_iteration_at_p
real*8 vi, dv
real*8 pv, gv
real*8 gi, gi_t
real*8 t, tt
real*8 fv, fv_d, fv_dd
real*8 fv_t, fv_td
real*8 a(2, 2), b(2), x(2)
real*8 ddv, dt, pv_old

tt = tp_temperature_si !TRIPLE POINT IN K

ice_vap_iteration_at_p = errorreturn

if((t_si <= 0d0) .or. &
  (d_vap_si <= 0d0)) then
  t_si = errorreturn
  d_ice_si = errorreturn
  d_vap_si = errorreturn
  g_ice_si = errorreturn
  g_vap_si = errorreturn
  p_vap_si = errorreturn
  return
endif

if(check_limits == 1) then
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

dv = d_vap_si
t = t_si

do it = 0, maxit

  !DERIVATIVES OF THE GIBBS & HELMHOLTZ FUNCTIONS FOR 2D-NEWTON ITERATION
  !WRITE(*,*)' '
  fv = flu_f_si(0, 0, t, dv)
  !WRITE(*,*)'IN ICE_VAP_ITERATION_AT_P: CALLED FLU_F_SI(0, 0, T, DV) WITH T,DV=',T,DV,' AND GOT ',FV
  if(fv == errorreturn) go to 999
  
  fv_d = flu_f_si(0, 1, t, dv)
  if(fv_d == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t, dv)
  if(fv_dd == errorreturn) go to 999
  
  fv_t = flu_f_si(1, 0, t, dv)
  if(fv_t == errorreturn) go to 999
  fv_td = flu_f_si(1, 1, t, dv)
  if(fv_td == errorreturn) go to 999

  pv_old = pv  !KEEP PREVIOUS PV TO CHECK ACCURACY
  !VAPOUR PRESSURE AND ICE + VAPOUR CHEMICAL POTENTIALS
  pv = dv ** 2 * fv_d    !PRESSURE OF VAPOUR
  gv = fv + dv * fv_d    !CHEM. POT. OF VAPOUR
  
  gi = ice_g_si(0, 0, t, pv)
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, pv)
  if(gi_t == errorreturn) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or.  &
       (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      p_vap_si = pv
      g_ice_si = gi
      g_vap_si = gv
      vi = ice_g_si(0, 1, t, pv)
      if(vi == errorreturn .or. vi <= 0) go to 999
      d_ice_si = 1d0 / vi
      d_vap_si = dv
      t_si = t
      ice_vap_iteration_at_p = isok
      go to 999
    end if
  end if
  
  if(it == maxit) go to 999
  
  !COEFFICIENT MATRIX
  a(1, 1) = dv * (2d0 * fv_d + dv * fv_dd)
  a(1, 2) = dv ** 2d0 * fv_td
  a(2, 1) = (2d0 * fv_d + dv * fv_dd)
  a(2, 2) = -gi_t + fv_t + dv * fv_td

  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = p_si - pv
  b(2) = gi - gv

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR
  ddv = x(1)
  dt = x(2)
  
  !update density & temperature
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
  d_vap_si <= flu_dmin .or. d_vap_si > flu_dmax .or. &
  p_vap_si < 0d0) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    p_vap_si = errorreturn
    ice_vap_iteration_at_p = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  d_ice_si <= 0d0) then
    t_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    p_vap_si = errorreturn
    ice_vap_iteration_at_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function ice_vap_iteration_at_t(t_si, maxit, eps, d_ice_si, d_vap_si,  &
                                g_ice_si, g_vap_si, p_vap_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ICE-VAPOUR PHASE EQUILIBRIUM FROM EQUAL PRESSURES
!AND CHEMICAL POTENTIALS OF THE TWO PHASES AT GIVEN TEMPERATURE, T_SI, IN K
!FROM AN INITIAL GUESS FOR THE VAPOUR DENSITY, D_VAP_SI, IN KG/M3
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN VAPOUR PRESSURE

!OUTPUT: 
!        ICE_VAP_ITERATION = ISOK IF SUCCESSFULLY DONE
!        ICE_VAP_ITERATION = ERRORRETURN IS RETURNED IF
!           - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!           - THE FUNCTION CALL TO FLU_F_SI OR ICE_G_SI HAS RETURNED AN ERROR
!           - DENSITIES HAVE TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION

!        D_VAP_SI: VAPOUR DENSITY IN KG/M3 AT THE ICE-VAPOUR EQUILIBRIUM
!        D_ICE_SI: ICE DENSITY IN KG/M3 AT THE THE ICE-VAPOUR EQUILIBRIUM
!        G_VAP_SI: VAPOUR CHEMICAL POTENTIAL IN J/KG AT THE ICE-VAPOUR EQUILIBRIUM
!        G_ICE_SI: ICE CHEMICAL POTENTIAL IN J/KG AT THE ICE-VAPOUR EQUILIBRIUM
!        P_VAP_SI: VAPOUR PRESSURE IN PA AT THE ICE-VAPOUR EQUILIBRIUM

!INPUT: 
!        T_SI: ABSOLUTE TEMPERATURE IN K
!        MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF VAPOUR PRESSURE
!             EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE PRESSURE TOLERANCE
!        D_VAP_SI: INITIAL GUESS OF VAPOUR DENSITY IN KG/M3

implicit none

integer it, maxit
real*8 ice_vap_iteration_at_t
real*8 vi, dv
real*8 pv
real*8 gi, gv
real*8 fv, fv_d, fv_dd
real*8 ddv, pv_old
real*8 t_si, eps, d_ice_si, d_vap_si,  &
       g_ice_si, g_vap_si, p_vap_si

ice_vap_iteration_at_t = errorreturn

if(d_vap_si <= 0d0) return
if(t_si <= 0d0) return

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    d_vap_si = errorreturn
    d_ice_si = errorreturn
    g_vap_si = errorreturn
    g_ice_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    d_vap_si = errorreturn
    d_ice_si = errorreturn
    g_vap_si = errorreturn
    g_ice_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
else
  if(t_si <= 0d0) then
    d_vap_si = errorreturn
    d_ice_si = errorreturn
    g_vap_si = errorreturn
    g_ice_si = errorreturn
    p_vap_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

dv = d_vap_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION FOR NEWTON ITERATION
  fv = flu_f_si(0, 0, t_si, dv)
  if(fv == errorreturn) go to 999
  fv_d = flu_f_si(0, 1, t_si, dv)
  if(fv_d == errorreturn) go to 999
  fv_dd = flu_f_si(0, 2, t_si, dv)
  if(fv_dd == errorreturn) go to 999
  
  pv_old = pv            !KEEP PREVIOUS PV TO CHECK ACCURACY
  !VAPOUR AND ICE PRESSURE  + VAPOUR CHEMICAL POTENTIAL
  pv = dv ** 2 * fv_d    !PRESSURE OF VAPOUR
  gv = fv + dv * fv_d    !CHEM. POT. OF VAPOUR
  
  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t_si, pv)
  if(gi == errorreturn) go to 999
  vi = ice_g_si(0, 1, t_si, pv)
  if(vi == errorreturn) go to 999
  if(vi <= 0d0) go to 999
  
  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(pv - pv_old) < eps) .or.  &
       (eps < 0d0 .and. abs(pv - pv_old) < -eps * pv)) then
      ice_vap_iteration_at_t = isok
      p_vap_si = pv
      g_ice_si = gi
      g_vap_si = gv
      d_ice_si = 1d0 / vi
      d_vap_si = dv
      go to 999
    end if
  end if

  if(it == maxit) go to 999
  
  !ITERATION STEP
  ddv = (2d0 * fv_d + dv * fv_dd) * (vi * dv - 1d0)
  if(ddv == 0d0) go to 999
  ddv = (gv - gi) / ddv
  
  !UPDATE VAPOUR DENSITY
  dv = dv + ddv
  if(dv <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_vap_si <= flu_dmin .or. d_vap_si > flu_dmax .or. &
  p_vap_si <= 0d0) then
    p_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    ice_vap_iteration_at_t = errorreturn
    return
  endif
  !ICE_LIMITS
  if(d_ice_si  <= 0d0) then
    p_vap_si = errorreturn
    g_ice_si = errorreturn
    g_vap_si = errorreturn
    d_ice_si = errorreturn
    d_vap_si = errorreturn
    ice_vap_iteration_at_t = errorreturn
  endif
endif

end function

!==========================================================================
function aux_vap_density_correlation(t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE VAPOUR DENSITY OF SUBLIMATION FROM THE
!CLAUSIUS-CLAPEYRON EQUATION AS A FIRST GUESS FOR THE SUBLIMATION ITERATION

implicit none

real*8 aux_vap_density_correlation
real*8 rt, t_si
real*8 r, q, pt, tt

r = gas_constant_h2o_si      !SPECIFIC GAS CONSTANT OF H2O IN J/(KG K)
q = tp_enthalpy_vap_si - tp_enthalpy_ice_si ! = 2834359.44543354   !SUBLIMATION ENTHALPY IN J/KG
pt = tp_pressure_iapws95_si  !IAPWS-95 TRIPLE POINT PRESSURE IN PA
tt = tp_temperature_si

rt = r * t_si

aux_vap_density_correlation = (pt / rt) * exp((t_si / tt - 1d0) * q / rt)

end function

!==========================================================================
function aux_temperature_correlation(p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE SUBLIMATION TEMPERATURE COMPUTED FROM THE PRESSURE BY A
!CLAUSIUS-CLAPEYRON EQUATION AS A FIRST GUESS FOR THE SUBLIMATION ITERATION

implicit none

real*8 aux_temperature_correlation
real*8 t, p_si
real*8 r, q, pt, tt

r = gas_constant_h2o_si      !SPECIFIC GAS CONSTANT OF H2O IN J/(KG K)
q = tp_enthalpy_vap_si - tp_enthalpy_ice_si ! = 2834359.44543354   !SUBLIMATION ENTHALPY IN J/KG
pt = tp_pressure_iapws95_si  !IAPWS-95 TRIPLE POINT PRESSURE IN PA
tt = tp_temperature_si

aux_temperature_correlation = errorreturn

if(p_si <= 0d0) return

t = 1d0 / tt - (r / q) * log(p_si / pt)
if(t <= 0d0) return

aux_temperature_correlation = 1d0 / t

end function

!==========================================================================
subroutine init_it_ctrl_ice_vap()
!==========================================================================

implicit none

real*8 dv_triple, t_triple

!TRIPLE POINT
dv_triple = tp_density_vap_iapws95_si
t_triple = tp_temperature_si

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR THE ITERATION
ctrl_mode_ice = 0
ctrl_mode_vapour = 0
ctrl_mode_temperature = 0
ctrl_loop_maximum = 100
ctrl_init_d_vap = dv_triple
ctrl_init_t = t_triple
ctrl_eps_exit_p_vap = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

!==========================================================================
subroutine set_it_ctrl_ice_vap(key, value)
!==========================================================================

implicit none

character*13 key
real*8 value

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!ICE-VAPOUR EQUILIBRIUM

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VAPOUR DENSITY TO COMPUTE VAPOUR PRESSURE

!'INIT_VAP_DENS'        0           USE DEFAULT VAPOUR DENSITY TO START =
!                                   VAPOUR DENSITY CORRELATION FUNCTION, D = AUX_VAP_DENSITY_CORRELATION(T)
!'INIT_VAP_DENS'       -1           USE TRIPLE POINT VAPOUR DENSITY TO START
!'INIT_VAP_DENS'        D > 0       USE VALUE D AS VAPOUR DENSITY TO START

!'INIT_TEMP    '        0           USE DEFAULT TEMPERATURE TO START ( = AUX_SUBL_TEMPERATURE(P))
!'INIT_TEMP    '       -1           USE TRIPLE POINT TEMPERATURE TO START
!'INIT_TEMP    '        T > 0       USE VALUE T AS TEMPERATURE TO START

!'TOL_VAP_PRESS'        0           USE DEFAULT EXIT ACCURACY FOR VAPOUR PRESSURE (0.1 PPM)
!'TOL_VAP_PRESS'        EPS         USE EPS AS EXIT ACCURACY FOR VAPOUR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)


call init_it_ctrl_ice_vap

call clear_ice_vap_state

select case (trim(key))

  case ('it_steps')                                      !ITERATION STEPS
     if(value == 0d0)then;     ctrl_loop_maximum = 100   !DEFAULT = 100
     elseif(value < 0d0)then; ctrl_loop_maximum = -1
     else;   ctrl_loop_maximum = nint(value)
     endif

  case ('init_vap_dens')                                 !START VAPOUR DENSITY
     if(value == 0d0)then;       ctrl_mode_vapour = 0    !DEFAULT = AUX POLYNOMIAL
     elseif(value < -1d0)then;                          !IGNORE IT
     elseif(value < 0d0)then;  ctrl_mode_vapour = nint(value)
     else; ctrl_mode_vapour = 1; ctrl_init_d_vap = value
     endif

  case ('init_temp')                                     !START TEMPERATURE
     if(value == 0d0)then; ctrl_mode_temperature = 0     !DEFAULT = AUX POLYNOMIAL
     elseif(value < (-1d0))then;                        !IGNORE IT
     elseif(value < 0d0)then; ctrl_mode_temperature = nint(value)
     else; ctrl_mode_temperature = 1; ctrl_init_t = value
     endif

  case ('tol_vap_press')                                 !REQUIRED VAPOUR DENSITY TOLERANCE
     if(value == 0d0)then; ctrl_eps_exit_p_vap = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
     else; ctrl_eps_exit_p_vap = value
     endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine set_ice_vap_state(t, pv, gi, gv, di, dv)
!==========================================================================

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

real*8 t, pv, gi, gv, di, dv

equi_ice_vap_done = isok

equi_ice_vap_t = t
equi_ice_vap_p_vap = pv
equi_ice_vap_g_ice = gi
equi_ice_vap_g_vap = gv
equi_ice_vap_d_ice = di
equi_ice_vap_d_vap = dv

end subroutine

!==========================================================================
subroutine clear_ice_vap_state
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

equi_ice_vap_done = 0

end subroutine

end module ice_vap_4
