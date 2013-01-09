module liq_ice_air_5

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

implicit none
private

character*16, private :: version = 'October 2012'

public :: liq_ice_air_pottemp_si, liq_ice_air_pottempequi_si, &
          liq_ice_air_pottempequisat_si, liq_ice_air_pottempequipseudo_si, &
          liq_ice_air_g_entropy_si, liq_ice_air_g_entropymoist_si, &
          liq_ice_air_g_temperatureequi_si, &
          liq_ice_air_h_temperature_si, liq_ice_air_massfraction_air_si

!CONTROL PARAMETERS OF THE EQUI TEMPERATURE ITERATION
integer ctrl_initialized

integer ctrl_mode_tempequi
integer ctrl_loop_maximum

real*8 ctrl_init_tempequi
real*8 ctrl_eps_exit_tempequi

!CONTROL PARAMETERS OF THE PSEUDO EQUIVALENT POTENTIAL TEMPERATURE ITERATION
integer ctrl_pottemppseudo_initialized

real*8 ctrl_delta_pres_pottemppseudo
real*8 ctrl_min_pres_pottemppseudo
real*8 ctrl_sat_pottemppseudo

contains

function liq_ice_air_pottemp_si(a_si, t_si, p_si, pr_si)
!THIS FUNCTION COMPUTES POTENTIAL TEMPERATURE AND SEAMLESSLY ACCOUNTS FOR ICE AND LIQUID WATER.
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!THETA(A,T,P,PR) ABSOLUTE POTENTIAL TEMPERATURE OF ICE AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

real*8 liq_ice_air_pottemp_si, a_si, t_si, p_si, pr_si
real*8 s

liq_ice_air_pottemp_si = errorreturn

if(a_si < 0d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return


if(p_si == pr_si) then
  liq_ice_air_pottemp_si = t_si
  return
end if

s=liq_ice_air_g_entropy_si(a_si,t_si,p_si)
if(s == errorreturn) return

liq_ice_air_pottemp_si = liq_ice_air_h_temperature_si(a_si, s, pr_si)

end function

function liq_ice_air_h_temperature_si(a_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE TEMPERATURE OF AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!T(A,ETA,P) ABSOLUTE TEMPERATURE OF MOIST AIR (LIQ-VAP-AIR OR ICE-VAP-AIR OR VAP-AIR) IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_h_temperature_si, a_si, eta_si, p_si
real*8 ifl, iml, t_freeze

liq_ice_air_h_temperature_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(p_si < 0d0) return

if(a_si == 1d0) then
    !No water:
    liq_ice_air_h_temperature_si = air_temperature_si(a_si, eta_si, p_si)
else
    !Try equilibrium without condensate:
    liq_ice_air_h_temperature_si = air_temperature_si(a_si, eta_si, p_si)
    if(liq_ice_air_h_temperature_si/=errorreturn) then
        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        if(t_freeze==errorreturn) return

        if(liq_ice_air_h_temperature_si<t_freeze) then
            if(ice_air_condensationpressure_si(a_si,liq_ice_air_h_temperature_si)>p_si) return
        else
            if(liq_air_condensationpressure_si(a_si,liq_ice_air_h_temperature_si)>p_si) return
        endif
    endif

    !Try mixed equilibrium using isentropic melting and freezing levels:
    iml=liq_ice_air_iml_si(a_si,eta_si)
    ifl=liq_ice_air_ifl_si(a_si,eta_si)
    if(iml==errorreturn.or.ifl==errorreturn) then
        !One the isentropic levels did not converge.
        !Try ice and then liquid.
        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        if(t_freeze==errorreturn) return

        liq_ice_air_h_temperature_si = ice_air_h_temperature_si(a_si, eta_si, p_si)
        if(liq_ice_air_h_temperature_si/=errorreturn.and.liq_ice_air_h_temperature_si<t_freeze) return

        liq_ice_air_h_temperature_si = liq_air_h_temperature_si(a_si, eta_si, p_si)
        if(liq_ice_air_h_temperature_si/=errorreturn.and.liq_ice_air_h_temperature_si>=t_freeze) return
    else
        !The isentropic levels converged.
        if(p_si>=ifl) then
            !Below freezing level -> pure liquid water
            liq_ice_air_h_temperature_si = liq_air_h_temperature_si(a_si, eta_si, p_si)
        elseif(p_si<=iml) then
            !Above melting level -> pure ice
            liq_ice_air_h_temperature_si = ice_air_h_temperature_si(a_si, eta_si, p_si)
        else !(p_si > iml .and. p_si < ifl)
            !Between the levels -> linear combination of the two
            !From comments in Liq_Ice_Air_4.F90:
            !set_liq_ice_air_eq_at_wa_eta_wt(a_si,eta_si,0d0)
            !iml=liq_ice_air_pressure_si()
            !liq_air_h_temperature_si(a_si, eta_si, iml) = liq_ice_air_temperature()
            !set_liq_ice_air_eq_at_wa_eta_wt(a_si,eta_si,1d0)
            !ifl=liq_ice_air_pressure_si()
            !ice_air_h_temperature_si(a_si, eta_si, ifl) = liq_ice_air_temperature()
            !
            !When p_si==iml --> liq_air_h_temperature
            !When p_si==ifl --> ice_air_h_temperature
            if(set_liq_ice_air_eq_at_wa_eta_wt(a_si,eta_si,(p_si-iml)/(ifl-iml)) == errorreturn) return
            liq_ice_air_h_temperature_si = liq_ice_air_temperature_si()
        endif
    endif
endif
end function

function liq_ice_air_massfraction_air_si(t_si, p_si)
!THIS FUNCTION COMPUTES THE MASSFRACTION OF SATURATED AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!A_SAT(T,P) MASSFRACTION OF SATURATED AIR (LIQ-VAP-AIR OR ICE-VAP-AIR OR VAP-AIR)

!INPUTS:
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_massfraction_air_si, t_si, p_si
real*8 t_freeze

liq_ice_air_massfraction_air_si = errorreturn

if(t_si < 0d0) return
if(p_si < 0d0) return

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
t_freeze=liq_ice_air_temperature_si()
if(t_freeze==errorreturn) return

if(t_si<t_freeze) then
    liq_ice_air_massfraction_air_si=ice_air_massfraction_air_si(t_si,p_si)
elseif(t_si>t_freeze) then
    liq_ice_air_massfraction_air_si=liq_air_massfraction_air_si(t_si,p_si)
else
    liq_ice_air_massfraction_air_si=liq_ice_air_airfraction_si()
endif

end function

function liq_ice_air_g_entropymoist_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE MOIST ENTROPY (EMANUEL 1994) OF WET AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!ETA_M(A,T,P) ABSOLUTE MOIST ENTROPY OF WET AIR (LIQ-VAP-AIR OR ICE-VAP-AIR OR VAP-AIR)

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_g_entropymoist_si, a_si, t_si, p_si
real*8 t_e

liq_ice_air_g_entropymoist_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

!Find the equivalent temperature
t_e=liq_ice_air_temperatureequi_si(a_si,t_si,p_si)

!Compute the entropy of the equivalent potential temperature:
liq_ice_air_g_entropymoist_si = liq_ice_air_g_entropy_si(a_si, t_e, p_si)

end function

function liq_ice_air_g_cond_entropy_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE CONDENSATION ENTROPY OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!ETA(A,T,P,PR) CONDENSATION ENTROPY OF ICE AIR OR LIQUID AIR 

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_cond_entropy_si, a_si, t_si, p_si
real*8 t_freeze

liq_ice_air_g_cond_entropy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
t_freeze=liq_ice_air_temperature_si()
if(t_freeze==errorreturn) return

if(t_si<t_freeze) then
    !condensation against ice
    liq_ice_air_g_cond_entropy_si=ice_air_g_cond_entropy_si(a_si, t_si, p_si,1d0)
elseif(t_si>t_freeze) then
    !condensation against liquid
    liq_ice_air_g_cond_entropy_si=liq_air_g_cond_entropy_si(a_si, t_si, p_si,1d0)
else
    !at freezing point: set entropy with only condensates, half in water, half in ice
    if(set_liq_ice_air_eq_at_wa_wl_wi(a_si,0.5*(1d0-a_si),0.5*(1-a_si)) == errorreturn) return
    liq_ice_air_g_cond_entropy_si=liq_ice_air_entropy_si()
endif

end function

function liq_ice_air_g_entropy_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE ENTROPY OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!ETA(A,T,P,PR) ENTROPY OF ICE AIR OR LIQUID AIR 

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_entropy_si, a_si, t_si, p_si
real*8 a_sat, t_freeze

liq_ice_air_g_entropy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

a_sat=liq_ice_air_massfraction_air_si(t_si,p_si)

if(a_sat<=a_si.and.a_sat/=errorreturn) then
    !no condensation
    liq_ice_air_g_entropy_si=air_g_entropy_si(a_si, t_si, p_si)
else
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return

    if(t_si<t_freeze) then
        if(ice_air_condensationpressure_si(a_si,t_si)>p_si) then
            !a_sat might have returned an error so we check if we are above
            !condensation level.
            liq_ice_air_g_entropy_si=air_g_entropy_si(a_si, t_si, p_si)
        else
            !condensation against ice
            liq_ice_air_g_entropy_si=ice_air_g_entropy_si(a_si, t_si, p_si)
        endif
    elseif(t_si>t_freeze) then
        if(liq_air_condensationpressure_si(a_si,t_si)>p_si) then
            !a_sat might have returned an error so we check if we are above
            !condensation level.
            liq_ice_air_g_entropy_si=air_g_entropy_si(a_si, t_si, p_si)
        else
            !condensation against liquid
            liq_ice_air_g_entropy_si=liq_air_g_entropy_si(a_si, t_si, p_si)
        endif
    else
        !at freeze point
        liq_ice_air_g_entropy_si=liq_ice_air_entropy_si()
    endif
endif

end function

function liq_ice_air_g_cond_cp_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE SPECIFIC HEAT OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!CP(A,T,P,PR) SPECIFIC HEAT OF ICE AIR OR LIQUID AIR 

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_cond_cp_si, a_si, t_si, p_si
real*8 t_freeze

liq_ice_air_g_cond_cp_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
t_freeze=liq_ice_air_temperature_si()
if(t_freeze==errorreturn) return

if(t_si<t_freeze) then
    !equilibrium against ice
    liq_ice_air_g_cond_cp_si=ice_air_g_cond_cp_si(a_si, t_si, p_si,1d0)
elseif(t_si>t_freeze) then
    liq_ice_air_g_cond_cp_si=liq_air_g_cond_cp_si(a_si, t_si, p_si,1d0)
else
    !at freeze point
    if(set_liq_ice_air_eq_at_wa_wl_wi(a_si,1d0-a_si,0d0) == errorreturn) return
    liq_ice_air_g_cond_cp_si=liq_ice_air_entropy_si()
endif

end function

function liq_ice_air_g_cp_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE SPECIFIC HEAT OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!CP(A,T,P,PR) SPECIFIC HEAT OF ICE AIR OR LIQUID AIR 

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_cp_si, a_si, t_si, p_si
real*8 a_sat, t_freeze

liq_ice_air_g_cp_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

a_sat=liq_ice_air_massfraction_air_si(t_si,p_si)
if(a_sat/=errorreturn.and.a_sat>a_si) then
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return

    if(t_si<t_freeze) then
        !equilibrium against ice
        liq_ice_air_g_cp_si=ice_air_g_cp_si(a_si, t_si, p_si)
    elseif(t_si>t_freeze) then
        liq_ice_air_g_cp_si=liq_air_g_cp_si(a_si, t_si, p_si)
    else
        !at freeze point
        liq_ice_air_g_cp_si=liq_ice_air_entropy_si()
    endif
else
    liq_ice_air_g_cp_si=air_g_cp_si(a_si, t_si, p_si)
endif

end function

function liq_ice_air_g_temperatureequi_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE EQUIVALENT TEMPERATURE OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!TEQ(A,T,P,PR) ABSOLUTE EQUIVALENT TEMPERATURE OF ICE AIR OR LIQUID AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_temperatureequi_si, a_si, t_si, p_si
real*8 t1, t2, eps
integer maxit

liq_ice_air_g_temperatureequi_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

if(a_si==1d0) then
    liq_ice_air_g_temperatureequi_si = t_si
    return
endif

call init_it_ctrl_tempequi

t1=min_liq_ice_air_temperatureequi_si(a_si,t_si,p_si)
t2=max_liq_ice_air_temperatureequi_si(a_si,t_si,p_si)

  if(ctrl_loop_maximum == 0) then
    maxit = 100
  elseif(ctrl_loop_maximum == (-1)) then
    liq_ice_air_g_temperatureequi_si = t1
    return
  elseif(ctrl_loop_maximum > 0) then
    maxit = ctrl_loop_maximum
  else
    return
  endif

eps = ctrl_eps_exit_tempequi
if(eps <= 0d0) return

!RUN THE ITERATION
liq_ice_air_g_temperatureequi_si = temperatureequiiteration_brent(a_si, t_si, p_si, t1, t2, maxit, eps)

end function

function liq_ice_air_pottempequi_si(a_si, t_si, p_si, pr_si)
!THIS FUNCTION COMPUTES THE EQUIVALENT POTENTIAL TEMPERATURE OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!THETA_E_STAR(A,T,P,PR) EQUIVALENT POTENTIAL TEMPERATURE OF ICE AIR OR LIQUID AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_pottempequi_si, a_si, t_si, p_si, pr_si

liq_ice_air_pottempequi_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return


liq_ice_air_pottempequi_si = liq_ice_air_pottemp_si(a_si,t_si,p_si,pr_si)

liq_ice_air_pottempequi_si = liq_ice_air_g_temperatureequi_si(a_si, &
                                                    liq_ice_air_pottempequi_si, pr_si)

end function

function liq_ice_air_pottempequipseudo_si(a_si, t_si, p_si, pr_si)
!THIS FUNCTION COMPUTES THE PSEUDO EQUIVALENT POTENTIAL TEMPERATURE OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!THETA_E_STAR(A,T,P,PR) PSEUDO EQUIVALENT POTENTIAL TEMPERATURE OF ICE AIR OR LIQUID AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_pottempequipseudo_si, a_si, t_si, p_si, pr_si
real*8 p_top, dp, t_prev, t_next, p_prev, p_next, a_sat, a_prev
integer i, nlev

liq_ice_air_pottempequipseudo_si=errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

!These parameters will influence the convergence of this algorithm
call init_it_ctrl_pottemppseudo
p_top=ctrl_min_pres_pottemppseudo
dp=ctrl_delta_pres_pottemppseudo

nlev = floor((pr_si-p_top)/dp)+1

if(pr_si<=p_top) return

!Move the parcel to the reference pressure:
t_prev=liq_ice_air_pottemp_si(a_si,t_si,p_si,pr_si)

!Remove all condensates
a_sat = liq_ice_air_massfraction_air_si(t_prev,p_si)
if(a_sat==errorreturn) return
if(a_sat>=a_si) then
    a_prev = a_sat
else
    a_prev = a_si
endif

p_prev=pr_si

dp = (pr_si-p_top)/nlev

do i=1,nlev
    p_next=p_prev-dp

    t_next=liq_ice_air_pottemp_si(a_prev,t_prev,p_prev,p_next)
    if(t_next==errorreturn.or.t_next<mix_air_tmin) exit

    a_sat = liq_ice_air_massfraction_air_si(t_next,p_next)
    if(a_sat==errorreturn) exit
    if(a_sat>=a_prev) a_prev = a_sat

    t_prev=t_next
    p_prev=p_next
enddo
a_sat=liq_ice_air_massfraction_air_si(t_prev,p_prev)
liq_ice_air_pottempequipseudo_si=liq_ice_air_pottemp_si(a_sat, t_prev, p_prev, pr_si)

end function

function liq_ice_air_pottempequisat_si(a_si, t_si, p_si, pr_si)
!THIS FUNCTION COMPUTES THE SATURATION EQUIVALENT POTENTIAL TEMPERATURE OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!THETA_E_STAR(A,T,P,PR) SATURATION EQUIVALENT POTENTIAL TEMPERATURE OF ICE AIR OR LIQUID AIR IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_pottempequisat_si, a_si, t_si, p_si, pr_si
real*8 a_sat

liq_ice_air_pottempequisat_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

!Find the saturation air mass fraction:
a_sat=liq_ice_air_massfraction_air_si(t_si,p_si)

if(a_si==1d0.or.a_si>=a_sat) then
    liq_ice_air_pottempequisat_si=liq_ice_air_pottempequi_si(a_sat, t_si, p_si, pr_si)
else
    liq_ice_air_pottempequisat_si=liq_ice_air_pottempequi_si(a_si, t_si, p_si, pr_si)
endif

end function
!=========================================================================
function min_liq_ice_air_temperatureequi_si(wa_si, t_si, p_si)
!==========================================================================
!THIS FUNCTION GIVES A LOWER BOUND ESTIMATE FOR THE EQUIVALENT TEMPERATURE OF WET ICE AIR 
!FROM ITS DRY-AIR FRACTION WA_SI IN KG/KG, ITS ABSOLUTE TEMPERATURE, T_SI IN K AND ITS
!PRESSURE P_SI IN PA

real*8 min_liq_ice_air_temperatureequi_si, wa_si, p_si, t_si

real*8 lv, t_temp, a, mixing_eta

min_liq_ice_air_temperatureequi_si = errorreturn

if(t_si <= 0d0) return
if(wa_si < 0d0 .or. wa_si > 1d0) return

if(set_liq_air_eq_at_t_p(tp_temperature_si, sealevel_pressure_si)==errorreturn) return
lv = liq_air_enthalpy_evap_si()

a = liq_ice_air_massfraction_air_si(t_si,p_si)
if(a/=errorreturn.and.a>=wa_si) then
    t_temp = t_si + (1d0-a)*lv/liq_ice_air_g_cond_cp_si(wa_si,t_si,p_si)
else
    t_temp = t_si + (1d0-wa_si)*lv/liq_ice_air_g_cond_cp_si(wa_si,t_si,p_si)
endif

if(p_si<(sealevel_pressure_si-tp_pressure_iapws95_si)) then 
    mixing_eta = -gas_constant_air_si/air_g_cp_si(1d0,t_si,p_si)*wa_si*log(p_si/(sealevel_pressure_si-tp_pressure_iapws95_si))
    t_temp = t_temp*exp(-mixing_eta)
endif

if(t_temp>=mix_air_tmax) then
    min_liq_ice_air_temperatureequi_si = mix_air_tmax
else
    min_liq_ice_air_temperatureequi_si = t_temp
endif

end function

!=========================================================================
function max_liq_ice_air_temperatureequi_si(wa_si, t_si, p_si)
!==========================================================================

!THIS FUNCTION GIVES AN UPPER BOUND ESTIMATE FOR THE EQUIVALENT TEMPERATURE OF WET ICE AIR 
!FROM ITS DRY-AIR FRACTION WA_SI IN KG/KG, ITS ABSOLUTE TEMPERATURE, T_SI IN K AND ITS
!PRESSURE P_SI IN PA

real*8 max_liq_ice_air_temperatureequi_si, wa_si, p_si, t_si

!real*8 tt, pt, cpa, cpi, sat, sit, ra, tmin
!real*8 p, t, eta

real*8 lv, t_freeze, t_temp, mixing_eta

max_liq_ice_air_temperatureequi_si = errorreturn

if(t_si <= 0d0) return
if(wa_si < 0d0 .or. wa_si > 1d0) return

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
t_freeze=liq_ice_air_temperature_si()
if(t_freeze==errorreturn) return

if(t_si < t_freeze) then
    if(set_ice_air_eq_at_t_p(t_si, p_si)==errorreturn) return
else
    if(set_ice_air_eq_at_t_p(tp_temperature_si, sealevel_pressure_si)==errorreturn) return
endif
lv = ice_air_enthalpy_subl_si()/liq_ice_air_g_cond_cp_si(1d0,t_si,p_si)

if(t_si<t_freeze) then
    t_temp = t_si*exp(lv*(1d0-wa_si)/t_si)
else
    t_temp = t_si*exp(lv*(1d0-wa_si)/t_freeze)
endif

if(p_si>(sealevel_pressure_si-tp_pressure_iapws95_si)) then 
    mixing_eta = -gas_constant_air_si/liq_ice_air_g_cond_cp_si(1d0,t_si,p_si)*&
                        wa_si*log(p_si/(sealevel_pressure_si-tp_pressure_iapws95_si))
    t_temp=t_temp*exp(-mixing_eta)
endif

if(t_temp>=mix_air_tmax) then
    max_liq_ice_air_temperatureequi_si = mix_air_tmax
else
    max_liq_ice_air_temperatureequi_si = t_temp
endif

end function

!==========================================================================
function temperatureequiiteration_brent(wa_si, t_si, p_si, t1, t2, &
                                    maxit, eps)
!=========================================================================

!THE FUNCTION COMPUTES THE EQUIVALENT TEMPERATURE AS A FUNCTION OF AIR FRACTION, ENTROPY AND PRESSURE,
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

real*8 temperatureequiiteration_brent, wa_si, t_si, p_si, t1, t2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
real*8 eta, mixing_eta
integer it, maxit

logical mflag

temperatureequiiteration_brent = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0 .or. &
  t1 <= 0d0 .or. &
  t2 <= 0d0 .or. &
  p_si <= 0d0 .or. &
  eps == 0d0 .or. &
  maxit <= 0) return

mixing_eta = -gas_constant_air_si*wa_si*log(p_si/(sealevel_pressure_si-tp_pressure_iapws95_si))

a = t1

eta = liq_ice_air_g_cond_entropy_si(wa_si, a, p_si)+mixing_eta
if(eta>liq_ice_air_g_entropy_si(wa_si, mix_air_tmax, p_si)) then
    fa = -mix_air_tmax
else
    fa = -liq_ice_air_h_temperature_si(wa_si,eta,p_si)
endif
fa = -fa - t_si

b = t2
eta = liq_ice_air_g_cond_entropy_si(wa_si, b, p_si)+mixing_eta
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
    temperatureequiiteration_brent = b
    go to 999
  end if

  if(eps > 0d0) then                !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      temperatureequiiteration_brent = b
      go to 999
    end if
  else                                 !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      temperatureequiiteration_brent = b
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

  eta = liq_ice_air_g_cond_entropy_si(wa_si, s, p_si) + mixing_eta
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
    temperatureequiiteration_brent = errorreturn
    return
  endif
  !AIR_LIMITS
  if(b < dry_air_tmin .or. b > dry_air_tmax) &
    temperatureequiiteration_brent = errorreturn
endif

end function

!==========================================================================
subroutine init_it_ctrl_pottemppseudo()
!=========================================================================

if(ctrl_pottemppseudo_initialized == -1) return

ctrl_pottemppseudo_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_min_pres_pottemppseudo = 1e4
ctrl_delta_pres_pottemppseudo = 1e3
ctrl_sat_pottemppseudo = 0d0
end subroutine

!==========================================================================
subroutine init_it_ctrl_tempequi()
!=========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_tempequi = 0           !default: minimum temperature is equal to absolute temperature
ctrl_init_tempequi = 273.15d0
ctrl_eps_exit_tempequi = 0.0001d0  !default = 0.1 mk

end subroutine

!==========================================================================
subroutine set_it_ctrl_liq_ice_air_tempequi(key, value)
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

call init_it_ctrl_tempequi

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
         ctrl_mode_tempequi = 0       !DEFAULT: TE=T_SI
      elseif(value > 0d0) then
         ctrl_mode_tempequi = 1
         ctrl_init_tempequi = nint(value)
      endif

  case ("tol_theta");                !REQUIRED THETA TOLERANCE
      if(value == 0d0) then
         ctrl_eps_exit_tempequi = 0.0001d0 !DEFAULT = 0.1 MK
      elseif(value > 0d0) then
         ctrl_eps_exit_tempequi = value
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

end module liq_ice_air_5
