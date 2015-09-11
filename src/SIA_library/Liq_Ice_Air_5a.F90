module liq_ice_air_5a

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
!     AIR_3A,       FILE AIR_3A.F90
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
use air_3a
use air_3b
use air_3c
use ice_liq_4

use flu_3a

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

public :: &
          liq_ice_air_h_temperature_si, &
          liq_ice_air_h_density_si, &
          liq_ice_air_h_enthalpy_si, &
          liq_ice_air_h_gc_gv_si, &
          liq_ice_air_h_gd_gc_si, &
          liq_ice_air_h_gd_si, &
          liq_ice_air_h_gc_si, &
          liq_ice_air_h_gh_gc_si


!CONTROL PARAMETERS OF THE EQUI TEMPERATURE ITERATION
integer ctrl_initialized

integer ctrl_mode_liquidfraction
integer ctrl_loop_maximum

real*8 ctrl_eps_exit_liquidfraction

contains

function liq_ice_air_h_temperature_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE TEMPERATURE OF AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!T(A,ETA,P) ABSOLUTE TEMPERATURE OF MOIST AIR (LIQ-VAP-AIR OR ICE-VAP-AIR OR VAP-AIR) IN K

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

real*8 liq_ice_air_h_temperature_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t

liq_ice_air_h_temperature_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    liq_ice_air_h_temperature_si = air_temperature_si(wa_si, eta_si, p_si)
else
    !Try equilibrium without condensate:
    t = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t/=errorreturn) then
        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        if(t_freeze==errorreturn) return

        if(t<t_freeze) then
            icl = ice_air_icl_si(wa_si,t,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                liq_ice_air_h_temperature_si = t
                return
            endif
        else
            icl = liq_air_icl_si(wa_si,t,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                liq_ice_air_h_temperature_si = t
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        liq_ice_air_h_temperature_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        liq_ice_air_h_temperature_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        liq_ice_air_h_temperature_si = liq_ice_air_temperature_si()
    endif
endif
end function

function liq_ice_air_h_density_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE DENSITY OF AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!D(A,ETA,P) = 1/(DG/DP)  DENSITY OF MOIST (LIQ-VAP-AIR OR ICE-VAP-AIR OR VAP-AIR) AIR IN KG/M3

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_density_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a_si

liq_ice_air_h_density_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    liq_ice_air_h_density_si = air_density_si(wa_si, t_si, p_si)
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then
        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_density_si = air_density_si(wa_si, t_si, p_si)
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_density_si = air_density_si(wa_si, t_si, p_si)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_density_si = liq_air_g_density_si(wa_si, t_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_density_si = ice_air_g_density_si(wa_si, t_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        liq_ice_air_h_density_si = liq_ice_air_density_si()
    endif

    if(liq_ice_air_h_density_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        a_si=1.0
        liq_ice_air_h_density_si = air_g_density_si(a_si, t_si, p_si)
    endif
endif
end function

function liq_ice_air_h_enthalpy_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE ENTHALPY OF AIR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!H - ENTHALPY

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_enthalpy_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a_si

liq_ice_air_h_enthalpy_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    liq_ice_air_h_enthalpy_si = air_g_enthalpy_si(wa_si, t_si, p_si)
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then
        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_enthalpy_si = air_g_enthalpy_si(wa_si, t_si, p_si)
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_enthalpy_si = air_g_enthalpy_si(wa_si, t_si, p_si)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_enthalpy_si = liq_air_g_enthalpy_si(wa_si, t_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_enthalpy_si = ice_air_g_enthalpy_si(wa_si, t_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        liq_ice_air_h_enthalpy_si = liq_ice_air_enthalpy_si()
    endif

    if(liq_ice_air_h_enthalpy_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        a_si=1.0
        liq_ice_air_h_enthalpy_si = air_g_enthalpy_si(a_si, t_si, p_si)
    endif

endif
end function

function liq_ice_air_h_gc_gv_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE DIFFERENCE BETWEEN THE GIBBS FUNCTION OF THE CONDENSATE AND
!THE GIBBS FUNCTION OF WATER VAPOR
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!GC-GV - DIFFERENCE BETWEEN CONDENSATE AND VAPOR GIBBS

!INPUTS:
!WA_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_gc_gv_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a, gv

liq_ice_air_h_gc_gv_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)

    !To find condensate, Find freezing temperature:
    if(set_liq_ice_air_eq_at_p(p_si) /= errorreturn) then
        !Above or at triple point
        t_freeze=liq_ice_air_temperature_si()
        if(t_freeze==errorreturn) return
        if(t_si<t_freeze) then
            !Freezing
            a = ice_air_massfraction_air_si(t_si, p_si) 
            gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
            liq_ice_air_h_gc_gv_si = (ice_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
        else
            !Not freezing
            a = liq_air_massfraction_air_si(t_si, p_si) 
            gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
            liq_ice_air_h_gc_gv_si = (liq_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
        endif
    else
        !Below triple point: pure Ice
        a = ice_air_massfraction_air_si(t_si, p_si) 
        gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
        liq_ice_air_h_gc_gv_si = (ice_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
    endif
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then

        if(set_liq_ice_air_eq_at_p(p_si) /= errorreturn) then
            !Aboce triple point
            t_freeze=liq_ice_air_temperature_si()
            !If the freezing temperature returns an error this function will return an error:
            if(t_freeze==errorreturn) return

            if(t_si<t_freeze) then
                !Freezing conditions, use ice_air:
                icl = ice_air_icl_si(wa_si,t_si,p_si)
                if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                    !Air is below the isentropic condensation level -> no condensation
                    a = ice_air_massfraction_air_si(t_si, p_si) 
                    gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
                    liq_ice_air_h_gc_gv_si = (ice_g_si(0,0, t_si, p_si) - gv)
                    return
                endif
            else
                !Non-Freezing conditions, use liq_air:
                icl = liq_air_icl_si(wa_si,t_si,p_si)
                if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                    !Air is below the isentropic condensation level -> no condensation
                    a = liq_air_massfraction_air_si(t_si, p_si) 
                    gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
                    liq_ice_air_h_gc_gv_si = (liq_g_si(0,0, t_si, p_si) - gv)
                    return
                endif
            endif
        else
            !Below triple point:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                a = ice_air_massfraction_air_si(t_si, p_si) 
                gv = air_g_si(0,0,0,wa_si, t_si, p_si)-wa_si*air_g_si(1,0,0,wa_si, t_si, p_si)
                liq_ice_air_h_gc_gv_si = (ice_g_si(0,0, t_si, p_si) - gv)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        a = liq_air_massfraction_air_si(t_si, p_si) 
        gv = air_g_si(0,0,0,a, t_si, p_si)-a*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gc_gv_si = (liq_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        a = ice_air_massfraction_air_si(t_si, p_si) 
        gv = air_g_si(0,0,0,a, t_si, p_si)-a*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gc_gv_si = (ice_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        t_si = liq_ice_air_temperature_si()
        a = liq_ice_air_airfraction_si()
        gv = air_g_si(0,0,0,a, t_si, p_si)-a*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gc_gv_si = (wt*liq_g_si(0,0, t_si, p_si)+(1.0-wt)*ice_g_si(0,0, t_si, p_si) - gv)*(1.0/a-1.0)
    endif
endif
end function

function liq_ice_air_h_gc_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES 
!THE GIBBS FUNCTION OF THE CONDENSATE
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!GD-GC - DIFFERENCE BETWEEN CONDENSATE AND VAPOR GIBBS

!INPUTS:
!WA_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_gc_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a_si

liq_ice_air_h_gc_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)

    !To find condensate, Find freezing temperature:
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return
    if(t_si<t_freeze) then
        !Freezing
        liq_ice_air_h_gc_si = ice_g_si(0,0, t_si, p_si)
    else
        !Not freezing
        liq_ice_air_h_gc_si = liq_g_si(0,0, t_si, p_si)
    endif
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then

        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_gc_si = ice_g_si(0,0, t_si, p_si)
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_gc_si = liq_g_si(0,0, t_si, p_si)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_gc_si = liq_g_si(0,0, t_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_gc_si = ice_g_si(0,0, t_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        t_si = liq_ice_air_temperature_si()
        liq_ice_air_h_gc_si = (wt*liq_g_si(0,0, t_si, p_si)+(1.0-wt)*ice_g_si(0,0, t_si, p_si))
    endif
    if(liq_ice_air_h_gc_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        !No water and ice:
        t_si = air_temperature_si(a_si, eta_si, p_si)
        liq_ice_air_h_gc_si = ice_g_si(0,0, t_si, p_si)
    endif
endif
end function

function liq_ice_air_h_gd_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE CHEMICAL POTENTIAL OF DRY AIR 
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!GD-GC - DIFFERENCE BETWEEN CONDENSATE AND VAPOR GIBBS

!INPUTS:
!WA_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_gd_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a, gd, a_si

liq_ice_air_h_gd_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)

    !To find condensate, Find freezing temperature:
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return
    if(t_si<t_freeze) then
        !Freezing
        gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    else
        !Not freezing
        gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    endif
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then

        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
                liq_ice_air_h_gd_si = gd
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
                liq_ice_air_h_gd_si = gd
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        a = liq_air_massfraction_air_si(t_si, p_si) 
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        a = ice_air_massfraction_air_si(t_si, p_si) 
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        t_si = liq_ice_air_temperature_si()
        a = liq_ice_air_airfraction_si()
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    endif
    if(liq_ice_air_h_gd_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        !No water and ice:
        a_si=1.0
        t_si = air_temperature_si(a_si, eta_si, p_si)
        gd = air_g_si(0,0,0,a_si, t_si, p_si)
        liq_ice_air_h_gd_si = gd
    endif
endif
end function



function liq_ice_air_h_gd_gc_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE DIFFERENCE BETWEEN THE GIBBS FUNCTION OF THE DRY AIR AND
!THE GIBBS FUNCTION OF THE CONDENSATE
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!GD-GC - DIFFERENCE BETWEEN CONDENSATE AND VAPOR GIBBS

!INPUTS:
!WA_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_gd_gc_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a, gd, a_si

liq_ice_air_h_gd_gc_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)

    !To find condensate, Find freezing temperature:
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return
    if(t_si<t_freeze) then
        !Freezing
        gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - ice_g_si(0,0, t_si, p_si)
    else
        !Not freezing
        gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - liq_g_si(0,0, t_si, p_si)
    endif
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then

        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
                liq_ice_air_h_gd_gc_si = gd - ice_g_si(0,0, t_si, p_si)
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                gd = air_g_si(0,0,0,wa_si, t_si, p_si)+(1.0-wa_si)*air_g_si(1,0,0,wa_si, t_si, p_si)
                liq_ice_air_h_gd_gc_si = gd - liq_g_si(0,0, t_si, p_si)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        a = liq_air_massfraction_air_si(t_si, p_si) 
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - liq_g_si(0,0, t_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        a = ice_air_massfraction_air_si(t_si, p_si) 
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - ice_g_si(0,0, t_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        t_si = liq_ice_air_temperature_si()
        a = liq_ice_air_airfraction_si()
        gd = air_g_si(0,0,0,a, t_si, p_si)+(1.0-a)*air_g_si(1,0,0,a, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - (wt*liq_g_si(0,0, t_si, p_si)+(1.0-wt)*ice_g_si(0,0, t_si, p_si))
    endif
    if(liq_ice_air_h_gd_gc_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        !No water and ice:
        a_si=1.0
        t_si = air_temperature_si(a_si, eta_si, p_si)
        gd = air_g_si(0,0,0,a_si, t_si, p_si)
        liq_ice_air_h_gd_gc_si = gd - ice_g_si(0,0, t_si, p_si)
    endif
endif
end function

function liq_ice_air_h_gh_gc_si(wa_si, eta_si, p_si)
!THIS FUNCTION COMPUTES THE DERIVATIVE OF THE GIBBS FUNCTION WRT TO WA_SI
!SHOULD BE VALID OVER THE WHOLE REGION OF VALIDILITY OF TEOS-10

!OUTPUT:
!GH-GC - DIFFERENCE BETWEEN CONDENSATE AND VAPOR GIBBS

!INPUTS:
!WA_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!ETA_SI      ABSOLUTE IN-SITU ENTROPY
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!
!Note that this is a good template for converting liq_air_g* and ice_air_g* functions to
!liq_ice_air_h*. The liq_ice_air_h* functions are more robust at the triple point since it is
!possible to obtain an exact expression for the mixed liquid/ice state.

real*8 liq_ice_air_h_gh_gc_si, wa_si, eta_si, p_si
real*8 t_freeze, icl, wt, t_si, a_si

liq_ice_air_h_gh_gc_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

if(wa_si == 1d0) then
    !No water:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    liq_ice_air_h_gh_gc_si = air_g_si(1,0,0,wa_si, t_si, p_si)
else
    !Try equilibrium without condensate:
    t_si = air_temperature_si(wa_si, eta_si, p_si)
    !Check if no condensation makes sense:
    if(t_si/=errorreturn) then

        if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
        t_freeze=liq_ice_air_temperature_si()
        !If the freezing temperature returns an error this function will return an error:
        if(t_freeze==errorreturn) return

        if(t_si<t_freeze) then
            !Freezing conditions, use ice_air:
            icl = ice_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_gh_gc_si = air_g_si(1,0,0,wa_si, t_si, p_si)
                return
            endif
        else
            !Non-Freezing conditions, use ice_air:
            icl = liq_air_icl_si(wa_si,t_si,p_si)
            if(icl==errorreturn.or.(icl/=errorreturn.and.icl<p_si)) then
                !Air is below the isentropic condensation level -> no condensation
                liq_ice_air_h_gh_gc_si = air_g_si(1,0,0,wa_si, t_si, p_si)
                return
            endif
        endif
    endif

    !If we have not returned, must check for condensation and mixed state:
    wt = liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
    if(wt==1d0) then
        !Pure water condensate:
        t_si = liq_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_gh_gc_si = liq_air_g_si(1,0,0,wa_si, t_si, p_si)
    elseif(wt==0d0) then
        !Pure ice condensate:
        t_si = ice_air_h_temperature_si(wa_si, eta_si, p_si)
        liq_ice_air_h_gh_gc_si = ice_air_g_si(1,0,0,wa_si, t_si, p_si)
    else
        !Mixed state water/ice:
        if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,wt) == errorreturn) return
        t_si = liq_ice_air_temperature_si()
        liq_ice_air_h_gh_gc_si = wt*liq_air_g_si(1,0,0,wa_si, t_si, p_si)+(1.0-wt)*ice_air_g_si(1,0,0,wa_si, t_si, p_si)
    endif
    if(liq_ice_air_h_gh_gc_si==errorreturn) then
        !If it fails, it is likely that we are outside the region of validity of moist air.
        !For lack of better alternative, use dry air values:
        !No water:
        a_si=1.0
        t_si = air_temperature_si(a_si, eta_si, p_si)
        liq_ice_air_h_gh_gc_si = air_g_si(1,0,0,a_si, t_si, p_si)
    endif
endif
end function

function liq_ice_air_h_liquidfraction_of_condensate_si(wa_si,eta_si,p_si)
!THIS FUNCTION RETURNS THE FRACTION OF LIQUID WATER IN WET ICE AIR IN KG/KG
!WHEN WATER + ICE ARE AT EQUILIBRIUM WITH HUMID AIR.
!
!RETURNS     
!LIQUIDFRACTION=1D0        ALL THE CONDENSATE IS WATER
!LIQUIDFRACTION=0D0        ALL THE CONDENSATE IS ICE
!0D0<LIQUIDFRACTION<1D0    MIXTURE OF THE TWO STATES          
!
!IF THERE IS NO CONDENSATE IT RETURNS ERRORRETURN
real*8 liq_ice_air_h_liquidfraction_of_condensate_si, wa_si, eta_si, p_si
real*8 ifl, iml, t_freeze, t
real*8 eps
integer maxit

liq_ice_air_h_liquidfraction_of_condensate_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(p_si < 0d0) return

!Try mixed equilibrium using isentropic melting and freezing levels:
iml=liq_ice_air_iml_si(wa_si,eta_si)
ifl=liq_ice_air_ifl_si(wa_si,eta_si)
if(iml==errorreturn.or.ifl==errorreturn) then
    !One the isentropic levels did not converge.
    !Try ice and then liquid.
    if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) return
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return

    t = ice_air_h_temperature_si(wa_si, eta_si, p_si)
    if(t/=errorreturn.and.t<t_freeze) then
        liq_ice_air_h_liquidfraction_of_condensate_si=0d0 
        return
    endif

    t = liq_air_h_temperature_si(wa_si, eta_si, p_si)
    if(t/=errorreturn.and.t>=t_freeze) then
        liq_ice_air_h_liquidfraction_of_condensate_si=1d0 
        return
    endif
else
    !The isentropic levels converged.
    if(p_si>=ifl) then
        !Below freezing level -> pure liquid water
        liq_ice_air_h_liquidfraction_of_condensate_si=1d0 
    elseif(p_si<=iml) then
        !Above melting level -> pure ice
        liq_ice_air_h_liquidfraction_of_condensate_si=0d0 
    else !(p_si > iml .and. p_si < ifl)
        !Between the levels -> linear combination of the two
        !From comments in Liq_Ice_Air_4.F90:
        !set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,0d0)
        !iml=liq_ice_air_pressure_si()
        !liq_air_h_temperature_si(wa_si, eta_si, iml) = liq_ice_air_temperature()
        !set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,1d0)
        !ifl=liq_ice_air_pressure_si()
        !ice_air_h_temperature_si(wa_si, eta_si, ifl) = liq_ice_air_temperature()
        !
        !When p_si==iml --> liq_air_h_temperature
        !When p_si==ifl --> ice_air_h_temperature
        !
        call init_it_ctrl_liquidfraction
        if(ctrl_loop_maximum == 0) then
            maxit = 100
        elseif(ctrl_loop_maximum > 0) then
            maxit = ctrl_loop_maximum
        else
            return
        endif
        eps = ctrl_eps_exit_liquidfraction
        if(eps <= 0d0) return
        liq_ice_air_h_liquidfraction_of_condensate_si=&
                liquidfractioniteration_brent(wa_si,eta_si,p_si,0d0,1d0,maxit,eps)
    endif
endif
end function

!==========================================================================
function liquidfractioniteration_brent(wa_si, eta_si, p_si, wt1, wt2, &
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

real*8 liquidfractioniteration_brent, wa_si, eta_si, p_si, wt1, wt2, eps
real*8 a, b, c, d, s
real*8 fa, fb, fc, fs
integer it, maxit

logical mflag

liquidfractioniteration_brent = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0 .or. &
  p_si <= 0d0 .or. &
  eps == 0d0 .or. &
  maxit <= 0) return

a = wt1
if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,a) == errorreturn) return
fa = -liq_ice_air_pressure_si()
fa = -fa - p_si

b = wt2
if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,b) == errorreturn) return
fb = -liq_ice_air_pressure_si()
fb = -fb - p_si

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
    liquidfractioniteration_brent = b
    go to 999
  end if

  if(eps > 0d0) then                !ABSOLUTE LIMIT
    if(abs(a - b) < eps) then
      liquidfractioniteration_brent = b
      go to 999
    end if
  else                                 !RELATIVE LIMIT
    if(abs(a - b) < -eps * b) then
      liquidfractioniteration_brent = b
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

  if(set_liq_ice_air_eq_at_wa_eta_wt(wa_si,eta_si,s) == errorreturn) return
  fs = -liq_ice_air_pressure_si()
  fs = -fs - p_si

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
  if(b < 0d0 .or. b > 1d0) then
    liquidfractioniteration_brent = errorreturn
    return
  endif
endif

end function

!==========================================================================
subroutine init_it_ctrl_liquidfraction()
!=========================================================================

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR DENSITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_liquidfraction = 0           !default: minimum temperature is equal to absolute temperature
ctrl_eps_exit_liquidfraction = 0.0001d0  !default = 0.1 mk

end subroutine

!==========================================================================
subroutine set_it_ctrl_liq_ice_air_liquidfraction(key, value)
!=========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!EQUIVALENT TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0    SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!TOL_THETA       0           USE DEFAULT EXIT ACCURACY FOR EQUIVALENT TEMPERATURE (0.1 MK)
!TOL_THETA       EPS > 0  USE EPS AS EXIT ACCURACY FOR EQUIVALENT TEMPERATURE

character*13 key
real*8 value

call init_it_ctrl_liquidfraction

select case (trim(key))

  case ("it_steps");                 !ITERATION STEPS
      if(value == 0d0) then
         ctrl_loop_maximum = 100     !DEFAULT = 100
      elseif(value < 0d0) then
         ctrl_loop_maximum = -1
      else
         ctrl_loop_maximum = nint(value)
      endif

  case ("tol_theta");                !REQUIRED THETA TOLERANCE
      if(value == 0d0) then
         ctrl_eps_exit_liquidfraction = 0.0001d0 !DEFAULT = 0.1 MK
      elseif(value > 0d0) then
         ctrl_eps_exit_liquidfraction = value
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

end module liq_ice_air_5a
