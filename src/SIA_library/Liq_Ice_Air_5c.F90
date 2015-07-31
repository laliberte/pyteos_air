module liq_ice_air_5c

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
use liq_ice_air_5b

implicit none
private

character*16, private :: version = 'October 2012'

public :: &
          liq_ice_air_g_cond_entropy_si,&
          liq_ice_air_g_wb_entropy_si,&
          liq_ice_air_g_cond_cp_si

contains

function liq_ice_air_g_wb_entropy_si(a_si, t_si, p_si)
!THIS FUNCTION COMPUTES THE CONDENSATION ENTROPY OF AIR AT FIXED PRESSURE AND FIXED DRY AIR MASS FRACTION.

!OUTPUT:
!ETA(A,T,P,PR) CONDENSATION ENTROPY OF ICE AIR OR LIQUID AIR 

!INPUTS:
!A_SI      ABSOLUTE DRY-AIR MASS FRACTION IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     REFERENCE PRESSURE IN PA

real*8 liq_ice_air_g_wb_entropy_si, a_si, t_si, p_si
real*8 t_freeze, a_sat

liq_ice_air_g_wb_entropy_si = errorreturn

if(a_si < 0d0 .or. a_si > 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return

a_sat=liq_ice_air_massfraction_air_si(t_si,p_si)
if(a_sat==errorreturn) return

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) then
    !Above or at triple point
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return

    if(t_si<t_freeze) then
        !condensation against ice
        liq_ice_air_g_wb_entropy_si=ice_air_g_cond_entropy_si(a_sat, t_si, p_si,a_si)
    elseif(t_si>t_freeze) then
        !condensation against liquid
        liq_ice_air_g_wb_entropy_si=liq_air_g_cond_entropy_si(a_sat, t_si, p_si,a_si)
    else
        liq_ice_air_g_wb_entropy_si=liq_ice_air_entropy_si()
    endif
else
    !Below triple point
    liq_ice_air_g_wb_entropy_si=ice_air_g_cond_entropy_si(a_sat, t_si, p_si,a_si)
endif

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

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) then
    !Above or at triple point
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
        liq_ice_air_g_cond_entropy_si=liq_ice_air_entropy_si()
    endif
else
    !Below triple point
    !condensation against ice
    liq_ice_air_g_cond_entropy_si=ice_air_g_cond_entropy_si(a_si, t_si, p_si,1d0)
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

if(set_liq_ice_air_eq_at_p(p_si) == errorreturn) then
    !Above or at triple point
    t_freeze=liq_ice_air_temperature_si()
    if(t_freeze==errorreturn) return

    if(t_si<t_freeze) then
        !equilibrium against ice
        liq_ice_air_g_cond_cp_si=ice_air_g_cond_cp_si(a_si, t_si, p_si,1d0)
    elseif(t_si>=t_freeze) then
    !elseif(t_si>t_freeze) then
        liq_ice_air_g_cond_cp_si=liq_air_g_cond_cp_si(a_si, t_si, p_si,1d0)
    !else
    !    !at freeze point
    !    wt=0.5d0
    !    if(set_liq_ice_air_eq_at_wa_wl_wi(a_si,wt*(1d0-a_si),(1d0-wt)*(1-a_si)) == errorreturn) return
    !    liq_ice_air_g_cond_cp_si=liq_ice_air_entropy_si()
    endif
else
    !Below triple point
        !equilibrium against ice
        liq_ice_air_g_cond_cp_si=ice_air_g_cond_cp_si(a_si, t_si, p_si,1d0)
endif

end function

end module liq_ice_air_5c
