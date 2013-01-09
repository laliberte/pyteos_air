module convert_0

!#########################################################################

!THIS MODULE IMPLEMENTS CONVERSION BETWEEN MASS FRACTIONS AND MOLE FRACTIONS
!AND PROVIDES A CONVERSION FROM PRACTICAL SALINITY TO ABSOLUTE SALINITY

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

!THIS MODULE REQUIRES THE LIBRARY MODULES
!     CONSTANTS_0, FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '28th June 2011'

public :: air_massfraction_air_si, air_massfraction_vap_si, &
          air_molar_mass_si, air_molfraction_air_si, &
          air_molfraction_vap_si, sa_from_sp, sp_from_sa

real*8 ::  mw = molar_mass_h2o_si      !MOLAR MASS OF H2O IN KG/MOL
real*8 ::  ma = molar_mass_air_L2000   !MOLAR MASS OF AIR IN KG/MOL USED BY LEMMON et al. 2000
integer :: flag_saar



contains

!==========================================================================
function air_molar_mass_si(a_si)
!==========================================================================

!OUTPUT:  
!air_molar_mass_si = MOLAR MASS OF HUMID AIR IN KG MOL-1

!INPUT:  
!A_SI = AIR MASS FRACTION, IN KG/KG

!CHECK VALUE: air_molar_mass_si(0.5) = 2.22122197774E-02
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: air_molar_mass_si, a_si

air_molar_mass_si = errorreturn

if (a_si < 0d0 .or. a_si > 1d0) return

air_molar_mass_si = 1d0 / ((1d0 - a_si) / mw + a_si / ma)

return
end function

!==========================================================================
function air_molfraction_vap_si(a_si)
!==========================================================================

!OUTPUT: 
!MOLE FRACTION OF VAPOUR IN HUMID AIR IN MOL MOL-1

!INPUT: 
!A_SI = AIR MASS FRACTION, IN KG/KG

!CHECK VALUE: air_molfraction_vap_si(0.5) = 0.616483190186

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: air_molfraction_vap_si, a_si

air_molfraction_vap_si = errorreturn

if (a_si < 0d0 .or. a_si > 1d0) return

air_molfraction_vap_si = (1d0 - a_si) / (1d0 - a_si * (1d0 - mw / ma))

return
end function

!==========================================================================
function air_molfraction_air_si(a_si)
!==========================================================================

!OUTPUT: 
!MOLE FRACTION OF DRY AIR IN HUMID AIR IN MOL MOL-1

!INPUT: 
!A_SI = AIR MASS FRACTION, IN KG/KG

!CHECK VALUE: air_molfraction_air_si(0.5) = 0.383516809814

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: air_molfraction_air_si, a_si

air_molfraction_air_si = errorreturn

if (a_si < 0d0 .or. a_si > 1d0) return

air_molfraction_air_si = (a_si * mw / ma) / (1d0 - a_si * (1d0 - mw / ma))

return
end function

!==========================================================================
function air_massfraction_air_si(x_si)
!==========================================================================

!OUTPUT: 
!A_SI = MASS FRACTION OF AIR IN HUMID AIR IN KG/KG

!INPUT: 
!X_SI = MOLE FRACTION OF AIR IN HUMID AIR IN MOL/MOL

!CHECK VALUE: air_massfraction_air_si(0.5) = 0.616483190186

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14):: air_massfraction_air_si, x_si

air_massfraction_air_si = errorreturn

if (x_si < 0d0 .or. x_si > 1d0) return

air_massfraction_air_si = x_si / (x_si + (1d0 - x_si) * mw / ma)

return
end function

!==========================================================================
function air_massfraction_vap_si(x_si)
!==========================================================================

!OUTPUT: 
!AIR_MASSFRACTION_VAP_SI = MASS FRACTION OF VAPOUR IN HUMID AIR IN KG/KG

!INPUT: 
!X_SI = DRY AIR MOLE FRACTION, IN MOL MOL-1

!CHECK VALUE: air_massfraction_vap_si(0.5) = 0.383516809814

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14):: air_massfraction_vap_si, x_si

air_massfraction_vap_si = errorreturn

if (x_si < 0d0 .or. x_si > 1d0) return

air_massfraction_vap_si = (1d0 - x_si) / (1d0 - x_si * (1d0 - ma / mw))

return
end function


!=====================================================
function sa_from_sp(sp_si,lon0,lat0,p_si)
!=====================================================

! CONVERT PRACTICAL SALINITY TO ABSOLUTE SALINITY 
!
! SP                  : PRACTICAL SALINITY                 [PSU]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! RESULT              : ABSOLUTE SALINITY                  [KG/KG]
!
!CHECK VALUE:
!ASAL_FROM_PSAL(0.0357, 201, -21, 101325+1023d4) = 0.035873322343341172
!ASAL_FROM_PSAL(0.035, 180, 40, 101325 + 2d7) =   0.035189093288995807  !NORTH PACIFIC SILICATE
!ASAL_FROM_PSAL(0.008, 20, 57, 101325) =     0.0081048377142857131      !CENTRAL BALTIC SEA

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa0, sp_si, sp0, lon0, lat0, p_si, p0
real (r14) :: sa_from_sp

p0 = (p_si - 101325d0)/1d4
sp0 = sp_si*1d3
sa0 = gsw_sa_from_sp(sp0,p0,lon0,lat0)
sa_from_sp = sa0*1d-3

return
end function


!=====================================================
function sp_from_sa(p_si,sa_si,lon0,lat0)
!=====================================================

! CONVERT ABSOLUTE SALINITY TO PRACTICAL SALINITY 
!
! SA_SI               : ABSOLUTE SALINITY                  [KG/KG]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! PSAL_FROM_ASAL  : PRACTICAL SALINITY                     [PSU]

!CHECK VALUE
!PSAL_FROM_ASAL(0.035873322343341172, 201, -21, 101325+1023d4) =   0.035700000000000010

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp0, sa_si, sa0, lon0, lat0, p_si, p0
real (r14) :: sp_from_sa

p0 = (p_si - 101325d0)/1d4
sa0 = sa_si*1d3
sp0 = gsw_sp_from_sa(sa0,p0,lon0,lat0)
sp_from_sa = sp0*1d-3

return
end function

!=====================================================
function sstar_from_sp(sp_si,p_si,lon0,lat0)
!=====================================================

! CONVERT PRACTICAL SALINITY TO PREFORMED SALINITY 
!
! SP_SI               : PRACTICAL SALINITY                 [PSU]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! SSTAR_FROM_SP       : PREFORMED SALINITY                 [PSU]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp_si, sp0, lon0, lat0, p_si, p0, sstar0, sstar_from_sp

p0 = (p_si - 101325d0)/1d4
sp0 = sp_si*1d3
sstar0 = gsw_sstar_from_sp(sp0,p0,lon0,lat0)
sstar_from_sp = sstar0*1d-3

return
end function

!=====================================================
function deltasa_from_sp(sp_si,p_si,lon0,lat0)
!=====================================================

! CONVERT PRACTICAL SALINITY TO ABSOLUTE SALINITY ANOMALY 
!
! SP_SI               : PRACTICAL SALINITY                 [PSU]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! DELTASA_FROM_SP     : ABSOLUTE SALINITY ANOMALY          [KG/KG]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp_si, sp0, lon0, lat0, p_si, p0, deltasa0, deltasa_from_sp

p0 = (p_si - 101325d0)/1d4
sp0 = sp_si*1d3
deltasa0 = gsw_deltasa_from_sp(sp0,p0,lon0,lat0)
deltasa_from_sp = deltasa0*1d-3

return
end function

!=====================================================
function sr_from_sp(sp_si)
!=====================================================

! CONVERT PRACTICAL SALINITY TO ABSOLUTE SALINITY ANOMALY 
!
! SP_SI               : PRACTICAL SALINITY                 [PSU]
!
! SR_SI               : PREFORMED SALINITY                 [KG/KG]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp_si, sp0, sr0, sr_from_sp

sp0 = sp_si*1d3
sr0 = gsw_sr_from_sp(sp0)
sr_from_sp = sr0*1d-3

return
end function

!=====================================================
function sstar_from_sa(sa_si,p_si,lon0,lat0)
!=====================================================

! CONVERT ABSOLUTE SALINITY TO PREFORMED SALINITY 
!
! SA_SI               : ABSOLUTE SALINITY                  [KG/KG]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! SSTAR_FROM_SA       : PREFORMED SALINITY                 [KG/KG]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa_si, sa0, lon0, lat0, p_si, p0, sstar0, sstar_from_sa

p0 = (p_si - 101325d0)/1d4
sa0 = sa_si*1d3
sstar0 = gsw_sstar_from_sa(sa0,p0,lon0,lat0)
sstar_from_sa = sstar0*1d-3

return
end function


!=====================================================
function sp_from_sstar(sstar_si,p_si,lon0,lat0)
!=====================================================

! CONVERT PREFORMED SALINITY TO PRACTICAL SALINITY
!
! SP_SI               : PREFORMED SALINITY                 [PSU]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! SP_FROM_SSTAR       : PRACTICAL SALINITY                 [PSU]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sstar_si, sp0, lon0, lat0, p_si, p0, sstar0, sp_from_sstar

p0 = (p_si - 101325d0)/1d4
sstar0 = sstar_si*1d3
sp0 = gsw_sp_from_sstar(sstar0,p0,lon0,lat0)
sp_from_sstar = sp0*1d-3

return
end function

!=====================================================
function sa_from_sstar(sstar_si,p_si,lon0,lat0)
!=====================================================

! CONVERT PREFORMED SALINITY TO PRACTICAL SALINITY
!
! SP_SI               : PREFORMED SALINITY                 [PSU]
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! SP_FROM_SSTAR       : ABSOLUTE SALINITY                  [KG/KG]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sstar_si, sa0, lon0, lat0, p_si, p0, sstar0, sa_from_sstar

p0 = (p_si - 101325d0)/1d4
sstar0 = sstar_si*1d3
sa0 = gsw_sa_from_sstar(sstar0,p0,lon0,lat0)
sa_from_sstar = sa0*1d-3

return
end function

!=====================================================
function saar(p_si,lon0,lat0)
!=====================================================

! CALCULATE THE ABSOLUTE SALINITY ANOMALY RATIO
!
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! SAAR                : ABSOLUTE SALINITY ANOMALY RATIO    [UNITLESS]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: lon0, lat0, p_si, p0, saar0, saar

p0 = (p_si - 101325d0)/1d4
saar0 = gsw_saar(p0,lon0,lat0)
saar = saar0*1d-3

return
end function

!=====================================================
function delta_sa_ref(p_si,lon0,lat0)
!=====================================================

! CALCULATE THE ABSOLUTE SALINITY ANOMALY REFERENCE VALUE
!
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! DELTA_SA_REF        : ABSOLUTE SALINITY ANOMALY REF. VALUE  [KG/KG]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: lon0, lat0, p_si, p0, delta_sa_ref0, delta_sa_ref

p0 = (p_si - 101325d0)/1d4
delta_sa_ref0 = gsw_delta_sa_ref(p0,lon0,lat0)
delta_sa_ref = delta_sa_ref0*1d-3

return
end function


!=====================================================
function fdelta(p_si,lon0,lat0)
!=====================================================

! CALCULATE THE RATIO OF ABSOLUTE SALINITY TO PREFORMED SALINITY MINUS 1
!
! P_SI                : ABSOLUTE PRESSURE                  [PA]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! DELTA_SA_REF        : RATIO OF SA TO SSTAR MINUS 1       [UNITLESS]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: lon0, lat0, p_si, p0, Fdelta0, Fdelta

p0 = (p_si - 101325d0)/1d4
fdelta0 = gsw_fdelta(p0,lon0,lat0)
fdelta = fdelta0*1d-3

return
end function


!==========================================================================
!==========================================================================
!==========================================================================
! THIS IS THE GSW ABSOLUTE SALINITY ALGORITHM IMPLEMENTED IN FORTRAN
!==========================================================================
!==========================================================================
!==========================================================================
!
! AS DEFINED IN 
!
! MCDOUGALL, T.J., JACKETT, D.R. AND MILLERO, F.J., 2009: AN ALGORITHM FOR ESTIMATING  
!            ABSOLUTE SALINITY IN THE GLOBAL OCEAN, OCEAN SCIENCE, SUBMITTED.  
!     
! 	
! DAVID JACKETT
! JANUARY 2009
!

!==========================================================================
! Gibbs SeaWater (GSW) Oceanographic Toolbox of TEOS-10 version 3.0 (Fortran)
!==========================================================================
!
! This is the salinity functions contained in the Gibbs SeaWater (GSW) 
! Oceanographic Toolbox of TEOS-10 (version 3.0).
! 
!
! salinity conversions
! gsw_sa_from_sp          - Absolute Salinity from Practical Salinity
! gsw_sstar_from_sp       - Preformed Salinity from Practical Salinity
!
! gsw_deltasa_from_sp     - Absolute Salinity Anomaly from Practical Salinity
! gsw_sr_from_sp          - Reference Salinity from Practical Salinity
! gsw_sp_from_sr          - Practical Salinity from Reference Salinity
! gsw_sp_from_sa          - Practical Salinity from Absolute Salinity
! gsw_sstar_from_sa       - Preformed Salinity from Absolute Salinity
! gsw_sp_from_sstar       - Practical Salinity from Preformed Salinity
! gsw_sa_from_sstar       - Absolute Salinity from Preformed Salinity
!
! Library functions of the GSW toolbox
! gsw_saar                - Absolute Salinity Anomaly Ratio (excluding the Baltic Sea)
! gsw_delta_sa_ref        - Absolute Salinity Anomaly ref. value (excluding the Baltic Sea)
! gsw_fdelta              - ratio of Absolute to Preformed Salinity, minus 1
! gsw_sa_from_sp_baltic   - Absolute Salinity Anomaly from Practical Salinity in the Baltic Sea
! gsw_sp_from_sa_baltic   - Practical Salinity from Absolute Salinity in the Baltic Sea
!
!
! Version 1.0 written by David Jackett
! Modified by Paul Barker (version 3.0)
!
! For help with this Oceanographic Toolbox email:- help_gsw@csiro.au
!
! This software is available from http://www.teos-10.org
!
!==========================================================================

!--------------------------------------------------------------------------
! salinity conversions
!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sp(sp,p,long,lat)
!==========================================================================

! Calculates Absolute Salinity, SA, from Practical Salinity, SP
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sa_from_sp   : Absolute Salinity                     [g/kg]

implicit none
integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p, gsw_sa_from_sp, saar
real (r14) :: gsw_sa_baltic

saar = gsw_saar(p,long,lat)

gsw_sa_from_sp = (35.16504d0/35.d0)*sp*(1.d0 + saar)

gsw_sa_baltic = gsw_sa_from_sp_baltic(sp,long,lat)

if (gsw_sa_baltic.lt.1d10) then
   gsw_sa_from_sp = gsw_sa_baltic
end if

if (saar.eq.9d15) then
   gsw_sa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sstar_from_sp(sp,p,long,lat)
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sp  : Preformed salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p
real (r14) :: saar, gsw_sstar_from_sp, sstar_baltic

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sp = (35.16504d0/35.d0)*sp*(1 - 0.35d0*saar);

!In the Baltic Sea, Sstar = sa.
sstar_baltic = gsw_sa_from_sp_baltic(sp,long,lat);

if (sstar_baltic.lt.1d10) then
    gsw_sstar_from_sp = sstar_baltic;
end if

if (saar.eq.9d15) then
    gsw_sstar_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_deltasa_from_sp(sp,p,long,lat)
!==========================================================================

! Calculates Absolute Salinity Anomaly, deltaSA, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_deltasa_from_sp : Absolute Salinty Anomaly           [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, long, lat, p
real (r14) :: gsw_deltasa_from_sp

gsw_deltasa_from_sp = gsw_sa_from_sp(sp,p,long,lat) - gsw_sr_from_sp(sp)

if (gsw_deltasa_from_sp.gt.1d10) then
    gsw_deltasa_from_sp = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sr_from_sp(sp)
!==========================================================================

! Calculates Reference Salinity, SR, from Practical Salinity, SP. 
!
! sp     : Practical Salinity                              [unitless]
!
! gsw_sr_from_sp : Reference Salinity                      [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sp, gsw_sr_from_sp

gsw_sr_from_sp = 1.004715428571429*sp;

if (gsw_sr_from_sp.ge.1.d10) then
    gsw_sr_from_sp = 9.d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sr(sr)
!==========================================================================

! Calculates Practical Salinity, sp, from Reference Salinity, SR. 
!
! SR     : Reference Salinity                              [g/kg]
!
! gsw_sp_from_sr  : Practical Salinity                     [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sr, gsw_sp_from_sr

gsw_sp_from_sr = 0.995306702338459*sr;

if (gsw_sp_from_sr.gt.1d10) then
    gsw_sp_from_sr = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sa(sa,p,long,lat)
!==========================================================================

! Calculates Practical salinity, sp, from Absolute salinity, sa  
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [DEG E]     
! lat    : latitude                                        [DEG N]
!
! gsw_sp_from_sa      : Practical Salinity                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p, gsw_sp_from_sa, saar
real (r14) :: gsw_sp_baltic

saar = gsw_saar(p,long,lat)

gsw_sp_from_sa = (35.d0/35.16504d0)*sa/(1d0 + saar)

gsw_sp_baltic = gsw_sp_from_sa_baltic(sa,long,lat);

if (gsw_sp_baltic.lt.1d10) then
   gsw_sp_from_sa = gsw_sp_baltic
end if

if (saar.eq.9d15) then
   gsw_sp_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sstar_from_sa(sa,p,long,lat)
!==========================================================================

! Calculates Preformed Salinity, Sstar, from Absolute Salinity, SA. 
!
! sa     : Absolute Salinity                               [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sstar_from_sa : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: sa, long, lat, p
real (r14) :: saar, gsw_sstar_from_sa

saar = gsw_saar(p,long,lat)

gsw_sstar_from_sa = sa*(1d0 - 0.35d0*saar)/(1d0 + saar)

! In the Baltic Sea, Sstar = sa, and note that gsw_saar returns zero
! for saar in the Baltic.

if (saar.eq.9d15) then
    gsw_sstar_from_sa = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sstar(sstar,p,long,lat)
!==========================================================================

! Calculates Absolute Salinity, SA, from Preformed Salinity, Sstar.
!
! Sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sa_from_sstar   : Absolute Salinity                  [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: long, lat, p
real (r14) :: saar, gsw_sa_from_sstar, sstar

saar = gsw_saar(p,long,lat)

gsw_sa_from_sstar = sstar*(1d0 + saar)/(1d0 - 0.35d0*saar)

! In the Baltic Sea, Sstar = SA, and note that gsw_saar returns zero
! for SAAR in the Baltic.

if (saar.eq.9d15) then
    gsw_sa_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sstar(sstar,p,long,lat)
!==========================================================================

! Calculates Practical Salinity, SP, from Preformed Salinity, Sstar. 
!
! sstar  : Preformed Salinity                              [g/kg]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_sp_from_Sstar : Preformed Salinity                   [g/kg]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: long, lat, p
real (r14) :: saar, gsw_sp_from_sstar, sp_baltic, Sstar

saar = gsw_saar(p,long,lat)

gsw_sp_from_sstar = (35.d0/35.16504d0)*Sstar/(1 - 0.35d0*saar);

!In the Baltic Sea, SA = Sstar.
sp_baltic = gsw_sp_from_sa_baltic(sstar,long,lat);

if (sp_Baltic.lt.1d10) then
    gsw_sp_from_sstar = sp_baltic;
end if

if (saar.eq.9d15) then
    gsw_sp_from_sstar = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! Library functions of the GSW toolbox
!--------------------------------------------------------------------------

!==========================================================================
function gsw_saar(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly Ratio, SAAR.
!
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_saar : Absolute Salinity Anomaly Ratio               [unitless]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled, k
integer :: flag_saar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: saar, saar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref
!real (r14), dimension(nz,ny,nx) :: delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_saar, p0_original, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, ndepth_max

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled/0/

save icalled, longs_ref, lats_ref, p_ref, ndepth_ref, saar_ref

flag_saar = 0
gsw_saar = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_saar = 9d15
 return
end if

if(long .lt. 0) then
 long = long + 360
end if

if(icalled.eq.0) then
   icalled = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_saar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   !read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  saar_ref = 9d15
   flag_saar = 0
2  continue
end if

if (flag_saar.eq.0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_saar"
end if

!Set gsw_saar = 9d15 and return if there is no data file present
if(flag_saar .eq. 0) then
 gsw_saar = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1 + (nx-1)*(long-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1
end if

indy0 = floor(1 + (ny-1)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if(indy0.eq.ny) then
   indy0 = ny-1
end if

ndepth_max = -1
do k = 1,4
   if(ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.d0) then
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
   end if
end do

if(ndepth_max.eq.-1.d0) then
  gsw_saar = 0d0 
   return
end if 

p0_original = p
if(p.gt.p_ref(int(ndepth_max))) then
 p = p_ref(int(ndepth_max))
end if
call indx(p_ref,nz,p,indz0)
    
r1 = (long-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0));
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0));
t1 = (p-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0));

do k = 1,4
   saar(k) = saar_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
  saar_old = saar
  call gsw_add_barrier(saar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)) .ge. 1d10) then 
   saar = gsw_add_mean(saar)
end if

sa_upper = (1.d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))

do k = 1,4
   saar(k) = saar_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   saar_old = saar
   call gsw_add_barrier(saar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,saar)
else if(abs(sum(saar)) .ge. 1d10) then 
   saar = gsw_add_mean(saar)
end if

sa_lower = (1.d0-s1)*(saar(1) + r1*(saar(2)-saar(1))) + s1*(saar(4) + r1*(saar(3)-saar(4)))
if(abs(sa_lower) .ge. 1d10) then
  sa_lower = sa_upper
end if
gsw_saar = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_saar).ge.1d10) then
   gsw_saar = 9d15
endif

p = p0_original
  
return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_delta_sa_ref(p,long,lat)
!==========================================================================

! Calculates the Absolute Salinity Anomaly reference value, delta_SA_ref.
!
! p      : sea pressure                                    [dbar]
! long   : longiture                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_delta_sa_ref : Absolute Salinity Anomaly reference value    [g/kg]

implicit none

!integer, parameter :: int9 = selected_int_kind(9) 
integer, parameter :: r14 = selected_real_kind(14,30)

integer, parameter :: nx=91, ny=45, nz=45

integer :: indx0, indy0, indz0, i, j, icalled2, k
integer :: flag_dsar
integer, dimension(4) :: deli, delj

real (r14), dimension(4) :: dsar, dsar_old
real (r14), dimension(nz) :: p_ref
real (r14), dimension(ny) :: lats_ref
real (r14), dimension(nx) :: longs_ref
real (r14), dimension(ny,nx) :: ndepth_ref 
real (r14), dimension(nz,ny,nx) :: saar_ref, delta_sa_ref
real (r14) :: p, lat, long, dlong, dlat
real (r14) :: gsw_delta_sa_ref, p0_original, sa_upper, sa_lower 
real (r14) :: r1, s1, t1, ndepth_max

data deli/0,1,1,0/, delj/0,0,1,1/

data icalled2/0/

save icalled2, longs_ref, lats_ref, p_ref, ndepth_ref, delta_sa_ref

flag_dsar = 0
gsw_delta_sa_ref = 9d15

if(lat .lt. -86d0 .or. lat .gt. 90d0) then
 gsw_delta_sa_ref = 9d15
 return
end if

if(long .lt. 0) then
 long = long + 360
end if

if(icalled2.eq.0) then
   icalled2 = 1
   open(10,file='gsw_data_v3_0.dat',status='old',err=1)
   flag_dsar = 1
   read(10,*) (longs_ref(i), i=1,nx)
   read(10,*) (lats_ref(i), i=1,ny)
   read(10,*) (p_ref(i), i=1,nz)
   read(10,*) ((ndepth_ref(j,i), j=1,ny), i=1,nx)
   read(10,*) (((saar_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   read(10,*) (((delta_sa_ref(k,j,i), k=1,nz), j=1,ny), i=1,nx)
   close(10)
   go to 2
1  delta_sa_ref = 9d15
   flag_dsar = 0
2  continue
end if

if (flag_dsar.eq.0) then
   write(*,*) "*** gsw_data_v3_0.dat is missing !!! ***"
   write(*,*) "Set the full path of gsw_data_v3_0.dat in gsw_delta_sa_ref"
end if

!Set gsw_delta_sa_ref = 9d15 and return if there is no data set present
if(flag_dsar.eq.0) then
 gsw_delta_sa_ref = 9d15
 return
endif

dlong = longs_ref(2)-longs_ref(1)
dlat = lats_ref(2)-lats_ref(1)

indx0 = floor(1 + (nx-1)*(long-longs_ref(1))/(longs_ref(nx)-longs_ref(1)))
if(indx0.eq.nx) then
   indx0 = nx-1
end if

indy0 = floor(1 + (ny-1)*(lat-lats_ref(1))/(lats_ref(ny)-lats_ref(1)))
if(indy0.eq.ny) then
   indy0 = ny-1
end if

ndepth_max = -1
do k = 1,4
   if(ndepth_ref(indy0+delj(k),indx0+deli(k)).gt.0.d0) then
      ndepth_max = max(ndepth_max,ndepth_ref(indy0+delj(k),indx0+deli(k)))
   end if
end do

if(ndepth_max.eq.-1.d0) then
  gsw_delta_sa_ref = 0d0 
   return
end if 

p0_original = p
if(p.gt.p_ref(int(ndepth_max))) then
 p = p_ref(int(ndepth_max))
end if
call indx(p_ref,nz,p,indz0)
    
r1 = (long-longs_ref(indx0))/(longs_ref(indx0+1)-longs_ref(indx0));
s1 = (lat-lats_ref(indy0))/(lats_ref(indy0+1)-lats_ref(indy0));
t1 = (p-p_ref(indz0))/(p_ref(indz0+1)-p_ref(indz0));

do k = 1,4
   dsar(k) = delta_sa_ref(indz0,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
  dsar_old = dsar
  call gsw_add_barrier(dsar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)) .ge. 1d10) then 
   dsar = gsw_add_mean(dsar)
end if

sa_upper = (1.d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))

do k = 1,4
   dsar(k) = delta_sa_ref(indz0+1,indy0+delj(k),indx0+deli(k))
end do

if(260.d0.le.long.and.long.le.291.999d0.and.3.4d0.le.lat.and.lat.le.19.55d0) then
   dsar_old = dsar
   call gsw_add_barrier(dsar_old,long,lat,longs_ref(indx0),lats_ref(indy0),dlong,dlat,dsar)
else if(abs(sum(dsar)) .ge. 1d10) then 
   dsar = gsw_add_mean(dsar)
end if

sa_lower = (1.d0-s1)*(dsar(1) + r1*(dsar(2)-dsar(1))) + s1*(dsar(4) + r1*(dsar(3)-dsar(4)))
if(abs(sa_lower) .ge. 1d10) then
  sa_lower = sa_upper
end if
gsw_delta_sa_ref = sa_upper + t1*(sa_lower-sa_upper)

if(abs(gsw_delta_sa_ref).ge.1d10) then
   gsw_delta_sa_ref = 9d15
endif

p = p0_original
  
return
end function

!--------------------------------------------------------------------------

!==========================================================================
subroutine gsw_add_barrier(input_data,long,lat,long_grid,lat_grid,dlong_grid,dlat_grid,output_data)
!==========================================================================

!  Adds a barrier through Central America (Panama) and then averages
!  over the appropriate side of the barrier
! 
!  data_in      :  data                                                     [unitless]
!  long         :  Longitudes of data in decimal degrees east               [ 0 ... +360 ]
!  lat          :  Latitudes of data in decimal degrees north               [ -90 ... +90 ]
!  longs_grid   :  Longitudes of regular grid in decimal degrees east       [ 0 ... +360 ]
!  lats_grid    :  Latitudes of regular grid in decimal degrees north       [ -90 ... +90 ]
!  dlongs_grid  :  Longitude difference of regular grid in decimal degrees  [ deg longitude ]
!  dlats_grid   :  Latitude difference of regular grid in decimal degrees   [ deg latitude ]
!
! gsw_add_barrier  : average of data depending on which side of the 
!                    Panama cannal it is on                                 [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer, dimension(4) :: above_line
integer k, nmean, above_line0, kk
real (r14), dimension(4) :: input_data, output_data
real (r14), dimension(6) :: longs_pan, lats_pan
real (r14) :: long, lat, r, lats_line, long_grid, lat_grid
real (r14) :: dlong_grid, dlat_grid, data_mean

data longs_pan/260.0000, 272.5900, 276.5000, 278.6500, 280.7300, 292.000/ 
data  lats_pan/ 19.5500,  13.9700,   9.6000,   8.1000,   9.3300,   3.400/ 

call indx(longs_pan,6,long,k)                            !   the long/lat point
r = (long-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat) then
   above_line0 = 1
else
   above_line0 = 0
end if

call indx(longs_pan,6,long_grid,k)                                     !  the 1 and 4 long/lat points 
r = (long_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(1) = 1
else
   above_line(1) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(4) = 1
else
   above_line(4) = 0
end if

call indx(longs_pan,6,long_grid+dlong_grid,k)                              !  the 2 and 3 long/lat points 
r = (long_grid+dlong_grid-longs_pan(k))/(longs_pan(k+1)-longs_pan(k))
lats_line = lats_pan(k) + r*(lats_pan(k+1)-lats_pan(k))

if(lats_line.le.lat_grid) then
   above_line(2) = 1
else
   above_line(2) = 0
end if

if(lats_line.le.lat_grid+dlat_grid) then
   above_line(3) = 1
else
   above_line(3) = 0
end if

nmean = 0 
data_mean = 0.d0

do kk = 1,4
   if ((abs(input_data(kk)).le.100d0).and.above_line0.eq.above_line(kk)) then
      nmean = nmean+1
      data_mean = data_mean+input_data(kk)
   end if
end do

if(nmean == 0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do kk = 1,4
   if((abs(input_data(kk)).ge.1d10).or.above_line0.ne.above_line(kk)) then
      output_data(kk) = data_mean
   else
      output_data(kk) = input_data(kk)
   end if
end do

return
end subroutine

!--------------------------------------------------------------------------

!==========================================================================
function gsw_add_mean(data_in)
!==========================================================================

! Replaces NaN's with non-nan mean of the 4 adjacent neighbours
!
! data   : data set of the 4 adjacent neighbours   
! p      : sea pressure                                        [dbar]
!
! gsw_add_mean : non-nan mean of the 4 adjacent neighbours     [unitless]

implicit none

integer, parameter :: int9 = selected_int_kind(9)
integer, parameter :: r14 = selected_real_kind(14,30)

integer :: k, nmean

real (r14), dimension(4) :: data_in, gsw_add_mean
real (r14) :: data_mean

nmean = 0
data_mean = 0.d0
gsw_add_mean = data_in

do k = 1,4
   if (abs(data_in(k)).le.100d0) then
      nmean = nmean+1
      data_mean = data_mean+data_in(k)
   end if
end do

if(nmean.eq.0)then
   data_mean = 0d0    !errorreturn
else
   data_mean = data_mean/nmean
endif

do k = 1,4
   if(abs(data_in(k)).ge.100d0) then
      gsw_add_mean(k) = data_mean
   end if
end do

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function xinterp1(x,y,n,x0)
!==========================================================================

! Linearly interpolate a real array   
!
! x      : y array (Must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! xinterp1 : Linearly interpolated value

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k

real (r14), dimension(n) :: x, y
real (r14) :: x0, r, xinterp1

call indx(x,n,x0,k)
r = (x0-x(k))/(x(k+1)-x(k))
xinterp1 = y(k) + r*(y(k+1)-y(k))

return
end function

!--------------------------------------------------------------------------

!==========================================================================
subroutine indx(x,n,z,k)
!==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  K      : index K - if X(K) <= Z < X(K+1), or
!  N-1     		    - if Z = X(N)
!

integer, parameter :: r14 = selected_real_kind(14,30)

integer :: n, k, ku, kl, km

real (r14), dimension(n) :: x
real (r14) :: z

if(z.gt.x(1).and.z.lt.x(n)) then
   kl=1
   ku=n
   do while (ku-kl.gt.1)
      km=(ku+kl)/2
      if(z.gt.x(km)) then
         kl=km
      else
         ku=km
      endif
   end do
   k=kl
   if(z.eq.x(k+1)) then 
     k = k+1
   end if
elseif (z.le.x(1)) then
      k = 1
elseif (z.ge.x(n)) then
      k = n-1
else
      write(*,*) 'ERROR in subroutine indx : out of range'
      write(*,*) 'z = ', z, 'n = ', n, 'x = ',x
end if

return
end subroutine

!--------------------------------------------------------------------------

!==========================================================================
function gsw_fdelta(p,long,lat)
!==========================================================================

! Calculates fdelta. 
!
! sp     : Practical Salinity                              [unitless]
! p      : sea pressure                                    [dbar]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
!
! gsw_fdelta : Absolute Salinty Anomaly                    [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14) :: long, lat, p, saar, gsw_fdelta

saar = gsw_saar(p,long,lat)

gsw_fdelta = ((1d0 + 0.35d0)*saar)/(1d0 - 0.35d0*saar);

if (saar.gt.1d10) then
    gsw_fdelta = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sa_from_sp_baltic(sp,long,lat)
!==========================================================================

! For the Baltic Sea, calculates Absolute Salinity with a value
! computed analytically from Practical Salinity
!
! sp     : Practical Salinity                              [unitless]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sa_from_sp_baltic : Absolute Salinity                [g/kg]                         

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), dimension(2) :: xb_right, yb_right
real (r14), dimension(3) :: xb_left, yb_left
real (r14) :: sp, long, lat, gsw_sa_from_sp_baltic, xx_left, xx_right

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

if(xb_left(2).lt.long .and. long.lt.xb_right(1) .and. yb_left(1).lt.lat .and. lat.lt.yb_left(3)) then
  
    xx_left = xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sa_from_sp_baltic =((35.16504d0 - 0.087d0)/35d0)*sp + 0.087d0
    else
        gsw_sa_from_sp_baltic = 9d15
    end if

else
    gsw_sa_from_sp_baltic = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================
function gsw_sp_from_sa_baltic(sa,long,lat)
!==========================================================================

! For the Baltic Sea, calculates Practical Salinity with a value
! computed analytically from Absolute Salinity
!
! sa     : Absolute Salinity                               [g/kg]
! long   : longitude                                       [deg E]     
! lat    : latitude                                        [deg N]
! p      : sea pressure                                    [dbar]
!
! gsw_sp_from_sa_baltic  : Practical Salinity              [unitless]

implicit none

integer, parameter :: r14 = selected_real_kind(14,30)

real (r14), dimension(2) :: xb_right, yb_right
real (r14), dimension(3) :: xb_left, yb_left
real (r14) :: sa, long, lat, gsw_sp_from_sa_baltic, xx_left, xx_right

data xb_left/12.6d0, 7.d0, 26.d0/, yb_left/50.d0, 59.d0, 69.d0/
data xb_right/45.d0, 26.d0/, yb_right/50.d0, 69.d0/

if(xb_left(2).lt.long .and. long.lt.xb_right(1) .and. yb_left(1).lt.lat .and. lat.lt.yb_left(3)) then
  
    xx_left = xinterp1(yb_left, xb_left, 3, lat)
    
    xx_right = xinterp1(yb_right, xb_right, 2, lat)
    
    if(xx_left.le.long .and. long.le.xx_right) then
        gsw_sp_from_sa_baltic = (35.d0/(35.16504d0 - 0.087d0))*(sa - 0.087d0)
    else
        gsw_sp_from_sa_baltic = 9d15
    end if
     
else
    gsw_sp_from_sa_baltic = 9d15
end if

return
end function

!--------------------------------------------------------------------------

!==========================================================================


end module convert_0
