module ice_flu_5

!#########################################################################

!THIS MODULE IMPLEMENTS THE CORRELATION FUNCTIONS FOR MELTING AND SUBLIMATION OF ICE,
!AS PROPOSED IN

!RAINER FEISTEL:
!DEVELOPMENT OF NEW EQUATIONS FOR MELTING PRESSURE AND SUBLIMATION PRESSURE.
!TASKGROUP REPORT,
!IAPWS THERMOPHYSICAL PROPERTIES OF WATER AND STEAM WG
!WITNEY, UK, 3 TO 8 SEPTEMBER 2006

!ANd

!REVISED RELEASE ON THE PRESSURE ALONG THE MELTING AND
!SUBLIMATION CURVES OF ORDINARY WATER SUBSTANCE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

!IMPLEMENTATION IN VB6 BY RAINER FEISTEL
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

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULE
!     CONSTANTS_0, file CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: fit_ice_liq_pressure_si, fit_ice_liq_temperature_si, &
          fit_ice_vap_pressure_si

contains

!==========================================================================
function fit_ice_liq_pressure_si(t_si)
!==========================================================================

!THIS IMPLEMENTS THE MELTING PRESSURE IN PA AS A FUNCTION OF TEMPERATURE IN K
!DETERMINED FROM A FIT OF THE FUNCTION ICE_LIQ_PRESSURE_SI,
!AS DEFINED IN EQ. (1) OF

!REVISED RELEASE ON THE PRESSURE ALONG THE MELTING AND
!SUBLIMATION CURVES OF ORDINARY WATER SUBSTANCE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

real*8 fit_ice_liq_pressure_si, t_si
real*8 pt, tt, a1, b1, a2, b2, a3, b3, tr, pr

pt = tp_pressure_exp_si
tt = tp_temperature_si

a1 = 1195393.37d0
b1 = 3d0
a2 = 80818.3159d0
b2 = 25.75d0
a3 = 3338.2686d0
b3 = 103.75d0

fit_ice_liq_pressure_si = errorreturn

if(t_si <= 0) return

tr = t_si / tt
pr = 1d0 + a1 * (1d0 - tr ** b1) + a2 * (1d0 - tr ** b2) + a3 * (1d0 - tr ** b3)

fit_ice_liq_pressure_si = pt * pr

end function

!==========================================================================
function fit_ice_liq_temperature_si(p_si)
!==========================================================================

!THIS IMPLEMENTS THE MELTING TEMPERATURE OF ICE IN K AS A FUNCTION OF PRESSURE IN PA
!DETERMINED FROM A FIT OF THE FUNCTION ICE_LIQ_TEMPERATURE_SI,
!AS DEFINED IN

!RAINER FEISTEL:
!DEVELOPMENT OF NEW EQUATIONS FOR MELTING PRESSURE AND SUBLIMATION PRESSURE.
!TASK GROUP REPORT, IAPWS THERMOPHYSICAL PROPERTIES OF WATER AND STEAM WG
!WITNEY, UK, 3 TO 8 SEPTEMBER 2006

real*8 fit_ice_liq_temperature_si, p_si
real*8 pt, tt, a1, a2, a3, a4, tr, pr

pt = tp_pressure_exp_si
tt = tp_temperature_si

a1 = -1.66356104484551d-07
a2 = -2.13519241979406d-13
a3 = 3.52967405341877d-20
a4 = -2.73184525236281d-26

fit_ice_liq_temperature_si = errorreturn

if(p_si <= 0) return

pr = p_si / pt

tr = 1d0 + a1 * (pr - 1d0) + a2 * (pr - 1d0) ** 2 + a3 * (pr - 1d0) ** 3 + a4 * (pr - 1d0) ** 4

fit_ice_liq_temperature_si = tt * tr

end function

!==========================================================================
function fit_ice_vap_pressure_si(t_si)
!==========================================================================

!THIS IMPLEMENTS THE SUBLIMATION PRESSURE IN PA AS A FUNCTION OF TEMPERATURE IN K
!DETERMINED FROM A FIT OF THE FUNCTION ICE_VAP_PRESSURE_SI,
!AS DEFINED IN EQ. (6) OF

!REVISED RELEASE ON THE PRESSURE ALONG THE MELTING AND
!SUBLIMATION CURVES OF ORDINARY WATER SUBSTANCE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

real*8 fit_ice_vap_pressure_si, t_si
real*8 pt, tt, a1, b1, a2, b2, a3, b3, tr, pr

pt = tp_pressure_exp_si
tt = tp_temperature_si

a1 = -21.2144006d0
b1 = 0.00333333333d0
a2 = 27.3203819d0
b2 = 1.20666667d0
a3 = -6.1059813d0
b3 = 1.70333333d0

fit_ice_vap_pressure_si = errorreturn

if(t_si <= 0) return

tr = t_si / tt
pr = a1 * tr ** b1 + a2 * tr ** b2 + a3 * tr ** b3

fit_ice_vap_pressure_si = pt * Exp(pr / tr)

end function

end module ice_flu_5

