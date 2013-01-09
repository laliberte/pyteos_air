module air_5

!#########################################################################

!THIS MODULE IMPLEMENTS PROPERTIES OF HUMID AIR THAT ARE EXPRESSED IN
!NON-BASIC SI UNITS OR REQUIRE ADDITIONAL EMPIRICAL CONSTANTS.

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

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0  FROM FILE CONSTANTS_0.F90 
!     AIR_3B       FROM FILE AIR_3B.F90 
!     LIQ_AIR_4A   FROM FILE LIQ_AIR_4A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0    FROM FILE CONVERT_0.F90 
!     MATHS_0      FROM FILE MATHS_0.F90 
!     AIR_1        FROM FILE AIR_1.F90 
!     FLU_1        FROM FILE FLU_1.F90 
!     AIR_2        FROM FILE AIR_2.F90 
!     FLU_2        FROM FILE FLU_2.F90 
!     AIR_3A       FROM FILE AIR_3A.F90 
!     FLU_3A       FROM FILE FLU_3A.F90 
!     LIQ_VAP_4    FROM FILE LIQ_VAP_4.F90 

!#########################################################################

use constants_0
use air_3b
use liq_air_4a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: air_lapserate_moist_c100m

contains

!==========================================================================
function air_lapserate_moist_c100m(rh_percent, t_c, p_hPa)
!==========================================================================

!THIS FUNCTION RETURNS THE ADIABATIC LAPSE RATE OF HUMID AIR IN DEG C PER 100 M ALTITUDE

!INPUT PARAMETERS
!RH_PERCENT:  RELATIVE HUMIDITY WMO DEFINITION IN % BETWEEN 0 AND 100
!T_C:         TEMPERATURE IN DEG C
!P_HPa:       PRESSURE IN HECTOPASCALIS (MILLIBARS)  !dbg

!CHECK VALUE, V1.0:
!AIR_LAPSERATE_MOIST_C100M(50, 20, 1000) = 0.971359680720

!CHECK VALUE, V1.1:
!AIR_LAPSERATE_MOIST_C100M(50, 20, 1000) = 0.971588085046

real*8 air_lapserate_moist_c100m, rh_percent, t_c, p_hPa
real*8 g, h, lr, d, a_si, t_si, p_si

g = 9.81d0  !m/s2   gravity
h = 100d0   !m      altitude difference

air_lapserate_moist_c100m = errorreturn

!convert to basic SI units:
t_si = 273.15d0 + t_c
p_si = p_hPa * 100d0

!compute dry-air fraction from RH
a_si = liq_air_a_from_rh_wmo_si(rh_percent * 0.01d0, t_si, p_si)
if (a_si == errorreturn) return

!compute lapse rate wrt pressure
lr = air_g_lapserate_si(a_si, t_si, p_si)
if (lr == errorreturn) return

!compute humid air density
d = air_g_density_si(a_si, t_si, p_si)
if (d == errorreturn) return

!get the lapse rate wrt 100 m sinking
air_lapserate_moist_c100m = lr * d * g * h

end function

end module air_5
