module ice_2

!#########################################################################

!THIS MODULE IMPLEMENTS THERMODYNAMIC PROPERTIES OF HEXAGONAL ICE I
!AS FUNCTIONS OF ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA,
!COMPUTED FROM THE GIBBS POTENTIAL ICE_G_SI(DRV_T, DRV_P, T_SI, P_SI)

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

!THIS MODULE REQUIRES THE LIBRARY MODULES 
!     CONSTANTS_0,       FILE CONSTANTS_0.F90
!     ICE_1 (IAPWS06),   FILE ICE_1.F90 

!#########################################################################

use constants_0
use ice_1

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: ice_chempot_si, ice_cp_si, ice_density_si, ice_enthalpy_si, &
          ice_entropy_si, ice_expansion_si, ice_helmholtz_energy_si, &
          ice_internal_energy_si, ice_kappa_s_si, ice_kappa_t_si, &
          ice_lapserate_si, ice_p_coefficient_si, ice_specific_volume_si

contains

!==========================================================================
function ice_chempot_si(t_si, p_si)
!==========================================================================

implicit none
real*8 t_si, p_si, ice_chempot_si

!OUTPUT:
!G(T,P) GIBBS ENERGY = CHEMICAL POTENTIAL IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_chempot_si(270,1E5) = -3786.74963128

ice_chempot_si = ice_g_si(0, 0, t_si, p_si)

end function

!==========================================================================
function ice_cp_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_cp_si
real*8 g_tt

!OUTPUT:
!CP(T,P) =  - T*(D2G/DT2)_P ISOBARIC HEAT CAPACITY IN J/(KG*K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_cp_si(270,1E5) = 2073.47946211

ice_cp_si = errorreturn

g_tt = ice_g_si(2, 0, t_si, p_si)
if (g_tt == errorreturn .or. g_tt >= 0d0) return

ice_cp_si = -t_si * g_tt

end function

!==========================================================================
function ice_density_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_density_si
real*8 g_p

!OUTPUT:
!1/(DG/DP)_T DENSITY IN KG/M3

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_density_si(270,1E5) = 917.181167192

ice_density_si = errorreturn

g_p = ice_g_si(0, 1, t_si, p_si)
if (g_p == errorreturn .or. g_p <= 0d0) return

ice_density_si = 1d0 / g_p

end function

!==========================================================================
function ice_enthalpy_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_enthalpy_si
real*8 g, g_t

!OUTPUT:
!H(T,P) = G - T*(DG/DT)_P ENTHALPY IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_enthalpy_si(270,1E5) = -339929.555499
 
ice_enthalpy_si = errorreturn

g = ice_g_si(0, 0, t_si, p_si)
g_t = ice_g_si(1, 0, t_si, p_si)
if (g == errorreturn .or. g_t == errorreturn) return

ice_enthalpy_si = g - t_si*g_t

end function

!==========================================================================
function ice_entropy_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_entropy_si
real*8 g_t

!OUTPUT:
!ETA(T,P) = -(DG/DT)_P ENTROPY IN J/(KG*K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_entropy_si(270,1E5) = -1244.97335506
 
ice_entropy_si = errorreturn

g_t = ice_g_si(1, 0, t_si, p_si)
if (g_t == errorreturn) return

ice_entropy_si = -g_t

end function

!==========================================================================
function ice_expansion_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_expansion_si
real*8 g_p, g_tp

!OUTPUT:
!ALPHA(T,P) = (D2G/DTDP)/(DG/DP)_T CUBIC THERMAL EXPANSION COEFFICIENT IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_expansion_si(270,1E5) = 1.58309329594E-04

ice_expansion_si = errorreturn

g_p = ice_g_si(0, 1, t_si, p_si)
g_tp = ice_g_si(1, 1, t_si, p_si)
if (g_tp == errorreturn .or. g_p == errorreturn .or. g_p <= 0d0) return

ice_expansion_si = g_tp / g_p

end function

!==========================================================================
function ice_helmholtz_energy_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_helmholtz_energy_si
real*8 g, g_p

!OUTPUT:
!F(T,P) = G - P*(DG/DP)_T HELMHOLTZ ENERGY (FREE ENERGY) IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_helmholtz_energy_si(270,1E5) = -3895.77934490

ice_helmholtz_energy_si = errorreturn

g = ice_g_si(0, 0, t_si, p_si)
g_p = ice_g_si(0, 1, t_si, p_si)
if (g == errorreturn .or. g_p == errorreturn .or. g_p <= 0d0) return

ice_helmholtz_energy_si = g - p_si * g_p

end function

!==========================================================================
function ice_internal_energy_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_internal_energy_si
real*8 g, g_t, g_p

!OUTPUT:
!U(T,P) = G - P*(DG/DP)_T - T*(DG/DT)_P INTERNAL ENERGY IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_internal_energy_si(270,1E5) = -340038.585212

ice_internal_energy_si = errorreturn

g = ice_g_si(0, 0, t_si, p_si)
g_t = ice_g_si(1, 0, t_si, p_si)
g_p = ice_g_si(0, 1, t_si, p_si)
if (g == errorreturn .or. g_t == errorreturn .or. g_p == errorreturn) return
if (g_p <= 0d0) return

ice_internal_energy_si = g - p_si * g_p - t_si * g_t

end function

!==========================================================================
function ice_kappa_s_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_kappa_s_si
real*8 g_p, g_tt, g_tp, g_pp

!OUTPUT:
!-(1/V)(DV/DP)_ETA ISENTROPIC COMPRESSIBILITY IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_kappa_s_si(270,1E5) = 1.13667916416E-10
 
ice_kappa_s_si = errorreturn

g_p = ice_g_si(0, 1, t_si, p_si)
if (g_p == errorreturn) return
if (g_p <= 0d0) return
g_tt = ice_g_si(2, 0, t_si, p_si)
if (g_tt == errorreturn) return
if (g_tt >= 0d0) return
g_tp = ice_g_si(1, 1, t_si, p_si)
if (g_tp == errorreturn) return
g_pp = ice_g_si(0, 2, t_si, p_si)
if (g_pp == errorreturn) return
if (g_pp >= 0d0) return

ice_kappa_s_si = (g_tp**2 - g_tt * g_pp) / (g_p * g_tt)

end function

!==========================================================================
function ice_kappa_t_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_kappa_t_si
real*8 g_p, g_pp

!OUTPUT:
!-(1/V)(DV/DP)_T ISOTHERMAL COMPRESSIBILITY IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_kappa_t_si(270,1E5) = 1.17226047281E-10

ice_kappa_t_si = errorreturn

g_p = ice_g_si(0, 1, t_si, p_si)
g_pp = ice_g_si(0, 2, t_si, p_si)
if (g_pp == errorreturn .or. g_pp >= 0d0 .or. &
  g_p == errorreturn .or. g_p <= 0d0) return

ice_kappa_t_si = -g_pp/g_p

end function

!==========================================================================
function ice_lapserate_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_lapserate_si
real*8 g_tp, g_tt

!OUTPUT:
!(DT/DP)_ETA ADIABATIC LAPSE RATE IN K/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_lapserate_si(270,1E5) = 2.24758128545E-08
 
ice_lapserate_si = errorreturn

g_tp = ice_g_si(1, 1, t_si, p_si)
g_tt = ice_g_si(2, 0, t_si, p_si)
if (g_tt == errorreturn .or. g_tt >= 0d0 .or. &
g_tp == errorreturn) return

ice_lapserate_si = -g_tp / g_tt

end function

!==========================================================================
function ice_p_coefficient_si(t_si, p_si)
!==========================================================================

implicit none

real*8 ice_p_coefficient_si
real*8 t_si, p_si
real*8 g_tp, g_pp

!OUTPUT:
!(DP/DT)_V  ISOCHORIC PRESSURE COEFFICIENT IN PA/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_p_coefficient_si(270,1E5) = 1350462.06254

ice_p_coefficient_si = errorreturn

g_tp = ice_g_si(1, 1, t_si, p_si)
if (g_tp == errorreturn) return
g_pp = ice_g_si(0, 2, t_si, p_si)
if (g_pp == errorreturn .or. g_pp >= 0d0) return

ice_p_coefficient_si = -g_tp / g_pp

end function

!==========================================================================
function ice_specific_volume_si(t_si, p_si)
!==========================================================================

implicit none

real*8 t_si, p_si, ice_specific_volume_si
real*8 g_p

!OUTPUT:
!V(T,P) = (DG/DP)_T SPECIFIC_VOLUME IN M3/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE: ice_specific_volume_si(270,1E5) = 1.09029713624E-03
 
ice_specific_volume_si = errorreturn

g_p = ice_g_si(0, 1, t_si, p_si)
if (g_p == errorreturn) return

ice_specific_volume_si = g_p

end function

end module ice_2
