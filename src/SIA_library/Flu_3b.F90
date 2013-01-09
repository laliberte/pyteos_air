module flu_3b

!#########################################################################

!THIS MODULE IMPLEMENTS THERMODYNAMIC PROPERTIES OF LIQUID WATER AND VAPOUR
!AS FUNCTIONS OF TEMPERATURE AND PRESSURE, COMPUTED NUMERICALLY FROM
!THE HELMHOLTZ FUNCTION OF FLUID WATER.

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
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     FLU_2,       FILE FLU_2.F90
!     FLU_3A,      FILE FLU_3A.F90
!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     FLU_1,       FILE FLU_1.F90

!#########################################################################

use constants_0
use flu_2
use flu_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: liq_cp_si, liq_cv_si, liq_enthalpy_si, liq_entropy_si, &
          liq_expansion_si, liq_gibbs_energy_si, liq_internal_energy_si, &
          liq_kappa_s_si, liq_kappa_t_si, liq_lapserate_si, liq_soundspeed_si, &
          vap_cp_si, vap_cv_si, vap_enthalpy_si, vap_entropy_si, &
          vap_expansion_si, vap_gibbs_energy_si, vap_internal_energy_si, &
          vap_kappa_s_si, vap_kappa_t_si, vap_lapserate_si, vap_soundspeed_si

contains

!=========================================================================
function liq_cp_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!CP(T,P) = T * (DETA/DT)_P  ISOBARIC HEAT CAPACITY OF LIQUID WATER IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_cp_si(300, 1e5) = 4180.63952202

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_cp_si, t_si, p_si
real*8 d

liq_cp_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_cp_si = flu_cp_si(t_si, d)

end function

!=========================================================================
function liq_cv_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!CV(T,P) = T * (DETA/DT)_D  ISOCHORIC HEAT CAPACITY OF LIQUID WATER IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_cv_si(300, 1e5) = 4130.17861503

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_cv_si, t_si, p_si
real*8 d

liq_cv_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_cv_si = flu_cv_si(t_si, d)

end function

!=========================================================================
function liq_enthalpy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!H(T,P) = F - T*(DF/DT)_D + D*(DF/DD)_T ENTHALPY OF LIQUID WATER IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_enthalpy_si(300, 1e5) = 112653.679689

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_enthalpy_si, t_si, p_si
real*8 d

liq_enthalpy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_enthalpy_si = flu_enthalpy_si(t_si, d)

end function

!=========================================================================
function liq_entropy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!ETA(T,P) = - (DF/DT)_D ENTROPY OF LIQUID WATER IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_entropy_si(300, 1e5) = 393.062433815

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_entropy_si, t_si, p_si
real*8 d

liq_entropy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_entropy_si = flu_entropy_si(t_si, d)

end function

!=========================================================================
function liq_expansion_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!ALPHA(T,P) = -(1/D) * (DD/DT)_P  THERMAL EXPANSION OF LIQUID WATER IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_expansion_si(300, 1e5) = 2.74803716256E-04


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_expansion_si, t_si, p_si
real*8 d

liq_expansion_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_expansion_si = flu_expansion_si(t_si, d)

end function

!=========================================================================
function liq_gibbs_energy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!G(T,P) = F + D*(DF/DD)_T   GIBBS ENERGY OF LIQUID WATER IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_gibbs_energy_si(300, 1e5) = -5265.05045577

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_gibbs_energy_si, t_si, p_si
real*8 d

liq_gibbs_energy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_gibbs_energy_si = flu_gibbs_energy_si(t_si, d)

end function

!=========================================================================
function liq_internal_energy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!U(T,P) = F - T*(DF/DT)_D INTERNAL ENERGY OF LIQUID WATER IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_internal_energy_si(300, 1e5) = 112553.334133


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_internal_energy_si, t_si, p_si
real*8 d

liq_internal_energy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_internal_energy_si = flu_internal_energy_si(t_si, d)

end function

!=========================================================================
function liq_kappa_s_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!KAPPA_S(T,P) = (1/D) * (DD/DP)_ETA  ISOTHERMAL COMPRESSIBILITY OF LIQUID WATER IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_kappa_s_si(300, 1e5) = 4.45077521253E-10

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_kappa_s_si, t_si, p_si
real*8 d

liq_kappa_s_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_kappa_s_si = flu_kappa_s_si(t_si, d)

end function

!=========================================================================
function liq_kappa_t_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!KAPPA_T(T,P) = (1/D) * (DD/DP)_ETA  ISOTHERMAL COMPRESSIBILITY OF LIQUID WATER IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_kappa_t_si(300, 1e5) = 4.50515304336E-10

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_kappa_t_si, t_si, p_si
real*8 d

liq_kappa_t_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_kappa_t_si = flu_kappa_t_si(t_si, d)

end function

!=========================================================================
function liq_lapserate_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!GAMMA(T,P) = (DT/DP)_ETA  ADIABATIC LAPSE RATE OF LIQUID WATER IN K/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_lapserate_si(300, 1e5) = 1.97878804448E-08

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_lapserate_si, t_si, p_si
real*8 d

liq_lapserate_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_lapserate_si = flu_lapserate_si(t_si, d)

end function

!=========================================================================
function liq_soundspeed_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!C(T,P) = SQRT[ (DP/dD)_ETA ] SOUND SPEED OF LIQUID WATER IN M/S

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: liq_soundspeed_si(300, 1e5) = 1501.52041506


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 liq_soundspeed_si, t_si, p_si
real*8 d

liq_soundspeed_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = liq_density_si(t_si, p_si)
if(d == errorreturn) return

liq_soundspeed_si = flu_soundspeed_si(t_si, d)

end function

!=========================================================================
function vap_cp_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!CP(T,P) = T * (DETA/DT)_P  ISOBARIC HEAT CAPACITY OF WATER VAPOUR IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_cp_si(300, 1e3) = 1874.12149028

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_cp_si, t_si, p_si
real*8 d

vap_cp_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_cp_si = flu_cp_si(t_si, d)

end function

!=========================================================================
function vap_cv_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!CV(T,P) = T * (DETA/DT)_D  ISOCHORIC HEAT CAPACITY OF WATER VAPOUR IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_cv_si(300, 1e3) = 1410.22845789

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_cv_si, t_si, p_si
real*8 d

vap_cv_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_cv_si = flu_cv_si(t_si, d)

end function

!=========================================================================
function vap_enthalpy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!H(T,P) = F - T*(DF/DT)_D + D*(DF/DD)_T ENTHALPY OF WATER VAPOUR IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_enthalpy_si(300, 1e3) = 2551013.47892

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_enthalpy_si, t_si, p_si
real*8 d

vap_enthalpy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_enthalpy_si = flu_enthalpy_si(t_si, d)

end function

!=========================================================================
function vap_entropy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!ETA(T,P) = - (DF/DT)_D ENTROPY OF WATER VAPOUR IN J/(KG K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_entropy_si(300, 1e3) = 9103.67940087

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_entropy_si, t_si, p_si
real*8 d

vap_entropy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_entropy_si = flu_entropy_si(t_si, d)

end function

!=========================================================================
function vap_expansion_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!ALPHA(T,P) = -(1/D) * (DD/DT)_P  THERMAL EXPANSION OF WATER VAPOUR IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_expansion_si(300, 1e3) = 3.34352010567E-03

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_expansion_si, t_si, p_si
real*8 d

vap_expansion_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_expansion_si = flu_expansion_si(t_si, d)

end function

!=========================================================================
function vap_gibbs_energy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!G(T,P) = F + D*(DF/DD)_T   GIBBS ENERGY OF WATER VAPOUR IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_gibbs_energy_si(300, 1e3) = -180090.341338

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_gibbs_energy_si, t_si, p_si
real*8 d

vap_gibbs_energy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_gibbs_energy_si = flu_gibbs_energy_si(t_si, d)

end function

!=========================================================================
function vap_internal_energy_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!U(T,P) = F - T*(DF/DT)_D INTERNAL ENERGY OF WATER VAPOUR IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_internal_energy_si(300, 1e3) = 2412625.00085

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_internal_energy_si, t_si, p_si
real*8 d

vap_internal_energy_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_internal_energy_si = flu_internal_energy_si(t_si, d)

end function

!=========================================================================
function vap_kappa_s_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!KAPPA_S(T,P) = (1/D) * (DD/DP)_ETA  ISENTROPIC COMPRESSIBILITY OF WATER VAPOUR IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_kappa_s_si(300, 1e3) = 7.52840457971E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_kappa_s_si, t_si, p_si
real*8 d

vap_kappa_s_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_kappa_s_si = flu_kappa_s_si(t_si, d)

end function

!=========================================================================
function vap_kappa_t_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!KAPPA_T(T,P) = (1/D) * (DD/DP)_ETA  ISOTHERMAL COMPRESSIBILITY OF WATER VAPOUR IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_kappa_t_si(300, 1e3) = 1.00048646242E-03

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_kappa_t_si, t_si, p_si
real*8 d

vap_kappa_t_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_kappa_t_si = flu_kappa_t_si(t_si, d)

end function

!=========================================================================
function vap_lapserate_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!GAMMA(T,P) = (DT/DP)_ETA  ADIABATIC LAPSE RATE OF WATER VAPOUR IN K/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_lapserate_si(300, 1e3) = 7.40674488635E-02


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_lapserate_si, t_si, p_si
real*8 d

vap_lapserate_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_lapserate_si = flu_lapserate_si(t_si, d)

end function

!=========================================================================
function vap_soundspeed_si(t_si, p_si)
!=========================================================================

!OUTPUT:
!C(T,P) = SQRT[ (DP/DD)_ETA ] SOUND SPEED OF WATER VAPOUR IN M/S

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K,
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: vap_soundspeed_si(300, 1e3) = 428.744430495

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A

implicit none
real*8 vap_soundspeed_si, t_si, p_si
real*8 d

vap_soundspeed_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = vap_density_si(t_si, p_si)
if(d == errorreturn) return

vap_soundspeed_si = flu_soundspeed_si(t_si, d)

end function

end module flu_3b
