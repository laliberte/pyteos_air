Module Sea_3c

!#########################################################################

!THIS MODULE IMPLEMENTS THE ENTROPY OF SEAWATER DEPENDING ON SALINITY, 
!ENTHALPY AND PRESSURE.

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
!     SEA_3A,      FILE SEA_3A.F90
!     SEA_3B,      FILE SEA_3B.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     FLU_1,       FILE FLU_1.F90
!     SAL_1,       FILE SAL_1.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_3A,      FILE FLU_3A.F90

!#########################################################################

use constants_0
use sea_3a
use sea_3b

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_eta_contraction_h_si, sea_eta_contraction_t_si, &
          sea_eta_contraction_theta_si, sea_eta_density_si, &
          sea_eta_entropy_si, sea_eta_expansion_h_si, &
          sea_eta_expansion_t_si, sea_eta_expansion_theta_si, &
          sea_eta_potdensity_si, sea_eta_pottemp_si, &
          sea_eta_temperature_si, set_it_ctrl_entropy

!Control parameters of the entropy iteration
Integer ctrl_initialized
Integer ctrl_mode_entropy
Integer ctrl_loop_maximum
real*8 ctrl_init_entropy
real*8 ctrl_eps_exit_entropy

contains

!==========================================================================
function sea_eta_entropy_si(sa_si, h_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES SPECIFIC ENTROPY OF SEAWATER FROM ENTHALPY
!(THE THERMODYNAMIC POTENTIAL IN TERMS OF SALINITY, ENTHALPY AND PRESSURE IS ENTROPY)

!RETURNS   ENTROPY OF SEAWATER IN J/(KG K),

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      IN-SITU ENTHALPY IN J/KG,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      POTENTIAL ENTHALPY IN J/KG,
!P_SI      ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ETA_ENTROPY_SI(0.035, 1E5, 1E5) = 350.310622663

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none

Integer maxit
real*8 sea_eta_entropy_si, sa_si, h_si, p_si
real*8 eta, eps

sea_eta_entropy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(p_si <= 0d0) return

sea_eta_entropy_si = errorreturn

call init_it_ctrl_entropy

select case (ctrl_mode_entropy)
  case (0);    eta = h_si / 273.15d0
  case Default; eta = ctrl_init_entropy
end select

if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  sea_eta_entropy_si = eta
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

eps = ctrl_eps_exit_entropy
if(eps <= 0d0) return

sea_eta_entropy_si = entropy_iteration(sa_si, h_si, p_si, eta, maxit, eps)

end function

!=========================================================================
function sea_eta_temperature_si(sa_si, h_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES IN-SITU TEMPERATURE OF SEAWATER FROM EITHER IN-SITU ENTHALPY OR POTENTIAL ENTHALPY

!RETURNS   TEMPERATURE OF SEAWATER IN K

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_TEMPERATURE_SI(0.035,1E5,1E7,1E5,"h") = 295.985682129
!SEA_ETA_TEMPERATURE_SI(0.035,1E5,1E7,1E5,"hpot") = 298.413424848

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none
character*13 key
real*8 sea_eta_temperature_si, sa_si, h_si, p_si, pref_si 
real*8 eta

sea_eta_temperature_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return
if (key /= "h" .and. key /= "hpot") return
if (key == "hpot" .and. pref_si <= 0d0) return

if (key == "hpot") then
  eta = sea_eta_entropy_si(sa_si, h_si, pref_si)  !ENTROPY FROM POT. ENTHALPY
else
  eta = sea_eta_entropy_si(sa_si, h_si, p_si)     !ENTROPY FROM IN-SITU ENTHALPY
end if
if (eta == errorreturn) return

sea_eta_temperature_si = sea_h_si(0, 1, 0, sa_si, eta, p_si)

end function

!=========================================================================
function sea_eta_pottemp_si(sa_si, h_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES  POTENTIAL TEMPERATURE OF SEAWATER FROM EITHER IN-SITU ENTHALPY OR POTENTIAL ENTHALPY

!RETURNS   ABSOLUTE POTENTIAL TEMPERATURE OF SEAWATER IN K,

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      NOT REQUIRED
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_POTTEMP_SI(0.035,1E5,1E7,1E5,"h") = 295.782199115
!SEA_ETA_POTTEMP_SI(0.035,1E5,1E7,1E5,"hpot") = 298.194955182

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none
character*13 key
real*8 sea_eta_pottemp_si, sa_si, h_si, p_si, pref_si
real*8 eta

sea_eta_pottemp_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return
if (key /= "h" .and. key /= "hpot") return
if (key == "hpot" .and. pref_si <= 0d0) return

if (key == "hpot") then
  eta = sea_eta_entropy_si(sa_si, h_si, pref_si)  !ENTROPY FROM POT. ENTHALPY
else
  eta = sea_eta_entropy_si(sa_si, h_si, p_si)     !ENTROPY FROM IN-SITU ENTHALPY
end if
if (eta == errorreturn) return

sea_eta_pottemp_si = sea_h_si(0, 1, 0, sa_si, eta, pref_si)

end function

!=========================================================================
function sea_eta_density_si(sa_si, h_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES IN-SITU DENSITY OF SEAWATER FROM EITHER IN-SITU ENTHALPY OR POTENTIAL ENTHALPY

!RETURNS   DENSITY OF SEAWATER IN KG/M3

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      POTENTIAL ENTHALPY IN J/KG,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_DENSITY_SI(0.035,1E5,1E7,1E5,"h") = 1028.10986556
!SEA_ETA_DENSITY_SI(0.035,1E5,1E7,1E5,"hpot") = 1027.36529798

implicit none
character*13 key
real*8 sea_eta_density_si, sa_si, h_si, p_si, pref_si
real*8 eta, v

sea_eta_density_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return
if (key /= "h" .and. key /= "hpot") return
if (key == "hpot" .and. pref_si <= 0d0) return

if (key == "hpot") then
  eta = sea_eta_entropy_si(sa_si, h_si, pref_si)  !ENTROPY FROM POT. ENTHALPY
else
  eta = sea_eta_entropy_si(sa_si, h_si, p_si)     !ENTROPY FROM IN-SITU ENTHALPY
end if
if (eta == errorreturn) return

v = sea_h_si(0, 0, 1, sa_si, eta, p_si)      !SPECIFIC VOLUME
if (v == errorreturn .or. v <= 0d0) return

sea_eta_density_si = 1d0 / v

end function

!=========================================================================
function sea_eta_potdensity_si(sa_si, h_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES POTENTIAL DENSITY OF SEAWATER FROM EITHER IN-SITU ENTHALPY OR POTENTIAL ENTHALPY

!RETURNS   POTENTIAL DENSITY OF SEAWATER IN KG/M3

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      IN-SITU ENTHALPY IN J/KG,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"


!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      POTENTIAL ENTHALPY IN J/KG,
!P_SI      NOT REQUIRED
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_POTDENSITY_SI(0.035,1E5,1E7,1E5,"h") = 1023.91737474
!SEA_ETA_POTDENSITY_SI(0.035,1E5,1E7,1E5,"hpot") = 1023.20527737

implicit none
character*13 key
real*8 sea_eta_potdensity_si, sa_si, h_si, p_si, pref_si
real*8 eta, v

sea_eta_potdensity_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

if (key == "hpot") then
  eta = sea_eta_entropy_si(sa_si, h_si, pref_si)  !ENTROPY FROM POT. ENTHALPY
else
  eta = sea_eta_entropy_si(sa_si, h_si, p_si)     !ENTROPY FROM IN-SITU ENTHALPY
end if
if (eta == errorreturn) return

v = sea_h_si(0, 0, 1, sa_si, eta, pref_si)      !POT. SPECIFIC VOLUME
if (v == errorreturn .or. v <= 0d0) return

sea_eta_potdensity_si = 1d0 / v

end function

!=========================================================================
function sea_eta_contraction_t_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE HALINE CONTRACTION COEFFICIENT OF SEAWATER AT CONSTANT TEMPERATURE
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   HALINE CONTRACTION -(1/V)*(DV/DS)_T_P OF SEAWATER IN KG/KG

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"
  
!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!sea_eta_contraction_t_SI(0.035,1E5,1E7,1E5,"h") = 0.728755239644
!SEA_ETA_CONTRACTION_T_SI(0.035,1E5,1E7,1E5,"hpot") = 0.726349317428
!SEA_ETA_CONTRACTION_T_SI(0.035,300,1E7,1E5,"t") = 0.724913833446
!SEA_ETA_CONTRACTION_T_SI(0.035,300,1E7,1E5,"tpot") = 0.724714253918

implicit none
character*13 key
real*8 sea_eta_contraction_t_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_contraction_t_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !ENTROPY
if (eta == errorreturn) return

sea_eta_contraction_t_si = sea_h_contraction_t_si(sa_si, eta, p_si)

end function

!=========================================================================
function sea_eta_contraction_theta_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE HALINE CONTRACTION COEFFICIENT OF SEAWATER AT CONSTANT POTENTIAL TEMPERATURE
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   HALINE CONTRACTION - (1/V)*(DV/DS) OF SEAWATER IN KG/KG

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG,
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_CONTRACTION_THETA_SI(0.035,1E5,1E7,1E5,"h") = 0.728499505688
!SEA_ETA_CONTRACTION_THETA_SI(0.035,1E5,1E7,1E5,"hpot") = 0.726099732703
!SEA_ETA_CONTRACTION_THETA_SI(0.035,300,1E7,1E5,"t") = 0.724667977117
!SEA_ETA_CONTRACTION_THETA_SI(0.035,300,1E7,1E5,"tpot") = 0.724468894946

implicit none
character*13 key
real*8 sea_eta_contraction_theta_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_contraction_theta_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !ENTROPY
if (eta == errorreturn) return

sea_eta_contraction_theta_si = sea_h_contraction_theta_si(sa_si, eta, p_si, pref_si)

end function

!=========================================================================
function sea_eta_contraction_h_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE HALINE CONTRACTION COEFFICIENT OF SEAWATER AT CONSTANT POTENTIAL ENTHALPY
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   HALINE CONTRACTION - (1/V)*(DV/DS) OF SEAWATER IN KG/KG

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_CONTRACTION_H_SI(0.035,1E5,1E7,1E5,"h") = 0.718452125957
!SEA_ETA_CONTRACTION_H_SI(0.035,1E5,1E7,1E5,"hpot") = 0.714531922616
!SEA_ETA_CONTRACTION_H_SI(0.035,300,1E7,1E5,"t") = 0.712069013013
!SEA_ETA_CONTRACTION_H_SI(0.035,300,1E7,1E5,"tpot") = 0.711718411190

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none
character*13 key
real*8 sea_eta_contraction_h_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_contraction_h_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !entropy
if (eta == errorreturn) return

sea_eta_contraction_h_si = sea_h_contraction_h_si(sa_si, eta, p_si, pref_si)

end function

!=========================================================================
function sea_eta_expansion_t_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT TEMPERATURE
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   THERMAL EXPANSION  (1/V)*(DV/DT) OF SEAWATER IN 1/K
!          AT CONSTANT SALINITY AND PRESSURE.

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT USED
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT USED
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_EXPANSION_T_SI(0.035,1E5,1E7,1E5,"h") = 2.89480851145E-04
!SEA_ETA_EXPANSION_T_SI(0.035,1E5,1E7,1E5,"hpot") = 3.07242256461E-04
!SEA_ETA_EXPANSION_T_SI(0.035,300,1E7,1E5,"t") = 3.18513471410E-04
!SEA_ETA_EXPANSION_T_SI(0.035,300,1E7,1E5,"tpot") = 3.20122324740E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none
character*13 key
real*8 sea_eta_expansion_t_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_expansion_t_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !entropy
if (eta == errorreturn) return

sea_eta_expansion_t_si = sea_h_expansion_t_si(sa_si, eta, p_si)

end function

!=========================================================================
function sea_eta_expansion_theta_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   EXPANSION  (1/V)*(DV/DTHETA) OF SEAWATER IN 1/K
!          AT CONSTANT SALINITY AND PRESSURE. THETA IS POTENTIAL TEMPERATURE

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_EXPANSION_THETA_SI(0.035,1E5,1E7,1E5,"h") = 2.91293979902E-04
!SEA_ETA_EXPANSION_THETA_SI(0.035,1E5,1E7,1E5,"hpot") = 3.09134848554E-04
!SEA_ETA_EXPANSION_THETA_SI(0.035,300,1E7,1E5,"t") = 3.20454167783E-04
!SEA_ETA_EXPANSION_THETA_SI(0.035,300,1E7,1E5,"tpot") = 3.22069710839E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE

implicit none
character*13 key
real*8 sea_eta_expansion_theta_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_expansion_theta_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !entropy
if (eta == errorreturn) return

sea_eta_expansion_theta_si = sea_h_expansion_theta_si(sa_si, eta, p_si, pref_si)

end function

!=========================================================================
function sea_eta_expansion_h_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT POTENTIAL ENTHALPY
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   THERMAL EXPANSION  (1/V)*(DV/DHPOT)_S_P OF SEAWATER IN KG / J
!          AT CONSTANT SALINITY AND PRESSURE. HPOT IS POTENTIAL ENTHALPY.

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      ENTROPY COMPUTED IN SEA_3B AND ON THOSE MADE IN THIS MODULE
!
!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ETA_EXPANSION_H_SI(0.035,1E5,1E7,1E5,"h") = 7.28514646021E-08
!SEA_ETA_EXPANSION_H_SI(0.035,1E5,1E7,1E5,"hpot") = 7.72876772245E-08
!SEA_ETA_EXPANSION_H_SI(0.035,300,1E7,1E5,"t") = 8.01009066333E-08
!SEA_ETA_EXPANSION_H_SI(0.035,300,1E7,1E5,"tpot") = 8.05023387611E-08

implicit none
character*13 key
real*8 sea_eta_expansion_h_si, sa_si, x_si, p_si, pref_si
real*8 eta

sea_eta_expansion_h_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

eta = sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)  !entropy
if (eta == errorreturn) return

sea_eta_expansion_h_si = sea_h_expansion_h_si(sa_si, eta, p_si, pref_si)

end function

!=========================================================================
function entropy_iteration(sa_si, h_si, p_si, eta_si, maxit, eps)
!==========================================================================

!RETURNS   ETA =  ENTROPY OF SEAWATER IN J/(KG K)
!          I.E. THE ENTROPY THAT SOLVES H_SI = SEA_ENTHALPY_SI(SA_SI, ETA_SI, P_SI)
!          FOR ETA_SI AT GIVEN H_SI

!SA_SI     ABSOLUTE SALINITY IN KG/KG
!H_SI      ENTHALPY IN J/KG
!P_SI      ABSOLUTE PRESSURE IN PA
!ETA_SI    ENTROPY IN J/(KG K), INITIAL VALUE
!MAXIT     MAX. NUMBER OF ITERATIONS
!EPS       REQUIRED TOLERANCE IN J/(KG K)

implicit none

Integer maxit, i
real*8 entropy_iteration, sa_si, h_si, p_si, eta_si, eps
real*8 h, eta, t, de

entropy_iteration = errorreturn

if((eps <= 0d0) .or. (maxit <= 0)) return

if(check_limits == 1) then
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) return
else
  if(sa_si < 0d0 .or. sa_si >= 1d0 .or. &
  p_si <= 0d0) return
endif

check_limits = check_limits - 1

eta = eta_si

do i = 1, maxit

  !GET ENTHALPY AND ITS FIRST DERIVATIVE FOR NEWTON ITERATION
  h = sea_h_si(0, 0, 0, sa_si, eta, p_si)
  if(h == errorreturn) go to 999
  t = sea_h_si(0, 1, 0, sa_si, eta, p_si)
  if(t == errorreturn) go to 999
  if(t <= 0d0) go to 999
  
  !NEXT ENTROPY IMPROVEMENT STEP
  de = (h_si - h) / t
  eta = eta + de

  if(abs(de) < eps) then
    entropy_iteration = eta
    go to 999
  end if
  
enddo

999 continue

check_limits = check_limits + 1

end function

!=========================================================================
function sea_eta_entropy_x_si(sa_si, x_si, p_si, pref_si, key)
!=========================================================================

!THIS FUNCTION COMPUTES THE ENTROPY OF SEAWATER
!FROM EITHER IN-SITU ENTHALPY (KEY="H") OR POTENTIAL ENTHALPY (KEY = "HPOT")
!OR FROM IN-SITU TEMPERATURE (KEY="T") OR POTENTIAL TEMPERATURE (KEY = "TPOT")

!RETURNS   ENTROPY OF SEAWATER IN J/(KG K)

!EITHER FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "H"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL ENTHALPY IN J/KG
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "HPOT"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   NOT REQUIRED
!KEY       "T"

!OR FROM
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!X_SI      POTENTIAL TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA
!KEY       "TPOT"

character*13 key
real*8 sea_eta_entropy_x_si, sa_si, x_si, p_si, pref_si

sea_eta_entropy_x_si = errorreturn

if (sa_si < 0d0 .or. sa_si >= 1d0) return
if (p_si <= 0d0) return
if (pref_si <= 0d0) return

select case (key)
  case ("hpot")
    if (pref_si <= 0) return
    sea_eta_entropy_x_si = sea_eta_entropy_si(sa_si, x_si, pref_si)  !ENTROPY FROM POT. ENTHALPY

  case ("h")
    sea_eta_entropy_x_si = sea_eta_entropy_si(sa_si, x_si, p_si)     !ENTROPY FROM IN-SITU ENTHALPY

  case ("tpot")
    if (pref_si <= 0d0) return
    if (x_si <= 0d0) return
    sea_eta_entropy_x_si = sea_entropy_si(sa_si, x_si, pref_si)      !ENTROPY FROM POT. TEMPERATURE

  case ("t")
    if (x_si <= 0d0) return
    sea_eta_entropy_x_si = sea_entropy_si(sa_si, x_si, p_si)         !ENTROPY FROM IN-SITU TEMPERATURE

   case default; write(*,*)key,' not supported'

end select

end function

!==========================================================================
subroutine set_it_ctrl_entropy(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!ENTROPY FROM ENTHALPY

!KEY            VALUE
!IT_STEPS       0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS       N > 0       SET ITERATION NUMBER TO N
!IT_STEPS      -1           DO NOT ITERATE, USE INITIAL VALUE

!MODE_ETA       0           USE DEFAULT ENTROPY ETA = H/273.15 TO START
!MODE_ETA       1           USE VALUE SET BY INIT_ETA AS ENTROPY TO START

!INIT_ETA       ETA         USE VALUE ETA AS ENTROPY TO START IF MODE_ETA = 1

!TOL_ETA        0           USE DEFAULT EXIT ACCURACY FOR ENTROPY (1E-4 J/(KG K))
!TOL_ETA        EPS > 0     USE EPS AS EXIT ACCURACY FOR ENTROPY

implicit none
character*13 key
real*8 value

call init_it_ctrl_entropy

select case (Trim(key))

  case ("it_steps")                   !ITERATION STEPS
      if(value == 0d0)then 
         ctrl_loop_maximum = 100      !DEFAULT = 100
      elseif(value < 0d0) then 
         ctrl_loop_maximum = -1
      else
         ctrl_loop_maximum = value
      endif

  case ("mode_eta")   !start eta
      if(value == 0d0)then 
         ctrl_mode_entropy = 0        !DEFAULT: ETA = H/T0
      elseif(value == 1d0) then
         ctrl_mode_entropy = 1        !DEFAULT: ETA = INIT_ETA
      endif

  case ("init_eta")                   !START ETA IF CTRL_MODE_ENTROPY = 1
      ctrl_init_entropy = value

  case ("tol_eta")                    !REQUIRED ETA TOLERANCE
      if(value == 0d0)then
         ctrl_eps_exit_entropy = 0.0001d0 !DEFAULT = 1e-4 J/(KG K)
      elseif(value > 0d0) then
         ctrl_eps_exit_entropy = value
      endif

   case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine init_it_ctrl_entropy()
!==========================================================================

implicit none

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!Set default values and modes for entropy iteration
ctrl_loop_maximum = 100
ctrl_mode_entropy = 0             !DEFAULT: ETA = H/T0
ctrl_eps_exit_entropy = 0.0001d0  !DEFAULT = 1E-4 J/(KG K)
ctrl_init_entropy = 0

end subroutine

end module sea_3c
