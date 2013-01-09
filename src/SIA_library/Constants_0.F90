module constants_0

!#########################################################################

!THIS MODULE INITIALIZES CONSTANTS AND PARAMETER VALUES USED IN THE LIBRARY
!IT IS SELF-CONTAINED AND DOES NOT REQUIRE OTHER LIBRARY MODULES

!THE VALUES OF THESE CONSTANTS ARE CONSISTENT WITH THE FUNCTIONS AND
!CONSTANTS USED IN THIS SOURCE CODE LIBRARY, RATHER THAN BEING INDEPENDENT 
!ENTITIES.  FUTURE UPDATES ARE ONLY POSSIBLE IF THIS CONSISTENCY IS 
!CAREFULLY MAINTAINED

!#########################################################################

implicit none

character*16, private :: version = '8 Jun 2010'

!TO EXTEND THE GIbbS FUNCTION AS DESCRIBED IN:

!R. FEISTEL
!EXTENDED EQUATION OF STATE FOR SEAWATER AT ELEVATED TEMPERATURE AND SALINITY
!DESALINATION, 250, 14-18, 2010.

!TO COVER THE MEASUREMENTS OF

!F. MILLERO AND F. HUANG
!OCEAN SCI. DISCUSS. 6 (2009) 153-169.
!HTTP://WWW.OCEAN-SCI-DISCUSS.NET/6/153/2009/

!SET THIS CONSTANT TO TRUE:
!logical, public :: IsExtension2010 = .true.
!NOTE THAT THIS OPTIONAL GIBBS FUNCTION IS NOT ENDORSED BY IAPWS

logical, public :: IsExtension2010 = .false.

!IMPOSE NECESSARY PARAMETER LIMITS
!TO INCLUDE RANGE CHECKS, INITIALIZE CHECK_LIMITS TO 1. IT IS USED INTERNALLY
!TO TURN ALL LIMIT CHECKING OFF WITHIN THE ITERATIVE SOLUTION ROUTINES.

!INITIALIZE CHECK_LIMITS TO 0 TO USE VERY LOOSE RANGE CHECKING.
!INITIALIZE CHECK_LIMITS TO 1 TO IMPOSE MORE RESTRICTIVE RANGE CHECKING. 
!IN THE LATTER CASE, CHECK_LIMITS IS USED INTERNALLY TO TURN ALL
!LIMIT CHECKING OFF WITHIN THE ITERATIVE SOLUTION ROUTINES.

integer, public :: check_limits = 1

real*8, public, parameter :: dry_air_tmin = 60d0
real*8, public, parameter :: dry_air_tmax = 873d0
real*8, public, parameter :: dry_air_dmin = 0d0
real*8, public, parameter :: dry_air_dmax = 1035d0

real*8, public, parameter :: mix_air_tmin = 193d0     !or perhaps 173d0
real*8, public, parameter :: mix_air_tmax = 473d0
real*8, public, parameter :: mix_air_dmin = 0d0
real*8, public, parameter :: mix_air_dmax = 100d0

real*8, public, parameter :: flu_tmin = 50d0          !130d0
real*8, public, parameter :: flu_tmax = 1273d0
real*8, public, parameter :: flu_dmin = 0d0
real*8, public, parameter :: flu_dmax = 1240d0        !1200

real*8, public, parameter :: ice_tmin = 0d0
real*8, public, parameter :: ice_tmax = 273.16d0
real*8, public, parameter :: ice_pmin = 0d0
real*8, public, parameter :: ice_pmax = 2d8

real*8, public, parameter :: sal_tmin = 262d0
real*8, public, parameter :: sal_tmax = 353d0
real*8, public, parameter :: sal_pmin = 100d0
real*8, public, parameter :: sal_pmax = 1d8
real*8, public, parameter :: sal_smin = 0d0
real*8, public, parameter :: sal_smax = 0.12d0

!PROCEDURE RETURN CODES
real*8, public, parameter :: errorreturn = 9.99999999d+98
integer, public, parameter :: isok = -1

!MATHEMATICAL CONSTANTS
real*8, public, parameter :: pi = 3.14159265358979d0

!FUNDAMENTAL PHYSICAL CONSTANTS
real*8, public, parameter :: gas_constant_molar_si = 8.314472d0      !MOLAR GAS CONSTANT J MOL-1 K-1,IAPWS 2005
real*8, public, parameter :: gas_constant_molar_L2000 = 8.31451d0    !MOLAR GAS CONSTANT J MOL-1 K-1, LEMMON ET AL. 2000

real*8, public, parameter :: molar_mass_h2o_si = 0.018015268d0       !MOLAR MASS OF H2O IN KG/MOL, IAPWS 2009
real*8, public, parameter :: molar_mass_seasalt_si = 0.0314038218d0  !MOLAR MASS OF SEA SALT IN KG/MOL, MILLERO ET AL. 2008
real*8, public, parameter :: molar_mass_air_si = 0.02896546d0        !MOLAR MASS OF DRY AIR IN KG MOL-1,  PICARD ET AL. 2008

!*****************************************************************************************************************************
!LEMMON ET AL. (2000) DECIDED TO USE THE LATEST MOLAR MASS WITH THEIR MOLAR FORMULATION.
!FOR CONSISTENCY WITH IAPWS-10, THE SIA LIBRARY VERSION 1.1 INTRODUCED IN May 2010 USES THIS VALUE

!VERSION 1.0 (V1.0)
!real*8, public, parameter :: molar_mass_air_L2000 = 0.0289586d0      !MOLAR MASS OF DRY AIR IN KG MOL-1,  LEMMON ET AL. 2000

!VERSION 1.1 (V1.1)
real*8, public, parameter :: molar_mass_air_L2000 = molar_mass_air_si !MOLAR MASS OF DRY AIR IN KG MOL-1,  LEMMON ET AL. 2000
!*****************************************************************************************************************************

real*8, public, parameter :: gas_constant_h2o_si = gas_constant_molar_si / molar_mass_h2o_si            !SPECIFIC GAS CONSTANT OF H2O IN J KG-1 K-1, IAPWS 2005


real*8, public, parameter :: gas_constant_h2o_iapws95 = 461.51805d0                                     !SPECIFIC GAS CONSTANT OF H2O IN J KG-1 K-1, IAPWS-95


real*8, public, parameter :: gas_constant_air_si = gas_constant_molar_si / molar_mass_air_si           !SPECIFIC GAS CONSTANT OF AIR IN J KG-1 K-1, PICARD ET AL. 2008
real*8, public, parameter :: gas_constant_air_L2000 = gas_constant_molar_L2000 / molar_mass_air_L2000  !SPECIFIC GAS CONSTANT OF AIR IN J KG-1 K-1, LEMMON ET AL. 2000

real*8, public, parameter :: sealevel_pressure_si = 101325d0                   !PA
real*8, public, parameter :: celsius_temperature_si = 273.15d0                 !K,          ITS-90 DEFINITION

!STANDARD OCEAN PROPERTIES,MILLERO ET AL. 2008
real*8, public, parameter :: so_salinity_si = 0.03516504d0                     !KG KG-1 OF KCL-NORMALISED SEAWATER, Millero et al. 2008
real*8, public, parameter :: so_temperature_si = celsius_temperature_si        !K
real*8, public, parameter :: so_pressure_si = sealevel_pressure_si             !PA

!CRITICAL POINT OF PURE WATER, IAPWS-95
real*8, public, parameter :: cp_density_si = 322d0                             !KG M-3
real*8, public, parameter :: cp_temperature_si = 647.096d0                     !K
real*8, public, parameter :: cp_pressure_si = 22064000d0                       !PA
real*8, public, parameter :: cp_chempot_si = -767471.156792841d0               !IN J KG-1

!TRIPLE POINT OF PURE WATER
real*8, public, parameter :: tp_temperature_si = 273.16d0                      !K   ITS-90 DEFINITION
real*8, public, parameter :: tp_pressure_exp_si = 611.657d0                    !PA  EXPERIMENTAL TRIPLE POINT PRESSURE
!REAL*8, PUBLIC, PARAMETER :: TP_PRESSURE_IAPWS95_SI =611.6547710078944264442598D-6    !FROM QUADRUPLE-PRECISION COMPUTATION
real*8, public, parameter :: tp_pressure_iapws95_si = 611.654771007894d0       !Pa
!REAL*8, PUBLIC, PARAMETER :: TP_DENSITY_VAP_IAPWS95_SI =4.854575724778588417176210D-3 !FROM QUADRUPLE-PRECISION COMPUTATION
real*8, public, parameter :: tp_density_vap_iapws95_si = 4.85457572477859d-03  !KG M-3
!REAL*8, PUBLIC, PARAMETER :: TP_DENSITY_LIQ_IAPWS95_SI =999.792520031620646603898354735d0 !FROM QUADRUPLE-PRECISION COMPUTATION
real*8, public, parameter :: tp_density_liq_iapws95_si = 999.792520031621d0    !KG M-3
real*8, public, parameter :: tp_density_ice_iapws95_si = 916.709492199488d0    !KG M-3
real*8, public, parameter :: tp_enthalpy_vap_si = 2500915.1914657d0            !J KG-1
real*8, public, parameter :: tp_enthalpy_ice_si = -333444.253967839d0          !J KG-1


!References:

!FEISTEL ET AL. 2008:
!FEISTEL, R., WRIGHT, D.G., MIYAGAWA, K., HARVEY, A.H., HRUBY, j., JACKETT, D.R., MCdOUGALL, T.J. AND WAGNER, W.:
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE, AND SEAWATER: a NEW STANDARD FOR OCEANOGRAPHY,
!OCEAN SCI., 4, 275-291, 2008, AVAILABLE AT: WWW.OCEAN-SCI.NET/4/275/2008/

!IAPWS-95:
!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!DOORWERTH, THE NETHERLANDS, SEPTEMBER 2009
!(ORIGINAL RELEASE: FREDERICIA, DENMARK, SEPTEMBER 1996)
!AVAILABLE AT:  HTTP://WWW.IAPWS.ORG

!IAPWS 2005:
!IAPWS: GUIDELINE ON THE USE OF FUNDAMENTAL PHYSICAL CONSTANTS AND BASIC CONSTANTS OF WATER.
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM. GAITHERSBURG, MAYLAND, USA, SEPTEMBER 2001.
!REVISION JULY 2005, AVAILABLE AT: HTTP://WWW.IAPWS.ORG

!LEMMON ET AL. 2000:
!LEMMON, E.W., JACOBSEN, R.T., PENONCELLO, S.G. AND FRIEND, D.G.:
!THERMODYNAMIC PROPERTIES OF AIR AND MIXTURES OF NITROGEN, ARGON AND OXYGEN FROM 60 TO 2000 K AT PRESSURES TO 2000 MPa,
!J. PHYS. CHEM. REF. DATA, 29, 331-362, 2000

!MILLERO ET AL. 2008:
!MILLERO, F.J., FEISTEL, R., WRIGHT, D.G. AND MCDOUGALL, T.J.:
!THE COMPOSITION OF STANDARD SEAWATER AND THE DEFINITION OF THE REFERENCE-COMPOSITION SALINITY SCALE,
!DEEP-SEA RES. I, 55, 50-72, 2008

!PICARD ET AL. 2008:
!PICARD, A., DAVIS, R.S., GLASER, M. AND FUJII, K.:
!REVISED FORMULA FOR THE DENSITY OF MOIST AIR (CIPM-2007),
!METROLOGIA, 45, 149-155, 2008

end module constants_0
