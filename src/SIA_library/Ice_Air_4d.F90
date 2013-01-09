module ice_air_4d

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS FUNCTION OF ICE AIR, I.E., OF THE COMPOSITE
!SYSTEM OF ICE AND HUMID AIR IN MUTUAL EQUILIBRIUM ("ICY AIR", E.G. CIRRUS CLOUDS)
!THEREFORE, THE AIR PROPERTIES COMPUTED HERE REFER TO SATURATED AIR.

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., WRIGHT, D.G., HERRMANN, S.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     CONVERT_0,    FILE CONVERT_0.F90
!     ICE_1,        FILE ICE_1.F90
!     AIR_3A,       FILE AIR_3A.F90
!     ICE_AIR_4A,   FILE ICE_AIR_4A.F90

!NESTED USE ASSOCIATIONS
!     MATHS_0,      FILE MATHS_0.F90
!     AIR_1,        FILE AIR_1.F90
!     FLU_1,        FILE FLU_1.F90
!     AIR_2,        FILE AIR_2.F90
!     FLU_2,        FILE FLU_2.F90
!     ICE_2,        FILE ICE_2.F90
!     AIR_3B,       FILE AIR_3B.F90
!     ICE_VAP_4,    FILE ICE_VAP_4.F90

!#########################################################################

use constants_0
use convert_0
use ice_1
use air_3a
use ice_air_4a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: ice_air_g_cond_si, ice_air_g_cond_cp_si, ice_air_g_cond_density_si, &
          ice_air_g_cond_enthalpy_si, ice_air_g_cond_entropy_si, &
          ice_air_g_cond_expansion_si, ice_air_g_cond_kappa_t_si,  &
          ice_air_g_cond_lapserate_si, ice_air_cond_solidfraction_si, &
          ice_air_cond_vapourfraction_si

contains

!==========================================================================
function ice_air_g_cond_si(drv_wa, drv_t, drv_p, wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!THIS IMPLEMENTS THE GIBBS FUNCTION OF ICE AIR COMPUTED FROM ITS HELMHOLTZ FUNCTION
!AND ITS FIRST AND SECOND DERIVATIVES WITH RESPECT TO THE AIR FRACTION WA_SI,
!THE ABSOLUTE TEMPERATURE T_SI AND THE ABSOLUTE PRESSURE P_SI

!NOTE: THIS GIBBS FUNCTION IS NOT DESIGNED FOR LIQUID OR CRITICAL AIR STATES
!NOTE: THE ACCURACY OF THIS FUNCTIONS DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!OUTPUT: ICE_AIR_G_SI = SPECIFIC GIBBS ENERGY IN J/KG OR ITS DERIVATIVE

!G(WA,T,P) = (1-W) * GI(T,P) + W * GH(A(T,P),T,P)

!W = MASS FRACTION OF HUMID AIR IN ICE AIR, W(WA, T, P) = WA/A(T,P)
!(1-W) = SOLID MASS FRACTION OF ICE AIR
!A = SATURATED AIR MASS FRACTION OF HUMID AIR
!(1-A) = SATURATED SPECIFIC HUMIDITY OF THE HUMID-AIR COMPONENT
!GI = GIBBS FUNCTION OF ICE
!GH = GIBBS FUNCTION OF HUMID AIR

!WA_SI = ABSOLUTE DRY-AIR FRACTION IN KG/KG, I.E. THE MASS FRACTION OF DRY AIR IN ICE AIR
!T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
!P_SI  = P = ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!ice_air_g_cond_si( 0, 0, 0, 0.5, 270, 1E5) =-2595.73785825
!ice_air_g_cond_si( 1, 0, 0, 0.5, 270, 1E5) = 2382.02354606
!ice_air_g_cond_si( 0, 1, 0, 0.5, 270, 1E5) = 610.261631357
!ice_air_g_cond_si( 0, 0, 1, 0.5, 270, 1E5) = 0.389737675174
!ice_air_g_cond_si( 2, 0, 0, 0.5, 270, 1E5) = 0
!ice_air_g_cond_si( 1, 1, 0, 0.5, 270, 1E5) =-1269.42344741
!ice_air_g_cond_si( 1, 0, 1, 0.5, 270, 1E5) = 0.777294756075
!ice_air_g_cond_si( 0, 2, 0, 0.5, 270, 1E5) =-7.00885819160
!ice_air_g_cond_si( 0, 1, 1, 0.5, 270, 1E5) = 1.60133877466E-03
!ice_air_g_cond_si( 0, 0, 2, 0.5, 270, 1E5) =-3.91271268873E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!ice_air_g_cond_si( 0, 0, 0, 0.5, 270, 1E5) =-2595.57166634
!ice_air_g_cond_si( 1, 0, 0, 0.5, 270, 1E5) = 2382.35592988
!ice_air_g_cond_si( 0, 1, 0, 0.5, 270, 1E5) = 610.264515318
!ice_air_g_cond_si( 0, 0, 1, 0.5, 270, 1E5) = 0.389645501224
!ice_air_g_cond_si( 2, 0, 0, 0.5, 270, 1E5) = 0
!ice_air_g_cond_si( 1, 1, 0, 0.5, 270, 1E5) =-1269.41767949
!ice_air_g_cond_si( 1, 0, 1, 0.5, 270, 1E5) = 0.777110408175
!ice_air_g_cond_si( 0, 2, 0, 0.5, 270, 1E5) =-7.00810930740
!ice_air_g_cond_si( 0, 1, 1, 0.5, 270, 1E5) = 1.60095965101E-03
!ice_air_g_cond_si( 0, 0, 2, 0.5, 270, 1E5) =-3.91178603885E-06


real*8 ice_air_g_cond_si, wa_si, t_si, p_si, a_sat_si
real*8 g, a, w, d
real*8 gh, gi
real*8 gh_p, gi_p
real*8 gh_t, gi_t
real*8 gh_pp, gi_pp, gc_pp
real*8 gh_tp, gi_tp, gc_tp
real*8 gh_tt, gi_tt, gc_tt

real*8 a_p, a_t

integer drv_wa, drv_t, drv_p

ice_air_g_cond_si = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_wa < 0) return
if(drv_t < 0 .or. drv_t > 2) return
if(drv_p < 0 .or. drv_p > 2) return

if(drv_wa > 1) then                        !G IS LINEAR IN WA_SI
  ice_air_g_cond_si = 0d0
  return
end if

a = ice_air_massfraction_air_si(t_si, p_si)    !SATURATED AIR MASS FRACTION
if(a == errorreturn) then
    if(ice_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif
if(a < 0d0 .or. a > 1d0) return
if(a < a_sat_si) a = a_sat_si

!NOTE: A <= WA_SI IS NECESSARY FOR VALIDITY, HOWEVER, THIS TEST MAY FAIL
!DUE TO INSUFFICIENT ACCURACY OF ESTIMATED INITIAL ITERATION VALUES
!IF THE DIFFERENCE BETWEEN A AND WA_SI IS VERY SMALL
!THEREFORE, THE FOLLOWING CHECK IS SUPPRESSED (22 JAN 2009):
!IF(A < WA_SI) RETURN

w = wa_si / a                                 !GASEOUS MASS FRACTION OF ICE AIR

if(drv_t + drv_p > 1) then
  d = a**2 * air_g_si(2, 0, 0, a, t_si, p_si) !AIR COEFFICIENT DA = A^2 * G_AA
end if

select case(drv_wa)

  case(0)

    select case(drv_t)

      case(0) 

        select case(drv_p)

          case(0)   
            gh = air_g_si(0, 0, 0, a, t_si, p_si)
            if(gh == errorreturn) return
            gi = ice_g_si(0, 0, t_si, p_si)
            if(gi == errorreturn) return
            g = (1d0 - w) * gi + w * gh                        !G

          case(1)
            gh_p = air_g_si(0, 0, 1, a, t_si, p_si)
            if(gh_p == errorreturn) return
            gi_p = ice_g_si(0, 1, t_si, p_si)
            if(gi_p == errorreturn) return
            g = (1d0 - w) * gi_p + w * gh_p                    !V = G_P

          case(2)
            a_p = ice_air_cond_a_si(0, 1, t_si, p_si, a_sat_si)
            if(a_p == errorreturn) return
            gc_pp = -d * a_p**2 / a                            !LATENT CONDENSATION DERIVATIVE
            gh_pp = air_g_si(0, 0, 2, a, t_si, p_si)
            if(gh_pp == errorreturn) return
            gi_pp = ice_g_si(0, 2, t_si, p_si)
            if(gi_pp == errorreturn) return
            g = (1d0 - w) * gi_pp + w * (gh_pp + gc_pp)        !G_PP

          case default; return

        end select

      case(1)     !D/DT

        select case(drv_p)

          case(0) 
            gh_t = air_g_si(0, 1, 0, a, t_si, p_si)
            if(gh_t == errorreturn) return
            gi_t = ice_g_si(1, 0, t_si, p_si)
            if(gi_t == errorreturn) return
            g = (1d0 - w) * gi_t + w * gh_t              !-ETA = G_T

          case(1) 
            a_t = ice_air_cond_a_si(1, 0, t_si, p_si, a_sat_si)
            if(a_t == errorreturn) return
            a_p = ice_air_cond_a_si(0, 1, t_si, p_si, a_sat_si)
            if(a_p == errorreturn) return
            gc_tp = -d * a_t * a_p / a                   !LATENT DERIVATIVE
            gh_tp = air_g_si(0, 1, 1, a, t_si, p_si)
            if(gh_tp == errorreturn) return
            gi_tp = ice_g_si(1, 1, t_si, p_si)
            if(gi_tp == errorreturn) return
            g = (1d0 - w) * gi_tp + w * (gh_tp + gc_tp)  !G_TP

          case default; return

        end select

      case(2)     !D2/DT2

        select case(drv_p)

          case(0)
            a_t = ice_air_cond_a_si(1, 0, t_si, p_si, a_sat_si)
            if(a_t == errorreturn) return
            gc_tt = -d * a_t**2 / a                      !LATENT DERIVATIVE
            gh_tt = air_g_si(0, 2, 0, a, t_si, p_si)
            if(gh_tt == errorreturn) return
            gi_tt = ice_g_si(2, 0, t_si, p_si)
            if(gi_tt == errorreturn) return
            g = (1d0 - w) * gi_tt + w * (gh_tt + gc_tt)  !G_TT

          case default; return

        end select

       case default; return

    end select

  case(1)       !g_w, g_wp, g_wt
    gh = air_g_si(0, drv_t, drv_p, a, t_si, p_si)
    if( gh == errorreturn) return
    gi = ice_g_si(drv_t, drv_p, t_si, p_si)
    if( gi == errorreturn) return
    g = (gh - gi) / a

  case(2); g = 0     !G_WW

  case default; return

end select

ice_air_g_cond_si = g

end function

!=========================================================================
function ice_air_g_cond_cp_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!CP(WA,T,P) = T * (D2G/DT2) HEAT CAPACITY OF ICE AIR IN J/(KG K)

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_cp_si(0.5, 270, 1E5) = 1892.39171173

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_cp_si(0.5, 270, 1E5) = 1892.18951300

real*8 ice_air_g_cond_cp_si, wa_si, t_si, p_si, a_sat_si
real*8 g_tt

ice_air_g_cond_cp_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tt = ice_air_g_cond_si(0, 2, 0, wa_si, t_si, p_si, a_sat_si)
if(g_tt == errorreturn) return

ice_air_g_cond_cp_si = -t_si * g_tt

end function

!=========================================================================
function ice_air_g_cond_density_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!D(WA,T,P) = 1/(DG/DP)_T  DENSITY OF ICE AIR IN KG/M3

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_density_si(0.5, 270, 1E5) = 2.56582841152

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_density_si(0.5, 270, 1E5) = 2.56643538000

real*8 ice_air_g_cond_density_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p

ice_air_g_cond_density_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = ice_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

ice_air_g_cond_density_si = 1d0 / g_p

end function

!=========================================================================
function ice_air_g_cond_enthalpy_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!H(WA,T,P) = G - T * (DG/DT)  ENTHALPY OF ICE AIR IN J/KG

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_enthalpy_si(0.5, 270, 1E5) = -167366.378325

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_enthalpy_si(0.5, 270, 1E5) = -167366.990802

real*8 ice_air_g_cond_enthalpy_si, wa_si, t_si, p_si, a_sat_si
real*8 g, g_t

ice_air_g_cond_enthalpy_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = ice_air_g_cond_si(0, 0, 0, wa_si, t_si, p_si, a_sat_si)
if(g == errorreturn) return
g_t = ice_air_g_cond_si(0, 1, 0, wa_si, t_si, p_si, a_sat_si)
if(g_t == errorreturn) return

ice_air_g_cond_enthalpy_si = g - t_si * g_t

end function

!=========================================================================
function ice_air_g_cond_entropy_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!ETA(WA,T,P) = - (DG/DT)  ENTROPY OF ICE AIR IN J/(KG K)

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_entropy_si(0.5, 270, 1E5) = -610.261631357

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_entropy_si(0.5, 270, 1E5) = -610.264515318

real*8 ice_air_g_cond_entropy_si, wa_si, t_si, p_si, a_sat_si
real*8 g_t

ice_air_g_cond_entropy_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0) return
if(p_si <= 0) return

g_t = ice_air_g_cond_si(0, 1, 0, wa_si, t_si, p_si, a_sat_si)
if(g_t == errorreturn) return

ice_air_g_cond_entropy_si = -g_t

end function

!=========================================================================
function ice_air_g_cond_expansion_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!ALPHAI(WA,T,P) = (D2G/DTDP)/(DG/DP)_T  THERMAL EXPANSION OF ICE AIR IN 1/K

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_expansion_si(0.5, 270, 1E5) = 4.10876052449E-03

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_expansion_si(0.5, 270, 1E5) = 4.10875949031E-03

real*8 ice_air_g_cond_expansion_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p, g_tp

ice_air_g_cond_expansion_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = ice_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_tp = ice_air_g_cond_si(0, 1, 1, wa_si, t_si, p_si, a_sat_si)
if(g_tp == errorreturn) return

ice_air_g_cond_expansion_si = g_tp / g_p

end function

!=========================================================================
function ice_air_g_cond_kappa_t_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!KAPPA_T(WA,T,P) = - (D2G/DP2)/(DG/DP)_T  ISOTHERMAL COMPRESSIBILITY OF ICE AIR IN 1/PA

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_kappa_t_si(0.5, 270, 1E5) = 1.00393493829E-05

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_kappa_t_si(0.5, 270, 1E5) = 1.00393460891E-05

real*8 ice_air_g_cond_kappa_t_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p, g_pp

ice_air_g_cond_kappa_t_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = ice_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_pp = ice_air_g_cond_si(0, 0, 2, wa_si, t_si, p_si, a_sat_si)
if(g_pp == errorreturn) return

ice_air_g_cond_kappa_t_si = -g_pp / g_p

end function

!=========================================================================
function ice_air_g_cond_lapserate_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!GAMMA(WA,T,P) = - (D2G/DTDP)/(D2G/DT2)_P "MOIST" ADIABATIC LAPSE RATE OF ICE AIR IN K/PA

!INPUTS:
!WA_SI     MASS FRACTION OF DRY AIR IN ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_g_cond_lapserate_si(0.5, 270, 1E5) = 2.28473558871E-04

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_g_cond_lapserate_si(0.5, 270, 1E5) = 2.28443875629E-04

real*8 ice_air_g_cond_lapserate_si, wa_si, t_si, p_si, a_sat_si
real*8 g_tp, g_tt

ice_air_g_cond_lapserate_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tp = ice_air_g_cond_si(0, 1, 1, wa_si, t_si, p_si, a_sat_si)
if(g_tp == errorreturn) return
if(g_tp <= 0d0) return
g_tt = ice_air_g_cond_si(0, 2, 0, wa_si, t_si, p_si, a_sat_si)
if(g_tt == errorreturn) return
if(g_tt >= 0d0) return

ice_air_g_cond_lapserate_si = -g_tp / g_tt

end function

!=========================================================================
function ice_air_cond_solidfraction_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!1 - W(WA, T, P) = 1 - WA/A(T,P) MASS FRACTION OF SOLID WATER IN ICE AIR (0 < W < 1)

!INPUTS:
!WA_SI     DRY AIR MASS FRACTION OF ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_solidfraction_si(0.5, 270, 1E5) = 0.498524741260

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_solidfraction_si(0.5, 270, 1E5) = 0.498525089434

real*8 ice_air_cond_solidfraction_si, wa_si, t_si, p_si, a_sat_si
real*8 a

ice_air_cond_solidfraction_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

a = ice_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    if(ice_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif
if(a < wa_si) return
if(a < a_sat_si) a = a_sat_si

ice_air_cond_solidfraction_si = 1d0 - wa_si / a

end function

!=========================================================================
function ice_air_cond_vapourfraction_si(wa_si, t_si, p_si, a_sat_si)
!=========================================================================

!OUTPUT:
!WV = WA*/(1/A(T,P)-1) MASS FRACTION OF VAPOUR IN ICE AIR (0 < W < 1)

!INPUTS:
!WA_SI     DRY-AIR MASS FRACTION OF ICE AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0: 
!ice_air_cond_vapourfraction_si(0.5, 270, 1E5) = 1.47525873976E-03

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1: 
!ice_air_cond_vapourfraction_si(0.5, 270, 1E5) = 1.47491056602E-03

real*8 ice_air_cond_vapourfraction_si, wa_si, t_si, p_si, a_sat_si
real*8 a

ice_air_cond_vapourfraction_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

a = ice_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    if(ice_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif
if(a < wa_si) return
if(a < a_sat_si) a = a_sat_si

ice_air_cond_vapourfraction_si = wa_si * (1d0 / a - 1d0)

end function

!=========================================================================
function ice_air_cond_a_si(drv_t, drv_p, t_si, p_si, a_sat_si)
!=========================================================================

!THIS FUNCTION RETURNS THE AIR FRACTION A OF SATURATED HUMID AIR IN KG/KG,
!AND ITS PARTIAL DERIVATIVES (D/DT)^DRV_T (D/DP)^DRV_P A(T,P)
!AS A FUNCTION OF ABSOLUTE TEMPERATURE, T_SI, IN K, AND ABSOLUTE PRESSURE, P_SI, IN PA

real*8 ice_air_cond_a_si, t_si, p_si, a_sat_si
real*8 d  !da coefficient
real*8 a  !air fraction

real*8 gh_p, gi_p, gh_ap
real*8 gh_t, gi_t, gh_at

integer drv_t, drv_p

ice_air_cond_a_si = errorreturn
d = errorreturn

if(drv_t < 0 .or. drv_t > 1) return
if(drv_p < 0 .or. drv_p > 1) return

if(t_si < 0) return
if(p_si < 0) return

a = ice_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    a=a_sat_si
endif
if(a < 0d0 .or. a > 1d0) return
if(a < a_sat_si) a = a_sat_si

if(drv_t > 0 .or. drv_p > 0) then
  if(a==a_sat_si) then
    ice_air_cond_a_si=0d0
    return
  end if
  d = a**2 * air_g_si(2, 0, 0, a, t_si, p_si) !air coefficient da = a^2 * g_aa
  if(d == errorreturn) return
  if(d <= 0) return

end if

select case(drv_t)

  case(0)

    select case(drv_p)

      case(0) 
        ice_air_cond_a_si = a

      case(1)
        gh_p = air_g_si(0, 0, 1, a, t_si, p_si)
        if(gh_p == errorreturn) return
        gh_ap = air_g_si(1, 0, 1, a, t_si, p_si)
        if(gh_ap == errorreturn) return
        gi_p = ice_g_si(0, 1, t_si, p_si)
        if(gi_p == errorreturn) return
        ice_air_cond_a_si = a * (gh_p - gi_p - a * gh_ap) / d

    end select

  case(1)

    select case(drv_p)

      case(0)
        gh_t = air_g_si(0, 1, 0, a, t_si, p_si)
        if(gh_t == errorreturn) return
        gh_at = air_g_si(1, 1, 0, a, t_si, p_si)
        if(gh_at == errorreturn) return
        gi_t = ice_g_si(1, 0, t_si, p_si)
        if(gi_t == errorreturn) return
        ice_air_cond_a_si = a * (gh_t - gi_t - a * gh_at) / d

    end select

end select

end function

end module ice_air_4d

