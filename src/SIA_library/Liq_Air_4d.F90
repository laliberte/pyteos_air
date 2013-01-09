module liq_air_4d

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS FUNCTION OF WET AIR, I.E., OF THE COMPOSITE
!SYSTEM OF LIQUID WATER AND HUMID AIR AT MUTUAL EQUILIBRIUM ("CLOUDY AIR")
!THEREFORE, THE HUMID-AIR PROPERTIES COMPUTED HERE REFER TO SATURATED AIR.

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
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     AIR_3A,        FILE AIR_3A.F90
!     FLU_3A,        FILE FLU_3A.F90
!     LIQ_AIR_4A,    FILE LIQ_AIR_4A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,     FILE CONVERT_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     AIR_1,         FILE AIR_1.F90
!     FLU_1,         FILE FLU_1.F90
!     FLU_2,         FILE FLU_2.F90
!     AIR_2,         FILE AIR_2.F90
!     AIR_3B,        FILE AIR_3A.F90
!     LIQ_VAP_4,     FILE LIQ_VAP_4.F90

!#########################################################################

use constants_0
use air_3a
use flu_3a
use liq_air_4a

implicit none
private

character*16, private :: version = 'Oct 2012'

public :: liq_air_g_cond_si, liq_air_g_cond_cp_si, liq_air_g_cond_density_si, &
          liq_air_g_cond_enthalpy_si,  liq_air_g_cond_entropy_si, &
          liq_air_g_cond_expansion_si, liq_air_g_cond_kappa_t_si,  &
          liq_air_g_cond_lapserate_si, liq_air_cond_liquidfraction_si, &
          liq_air_cond_vapourfraction_si

contains

!==========================================================================
function liq_air_g_cond_si(drv_wa, drv_t, drv_p, wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!THIS IMPLEMENTS THE GIBBS FUNCTION OF WET AIR COMPUTED FROM ITS HELMHOLTZ FUNCTION
!AND ITS FIRST AND SECOND DERIVATIVES WITH RESPECT TO THE AIR FRACTION WA_SI,
!THE ABSOLUTE TEMPERATURE T_SI AND THE ABSOLUTE PRESSURE P_SI

!NOTE: THIS GIBBS FUNCTION IS NOT DESIGNED FOR LIQUID AIR OR CRITICAL AIR STATES
!NOTE: THE ACCURACY OF THIS FUNCTIONS DEPENDS ON THE ITERATION SETTINGS
!IN AIR_3A TO COMPUTE THE GIBBS FUNCTION OF HUMID AIR AND ON THE SETTINGS
!IN FLU_3A TO COMPUTE THE GIBBS FUNCTION OF LIQUID WATER

!OUTPUT: 
!LIQ_AIR_G = SPECIFIC GIBBS ENERGY IN J/KG OR ITS DERIVATIVE

!G(WA,T,P) = (1-W) * GL(T,P) + W * GH(A(T,P),T,P)

!W =     MASS FRACTION OF HUMID AIR IN WET AIR, W(WA, T, P) = WA/A(T,P)
!(1-W) = LIQUID MASS FRACTION OF WET AIR
!A =     SATURATED AIR MASS FRACTION OF HUMID AIR
!(1-A) = SATURATED SPECIFIC HUMIDITY OF THE HUMID-AIR COMPONENT
!GL =    GIBBS FUNCTION OF LIQUID WATER
!GH =    GIBBS FUNCTION OF HUMID AIR

!INPUTS:
!DRV_WA = ORDER OF THE PARTIAL WA-DERIVATIVE
!DRV_T = ORDER OF THE PARTIAL T-DERIVATIVE
!DRV_P = ORDER OF THE PARTIAL P-DERIVATIVE
!WA_SI = ABSOLUTE AIR FRACTION IN KG/KG, I.E. THE MASS FRACTION OF DRY AIR IN WET AIR
!T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
!P_SI  = P = ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_SI( 0, 0, 0, 0.5, 300, 1E5) =-5397.43216342
!LIQ_AIR_G_SI( 1, 0, 0, 0.5, 300, 1E5) =-264.763415292
!LIQ_AIR_G_SI( 0, 1, 0, 0.5, 300, 1E5) =-343.818233768
!LIQ_AIR_G_SI( 0, 0, 1, 0.5, 300, 1E5) = 0.446835172414
!LIQ_AIR_G_SI( 2, 0, 0, 0.5, 300, 1E5) = 0
!LIQ_AIR_G_SI( 1, 1, 0, 0.5, 300, 1E5) = 98.4884000933
!LIQ_AIR_G_SI( 1, 0, 1, 0.5, 300, 1E5) = 0.891663433710
!LIQ_AIR_G_SI( 0, 2, 0, 0.5, 300, 1E5) =-14.1012539816
!LIQ_AIR_G_SI( 0, 1, 1, 0.5, 300, 1E5) = 2.43241086432E-03
!LIQ_AIR_G_SI( 0, 0, 2, 0.5, 300, 1E5) =-4.62469737987E-06

!CHECK VALUES WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_SI( 0, 0, 0, 0.5, 300, 1E5) =-5396.77820137
!LIQ_AIR_G_SI( 1, 0, 0, 0.5, 300, 1E5) =-263.45549120|3
!LIQ_AIR_G_SI( 0, 1, 0, 0.5, 300, 1E5) =-343.783393872
!LIQ_AIR_G_SI( 0, 0, 1, 0.5, 300, 1E5) = 0.446729465555
!LIQ_AIR_G_SI( 2, 0, 0, 0.5, 300, 1E5) = 0
!LIQ_AIR_G_SI( 1, 1, 0, 0.5, 300, 1E5) = 98.5580798842
!LIQ_AIR_G_SI( 1, 0, 1, 0.5, 300, 1E5) = 0.891452019991
!LIQ_AIR_G_SI( 0, 2, 0, 0.5, 300, 1E5) =-14.0995955397
!LIQ_AIR_G_SI( 0, 1, 1, 0.5, 300, 1E5) = 2.43183979422E-03
!LIQ_AIR_G_SI( 0, 0, 2, 0.5, 300, 1E5) =-4.62360294023E-06

real*8 liq_air_g_cond_si, wa_si, t_si, p_si, a_sat_si
real*8 g, a, w, d
real*8 gh, gl
real*8 gh_p, gl_p
real*8 gh_t, gl_t
real*8 gh_pp, gl_pp, gc_pp
real*8 gh_tp, gl_tp, gc_tp
real*8 gh_tt, gl_tt, gc_tt
real*8 a_p, a_t

integer drv_wa, drv_t, drv_p

liq_air_g_cond_si = errorreturn
d = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_wa < 0) return
if(drv_t < 0) return
if(drv_p < 0) return
if(drv_wa + drv_t + drv_p > 2) return

if(drv_wa > 1) then                        !G IS LINEAR IN WA_SI
  liq_air_g_cond_si = 0d0
  return
end if

a = liq_air_massfraction_air_si(t_si, p_si)    !SATURATED AIR MASS FRACTION
if(a == errorreturn) then
    if(liq_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif

if(a < 0d0 .or. a > 1d0) return
if(a<=a_sat_si) a=a_sat_si

if(a < wa_si) then                         !NO LIQUID FRACTION AT THIS T AND P
                                           !INSTEAD, INVOKE THE GIBBS FUNCTION 
                                           !OF UNSATURATED HUMID AIR
  !liq_air_g_cond_si = air_g_si(drv_wa, drv_t, drv_p, wa_si, t_si, p_si)  !dbg
  return
end if

w = wa_si / a                                 !GASEOUS MASS FRACTION OF WET AIR

if(drv_t + drv_p > 1) then
  d = a**2 * air_g_si(2, 0, 0, a, t_si, p_si) !AIR COEFFICIENT DA = A^2 * G_AA
end if

select case(drv_wa)

  case(0)

    select case(drv_t)

      case(0)     !D/DT

        select case(drv_p)

          case(0) 
            gh = air_g_si(0, 0, 0, a, t_si, p_si)
            if(gh == errorreturn) return
            gl = liq_g_si(0, 0, t_si, p_si)
            if(gl == errorreturn) return
            g = (1d0 - w) * gl + w * gh                      !g

          case(1)
            gh_p = air_g_si(0, 0, 1, a, t_si, p_si)
            if(gh_p == errorreturn) return
            gl_p = liq_g_si(0, 1, t_si, p_si)
            if(gl_p == errorreturn) return
            g = (1d0 - w) * gl_p + w * gh_p                  !V = G_P

          case(2) 
            a_p = liq_air_cond_a_si(0, 1, t_si, p_si, a_sat_si)
            if(a_p == errorreturn) return
            gc_pp = -d * a_p**2 / a                          !LATENT CONDENSATION DERIVATIVE
            gh_pp = air_g_si(0, 0, 2, a, t_si, p_si)
            if(gh_pp == errorreturn) return
            gl_pp = liq_g_si(0, 2, t_si, p_si)
            if(gl_pp == errorreturn) return
            g = (1d0 - w) * gl_pp + w * (gh_pp + gc_pp)      !G_PP

          case default; return
        end select

      case(1)     !D/DT

        select case(drv_p)

          case(0)
            gh_t = air_g_si(0, 1, 0, a, t_si, p_si)
            if(gh_t == errorreturn) return
            gl_t = liq_g_si(1, 0, t_si, p_si)
            if(gl_t == errorreturn) return
            g = (1d0 - w) * gl_t + w * gh_t                  !-ETA = G_T

          case(1)
            a_t = liq_air_cond_a_si(1, 0, t_si, p_si, a_sat_si)
            if(a_t == errorreturn) return
            a_p = liq_air_cond_a_si(0, 1, t_si, p_si, a_sat_si)
            if(a_p == errorreturn) return
            gc_tp = -d * a_t * a_p / a                       !LATENT DERIVATIVE
            gh_tp = air_g_si(0, 1, 1, a, t_si, p_si)
            if(gh_tp == errorreturn) return
            gl_tp = liq_g_si(1, 1, t_si, p_si)
            if(gl_tp == errorreturn) return
            g = (1d0 - w) * gl_tp + w * (gh_tp + gc_tp)      !G_TP

          case default; return
        end select

      case(2)     !D2/DT2

        select case(drv_p)

          case(0)
            a_t = liq_air_cond_a_si(1, 0, t_si, p_si, a_sat_si)
            if(a_t == errorreturn) return
            gc_tt = -d * a_t**2 / a                          !LATENT DERIVATIVE
            gh_tt = air_g_si(0, 2, 0, a, t_si, p_si)
            if(gh_tt == errorreturn) return
            gl_tt = liq_g_si(2, 0, t_si, p_si)
            if(gl_tt == errorreturn) return
            g = (1d0 - w) * gl_tt + w * (gh_tt + gc_tt)      !G_TT

          case default; return
        end select

       case default; return
    end select

  case(1)       !g_w, g_wp, g_wt
    gh = air_g_si(0, drv_t, drv_p, a, t_si, p_si)
    if (gh == errorreturn) return
    gl = liq_g_si(drv_t, drv_p, t_si, p_si)
    if (gl == errorreturn) return
    g = (gh - gl) / a

  case(2); g = 0d0
    
  case default; return

end select

liq_air_g_cond_si = g

end function

!=========================================================================
function liq_air_g_cond_cp_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!CP(WA,T,P) = T * (D2G/DT2) HEAT CAPACITY OF WET AIR IN J/(KG K)

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_CP_SI(0.5, 300, 1E5) = 4230.37619447

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_CP_SI(0.5, 300, 1E5) = 4229.87866191

real*8 liq_air_g_cond_cp_si, wa_si, t_si, p_si, a_sat_si
real*8 g_tt

liq_air_g_cond_cp_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tt = liq_air_g_cond_si(0, 2, 0, wa_si, t_si, p_si, a_sat_si)
if(g_tt == errorreturn) return

liq_air_g_cond_cp_si = -t_si * g_tt

end function

!=========================================================================
function liq_air_g_cond_density_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!D(WA,T,P) = 1/(DG/DP)  DENSITY OF WET AIR IN KG/M3

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_DENSITY_SI(0.5, 300, 1E5) = 2.23796169535

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_DENSITY_SI(0.5, 300, 1E5) = 2.23849125053

real*8 liq_air_g_cond_density_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p

liq_air_g_cond_density_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = liq_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

liq_air_g_cond_density_si = 1d0 / g_p

end function

!=========================================================================
function liq_air_g_cond_enthalpy_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!H(WA,T,P) = G - T * (DG/DT)  ENTHALPY OF WET AIR IN J/KG

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_ENTHALPY_SI(0.5, 300, 1E5) = 97748.0379670

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_ENTHALPY_SI(0.5, 300, 1E5) = 97738.2399604

real*8 liq_air_g_cond_enthalpy_si, wa_si, t_si, p_si, a_sat_si
real*8 g, g_t

liq_air_g_cond_enthalpy_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = liq_air_g_cond_si(0, 0, 0, wa_si, t_si, p_si, a_sat_si)
if(g == errorreturn) return
g_t = liq_air_g_cond_si(0, 1, 0, wa_si, t_si, p_si, a_sat_si)
if(g_t == errorreturn) return

liq_air_g_cond_enthalpy_si = g - t_si * g_t

end function

!=========================================================================
function liq_air_g_cond_entropy_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!ETA(WA,T,P) = - (DG/DT)  ENTROPY OF WET AIR IN J/(KG K)

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_ENTROPY_SI(0.5, 300, 1E5) = 343.818233768

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_ENTROPY_SI(0.5, 300, 1E5) = 343.783393872

real*8 liq_air_g_cond_entropy_si, wa_si, t_si, p_si, a_sat_si
real*8 g_t

liq_air_g_cond_entropy_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_t = liq_air_g_cond_si(0, 1, 0, wa_si, t_si, p_si, a_sat_si)
if(g_t == errorreturn) return

liq_air_g_cond_entropy_si = -g_t

end function

!=========================================================================
function liq_air_g_cond_expansion_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!ALPHA(WA,T,P) = (D2G/DTDP)/(DG/DP)  THERMAL EXPANSION OF WET AIR IN 1/K

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_EXPANSION_SI(0.5, 300, 1E5) = 5.44364234170E-03

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_EXPANSION_SI(0.5, 300, 1E5) = 5.44365210207E-03

real*8 liq_air_g_cond_expansion_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p, g_tp

liq_air_g_cond_expansion_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = liq_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_tp = liq_air_g_cond_si(0, 1, 1, wa_si, t_si, p_si, a_sat_si)
if(g_tp == errorreturn) return

liq_air_g_cond_expansion_si = g_tp / g_p

end function

!=========================================================================
function liq_air_g_cond_kappa_t_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!KAPPA_T(WA,T,P) = - (D2G/DP2)/(DG/DP)  ISOTHERMAL COMPRESSIBILITY OF WET AIR IN 1/PA

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_KAPPA_T_SI(0.5, 300, 1E5) = 1.03498955887E-05

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_KAPPA_T_SI(0.5, 300, 1E5) = 1.03498947276E-05

real*8 liq_air_g_cond_kappa_t_si, wa_si, t_si, p_si, a_sat_si
real*8 g_p, g_pp

liq_air_g_cond_kappa_t_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = liq_air_g_cond_si(0, 0, 1, wa_si, t_si, p_si, a_sat_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_pp = liq_air_g_cond_si(0, 0, 2, wa_si, t_si, p_si, a_sat_si)
if(g_pp == errorreturn) return

liq_air_g_cond_kappa_t_si = -g_pp / g_p

end function

!=========================================================================
function liq_air_g_cond_lapserate_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!GAMMA(WA,T,P) = - (D2G/DTDP)/(D2G/DT2) "MOIST" ADIABATIC LAPSE RATE OF WET AIR IN K/PA

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_G_LAPSERATE_SI(0.5, 300, 1E5) = 1.72496067902E-04

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_G_LAPSERATE_SI(0.5, 300, 1E5) = 1.72475854884E-04

real*8 liq_air_g_cond_lapserate_si, wa_si, t_si, p_si, a_sat_si
real*8 g_tp, g_tt

liq_air_g_cond_lapserate_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tp = liq_air_g_cond_si(0, 1, 1, wa_si, t_si, p_si, a_sat_si)
if(g_tp == errorreturn) return
if(g_tp <= 0d0) return
g_tt = liq_air_g_cond_si(0, 2, 0, wa_si, t_si, p_si, a_sat_si)
if(g_tt == errorreturn) return
if(g_tt >= 0d0) return

liq_air_g_cond_lapserate_si = -g_tp / g_tt

end function

!=========================================================================
function liq_air_cond_liquidfraction_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!WW = 1 - WA/A(T,P) MASS FRACTION OF LIQUID WATER IN WET AIR (0 < WW < 1)

!INPUT:
!WA_SI     MASS FRACTION OF DRY AIR IN WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_LIQUIDFRACTION_SI(0.5, 300, 1E5) = 0.488543697765

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_LIQUIDFRACTION_SI(0.5, 300, 1E5) = 0.488546404734

real*8 liq_air_cond_liquidfraction_si, wa_si, t_si, p_si, a_sat_si
real*8 a

liq_air_cond_liquidfraction_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

a = liq_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    if(liq_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif
if(a < wa_si) return
if(a<a_sat_si) a=a_sat_si

liq_air_cond_liquidfraction_si = 1d0 - wa_si / a

end function

!=========================================================================
function liq_air_cond_vapourfraction_si(wa_si, t_si, p_si, a_sat_si)
!==========================================================================

!OUTPUT:
!WV = WA*(1/A(T,P)-1) MASS FRACTION OF VAPOUR IN WET AIR (0 < WV < 1)

!INPUT:
!WA_SI     DRY-AIR MASS FRACTION OF WET AIR IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS, V1.0:
!LIQ_AIR_VAPOURFRACTION_SI(0.5, 300, 1E5) = 1.14563022353E-02

!CHECK VALUE WITH DEFAULT SETTINGS, V1.1:
!LIQ_AIR_VAPOURFRACTION_SI(0.5, 300, 1E5) = 1.14535952655E-02

real*8 liq_air_cond_vapourfraction_si, wa_si, t_si, p_si, a_sat_si
real*8 a

liq_air_cond_vapourfraction_si = errorreturn

if(wa_si < 0d0 .or. wa_si > 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

a = liq_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    if(liq_air_condensationpressure_si(wa_si,t_si)>p_si) then
        !no condensation but massfraction fails. We impose
        !our saturation pressure
        a=a_sat_si
    else
        return
    endif
endif
if(a < wa_si) return
if(a<a_sat_si) a=a_sat_si

liq_air_cond_vapourfraction_si = wa_si * (1d0 / a - 1d0)

end function

!=========================================================================
function liq_air_cond_a_si(drv_t, drv_p, t_si, p_si, a_sat_si)
!==========================================================================

!THIS FUNCTION RETURNS THE AIR FRACTION A OF SATURATED HUMID AIR IN KG/KG,
!AND ITS PARTIAL DERIVATIVES (D/DT)^DRV_T (D/DP)^DRV_P A(T,P)
!AS A FUNCTION OF ABSOLUTE TEMPERATURE, T_SI, IN K, AND ABSOLUTE PRESSURE, P_SI, IN PA

real*8 liq_air_cond_a_si, t_si, p_si, a_sat_si
real*8 d  !da coefficient
real*8 a  !air fraction

real*8 gh_p, gl_p, gh_ap
real*8 gh_t, gl_t, gh_at

integer drv_t, drv_p

liq_air_cond_a_si = errorreturn
d = errorreturn

if(drv_t < 0 .or. drv_t > 1) return
if(drv_p < 0 .or. drv_p > 1) return

if(t_si < 0d0) return
if(p_si < 0d0) return

a = liq_air_massfraction_air_si(t_si, p_si)
if(a == errorreturn) then
    a=a_sat_si
endif
if(a < 0d0 .or. a > 1d0) return
if(a<a_sat_si) a=a_sat_si

if(drv_t > 0 .or. drv_p > 0) then
  if(a==a_sat_si) then
    liq_air_cond_a_si=0d0
    return
  end if
  d = a**2 * air_g_si(2, 0, 0, a, t_si, p_si) !air coefficient da = a^2 * g_aa
  if(d == errorreturn) return
  if(d <= 0d0) return
end if

select case(drv_t)

  case(0)

    select case(drv_p)

      case(0)
        liq_air_cond_a_si = a

      case(1)
        gh_p = air_g_si(0, 0, 1, a, t_si, p_si)
        if(gh_p == errorreturn) return
        gh_ap = air_g_si(1, 0, 1, a, t_si, p_si)
        if(gh_ap == errorreturn) return
        gl_p = liq_g_si(0, 1, t_si, p_si)
        if(gl_p == errorreturn) return
        liq_air_cond_a_si = a * (gh_p - gl_p - a * gh_ap) / d

    end select

  case(1)

    select case(drv_p)

      case(0)
        gh_t = air_g_si(0, 1, 0, a, t_si, p_si)
        if(gh_t == errorreturn) return
        gh_at = air_g_si(1, 1, 0, a, t_si, p_si)
        if(gh_at == errorreturn) return
        gl_t = liq_g_si(1, 0, t_si, p_si)
        if(gl_t == errorreturn) return
        liq_air_cond_a_si = a * (gh_t - gl_t - a * gh_at) / d
    
    end select

end select

end function

end module liq_air_4d

