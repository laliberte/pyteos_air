module sea_ice_4

!#########################################################################

!THIS MODULE IMPLEMENTS THE PHASE EQUILBRIA PROPERTIES OF ICE WITH SEAWATER
!COMPUTED FROM IAPWS-95, IAPWS-06 AND IAPWS-08:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996

!RELEASE ON AN EQUATION OF STATE FOR H2O ICE IH
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!WITNEY, UK, SEPTEMBER 2006

!RELEASE ON THE IAPWS FORMULATION 2008 FOR THE THERMODYNAMIC PROPERTIES OF SEAWATER
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!BERLIN, GERMANY, SEPTEMBER 2008

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
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     CONVERT_0,    FILE CONVERT_0.F90
!     MATHS_0,      FILE MATHS_0.F90
!     FLU_1,        FILE FLU_1.F90
!     ICE_1,        FILE ICE_1.F90
!     SAL_1,        FILE SAL_1.F90
!     ICE_2,        FILE ICE_2.F90
!     SAL_2,        FILE SAL_2.F90
!     FLU_3A,       FILE FLU_3A.F90
!     FLU_3B,       FILE FLU_3B.F90
!     SEA_3A,       FILE SEA_3A.F90

!NESTED USE ASSOCIATIONS
!     FLU_2,        FILE FLU_2.F90

!#########################################################################

use constants_0
use convert_0
use maths_0
use flu_1
use ice_1
use sal_1
use ice_2
use sal_2
use flu_3a
use flu_3b
use sea_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_ice_brinefraction_seaice_si, sea_ice_brinesalinity_si, &
          sea_ice_cp_seaice_si, sea_ice_density_ice_si, &
          sea_ice_density_sea_si, sea_ice_density_seaice_si, &
          sea_ice_dtfdp_si, sea_ice_dtfds_si, &
          sea_ice_enthalpy_ice_si, sea_ice_enthalpy_melt_si, &
          sea_ice_enthalpy_sea_si, sea_ice_enthalpy_seaice_si, &
          sea_ice_entropy_ice_si, sea_ice_entropy_sea_si, &
          sea_ice_entropy_seaice_si, sea_ice_expansion_seaice_si, &
          sea_ice_freezingtemperature_si, sea_ice_g_si, &
          sea_ice_kappa_t_seaice_si, sea_ice_meltingpressure_si, &
          sea_ice_pressure_si, sea_ice_salinity_si, &
          sea_ice_temperature_si, sea_ice_volume_melt_si, &
          set_it_ctrl_sea_ice, set_sea_ice_eq_at_s_p, &
          set_sea_ice_eq_at_s_t, set_sea_ice_eq_at_t_p

integer, private :: equi_sea_ice_done
integer, private :: ctrl_loop_maximum, ctrl_mode_sea, ctrl_mode_temperature, &
         ctrl_mode_liquid, ctrl_initialized
real*8, private :: ctrl_init_s_sea, ctrl_init_t, ctrl_init_d_liq, &
        ctrl_eps_exit_s, ctrl_eps_exit_t, ctrl_eps_exit_p
real*8, private :: equi_sea_ice_s, equi_sea_ice_t, equi_sea_ice_p, &
        equi_sea_ice_d_liq

contains

!==========================================================================
function sea_ice_brinesalinity_si(t_si, p_si)
!==========================================================================

implicit none
real*8 sea_ice_brinesalinity_si, t_si, p_si

!THIS FUNCTION RETURNS THE ABSOLUTE SALINITY IN KG/KG OF BRINE POCKETS AT THE
!ABSOLUTE TEMPERATURE T_SI IN K AND THE ABSOLUTE PRESSURE P_SI IN PA

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_BRINESALINITY_SI(270, 1E5) = 0.0560264150322

sea_ice_brinesalinity_si = errorreturn

if(set_sea_ice_eq_at_t_p(t_si, p_si) == errorreturn) return

sea_ice_brinesalinity_si = equi_sea_ice_s

end function

!=========================================================================
function sea_ice_brinefraction_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!W(SSI,T,P) = SSI/S(T,P)  MASS FRACTION OF BRINE IN SEA ICE (0 < W < 1)

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!p_si      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_BRINEFRACTION_SEAICE_SI(0.035, 270, 1E5) = 0.624705328368

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_brinefraction_seaice_si, ssi_si, t_si, p_si
real*8 s

sea_ice_brinefraction_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

s = sea_ice_brinesalinity_si(t_si, p_si)
if(s == errorreturn) return
if(s <= ssi_si) return

sea_ice_brinefraction_seaice_si = ssi_si / s

end function

!==========================================================================
function sea_ice_meltingpressure_si(sa_si, t_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE PRESSURE IN PA AT WHICH ICE MELTS AT THE
!ABSOLUTE TEMPERATURE T_SI IN K AND THE AMBIENT SEAWATER SALINITY SA_SI IN KG/KG

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_MELTINGPRESSURE_SI(0.035, 270) = 16132047.4385

implicit none
real*8 sea_ice_meltingpressure_si, sa_si, t_si

sea_ice_meltingpressure_si = errorreturn

if(set_sea_ice_eq_at_s_t(sa_si, t_si) == errorreturn) return

sea_ice_meltingpressure_si = equi_sea_ice_p

end function

!==========================================================================
function sea_ice_freezingtemperature_si(sa_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE TEMPERATURE IN K AT WHICH SEAWATER WITH
!ABSOLUTE SALINITY SA_SI IN KG/KG FREEZES AT THE ABSOLUTE PRESSURE P IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_FREEZINGTEMPERATURE_SI(0.035, 1E5) = 271.240373585

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_freezingtemperature_si, sa_si, p_si

sea_ice_freezingtemperature_si = errorreturn

if(set_sea_ice_eq_at_s_p(sa_si, p_si) == errorreturn) return

sea_ice_freezingtemperature_si = equi_sea_ice_t

end function

!==========================================================================
function sea_ice_dtfdp_si(sa_si, p_si)
!==========================================================================

!CLAUSIUS-CLAPEYRON: FREEZING POINT LOWERING DUE TO PRESSURE, DT/DP

!CHECKVALUE: SEA_ICE_DTDP_SI(0.035,1E5) = 7.48210942879E-08

real*8 sea_ice_dtfdp_si, sa_si, p_si
real*8 tf, s_t, s_p

sea_ice_dtfdp_si = ErrorReturn

tf = sea_ice_freezingtemperature_si(sa_si, p_si)
if(tf == errorreturn) return
s_t = sea_ice_s_si(1, 0, tf, p_si)
if(s_t == errorreturn .or. s_t == 0d0) return
s_p = sea_ice_s_si(0, 1, tf, p_si)
if(s_p == errorreturn) return

sea_ice_dtfdp_si = s_p / s_t

end function

!==========================================================================
function sea_ice_dtfds_si(sa_si, p_si)
!==========================================================================

!RAOULT'S LAW: FREEZING POINT LOWERING DUE TO SALT, DT/DS

!CHECKVALUE: SEA_ICE_DTDS_SI(0.035,1E5) = -56.8751336296

real*8 sea_ice_dtfds_si, sa_si, p_si
real*8 tf, s_t

sea_ice_dtfds_si = errorreturn

tf = sea_ice_freezingtemperature_si(sa_si, p_si)
if(tf == errorreturn) return
s_t = sea_ice_s_si(1, 0, tf, p_si)
if(s_t == errorreturn .or. s_t == 0) return

sea_ice_dtfds_si = 1d0 / s_t

End Function

!==========================================================================
function sea_ice_g_si(drv_s, drv_t, drv_p, s_si, t_si, p_si)
!==========================================================================

!OUTPUT: 
!THE S-T-P DERIVATIVE (D/DS)^DRV_S (D/DT)^DRV_T (D/DP)^DRV_P GSI(S,T,P)
!OF THE GIBBS FUNCTION GSI OF SEA ICE, AS A FUNCTION OF SEA-ICE SALINITY IN KG/KG,
!ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA

! GSI(S,T,P) = (1-W) * GIH(T,P) + W * GSW(SB(T,P),T,P)

! W = MASS FRACTION OF BRINE, W(S, T, P) = S/SB(T,P)
! SB = BRINE SALINITY

!INPUT:
!S_SI = S = ABSOLUTE SEA-ICE SALINITY IN KG/KG, I.E. THE MASS FRACTION OF SALT IN SEA ICE
!T_SI  = T = ABSOLUTE TEMPERATURE IN K, ITS-90
!P_SI  = P = ABSOLUTE PRESSURE IN PA
!SAL_G_SI  = SPECIFIC GIBBS ENERGY IN J/KG

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_ICE_G_SI( 0, 0, 0, 0.035, 270, 1E5) =-414.017574547
!SEA_ICE_G_SI( 1, 0, 0, 0.035, 270, 1E5) = 96363.7730495
!SEA_ICE_G_SI( 0, 1, 0, 0.035, 270, 1E5) = 500.44544418|1
!SEA_ICE_G_SI( 0, 0, 1, 0.035, 270, 1E5) = 1.00689072300E-03
!SEA_ICE_G_SI( 2, 0, 0, 0.035, 270, 1E5) = 0
!SEA_ICE_G_SI( 1, 1, 0, 0.035, 270, 1E5) =-1144.02883419
!SEA_ICE_G_SI( 1, 0, 1, 0.035, 270, 1E5) =-8.62856321467E-04
!SEA_ICE_G_SI( 0, 2, 0, 0.035, 270, 1E5) =-232.847783380
!SEA_ICE_G_SI( 0, 1, 1, 0.035, 270, 1E5) =-1.6586644669|4E-05
!SEA_ICE_G_SI( 0, 0, 2, 0.035, 270, 1E5) =-1.5759193211|8E-12

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
integer drv_t, drv_s, drv_p
real*8 sea_ice_g_si, s_si, t_si, p_si
real*8 g, s, w, d
real*8 gb, gi
real*8 gb_p, gi_p
real*8 gb_t, gi_t
real*8 gb_pp, gi_pp, gm_pp
real*8 gb_tp, gi_tp, gm_tp
real*8 gb_tt, gi_tt, gm_tt
real*8 gb_s, gb_sp, gb_st
real*8 s_p, s_t

sea_ice_g_si = errorreturn
g = errorreturn
d = errorreturn

if(s_si < 0d0 .or. s_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_s < 0) return
if(drv_t < 0) return
if(drv_p < 0) return
if(drv_s + drv_t +drv_p > 2) return

if(drv_s > 1) then       !G IS LINEAR IN S_SI
  sea_ice_g_si = 0d0
  return
end if

if(set_sea_ice_eq_at_t_p(t_si, p_si) == errorreturn) return

s = sea_ice_salinity_si()   !BRINE SALINITY
if(s == errorreturn) return
if(s <= 0d0 .or. s >= 1d0) return
if(s < s_si) return

w = s_si / s                !LIQUID MASS FRACTION OF SEA ICE

if(drv_t + drv_p > 1) then
  d = sal_dilution_si(s, t_si, p_si)
end if

select case (drv_s)

  case (0)
    select case (drv_t)

      case (0)     !D/DT
        select case (drv_p)

          case (0)
            gb = sea_g_si(0, 0, 0, s, t_si, p_si)
            if(gb == errorreturn) return
            gi = ice_g_si(0, 0, t_si, p_si)
            if(gi == errorreturn) return
            g = (1d0 - w) * gi + w * gb                    !G

          case (1) 
            gb_p = sea_g_si(0, 0, 1, s, t_si, p_si)
            if(gb_p == errorreturn) return
            gi_p = ice_g_si(0, 1, t_si, p_si)
            if(gi_p == errorreturn) return
            g = (1d0 - w) * gi_p + w * gb_p                !V = G_P

          case (2) 
            s_p = sea_ice_s_si(0, 1, t_si, p_si)
            if(s_p == errorreturn) return
            gm_pp = -d * s_p ** 2 / s                      !LATENT DERIVATIVE
            gb_pp = sea_g_si(0, 0, 2, s, t_si, p_si)
            if(gb_pp == errorreturn) return
            gi_pp = ice_g_si(0, 2, t_si, p_si)
            if(gi_pp == errorreturn) return
            g = (1d0 - w) * gi_pp + w * (gb_pp + gm_pp)    !G_PP

          case default; return
        end select

      case (1)     !D/DT

        select case (drv_p)

          case (0) 
            gb_t = sea_g_si(0, 1, 0, s, t_si, p_si)
            if(gb_t == errorreturn) return
            gi_t = ice_g_si(1, 0, t_si, p_si)
            if(gi_t == errorreturn) return
            g = (1d0 - w) * gi_t + w * gb_t                !-ETA = G_T

          case (1) 
            s_t = sea_ice_s_si(1, 0, t_si, p_si)
            if(s_t == errorreturn) return
            s_p = sea_ice_s_si(0, 1, t_si, p_si)
            if(s_p == errorreturn) return
            gm_tp = -d * s_t * s_p / s                     !LATENT DERIVATIVE
            gb_tp = sea_g_si(0, 1, 1, s, t_si, p_si)
            if(gb_tp == errorreturn) return
            gi_tp = ice_g_si(1, 1, t_si, p_si)
            if(gi_tp == errorreturn) return
            g = (1d0 - w) * gi_tp + w * (gb_tp + gm_tp)    !G_TP

          case default; return
        end select

      case (2)     !D2/DT2
        select case (drv_p)

          case (0) 
            s_t = sea_ice_s_si(1, 0, t_si, p_si)
            if(s_t == errorreturn) return
            gm_tt = -d * s_t ** 2 / s                      !LATENT DERIVATIVE
            gb_tt = sea_g_si(0, 2, 0, s, t_si, p_si)
            if(gb_tt == errorreturn) return
            gi_tt = ice_g_si(2, 0, t_si, p_si)
            if(gi_tt == errorreturn) return
            g = (1d0 - w) * gi_tt + w * (gb_tt + gm_tt)    !G_TT

          case default; return
        end select

       case default; return
    end select

  case (1)       !D/DS
  
    select case (drv_t)
    
      case (0)
        select case (drv_p)

          case (0) 
            g = sal_g_si(1, 0, 0, s, t_si, p_si)           !G_S

          case (1) 
            gb_sp = sal_g_si(1, 0, 1, s, t_si, p_si)       !G_SP
            if(gb_sp == errorreturn) return
            gb_s = sal_g_si(1, 0, 0, s, t_si, p_si)
            if(gb_s == errorreturn) return
            s_p = sea_ice_s_si(0, 1, t_si, p_si)
            if(s_p == errorreturn) return
            g = gb_sp + gb_s * s_p

          case default; return
        end select

      case (1)   !D2/DSDT

        select case (drv_p)

          case (0) 
            gb_st = sal_g_si(1, 1, 0, s, t_si, p_si)       !G_ST
            if(gb_st == errorreturn) return
            gb_s = sal_g_si(1, 0, 0, s, t_si, p_si)
            if(gb_s == errorreturn) return
            s_t = sea_ice_s_si(1, 0, t_si, p_si)
            if(s_t == errorreturn) return
            g = gb_st + gb_s * s_t

          case default; return
        end select

       case default; return
    end select
end select

sea_ice_g_si = g

end function

!=========================================================================
function sea_ice_cp_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!CP(SSI,T,P) = T * (D2G/DT2)_SSI_P ISOBARIC HEAT CAPACITY OF SEA ICE IN J/(KG K),

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_CP_SEAICE_SI(0.035, 270, 1E5) = 62868.9015126

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_cp_seaice_si, ssi_si, t_si, p_si
real*8 g_tt

sea_ice_cp_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_tt = sea_ice_g_si(0, 2, 0, ssi_si, t_si, p_si)
if(g_tt == errorreturn) return

sea_ice_cp_seaice_si = -t_si * g_tt

end function

!=========================================================================
function sea_ice_density_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!D(SSI,T,P) = 1/(DG/DP)_SSI_T  DENSITY OF SEA ICE IN KG/M3

!INPUT:
!ssi_si    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!p_si      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_DENSITY_SEAICE_SI(0.035, 270, 1E5) = 993.156434117

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_density_seaice_si, ssi_si, t_si, p_si
real*8 g_p

sea_ice_density_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_ice_g_si(0, 0, 1, ssi_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return

sea_ice_density_seaice_si = 1d0 / g_p

end function

!=========================================================================
function sea_ice_enthalpy_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!H(SSI,T,P) = G - T * (DG/DT)_SSI_P  ENTHALPY OF SEA ICE IN J/KG

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_ENTHALPY_SEAICE_SI(0.035, 270, 1E5) = -135534.287503

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_enthalpy_seaice_si, ssi_si, t_si, p_si
real*8 g, g_t

sea_ice_enthalpy_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g = sea_ice_g_si(0, 0, 0, ssi_si, t_si, p_si)
if(g == errorreturn) return
g_t = sea_ice_g_si(0, 1, 0, ssi_si, t_si, p_si)
if(g_t == errorreturn) return

sea_ice_enthalpy_seaice_si = g - t_si * g_t

end function

!=========================================================================
function sea_ice_entropy_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!ETA(SSI,T,P) = - (DGSI/DT)_SSI_P  ENTROPY OF SEA ICE IN J/(KG K)

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_ENTROPY_SEAICE_SI(0.035, 270, 1E5) = -500.44544418|1

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_entropy_seaice_si, ssi_si, t_si, p_si
real*8 g_t

sea_ice_entropy_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_t = sea_ice_g_si(0, 1, 0, ssi_si, t_si, p_si)
if(g_t == errorreturn) return

sea_ice_entropy_seaice_si = -g_t

end function

!=========================================================================
function sea_ice_expansion_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!ALPHA(SSI,T,P) = (D2G/DTDP)_SSA/(DG/DP)_SSI_T  THERMAL EXPANSION OF SEA ICE IN 1/K

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_EXPANSION_SEAICE_SI(0.035, 270, 1E5) = -1.64731328738

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_expansion_seaice_si, ssi_si, t_si, p_si
real*8 g_p, g_tp

sea_ice_expansion_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_ice_g_si(0, 0, 1, ssi_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_tp = sea_ice_g_si(0, 1, 1, ssi_si, t_si, p_si)
if(g_tp == errorreturn) return

sea_ice_expansion_seaice_si = g_tp / g_p

end function

!=========================================================================
function sea_ice_kappa_t_seaice_si(ssi_si, t_si, p_si)
!==========================================================================

!OUTPUT:   
!KAPPA_T(SSI,T,P) = - (D2G/DP2)_SSI_T/(DG/DP)_SSI_T  ISOTHERMAL COMPRESSIBILITY OF SEA ICE IN 1/PA,

!INPUT:
!SSI_SI    ABSOLUTE SEA-ICE SALINITY IN KG/KG
!T_SI      ABSOLUTE TEMPERATURE IN K
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_ICE_KAPPA_T_SEAICE_SI(0.035, 270, 1E5) = 1.56513441348E-09

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_kappa_t_seaice_si, t_si, p_si, ssi_si
real*8 g_p, g_pp

sea_ice_kappa_t_seaice_si = errorreturn

if(ssi_si < 0d0 .or. ssi_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

g_p = sea_ice_g_si(0, 0, 1, ssi_si, t_si, p_si)
if(g_p == errorreturn) return
if(g_p <= 0d0) return
g_pp = sea_ice_g_si(0, 0, 2, ssi_si, t_si, p_si)
if(g_pp == errorreturn) return

sea_ice_kappa_t_seaice_si = -g_pp / g_p

end function

!=========================================================================
function sea_ice_s_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS ABSOLUTE BRINE SALINITY S IN KG/KG, AND ITS PARTIAL DERIVATIVES
!(D/DT)^DRV_T (D/DP)^DRV_P S(T,P)
!AS A FUNCTION OF ABSOLUTE TEMPERATURE, T_SI, IN K, AND ABSOLUTE PRESSURE, P_SI, IN PA

!AS DEFINED IN EQ. (4.10) OF
!R. FEISTEL, E. HAGEN
!A GIBBS THERMODYNAMIC POTENTIAL OF SEA ICE
!COLD REGIONS SCIENCE AND TECHNOLOGY 28 (1998) 83, 142

implicit none
integer drv_t, drv_p
real*8 sea_ice_s_si, t_si, p_si
real*8 d  !DILUTION COEFF
real*8 s  !BRINE SALINITY

real*8 gb_p, gi_p, gb_sp
real*8 gb_t, gi_t, gb_st

sea_ice_s_si = errorreturn
d = errorreturn

if(drv_t < 0 .or. drv_t > 1) return
if(drv_p < 0 .or. drv_p > 1) return

if(t_si < 0d0) return
if(p_si < 0d0) return

!THIS MUST BE CALLED OUTSIDE BEFORE  !dbg comment is unclear
!IF(SET_SEA_ICE_EQUILIBRIUM_TO_T_P(T_SI, P_SI) == ERRORRETURN) RETURN

!IF NOT, EXIT HERE
if(equi_sea_ice_done /= isok) return

s = sea_ice_salinity_si()  !BRINE SALINITY
if(s == errorreturn) return
if(s < 0d0 .or. s >= 1d0) return

if(drv_t > 0 .or. drv_p > 0) then
  d = sal_dilution_si(s, t_si, p_si)
  if(d == errorreturn) return
  if(d <= 0d0) return
end if

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0) 
        sea_ice_s_si = s

      case (1) 
        gb_p = sea_g_si(0, 0, 1, s, t_si, p_si)
        if(gb_p == errorreturn) return
        gb_sp = sea_g_si(1, 0, 1, s, t_si, p_si)
        if(gb_sp == errorreturn) return
        gi_p = ice_g_si(0, 1, t_si, p_si)
        if(gi_p == errorreturn) return
        sea_ice_s_si = (gb_p - gi_p - s * gb_sp) / d
    end select

  case (1)

    select case (drv_p)

      case (0) 
        gb_t = sea_g_si(0, 1, 0, s, t_si, p_si)
        if(gb_t == errorreturn) return
        gb_st = sea_g_si(1, 1, 0, s, t_si, p_si)
        if(gb_st == errorreturn) return
        gi_t = ice_g_si(1, 0, t_si, p_si)
        if(gi_t == errorreturn) return
        sea_ice_s_si = (gb_t - gi_t - s * gb_st) / d
    end select

end select

end function

!==========================================================================
function set_sea_ice_eq_at_t_p(t_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE EQUILIBRIUM AT GIVEN ABSOLUTE TEMPERATURE T_SI 
!IN K AND ABSOLUTE PRESSURE P_SI IN PA BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_ice_eq_at_t_p, t_si, p_si
real*8 s_norm
real*8 d
real*8 s, eps

s_norm = so_salinity_si  !KG/KG OF KCL NORMALISED SEAWATER

if(equi_sea_ice_done == isok .and. &
  t_si == equi_sea_ice_t .and. p_si == equi_sea_ice_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_ice_eq_at_t_p = isok
  return
end if

call clear_sea_ice_state !WASTE ANY PREVIOUS STATE

set_sea_ice_eq_at_t_p = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_ice

!SET INITIAL BRINE SALINITY GUESS
select case (ctrl_mode_sea)
  case (0);  s = aux_brinesalinity_si(t_si, p_si)
  case (-1); s = s_norm
  case (1);  s = ctrl_init_s_sea
  case default; return
end select

if(s <= 0d0 .or. s >= 1d0) return
if(s == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  d = liq_density_si(t_si, p_si)
  call set_sea_ice_state (s, t_si, p_si, d)
  set_sea_ice_eq_at_t_p = isok
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else  
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_s
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_ice_iteration_at_t_p(t_si, p_si, maxit, eps, s) == errorreturn) then
  return
end if

d = liq_density_si(t_si, p_si)     !ONLY FOR COMPLETENESS OF THE STATE DESCRIPTOR
call set_sea_ice_state (s, t_si, p_si, d)

set_sea_ice_eq_at_t_p = isok

end function

!==========================================================================
function set_sea_ice_eq_at_s_p(sa_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE EQUILIBRIUM AT GIVEN ABSOLUTE SALINITY SA_SI
!IN KG/KG AND ABSOLUTE PRESSURE P_SI IN PA BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_ice_eq_at_s_p, sa_si, p_si
real*8 t, d, eps
real*8 tt, dt

tt = tp_temperature_si
dt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER

if(equi_sea_ice_done == isok .and. &
  sa_si == equi_sea_ice_s .and. p_si == equi_sea_ice_p) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_ice_eq_at_s_p = isok
  return
end if

call clear_sea_ice_state !WASTE ANY PREVIOUS STATE

set_sea_ice_eq_at_s_p = errorreturn

if(sa_si <= 0d0) return
if(sa_si >= 1d0) return
if(p_si <= 0d0) return

call init_it_ctrl_sea_ice

!SET INITIAL TEMPERATURE GUESS
select case (ctrl_mode_temperature)
  case (0);  t = aux_freezingtemperature_si(sa_si, p_si)
  case (-1); t = tt
  case (1);  t = ctrl_init_t
  case default; return
end select
if(t <= 0d0) return
if(t == errorreturn) return

!SET INITIAL LIQUID DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0);  d = aux_density_eos80_si(t, p_si)
  case (-1); d = dt
  case (1);  d = ctrl_init_d_liq
  case default; return
end select
if(d <= 0d0) return
if(d == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  call set_sea_ice_state (sa_si, t, p_si, d)
  set_sea_ice_eq_at_s_p = isok
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_t
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_ice_iteration_at_s_p(sa_si, p_si, maxit, eps, d, t) == errorreturn) then
  return
end if

call set_sea_ice_state (sa_si, t, p_si, d)

set_sea_ice_eq_at_s_p = isok

end function

!==========================================================================
function set_sea_ice_eq_at_s_t(sa_si, t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE EQUILIBRIUM AT GIVEN ABSOLUTE TEMPERATURE T_SI
!IN K AND ABSOLUTE SALINITY SA_SI IN KG/KG BY NEWTON ITERATION

implicit none
integer maxit
real*8 set_sea_ice_eq_at_s_t, sa_si, t_si
real*8 eps
real*8 d_si, p_si
real*8 dt

dt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER

if(equi_sea_ice_done == isok .and. &
  sa_si == equi_sea_ice_s .and. t_si == equi_sea_ice_t) then
  !THE REQUESTED STATE HAS ALREADY BEEN COMPUTED EARLIER
  set_sea_ice_eq_at_s_t = isok
  return
end if

call clear_sea_ice_state !WASTE ANY PREVIOUS STATE

set_sea_ice_eq_at_s_t = errorreturn

if(sa_si <= 0d0) return
if(sa_si >= 1d0) return
if(t_si <= 0d0) return

call init_it_ctrl_sea_ice

!SET INITIAL LIQUID WATER DENSITY GUESS
select case (ctrl_mode_liquid)
  case (0)
    d_si = aux_liq_density_correlation_si(t_si)
    if(d_si == errorreturn) return
    if(d_si <= 0d0) return
  case (-1); d_si = dt
  case (1);  d_si = ctrl_init_d_liq
  case default; return
end select
if(d_si <= 0d0) return
if(d_si == errorreturn) return

!SET MAX. ITERATION NUMBER
if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == (-1))then
  p_si = d_si ** 2 * flu_f_si(0, 1, t_si, d_si)
  call set_sea_ice_state (sa_si, t_si, p_si, d_si)
  set_sea_ice_eq_at_s_t = isok
  return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

!SET ITERATION ACCURACY LIMIT
eps = ctrl_eps_exit_p
if(eps == 0d0) return

!RUN ITERATION LOOP
if(sea_ice_iteration_at_s_t(sa_si, t_si, maxit, eps, d_si, p_si) == errorreturn) then
  return
end if

call set_sea_ice_state (sa_si, t_si, p_si, d_si)

set_sea_ice_eq_at_s_t = isok

end function

!==========================================================================
function sea_ice_density_ice_si()
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY IN KG/M3 OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_DENSITY_ICE_SI = 917.000739687

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_DENSITY_ICE_SI = 918.898527655

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_DENSITY_ICE_SI = 917.181167192

implicit none
real*8 sea_ice_density_ice_si

sea_ice_density_ice_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_density_ice_si = ice_density_si(equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_density_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY IN KG/M3 OF SEAWATER AT THE EQUILIBRIUM WITH ICE,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_DENSITY_SEA_SI = 1028.05199645

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_DENSITY_SEA_SI = 1035.73670169

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_DENSITY_SEA_SI = 1045.16805918

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_density_sea_si

sea_ice_density_sea_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_density_sea_si = sea_density_si(equi_sea_ice_s, equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_enthalpy_melt_si()
!==========================================================================

!THIS FUNCTION RETURNS THE MELTING ENTHALPY IN J/KG OF ICE IN SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_ENTHALPY_MELT_si = 329942.976285

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_ENTHALPY_MELT_SI = 326829.393605

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_ENTHALPY_MELT_SI = 328249.119579

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_enthalpy_melt_si
real*8 s, t, p, d

sea_ice_enthalpy_melt_si = errorreturn

if(equi_sea_ice_done /= isok) return

s = equi_sea_ice_s
t = equi_sea_ice_t
p = equi_sea_ice_p

d = sal_dilution_si(s, t, p)
if(d == errorreturn) return

sea_ice_enthalpy_melt_si = -t * d * sea_ice_s_si(1, 0, t, p)

end function

!==========================================================================
function sea_ice_volume_melt_si()
!==========================================================================

!THIS FUNCTION RETURNS THE MELTING VOLUME IN M3/KG OF ICE IN SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_VOLUME_MELT_SI = -9.10140854473E-05

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_VOLUME_MELT_SI = -9.67135426848E-05

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_VOLUME_MELT_SI = -9.18186917900E-05

implicit none
real*8 sea_ice_volume_melt_si
real*8 s, t, p, d

sea_ice_volume_melt_si = errorreturn

if(equi_sea_ice_done /= isok) return

s = equi_sea_ice_s
t = equi_sea_ice_t
p = equi_sea_ice_p

d = sal_dilution_si(s, t, p)
if(d == errorreturn) return

sea_ice_volume_melt_si = d * sea_ice_s_si(0, 1, t, p)

end function

!==========================================================================
function sea_ice_enthalpy_ice_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTHALPY IN J/KG OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_ENTHALPY_ICE_SI = -337351.999358

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_ENTHALPY_ICE_SI = -323205.968289

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_ENTHALPY_ICE_SI = -339929.555499

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_enthalpy_ice_si

sea_ice_enthalpy_ice_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_enthalpy_ice_si = ice_enthalpy_si(equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_enthalpy_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTHALPY IN J/KG OF SEAWATER AT THE EQUILIBRIUM WITH ICE,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_ENTHALPY_SEA_SI = -7613.1933791|9

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_ENTHALPY_SEA_SI = 2832.9491040|7

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_ENTHALPY_SEA_SI = -12742.8664892

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_enthalpy_sea_si

sea_ice_enthalpy_sea_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_enthalpy_sea_si = sea_enthalpy_si(equi_sea_ice_s, equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_entropy_ice_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTROPY IN J/(KG K) OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_ENTROPY_ICE_SI = -1235.44872812

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_ENTROPY_ICE_SI = -1247.71314646

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_ENTROPY_ICE_SI = -1244.97335506

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_entropy_ice_si 

sea_ice_entropy_ice_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_entropy_ice_si = ice_entropy_si(equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_entropy_sea_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ENTHALPY IN J/KG OF SEAWATER AT THE EQUILIBRIUM WITH ICE,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_ENTROPY_SEA_SI = -27.926459810|3

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_ENTROPY_SEA_SI = -46.7361169560

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_ENTROPY_SEA_SI = -53.1667911144

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_entropy_sea_si

sea_ice_entropy_sea_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_entropy_sea_si = sea_entropy_si(equi_sea_ice_s, equi_sea_ice_t, equi_sea_ice_p)

end function

!==========================================================================
function sea_ice_pressure_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE PRESSURE IN PA OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_PRESSURE_SI = 100000

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_PRESSURE_SI = 16132047.4385

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_PRESSURE_SI = 100000

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_pressure_si

sea_ice_pressure_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_pressure_si = equi_sea_ice_p

end function

!==========================================================================
function sea_ice_temperature_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE TEMPERATURE IN K OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_TEMPERATURE_SI = 271.240373585159

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_TEMPERATURE_SI = 270

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_TEMPERATURE_SI = 270


!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_temperature_si

sea_ice_temperature_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_temperature_si = equi_sea_ice_t

end function

!==========================================================================
function sea_ice_salinity_si()
!==========================================================================

!THIS FUNCTION RETURNS THE ABSOLUTE BRINE SALINITY IN KG/KG OF ICE AT THE EQUILIBRIUM WITH SEAWATER,
!SET BY A PREVIOUS CALL OF EITHER SET_SEA_ICE_EQUILIBRIUM_TO_S_P,
!SET_SEA_ICE_EQUILIBRIUM_TO_S_T OR SET_SEA_ICE_EQUILIBRIUM_TO_T_P

!CHECK VALUES WITH DEFAULT SETTINGS:
!SET_SEA_ICE_EQ_AT_S_P 0.035, 1E5
!SEA_ICE_SALINITY_SI = 0.035

!SET_SEA_ICE_EQ_AT_S_T 0.035, 270
!SEA_ICE_SALINITY_SI = 0.035

!SET_SEA_ICE_EQ_AT_T_P 270, 1E5
!SEA_ICE_SALINITY_SI = 0.0560264150322

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS OF THIS MODULE

implicit none
real*8 sea_ice_salinity_si

sea_ice_salinity_si = errorreturn

if(equi_sea_ice_done /= isok) return

sea_ice_salinity_si = equi_sea_ice_s

end function

!==========================================================================
function sea_ice_iteration_at_t_p(t_si, p_si, maxit, eps, s_sea_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER IN THE TWO PHASES AT GIVEN ABSOLUTE TEMPERATURE T_SI IN K
!AND ABSOLUTE PRESSURE P_SI IN PA FROM AN INITIAL GUESS FOR ABSOLUTE BRINE SALINITY S_SEA_SI 
!IN KG/KG.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN BRINE SALINITY.

!OUTPUT: 
!        SEA_ICE_ITERATION_AT_T_P = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_ITERATION_AT_T_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS FUNCTION HAS RETURNED AN ERROR
!        - SALINITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - T-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
! S_SEA_SI: BRINE SALINITY IN KG/KG IN THE SEAWATER-ICE EQUILIBRIUM

!INPUT: 
!       T_SI: ABSOLUTE TEMPERATURE IN K
!       P_SI: ABSOLUTE PRESSURE IN PA
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF SALINITY
!             EPS > 0: ABSOLUTE SALINITY TOLERANCE IN KG/KG
!             EPS < 0: RELATIVE SALINITY TOLERANCE
!   S_SEA_SI: INITIAL GUESS OF BRINE SALINITY IN KG/KG

implicit none
integer it, maxit
real*8 sea_ice_iteration_at_t_p, t_si, p_si, eps, s_sea_si
real*8 s, ds
real*8 mus, sgss
real*8 gl, gi

sea_ice_iteration_at_t_p = errorreturn

if(s_sea_si <= 0d0 .or. s_sea_si >= 1d0) then
  s_sea_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    s_sea_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax  &
  .or. p_si < sal_pmin .or. p_si > sal_pmax) then
    s_sea_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) then
    s_sea_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0 .or. &
    t_si <= 0d0) then
    s_sea_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

!GET THE CHEMICAL POTENTIAL OF LIQUID WATER GL(T,P) (THIS IMPLIES ITERATION)
gl = liq_gibbs_energy_si(t_si, p_si)
if(gl == errorreturn) go to 999

!GET THE CHEMICAL POTENTIAL OF ICE GI(T,P)
gi = ice_g_si(0, 0, t_si, p_si)
if(gi == errorreturn) go to 999

if(gl <= gi) go to 999  !NO ICE IN SEAWATER

s = s_sea_si

do it = 1, maxit

  !DERIVATIVES OF THE GIBBS FUNCTIONS FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(s, t_si, p_si) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) goto 999
  sgss = sal_dilution_si(s, t_si, p_si)  !S* D2G/DS2
  if(sgss == errorreturn) goto 999
  if(sgss <= 0d0) goto 999  !PHYSICALLY, THIS IS IMPOSSIBLE

  !BRINE SALINITY ITERATION STEP
  ds = (gl + mus - gi) / sgss

  !UPDATE BRINE SALINITY
  s = s + ds
  if(s <= 0d0) goto 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(ds) < eps) .or. &
     (eps < 0 .and. abs(ds) < -eps * s)) then
    if(s >= 1d0) goto 999
    sea_ice_iteration_at_t_p = isok
    s_sea_si = s
    goto 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !SAL_LIMITS
  if(s_sea_si < sal_smin .or. s_sea_si > sal_smax) then
    s_sea_si = errorreturn
    sea_ice_iteration_at_t_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_ice_iteration_at_s_t(sa_si, t_si, maxit, eps, d_si, p_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE PHASE EQUILIBRIUM FROM
!EQUAL CHEMICAL POTENTIALS OF WATER IN THE TWO PHASES AT GIVEN ABSOLUTE SALINITY SA_SI IN KG/KG 
!AND ABSOLUTE TEMPERATURE T_SI IN K, FROM AN INITIAL GUESS FOR THE LIQUID DENSITY D_SI IN KG/KG.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN PRESSURE.

!OUTPUT: 
!        SEA_ICE_ITERATION_AT_S_T = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_ITERATION_AT_S_T = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - S-T VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
! D_SI: DENSITY OF WATER AT (T_SI, P_SI) IN KG/M3
! P_SI: ABSOLUTE PRESSURE IN PA AT THE SEAWATER-ICE EQUILIBRIUM

!INPUT: 
!      SA_SI: ABSOLUTE SALINITY IN KG/KG
!       T_SI: ABSOLUTE TEMPERATURE IN K
!      MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!        EPS: REQUIRED ACCURACY OF PRESSURE
!             EPS > 0: ABSOLUTE PRESSURE TOLERANCE IN PA
!             EPS < 0: RELATIVE PRESSURE TOLERANCE
!       D_SI: INITIAL GUESS FOR DENSITY OF WATER IN KG/M3 AT (T_SI, P_SI)

implicit none
integer it, maxit
real*8 sea_ice_iteration_at_s_t, sa_si, t_si, eps, d_si, p_si
real*8 p, d, dd
real*8 mus, p_old
real*8 gl, gi, gi_p
real*8 gs_p, gs_sp
real*8 f, f_d, f_dd

sea_ice_iteration_at_s_t = errorreturn

if(d_si <= 0d0) then
  p_si = errorreturn
  d_si = errorreturn
  return
endif

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) then
    p_si = errorreturn
    d_si = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si < sal_smin .or. sa_si > sal_smax) then
    p_si = errorreturn
    d_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    p_si = errorreturn
    d_si = errorreturn
    return
  endif
else
  if(t_si <= 0d0 .or. &
    sa_si <= 0d0 .or. sa_si >= 1d0) then
    p_si = errorreturn
    d_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

d = d_si

do it = 0, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF FLUID WATER FOR NEWTON ITERATION
  f = flu_f_si(0, 0, t_si, d)
  if(f == errorreturn) go to 999
  f_d = flu_f_si(0, 1, t_si, d)
  if(f_d == errorreturn) go to 999
  f_dd = flu_f_si(0, 2, t_si, d)
  if(f_dd == errorreturn) go to 999
  p_old = p
  p = d ** 2 * f_d
  gl = f + d * f_d

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t_si, p)
  if(gi == errorreturn) go to 999
  gi_p = ice_g_si(0, 1, t_si, p)
  if(gi_p == errorreturn) go to 999

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(sa_si, t_si, p) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  gs_p = sal_g_si(0, 0, 1, sa_si, t_si, p)
  if(gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, sa_si, t_si, p)
  if(gs_sp == errorreturn) go to 999

  if(it > 0) then
    !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
    if((eps > 0d0 .and. abs(p - p_old) < eps) .or. &
      (eps < 0 .and. abs(p - p_old) < -eps * abs(p))) then
      if(p <= 0d0) go to 999
      sea_ice_iteration_at_s_t = isok
      p_si = p
      d_si = d
      go to 999
    end if
  end if

  if(it == maxit) go to 999

  !DENSITY ITERATION STEP
  dd = (2d0 * f_d + d * f_dd) * (1d0 + d * (-gi_p + gs_p - sa_si * gs_sp))
  if(dd == 0d0) go to 999
  dd = (gi - gl - mus) / dd

  !UPDATE LIQUID DENSITY
  d = d + dd
  if(d <= 0d0) go to 999

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(d_si <= flu_dmin .or. d_si > flu_dmax) then
    p_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_t = errorreturn
    return
  endif
  !SAL_LIMITS
  if(p_si < sal_pmin .or. p_si > sal_pmax) then
    p_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_t = errorreturn
    return
  endif
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    p_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_t = errorreturn
    return
  endif
endif

end function

!==========================================================================
function sea_ice_iteration_at_s_p(sa_si, p_si, maxit, eps, d_si, t_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SEAWATER-ICE PHASE EQUILIBRIUM FROM EQUAL CHEMICAL POTENTIALS 
!OF WATER AT GIVEN ABSOLUTE SALINITY SA_SI IN KG/KG AND ABSOLUTE PRESSURE P_SI IN PA,
!FROM AN INITIAL GUESSES FOR THE ABSOLUTE TEMPERATURE T_SI IN K AND THE LIQUID DENSITY D_SI
!IN KG/M3.
!THE ITERATION LIMIT EPS REFERS TO THE ERROR IN TEMPERATURE.

!OUTPUT: 
!        SEA_ICE_ITERATION_AT_S_P = ISOK IF SUCCESSFULLY DONE
!        SEA_ICE_ITERATION_AT_S_P = ERRORRETURN IS RETURNED IF
!        - THE MAXIMUM NUMBER OF ITERATIONS IS EXCEEDED WITHOUT MEETING THE EXIT CRITERION
!        - THE FUNCTION CALL TO A GIBBS/HELMHOLTZ FUNCTION HAS RETURNED AN ERROR
!        - DENSITY OR TEMPERATURE HAS TAKEN A ZERO OR NEGATIVE VALUE DURING THE ITERATION
!        - S-P VALUES DO NOT PERMIT THE STABLE EXISTENCE OF ICE
!  D_SI: DENSITY OF WATER AT (T_SI, P_SI)
!  T_SI: ABSOLUTE TEMPERATURE IN K OF THE SEAWATER-ICE EQUILIBRIUM

!INPUT: 
!       SA_SI: ABSOLUTE SALINITY IN KG/KG
!        P_SI: ABSOLUTE PRESSURE IN PA
!       MAXIT: MAXIMUM NUMBER OF ITERATION STEPS TO BE DONE
!         EPS: REQUIRED ACCURACY OF PRESSURE
!              EPS > 0: ABSOLUTE TEMPERATURE TOLERANCE IN K
!              EPS < 0: RELATIVE TEMPERATURE TOLERANCE
!        D_SI: INITIAL GUESS FOR DENSITY OF WATER IN KG/M3 AT (T_SI, P_SI)
!        T_SI: INITIAL GUESS FOR ABSOLUTE FREEZING TEMPERATURE IN K

implicit none
integer it, maxit
real*8 sea_ice_iteration_at_s_p, sa_si, p_si, eps, d_si, t_si
real*8 d, dd, t, dt
real*8 mus, p
real*8 gl, gi, gi_t
real*8 gs_t, gs_st
real*8 f, f_t, f_td
real*8 f_d, f_dd

real*8 a(2, 2), b(2), x(2)

sea_ice_iteration_at_s_p = errorreturn

if((d_si <= 0d0) .or. &
  (t_si <= 0d0)) then
  d_si = errorreturn
  t_si = errorreturn
  return
endif

if(check_limits == 1) then
  !SAL_LIMITS
  if(sa_si < sal_smin .or. sa_si > sal_smax &
  .or. p_si < sal_pmin .or. p_si > sal_pmax) then
    d_si = errorreturn
    t_si = errorreturn
    return
  endif
  !ICE_LIMITS
  if(p_si <= ice_pmin .or. p_si > ice_pmax) then
    d_si = errorreturn
    t_si = errorreturn
    return
  endif
else
  if(p_si <= 0d0 .or. &
    sa_si <= 0d0 .or. sa_si >= 1d0) then
    d_si = errorreturn
    t_si = errorreturn
    return
  endif
endif

check_limits = check_limits - 1

t = t_si
d = d_si

do it = 1, maxit

  !DERIVATIVES OF THE HELMHOLTZ FUNCTION OF FLUID WATER FOR NEWTON ITERATION
  f = flu_f_si(0, 0, t, d)
  if(f == errorreturn) go to 999
  f_t = flu_f_si(1, 0, t, d)
  if(f_t == errorreturn) go to 999
  f_d = flu_f_si(0, 1, t, d)
  if(f_d == errorreturn) go to 999
  f_td = flu_f_si(1, 1, t, d)
  if(f_td == errorreturn) go to 999
  f_dd = flu_f_si(0, 2, t, d)
  if(f_dd == errorreturn) go to 999
  p = d ** 2 * f_d
  gl = f + d * f_d

  !DERIVATIVES OF THE GIBBS FUNCTION OF ICE FOR NEWTON ITERATION
  gi = ice_g_si(0, 0, t, p_si)  !dbg DW chnaged p to p_si
  if(gi == errorreturn) go to 999
  gi_t = ice_g_si(1, 0, t, p_si)
  if(gi_t == errorreturn) go to 999

  !DERIVATIVES OF THE SALINE GIBBS FUNCTION FOR NEWTON ITERATION
  mus = sal_chempot_h2o_si(sa_si, t, p_si) !SALINE CHEMICAL POTENTIAL GS + S*DG/DS
  if(mus == errorreturn) go to 999
  gs_t = sal_g_si(0, 1, 0, sa_si, t, p_si)
  if(gs_t == errorreturn) go to 999
  gs_st = sal_g_si(1, 1, 0, sa_si, t, p_si)
  if(gs_st == errorreturn) go to 999

  !COEFFICIENT MATRIX
  a(1, 1) = f_t + d * f_td + gs_t - sa_si * gs_st - gi_t
  a(1, 2) = 2d0 * f_d + d * f_dd
  a(2, 1) = d ** 2 * f_td
  a(2, 2) = d * (2d0 * f_d + d * f_dd)
  !RIGHT-HAND SIDES, MUST VANISH IN EQUILIBRIUM
  b(1) = gi - gl - mus
  b(2) = p_si - p

  !SOLVE EQUATIONS
  if(matrix_solve(a, b, x, 2) /= 0) go to 999 !MATRIX SINGULAR
  dt = x(1)
  dd = x(2)

  !UPDATE DENSITY AND TEMPERATURE
  d = d + dd
  if(d <= 0d0) go to 999
  t = t + dt
  if(t <= 0d0) go to 999

  !CHECK ABSOLUTE OR RELATIVE ERROR LIMIT
  if((eps > 0d0 .and. abs(dt) < eps) .or. &
    (eps < 0 .and. abs(dt) < -eps * abs(t))) then
    if(t <= 0d0) go to 999
    sea_ice_iteration_at_s_p = isok
    t_si = t
    d_si = d
    go to 999
  end if

enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_si <= flu_dmin .or. d_si > flu_dmax) then
    t_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_p = errorreturn
    return
  endif
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax) then
    t_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_p = errorreturn
    return
  endif
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax) then
    t_si = errorreturn
    d_si = errorreturn
    sea_ice_iteration_at_s_p = errorreturn
    return
  endif
endif

end function

!==========================================================================
function aux_brinesalinity_si(t_si, p_si)
!==========================================================================

!COMPUTES A BRINE SALINITY ESTIMATE IN KG/KG FOR SEAWATER-ICE EQUILIBRIUM AT GIVEN
!ABSOLUTE TEMPERATURE IN K AND ABSOLUTE PRESSURE IN PA FROM CLAUSIUS-CLAPEYRON
!AND RAOULT LAWS

implicit none
real*8 aux_brinesalinity_si, t_si, p_si
real*8 tt, pt, alpha, chi

tt = tp_temperature_si
pt = tp_pressure_iapws95_si !IAPWS-95 TRIPLE POINT PRESSURE IN Pa

alpha = -0.217d0            !RAOULT COFFICIENT, (T-TT)/T = ALPHA * S,  IN KG/KG
chi = -0.0000000743d0       !CLAUSIUS-CLAPEYRON COEFFICIENT,  (T-TT) = CHI * (P-PT), IN K/PA !dbg [CHI] in K/PA?

aux_brinesalinity_si = errorreturn
if(p_si <= 0d0) return
if(t_si <= 0d0) return

aux_brinesalinity_si = (t_si - tt - chi * (p_si - pt)) / (alpha * t_si)

end function

!==========================================================================
function aux_freezingtemperature_si(sa_si, p_si)
!==========================================================================

!COMPUTES A FREEZING TEMPERATURE ESTIMATE IN K FOR SEAWATER-ICE EQUILIBRIUM AT GIVEN
!BRINE SALINITY IN KG/KG AND ABSOLUTE PRESSURE IN PA FROM CLAUSIUS-CLAPEYRON
!AND RAOULT LAWS

implicit none
real*8 aux_freezingtemperature_si, sa_si, p_si
real*8 tt, pt, alpha, chi

tt = tp_temperature_si
pt = tp_pressure_iapws95_si !IAPWS-95 TRIPLE POINT PRESSURE IN Pa

alpha = -0.217d0            !RAOULT COFFICIENT, (T-TT)/T = ALPHA * S,  IN KG/KG
chi = -0.0000000743d0       !CLAUSIUS-CLAPEYRON COEFFICIENT,  (T-TT) = CHI * (P-PT), IN K/PA !dbg [CHI] in K/PA

aux_freezingtemperature_si = errorreturn
if(p_si <= 0d0) return
if(sa_si < 0d0) return
if(sa_si >= 1d0) return

aux_freezingtemperature_si = (tt + chi * (p_si - pt)) / (1 - alpha * sa_si)

end function

!==========================================================================
function aux_liq_density_correlation_si(t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE DENSITY IN KG/M3 OF THE FREEZING LIQUID
!AS A FUNCTION OF ABSOLUTE TEMPERATURE T_SI IN K FROM AN APPROXIMATE
!CORRELATION FIT BETWEEN 252 AND 273 K WITH RMS = 1.2 E-4 IN D/DT

implicit none
real*8 aux_liq_density_correlation_si, t_si
real*8 tt, dt, a1, a2, a3
real*8 tr, dr

tt = tp_temperature_si
dt = tp_density_liq_iapws95_si !TRIPLE POINT DENSITY OF LIQUID WATER

a1 = -1.78582981492113d0
a2 = -12.2325084306734d0
a3 = -52.8236936433529d0

aux_liq_density_correlation_si = errorreturn

if(t_si <= 0d0) return

tr = t_si / tt - 1d0
dr = a1 * tr + a2 * tr ** 2 + a3 * tr ** 3

aux_liq_density_correlation_si = dt * (dr + 1d0)

end function

!==========================================================================
function aux_meltingpressure_si(sa_si, t_si)
!==========================================================================

!COMPUTES A MELTING PRESSURE ESTIMATE IN PA FOR SEAWATER-ICE EQUILIBRIUM AT GIVEN
!ABSOLUTE TEMPERATURE IN K AND BRINE SALINITY IN KG/KG FROM CLAUSIUS-CLAPEYRON
!AND RAOULT LAWS

implicit none
real*8 aux_meltingpressure_si, sa_si, t_si
real*8 tt, pt, alpha, chi

tt = tp_temperature_si
pt = tp_pressure_iapws95_si !IAPWS-95 TRIPLE POINT PRESSURE IN Pa

alpha = -0.217d0            !RAOULT COFFICIENT, (T-TT)/T = ALPHA * S,  S IN KG/KG
chi = -0.0000000743d0       !CLAUSIUS-CLAPEYRON COEFFICIENT,  (T-TT) = CHI * (P-PT)

aux_meltingpressure_si = errorreturn
if(sa_si < 0d0) return
if(sa_si >= 1d0) return
if(t_si <= 0d0) return

aux_meltingpressure_si = pt + (t_si - tt - alpha * t_si * sa_si) / chi

end function

!==========================================================================
function aux_density_eos80_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER COMPUTED FROM THE
!INTERNATIONAL EQUATION OF STATE OF SEAWATER 1980, EOS-80, AS A FUNCTION OF TEMPERATURE
!AND PRESSURE

!OUTPUT:  
!aux_density_eos80_si: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!p_si: ABSOLUTE PRESSURE IN PA

implicit none
real*8 aux_density_eos80_si, t_si, p_si
real*8 aw, bw, kw, rw
real*8 t, p

aux_density_eos80_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

!R.C.MILLARD JR.: INTERNATIONAL OCEANOGRAPHIC TABLES VOL.4
!UNESCO TECHNICAL PAPERS IN MARINE SCIENCE 40
!UNESCO 1987

t = 1.00024d0*(t_si - 273.15d0) !CONVERT TO DEGC(T90) THEN USE SAUNDER'S FORMULA
p = (p_si - 101325d0)/1d5

bw = 0.0000850935d0 + t * (-0.00000612293d0 + t * 0.000000052787d0)
aw = 3.239908d0 + t * (0.00143713d0 + t * (0.000116092d0 - t * 0.000000577905d0))
kw = 0.01360477d0 - t * 0.00005155288d0
kw = 19652.21d0 + t * (148.4206d0 + t * (-2.327105d0 + t * kw))
kw = kw + p * (aw + p * bw)
if(kw == 0) return

kw = 1d0 - p / kw
if(kw == 0d0) return

rw = 0.0001001685d0 + t * (-0.000001120083d0 + t * 0.000000006536332d0)
rw = 999.842594d0 + t * (0.06793952d0 + t * (-0.00909529d0 + t * rw))

aux_density_eos80_si = rw / kw

end function

!==========================================================================
subroutine set_it_ctrl_sea_ice(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE ITERATION USED TO COMPUTE
!SEAWATER-ICE EQUILIBRIUM

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VALUES TO COMPUTE EQUILIBRIUM

!'INIT_LIQ_DENS'        0           USE DEFAULT LIQUID DENSITY TO START =
!                                   LIQUID DENSITY CORRELATION FUNCTION, D = AUX_LIQ_DENSITY_CORRELATION(T)
!'INIT_LIQ_DENS'       -1           USE TRIPLE POINT LIQUID DENSITY TO START
!'INIT_LIQ_DENS'        D > 0       USE VALUE D AS LIQUID DENSITY TO START

!'INIT_BRINE_SA'        0           USE DEFAULT BRINE SALINITY TO START ( = AUX_BRINESALINITY_SI(T, P))
!'INIT_BRINE_SA'       -1           USE NORMAL SALINITY TO START
!'INIT_BRINE_SA'        S > 0       USE VALUE S AS TEMPERATURE TO START

!'INIT_TEMP    '        0           USE DEFAULT TEMPERATURE TO START ( = AUX_FREEZINGTEMPERATURE_SI(S, P))
!'INIT_TEMP    '       -1           USE 0 C TEMPERATURE TO START
!'INIT_TEMP    '        T > 0       USE VALUE T AS TEMPERATURE TO START

!'TOL_BRINE_SA '        0           USE DEFAULT EXIT ACCURACY FOR BRINE SALINITY (0.1 PPM)
!'TOL_BRINE_SA '        EPS         USE EPS AS EXIT ACCURACY FOR BRINE SALINITY (EPS < 0 MEANS RELATIVE ERROR)

!'TOL_TEMP     '        0           USE DEFAULT EXIT ACCURACY FOR FREEZING TEMPERATURE (0.1 PPM)
!'TOL_TEMP     '        EPS         USE EPS AS EXIT ACCURACY FOR TEMPERATURE (EPS < 0 MEANS RELATIVE ERROR)

!'TOL_PRESS    '        0           USE DEFAULT EXIT ACCURACY FOR MELTING PRESSURE (0.1 PPM)
!'TOL_PRESS    '        EPS         USE EPS AS EXIT ACCURACY FOR PRESSURE (EPS < 0 MEANS RELATIVE ERROR)

implicit none
character*13 key
real*8 value

call init_it_ctrl_sea_ice

call clear_sea_ice_state

select case (trim(key))

  case ('it_steps')              !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100    !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ('init_liq_dens')         !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_liquid = 0       !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_liquid = nint(value)
    else
      ctrl_mode_liquid = 1
      ctrl_init_d_liq = value
    endif

  case ('init_brine_sa')         !START LIQUID DENSITY
    if(value == 0d0)then
      ctrl_mode_sea = 0          !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_sea = nint(value)
    else
      ctrl_mode_sea = 1
      ctrl_init_s_sea = value
    endif

  case ('init_temp')             !START TEMPERATURE
    if(value == 0d0)then
      ctrl_mode_temperature = 0  !DEFAULT = AUX POLYNOMIAL
    elseif(value < (-1d0))then
      !ignore it
    elseif(value < 0d0)then
      ctrl_mode_temperature = nint(value)
    else
      ctrl_mode_temperature = 1
      ctrl_init_t = value
    endif

  case ('tol_brine_sa')          !REQUIRED TOLERANCE FOR BRINE SALINITY
    if(value == 0d0)then
      ctrl_eps_exit_s = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_s = value
    endif

  case ('tol_temp')              !REQUIRED TOLERANCE FOR TEMPERATURE
    if(value == 0d0)then
      ctrl_eps_exit_t = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_t = value
    endif

  case ('tol_press')             !REQUIRED TOLERANCE FOR PRESSURE
    if(value == 0d0)then
      ctrl_eps_exit_p = -0.0000001d0   !DEFAULT = 0.1 PPM RELATIVE
    else
      ctrl_eps_exit_p = value
    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

!==========================================================================
subroutine set_sea_ice_state(sa_si, t_si, p_si, d_si)
!==========================================================================

!STORES THE ACTUAL PROPERTIES AS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none
real*8 sa_si, t_si, p_si, d_si

equi_sea_ice_done = isok

equi_sea_ice_s = sa_si
equi_sea_ice_t = t_si
equi_sea_ice_p = p_si
equi_sea_ice_d_liq = d_si   !DENSITY OF PURE WATER

end subroutine

!==========================================================================
subroutine clear_sea_ice_state
!==========================================================================

!CLEARS THE CURRENT EQUILIBRIUM STATE DESCRIPTOR

implicit none

equi_sea_ice_done = 0

end subroutine

!==========================================================================
subroutine init_it_ctrl_sea_ice
!==========================================================================

implicit none
real*8 s_norm, tt, dt

s_norm = so_salinity_si  !KG/KG OF KCL NORMALISED SEAWATER
tt = tp_temperature_si
dt = tp_density_liq_iapws95_si  !TRIPLE POINT DENSITY OF LIQUID WATER

if(ctrl_initialized == isok) return

ctrl_initialized = isok

!SET DEFAULT VALUES AND MODES FOR THE ITERATION
ctrl_mode_liquid = 0
ctrl_mode_sea = 0
ctrl_mode_temperature = 0
ctrl_loop_maximum = 100
ctrl_init_s_sea = s_norm
ctrl_init_t = tt
ctrl_init_d_liq = dt

ctrl_eps_exit_s = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_t = -0.0000001d0 !RELATIVE, 0.1 PPM
ctrl_eps_exit_p = -0.0000001d0 !RELATIVE, 0.1 PPM

end subroutine

end module sea_ice_4
