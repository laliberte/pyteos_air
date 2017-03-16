module sea_3b

!#########################################################################

!THIS MODULE IMPLEMENTS THE ENTHALPY OF SEAWATER DEPENDING ON SALINITY, ENTROPY
!AND PRESSURE, AS WELL AS ITS PARTIAL DERIVATIVES, AS DEFINED IN:

!FEISTEL, R.
!A GIBBS FUNCTION FOR SEAWATER THERMODYNAMICS
!FOR -6 DEG C TO 80 DEG C AND SALINITY UP TO 120 G KG-1
!DEEP-SEA RESEARCH I, SUBMITTED 2007

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
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_3A,      FILE FLU_3A.F90
!     SEA_3A,      FILE SEA_3A.F90

!NESTED USE ASSOCIATIONS
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90
!     FLU_1,       FILE FLU_1.F90
!     SAL_1,       FILE SAL_1.F90

!#########################################################################

use constants_0
use sal_2
use flu_3a
use sea_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_h_si, sea_h_contraction_h_si, sea_h_contraction_t_si, &
          sea_h_contraction_theta_si, sea_h_expansion_h_si, &
          sea_h_expansion_t_si, sea_h_expansion_theta_si, &
          sea_potdensity_si, sea_potenthalpy_si,  &
          sea_pottemp_si, sea_temperature_si, set_it_ctrl_pottemp

integer, private :: ctrl_initialized, ctrl_loop_maximum, ctrl_mode_pottemp
real*8, private :: ctrl_init_pottemp, ctrl_eps_exit_pottemp

contains

!==========================================================================
function sea_h_si(drv_s, drv_eta, drv_p, sa_si, eta_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS ENTHALPY AS A THERMODYNAMIC POTENTIAL, DEPENDING
!ON SALINITY, ENTROPY AND PRESSURE

!OUTPUT:
!SEA_H_SI AS THE S-ETA-P DERIVATIVE

!(D/DS)^DRV_S (D/DETA)^DRV_ETA (D/DP)^DRV_P H(S,ETA,P)

!OF THE SPECIFIC ENTHALPY OF SEAWATER, H(S,ETA,P), IN J/KG

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    SPECIFIC ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE PRESSURE IN PA

!CHECK VALUES WITH DEFAULT SETTINGS:
!SEA_H_SI( 0, 0, 0, 0.035, 500, 1E5) = 145481.970750
!SEA_H_SI( 1, 0, 0, 0.035, 500, 1E5) = 86860.7980486
!SEA_H_SI( 0, 1, 0, 0.035, 500, 1E5) = 309.557955853
!SEA_H_SI( 0, 0, 1, 0.035, 500, 1E5) = 9.81092930969E-04
!SEA_H_SI( 2, 0, 0, 0.035, 500, 1E5) = 2393730.16716
!SEA_H_SI( 1, 1, 0, 0.035, 500, 1E5) = 72.5298236488
!SEA_H_SI( 1, 0, 1, 0.035, 500, 1E5) =-6.84629317367E-04
!SEA_H_SI( 0, 2, 0, 0.035, 500, 1E5) = 7.72873234085E-02
!SEA_H_SI( 0, 1, 1, 0.035, 500, 1E5) = 2.86305358702E-08
!SEA_H_SI( 0, 0, 2, 0.035, 500, 1E5) =-3.96880481108E-13

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

implicit none
integer drv_s, drv_eta, drv_p
real*8 sea_h_si, sa_si, eta_si, p_si
real*8 t

sea_h_si = errorreturn

if(drv_s < 0 .or. drv_s > 2) return
if(drv_eta < 0 .or. drv_eta > 2) return
if(drv_p < 0 .or. drv_p > 2) return

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(p_si < 0d0) return

!COMPUTE TEMPERATURE FROM ENTROPY.
!IF P_SI = IN-SITU PRESSURE, THIS IS IN-SITU TEMPERATURE
!IF P_SI = REFERENCE PRESSURE, THIS IS POTENTIAL TEMPERATURE
t = sea_temperature_si(sa_si, eta_si, p_si)
if(t == errorreturn) return

sea_h_si = sea_s_eta_p_derivatives_si(drv_s, drv_eta, drv_p, sa_si, t, p_si)

end function

!==========================================================================
function sea_pottemp_si(sa_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL TEMPERATURE OF SEAWATER

!OUTPUT:
!THETA(S,T,P,PR) ABSOLUTE POTENTIAL TEMPERATURE OF SEAWATER IN K

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_POTTEMP_SI(0.035,300,1e7,1e5) = 299.771869405

implicit none
real*8 sea_pottemp_si, sa_si, t_si, p_si, pr_si
real*8 s

sea_pottemp_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

if(p_si == pr_si) then
  sea_pottemp_si = t_si
  return
end if

s = sea_entropy_si(sa_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

sea_pottemp_si = sea_temperature_si(sa_si, s, pr_si)

end function

!==========================================================================
function sea_potdensity_si(sa_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL DENSITY OF SEAWATER

!OUTPUT:
!RHO_THETA(S,T,P,PR) POTENTIAL DENSITY OF SEAWATER IN KG/M^3,

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_POTDENSITY_SI(0.035,300,1e7,1e5) = 1022.71520130

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

implicit none
real*8 sea_potdensity_si, sa_si, t_si, p_si, pr_si
real*8 s, v

sea_potdensity_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = sea_entropy_si(sa_si, t_si, p_si)  !specific entropy in-situ
if(s == errorreturn) return

v = sea_h_si(0, 0, 1, sa_si, s, pr_si) !specific volume at reference pressure
if(v == errorreturn) return
if(v <= 0d0) return

sea_potdensity_si = 1d0 / v

end function

!==========================================================================
function sea_potenthalpy_si(sa_si, t_si, p_si, pr_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL ENTHALPY OF SEAWATER

!OUTPUT:
!H_THETA(S,T,P,PR) POTENTIAL ENTHALPY OF SEAWATER IN J/KG

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!T_SI      ABSOLUTE IN-SITU TEMPERATURE IN K
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PR_SI     ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_POTENTHALPY_SI(0.035,300,1e7,1e5) = 106307.996083

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

implicit none
real*8 sea_potenthalpy_si, sa_si, t_si, p_si, pr_si
real*8 s

sea_potenthalpy_si = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si < 0d0) return
if(p_si < 0d0) return
if(pr_si < 0d0) return

s = sea_entropy_si(sa_si, t_si, p_si)  !SPECIFIC ENTROPY IN-SITU
if(s == errorreturn) return

sea_potenthalpy_si = sea_h_si(0, 0, 0, sa_si, s, pr_si)

end function

!==========================================================================
Function sea_h_expansion_t_si(sa_si, eta_si, p_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THERMAL EXPANSION COEFFICIENT OF SEAWATER AS FUNCTION OF ENTROPY

!RETURNS   SEA_H_EXPANSION_T_SI(S,ETA,P) = (1/V)*(DV/DT)_S_P THERMAL EXPANSION COEFFICIENT IN 1/K,
!                                          AT CONSTANT SALINITY AND PRESSURE
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_EXPANSION_T_SI(0.035, 500, 1e5) = 3.77581809091E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_expansion_t_si, sa_si, eta_si, p_si
Real*8 hp, hep, hee

sea_h_expansion_t_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
hee = sea_h_si(0, 2, 0, sa_si, eta_si, p_si)
If (hee == ErrorReturn .or. hee == 0d0) Return

sea_h_expansion_t_si = hep / (hp * hee)

End Function

!==========================================================================
Function sea_h_expansion_theta_si(sa_si, eta_si, p_si, pref_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THERMAL EXPANSION COEFFICIENT WRT POTENTIAL TEMPERATURE 
!OF SEAWATER AS A FUNCTION OF ENTROPY

!RETURNS   SEA_H_EXPANSION_THETA_SI(S,ETA,P) = (1/V)*(DV/DTHETA)_S_P THERMAL EXPANSION COEFFICIENT 
!                                              IN 1/K, AT CONSTANT SALINITY AND PRESSURE
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_EXPANSION_THETA_SI(0.035, 500, 1e7, 1e5) = 3.84755380181E-04

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_expansion_theta_si, sa_si, eta_si, p_si, pref_si
Real*8 hp, hep, hee

sea_h_expansion_theta_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return
If (pref_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
hee = sea_h_si(0, 2, 0, sa_si, eta_si, pref_si)
If (hee == ErrorReturn .or. hee == 0d0) Return

sea_h_expansion_theta_si = hep / (hp * hee)

End Function

!==========================================================================
Function sea_h_expansion_h_si(sa_si, eta_si, p_si, pref_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT WRT POTENTIAL ENTHALPY
!OF SEAWATER AS FUNCTION OF ENTROPY

!RETURNS   SEA_H_EXPANSION_H_SI(S,ETA,P) = (1/V)*(DV/DH)_S_P THERMAL EXPANSION COEFFICIENT IN (KG K)/J,
!                                          AT CONSTANT SALINITY AND PRESSURE. HPOT IS POTENTIAL ENTHALPY
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_EXPANSION_H_SI(0.035, 500, 1e7, 1e5) = 9.60618615640E-08

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_expansion_h_si, sa_si, eta_si, p_si, pref_si
Real*8 hp, hep, he

sea_h_expansion_h_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return
If (pref_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
he = sea_h_si(0, 1, 0, sa_si, eta_si, pref_si)
If (he == ErrorReturn .or. he == 0d0) Return

sea_h_expansion_h_si = hep / (hp * he)

End Function

!==========================================================================
Function sea_h_contraction_t_si(sa_si, eta_si, p_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE ADIABATIC HALINE CONTRACTION COEFFICIENT OF SEAWATER

!RETURNS   SEA_H_CONTRACTION_T_SI(S,ETA,P) = - (1/V)*(DV/DS)_T_P HALINE CONTRACTION COEFFICIENT IN KG/KG,
!                                            (KG SEAWATER / KG SALT) AT CONSTANT TEMPERATURE AND PRESSURE
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_CONTRACTION_T_SI(0.035, 500, 1e5) = 0.717342103505

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_contraction_t_si, sa_si, eta_si, p_si
Real*8 hp, hep, hse, hee, hsp

sea_h_contraction_t_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hse = sea_h_si(1, 1, 0, sa_si, eta_si, p_si)
If (hse == ErrorReturn) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
hsp = sea_h_si(1, 0, 1, sa_si, eta_si, p_si)
If (hsp == ErrorReturn) Return
hee = sea_h_si(0, 2, 0, sa_si, eta_si, p_si)
If (hee == ErrorReturn .or. hee == 0d0) Return

sea_h_contraction_t_si = (hse * hep - hsp * hee) / (hp * hee)

End Function

!==========================================================================
Function sea_h_contraction_theta_si(sa_si, eta_si, p_si, pref_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE ADIABATIC HALINE CONTRACTION COEFFICIENT OF SEAWATER
!AT CONSTANT POTENTIAL TEMPERATURE

!RETURNS   SEA_H_CONTRACTION_THETA_SI(S,ETA,P) = - (1/V)*(DV/DS)_THETA_P HALINE CONTRACTION COEFFICIENT 
!                                                IN KG/KG, (KG SEAWATER / KG SALT) AT CONSTANT 
!                                                POTENTIAL TEMPERATURE AND PRESSURE

!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_CONTRACTION_THETA_SI(0.035, 500, 1e7, 1e5) = 0.697779873590

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_contraction_theta_si, sa_si, eta_si, p_si, pref_si
Real*8 hp, hep, hse, hee, hsp

sea_h_contraction_theta_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return
If (pref_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hse = sea_h_si(1, 1, 0, sa_si, eta_si, pref_si)
If (hse == ErrorReturn) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
hsp = sea_h_si(1, 0, 1, sa_si, eta_si, p_si)
If (hsp == ErrorReturn) Return
hee = sea_h_si(0, 2, 0, sa_si, eta_si, pref_si)
If (hee == ErrorReturn .or. hee == 0d0) Return

sea_h_contraction_theta_si = (hse * hep - hsp * hee) / (hp * hee)

End Function

!==========================================================================
Function sea_h_contraction_h_si(sa_si, eta_si, p_si, pref_si)
!==========================================================================
                               
!THIS FUNCTION COMPUTES THE HALINE CONTRACTION COEFFICIENT OF SEAWATER AT CONSTANT POTENTIAL ENTHALPY

!RETURNS   SEA_H_CONTRACTION_H_SI(S,ETA,P) = - (1/V)*(DV/DS)_H_P HALINE CONTRACTION COEFFICIENT IN KG/KG,
!                                              (KG SEAWATER / KG SALT) AT CONSTANT POTENTIAL ENTHALPY 
!                                              AND PRESSURE
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE IN-SITU PRESSURE IN PA
!PREF_SI   ABSOLUTE REFERENCE PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_H_CONTRACTION_H_SI(0.035, 500, 1e7, 1e5) = 0.697779873590

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A
!      AND ON THE ITERATION SETTINGS FOR TEMPERATURE BY SET_IT_CTRL_POTTEMP OF THIS MODULE

Real*8 sea_h_contraction_h_si, sa_si, eta_si, p_si, pref_si
Real*8 hp, hep, hs, he, hsp

sea_h_contraction_h_si = ErrorReturn

If (sa_si < 0d0 .or. sa_si >= 1d0) Return
If (p_si <= 0d0) Return
If (pref_si <= 0d0) Return

hp = sea_h_si(0, 0, 1, sa_si, eta_si, p_si)
If (hp == ErrorReturn .or. hp == 0d0) Return
hs = sea_h_si(1, 0, 0, sa_si, eta_si, pref_si)
If (hs == ErrorReturn) Return
hep = sea_h_si(0, 1, 1, sa_si, eta_si, p_si)
If (hep == ErrorReturn) Return
hsp = sea_h_si(1, 0, 1, sa_si, eta_si, p_si)
If (hsp == ErrorReturn) Return
he = sea_h_si(0, 1, 0, sa_si, eta_si, pref_si)
If (he == ErrorReturn .or. he == 0d0) Return

sea_h_contraction_h_si = (hs * hep - hsp * he) / (hp * he)

End Function

!=========================================================================
function sea_temperature_si(sa_si, eta_si, p_si)
!==========================================================================
                                 
!OUTPUT:
!T(S,ETA,P) (POTENTIAL) TEMPERATURE OF SEAWATER IN K

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA

!CHECK VALUE WITH DEFAULT SETTINGS: 
!SEA_TEMPERATURE_SI(0.035, 500, 1e5) = 309.557955853

!THIS IS THE INVERSE FUNCTION TO SEA_ENTROPY_SI(SA_SI, T_SI, P_SI) IN SEA_3A

!NOTE: THE ACCURACY OF THIS FUNCTION DEPENDS ON THE ITERATION SETTINGS FOR
!      DENSITY COMPUTED IN FLU_3A AND ON THOSE MADE IN THIS MODULE

implicit none
integer maxit
real*8 sea_temperature_si, sa_si, eta_si, p_si
real*8 t, eps

sea_temperature_si = errorreturn

call init_it_ctrl_pottemp

select case (ctrl_mode_pottemp)
  case (0);    t = 273.15d0 + eta_si / 4000d0
  case default; t = ctrl_init_pottemp
end select

if(ctrl_loop_maximum == 0)then
   maxit = 100
elseif(ctrl_loop_maximum == (-1))then
   sea_temperature_si = t
   return
elseif(ctrl_loop_maximum > 0)then
   maxit = ctrl_loop_maximum
else
   return
endif

eps = ctrl_eps_exit_pottemp
if(eps <= 0d0) return

!RUN THE ITERATION
sea_temperature_si = pottemp_iteration(sa_si, eta_si, p_si, t, maxit, eps)

end function

!==========================================================================
function sea_s_eta_p_derivatives_si(drv_s, drv_eta, drv_p, sa_si, t_si, p_si)
!==========================================================================
                                              
!THIS FUNCTION COMPUTES SEAWATER S-ETA-P DERIVATIVES OF H FROM S-T-P DERIVATIVES OF G

implicit none
integer drv_s, drv_eta, drv_p
real*8 sea_s_eta_p_derivatives_si, sa_si, t_si, p_si
real*8 g, g_t, g_p, g_s
real*8 g_ss, g_st, g_sp
real*8 g_tt, g_pp, g_tp
real*8 h

!IN ONE CASE WE DO NOT AT ALL NEED TO COMPUTE THE GIBBS FUNCTION:
if(drv_s == 0 .and. drv_eta == 1 .and. drv_p == 0) then
  sea_s_eta_p_derivatives_si = t_si
  return
end if

sea_s_eta_p_derivatives_si = errorreturn

select case (drv_s)

  case (0)
    select case (drv_eta)
    
      case (0)
        select case (drv_p)
        
          case (0)
            g = sea_g_si(0, 0, 0, sa_si, t_si, p_si)
            if(g == errorreturn) go to 999
            g_t = sea_g_si(0, 1, 0, sa_si, t_si, p_si)
            if(g_t == errorreturn) go to 999
            h = g - t_si * g_t                           !H
          
          case (1)
            g_p = sea_g_si(0, 0, 1, sa_si, t_si, p_si)
            if(g_p == errorreturn) go to 999
            h = g_p                                      !DH/DP
                  
          case (2)
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
            if(g_tp == errorreturn) go to 999
            g_pp = sea_g_si(0, 0, 2, sa_si, t_si, p_si)
            if(g_pp == errorreturn) go to 999
            h = (g_tt * g_pp - g_tp ** 2) / g_tt         !D2H/DP2
                  
          case default; go to 999
        end select
        
      case (1)
        select case (drv_p)
          case (0); h = t_si                             !DH/DETA, HAS ALREADY BEEN HANDLED INITIALLY
          
          case (1)
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
            if(g_tp == errorreturn) go to 999
            h = -g_tp / g_tt                             !D2H/DETADP
             
          case default; go to 999
        end select

      case (2)
        select case (drv_p)
          case (0) 
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            h = -1d0 / g_tt                              !D2H/DETA2
                  
          case default; go to 999
        end select

      case default; go to 999
    end select
    
  case (1)
    select case (drv_eta)
    
      case (0)
        select case (drv_p)
        
          case (0)
            g_s = sea_g_si(1, 0, 0, sa_si, t_si, p_si)
            if(g_s == errorreturn) go to 999
            h = g_s                                      !DH/DS
          
          case (1)
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            g_st = sea_g_si(1, 1, 0, sa_si, t_si, p_si)
            if(g_st == errorreturn) go to 999
            g_sp = sea_g_si(1, 0, 1, sa_si, t_si, p_si)
            if(g_sp == errorreturn) go to 999
            g_tp = sea_g_si(0, 1, 1, sa_si, t_si, p_si)
            if(g_tp == errorreturn) go to 999
            h = (g_tt * g_sp - g_st * g_tp) / g_tt       !D2H/DSDP
          
          case default; go to 999
        end select

      case (1)
        select case (drv_p)
        
          case (0)
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            g_st = sea_g_si(1, 1, 0, sa_si, t_si, p_si)
            if(g_st == errorreturn) go to 999
            h = -g_st / g_tt                             !D2H/DSDETA
          
          case default; go to 999
        end select

      case default; go to 999
    end select

  
  case (2)
    select case (drv_eta)
    
      case (0)
        select case (drv_p)
        
          case (0)
            g_tt = sea_g_si(0, 2, 0, sa_si, t_si, p_si)
            if(g_tt == errorreturn) go to 999
            if(g_tt == 0d0) go to 999
            g_st = sea_g_si(1, 1, 0, sa_si, t_si, p_si)
            if(g_st == errorreturn) go to 999
            g_ss = sea_g_si(2, 0, 0, sa_si, t_si, p_si)
            if(g_ss == errorreturn) go to 999
            h = (g_tt * g_ss - g_st ** 2) / g_tt         !D2H/DS2
          
          case default; go to 999
        end select

      case default; go to 999
    end select

  case default; go to 999
end select

sea_s_eta_p_derivatives_si = h

999 continue

end function

!=========================================================================
function pottemp_iteration(sa_si, eta_si, p_si, t_si, maxit, eps)
!==========================================================================

!OUTPUT:
!THETA =   POTENTIAL TEMPERATURE OF SEAWATER IN K
!          I.E. THE TEMPERATURE THAT SOLVES ETA_SI = SEA_ENTROPY_SI(SA_SI, T_SI, P_SI)

!INPUT:
!SA_SI     ABSOLUTE SALINITY IN KG/KG
!ETA_SI    ENTROPY IN J/(KG K)
!P_SI      ABSOLUTE (REFERENCE) PRESSURE IN PA
!T_SI      ABSOLUTE (POTENTIAL) TEMPERATURE IN K, INITIAL VALUE
!MAXIT     MAX. NUMBER OF ITERATIONS
!EPS       REQUIRED TOLERANCE IN K

implicit none
integer i, maxit
real*8 pottemp_iteration, sa_si, eta_si, p_si, t_si, eps
real*8 s, theta, cp, dt

pottemp_iteration = errorreturn

if(eps <= 0d0) return
if(maxit <= 0d0) return

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax) return
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  sa_si > sal_smax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) return
else
if(sa_si < 0d0 .or. sa_si > 1d0 .or. &
  t_si <= 0d0 .or. &
  p_si <= 0d0) return
endif

check_limits = check_limits - 1

theta = t_si

do i = 1, maxit

  !get entropy and its first derivative for newton iteration
  s = sea_entropy_si(sa_si, theta, p_si)
  cp = sea_cp_si(sa_si, theta, p_si)
  
  if(s == errorreturn) go to 999
  if(cp == errorreturn) go to 999
  if(cp <= 0d0) go to 999
  
  !next temperature improvement step
  dt = theta * (eta_si - s) / cp
  theta = theta + dt
  if(theta <= 0d0) go to 999
  
  if(abs(dt) < eps) then
    pottemp_iteration = theta
    go to 999
  end if
  
enddo

999 continue

check_limits = check_limits + 1

end function

!==========================================================================
subroutine init_it_ctrl_pottemp()
!==========================================================================

implicit none

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR TEMPERATURE ITERATION
ctrl_loop_maximum = 100
ctrl_mode_pottemp = 0           !default: theta = t0 + 273*eta/4000
ctrl_init_pottemp = 273.15d0
ctrl_eps_exit_pottemp = 0.0001d0  !default = 0.1 mk

end subroutine

!==========================================================================
subroutine set_it_ctrl_pottemp(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!POTENTIAL TEMPERATURE FROM REFERENCE PRESSURE

!KEY                    VALUE
!'IT_STEPS     '        0           SET ITERATION NUMBER TO DEFAULT (100)
!'IT_STEPS     '        N > 0       SET ITERATION NUMBER TO N
!'IT_STEPS     '       -1           DO NOT ITERATE, USE INITIAL VALUE

!'INIT_THETA   '        0           USE DEFAULT POTENTIAL TEMPERATURE THETA = T0 + 273*ETA/4000 TO START
!'INIT_THETA   '        T > 0       USE VALUE T AS POTENTIAL TEMPERATURE TO START

!'TOL_THETA    '        0           USE DEFAULT EXIT ACCURACY FOR POTENTIAL TEMPERATURE (0.1 MK)
!'TOL_THETA    '        EPS > 0     USE EPS AS EXIT ACCURACY FOR POTENTIAL TEMPERATURE

implicit none
character*13 key
real*8 value

call init_it_ctrl_pottemp

select case (trim(key))

  case ('it_steps')               !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100     !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = nint(value)
    endif

  case ('init_theta')             !START THETA
    if(value == 0d0)then
      ctrl_mode_pottemp = 0       !DEFAULT: THETA = T0 + ETA/4000
    elseif(value > 0d0)then
      ctrl_mode_pottemp = 1
      ctrl_init_pottemp = value
    endif

  case ('tol_theta')              !REQUIRED THETA TOLERANCE
    if(value == 0d0)then
      ctrl_eps_exit_pottemp = 0.0001d0 !DEFAULT = 0.1 MK
    elseif(value > 0d0)then
      ctrl_eps_exit_pottemp = value
    endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

end module sea_3b

