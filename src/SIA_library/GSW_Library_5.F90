module gsw_library_5

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS SEAWATER (GSW) LIBRARY, A SUBSET OF
!TAYLORED FUNCTIONS BASED ON THE GENERAL SEA-ICE-AIR (SIA) LIBRARY

!#########################################################################

!IMPLEMENTATION IN FORTRAN BY D.G. WRIGHT AND D. JACKETT

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES
!    CONSTANTS_0,  FILE CONSTANTS_0.F90  
!    CONVERT_0,    FILE CONVERT_0.F90
!    MATHS_0,      FILE MATHS_0.F90  
!    FLU_1,        FILE FLU_1.F90  
!    SAL_1,        FILE SAL_1.F90  
!    SAL_2,        FILE SAL_2.F90  
!    FLU_3A,       FILE FLU_3A.F90  
!    SEA_3A,       FILE SEA_3A.F90  
!    SEA_3B,       FILE SEA_3B.F90  
!    LIQ_F03_5,    FILE LIQ_F03_5.F90  
!    SEA_5A,       FILE SEA_5A.F90  

!#########################################################################

use constants_0
use convert_0
use maths_0
use flu_1
use sal_1
use sal_2
use flu_3a
use sea_3a
use sea_3b
use liq_f03_5
use sea_5a

implicit none
private

character*16, private :: version = '1 Feb 2011'

public :: gsw_alpha_ct, gsw_alpha_pt0, &
          gsw_alpha_t, gsw_asal_from_psal, &
          gsw_beta_ct, gsw_beta_pt0, gsw_beta_t, &
          gsw_cabb_ct, gsw_cabb_pt0, &
          gsw_cp, gsw_ctmp_from_ptmp0, gsw_dens, &
          gsw_enthalpy, gsw_entropy, gsw_g, &
          gsw_kappa, gsw_kappa_t, gsw_pden, &
          gsw_psal_from_asal, gsw_ptmp, &
          gsw_ptmp0_from_ctmp, gsw_specvol, gsw_svel, &
          gsw_thrmb_ct, gsw_thrmb_pt0



!============================================================================
! This is the gsw thermodynamic Gibbs function library implemented in Fortran
!============================================================================
!
! consisting of the sum of pure water and saline components defined by 
!
!	Feistel, R., 2003: A new extended Gibbs thermodynamic potential of seawater,  
!                      Progr. Oceanogr., 58, 43-114.
!
! for the pure water component and by 
!
!	Feistel, R., 2008: A Gibbs function for seawater thermodynamics for -6 to 80 C  
!                      and salinity up to 120 G/KG, Deep-Sea Res. I, 55, 1639-1671, and  
!     
!  	IAPWS 2008: Release on the IAPWS Formulation 2008 for the Thermodynamic Properties of Seawater, 
!               The International Association for the Properties of Water and Steam, 
!               Berlin, Germany, September 2008, available at http://www.iapws.org
!  	
! for the saline componentr., 
!
! The basic routines are copies of those written by David Jackett for the GSW library
! on the TEOS-10 site.  However, they been modified to use the routines and approaches 
! used in the SIA library in order to provide a cross-check between the two libraries.
!
! September 2009
!============================================================================

contains

!============================================================================
function gsw_g(drv_s, drv_t, drv_p, sa, t, p)
!============================================================================

! SEAWATER SPECIFIC GIBBS FREE ENERGY AND DERIVATIVES UP TO ORDER 2
!
! DRV_S               : ORDER OF S DERIVATIVE
! DRV_T               : ORDER OF T DERIVATIVE
! DRV_P               : ORDER OF P DERIVATIVE
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
!OUTPUT: gsw_g        : SPECIFIC GIBBS ENERGY OR ITS DERIVATIVE
!GSW_G IS THE S-T-P DERIVATIVE
!(D/DS)^DRV_S (D/DT)^DRV_T (D/DP)^DRV_P G(S,T,P)
!OF THE SPECIFIC GIBBS ENERGY OF SEAWATER, G(S,T,P), IN J/KG
!T AND P DERIVATIVES ARE IN UNITS OF 1/K AND 1/PA
!S DERIVATIVES ARE IN UNITS OF 1/(G/KG)

!CHECK VALUES:
!GSW_G(0,0,0,35,26.85,0d0) = -5113.70064124
!GSW_G(1,0,0,35,26.85,0d0) =  78.5928261339
!GSW_G(0,1,0,35,26.85,0d0) = -374.452000830
!GSW_G(0,0,1,35,26.85,0d0) =  0.977858058750E-03
!GSW_G(2,0,0,35,26.85,0d0) =  2.24755137017
!GSW_G(1,1,0,35,26.85,0d0) =  0.789935187192
!GSW_G(1,0,1,35,26.85,0d0) = -0.716680931996E-06
!GSW_G(0,2,0,35,26.85,0d0) = -13.3358337534
!GSW_G(0,1,1,35,26.85,0d0) =  0.304607508052E-06
!GSW_G(0,0,2,35,26.85,0d0) = -0.410939723950E-12

implicit none
integer drv_s, drv_t, drv_p
real*8 gsw_g, sa, t, p, sa_si, t_si, p_si
real*8 gw, gs

sa_si = sa * 1d-3
t_si = t + 273.15d0
p_si = 1d4 * p + 101325d0

gsw_g = errorreturn

if(sa_si < 0d0 .or. sa_si >= 1d0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

if(drv_s < 0 .or. drv_s > 3) return
if(drv_t < 0) return
if(drv_p < 0) return
if(drv_s + drv_t + drv_p > 2) return

!WATER PART
gw = 0d0
if(drv_s == 0) then
  gw = fit_liq_g_f03_si(drv_t, drv_p, t_si, p_si)
  if(gw == errorreturn) return
end if

!SALINE PART
gs = sal_g_si(drv_s, drv_t, drv_p, sa_si, t_si, p_si)
if(gs == errorreturn) return

gsw_g = gw + gs
    
!It appears that the gsw library works with gsw_g in G/KG, DEGC, Pa
gsw_g = gsw_g * (1d-3)**drv_s

return
end function


!============================================================================
function gsw_alpha_t(sa, t, p)
!============================================================================

! THERMAL EXPANSION COEFFICIENT OF SEAWATER 
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
!
! GSW_ALPHA_T         : THERMAL EXPANSION COEFFICIENT      [1/K]
!                       WRT (IN SITU) TEMPERATURE

!CHECKVALUE: gsw_alpha_t(35.7d0, 25.5d0, 1023d0) = 0.000309837839319d0

implicit none
real*8 sa, t, p, gsw_alpha_t, gtp, gp

gsw_alpha_t = errorreturn

gtp = gsw_g(0, 1, 1, sa, t, p)
if(gtp == 0d0) return
gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return

gsw_alpha_t = gtp / gp

end function


!============================================================================
function gsw_beta_t(sa, t, p)
!============================================================================

! HALINE CONTRACTION COEFFICIENT OF SEAWATER
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : HALINE CONTRACTION COEFFICIENT     [KG/G]

!CHECKVALUE: gsw_beta_t(35.7d0, 25.5d0, 1023d0) = 0.000725729797839d0

implicit none
real*8 sa, t, p, gsw_beta_t, gsp, gp

gsw_beta_t = errorreturn

gsp = gsw_g(1, 0, 1, sa, t, p)
if(gsp == errorreturn) return
gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return

gsw_beta_t = -gsp / gp

return
end function


!============================================================================
function gsw_cp(sa, t, p)
!============================================================================

! ISOBARIC HEAT CAPACITY OF SEAWATER
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : HEAT CAPACITY                      [J/(KG K)]

!CHECKVALUE: gsw_cp(35.7d0, 25.5d0, 1023d0) = 3974.42541260

implicit none
real*8 sa, t, p, gsw_cp, gtt

gsw_cp = errorreturn

gtt = gsw_g(0, 2, 0, sa, t, p)
if(gtt == errorreturn)return

gsw_cp = -(T + 273.15d0) * gtt

return
end function


!============================================================================
function gsw_dens(sa, t, p)
!============================================================================

! DENSITY OF SEAWATER
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : DENSITY                            [KG/M^3]

!CHECKVALUE: gsw_dens(35.7d0, 25.5d0, 1023d0) = 1027.95249316d0

implicit none
real*8 sa, t, p, gsw_dens, gp

gsw_dens = errorreturn

gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return

gsw_dens = 1d0 / gp

return
end function


!============================================================================
function gsw_enthalpy(sa, t, p)
!============================================================================

! SPECIFIC ENTHALPY OF SEAWATER
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEG C]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : SPECIFIC ENTHALPY                  [J/KG]

!CHECKVALUE: gsw_enthalpy(35.7d0, 25.5d0, 1023d0) = 110776.712409d0

implicit none
real*8 sa, t, p, gsw_enthalpy, g, gt

gsw_enthalpy = errorreturn

g = gsw_g(0, 0, 0, sa, t, p)
if(g == errorreturn) return
gt = gsw_g(0, 1, 0, sa, T, P)
if(gt == errorreturn) return

gsw_enthalpy = g - (t + 273.15d0) * gt


return
end function


!============================================================================
function gsw_entropy(sa, t, p)
!============================================================================

! SPECIFIC ENTROPY OF SEAWATER
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : SPECIFIC ENTROPY                   [J/(KG K)]

!CHECKVALUE: gsw_entropy(35.7d0, 25.5d0, 1023d0) = 352.818797715d0

implicit none
real*8 sa, t, p, gsw_entropy, gt

gsw_entropy = errorreturn

gt = gsw_g(0, 1, 0, sa, t, p)
if(gt == errorreturn) return

gsw_entropy = -gt


return
end function


!============================================================================
function gsw_kappa(sa, t, p)
!============================================================================

! ISENTROPIC COMPRESSIBILITY OF SEAWATER
!
! SA                  : ASOLUTE SALINITY                   [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : ISENTROPIC COMPRESSIBILITY         [1/DBAR]

!CHECKVALUE: gsw_kappa(35.7d0, 25.5d0, 1023d0) = 4.03386268546d-6

implicit none
real*8 sa, t, p, gp, gpp, gtt, gtp, gsw_kappa

gsw_kappa = errorreturn

gtt = gsw_g(0, 2, 0, sa, t, p)
if(gtt == errorreturn .or. gtt >= 0) return
gtp = gsw_g(0, 1, 1, sa, t, p)
if(gtp == errorreturn) return
gpp = gsw_g(0, 0, 2, sa, t, p)
if(gpp == errorreturn) return
gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return

gsw_kappa = 10000d0 * (gtp * gtp - gtt * gpp) / (gp * gtt)

return
end function


!============================================================================
function gsw_kappa_t(sa, t, p)
!============================================================================

! ISOTHERMAL COMPRESSIBILITY OF SEAWATER
!
! SA                  : ASOLUTE SALINITY                   [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : ISOTHERMAL COMPRESSIBILITY         [1/DBAR]

!CHECKVALUE: gsw_kappa_t(35.7d0, 25.5d0, 1023d0) = 4.10403794615d-6

implicit none
real*8 sa, t, p, gsw_kappa_t, gpp, gp

gsw_kappa_t = errorreturn

gpp = gsw_g(0, 0, 2, sa, t, p)
if(gpp == errorreturn) return
gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return

gsw_kappa_t = -10000d0 * gpp / gp

return
end function



!============================================================================
function gsw_specvol(sa, t, p)
!============================================================================

! SPECIFIC VOLUME OF SEAWATER
!
! SA                  : ASOLUTE SALINITY                   [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : SPECIFIC VOLUME                    [M^3/KG]  !dbg invert units

!CHECKVALUE: gsw_specvol(35.7d0, 25.5d0, 1023d0) = 0.000972807602158d0

implicit none
real*8 sa, t, p, gsw_specvol

gsw_specvol = gsw_g(0,0,1,sa,t,p)

return
end function


!============================================================================
function gsw_svel(sa, t, p)
!============================================================================

! SOUND SPEED OF SEAWATER
!
! sa                  : ASOLUTE SALINITY                   [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! 
! RESULT              : SOUND SPEED                        [M/S]

!CHECKVALUE: gsw_svel(35.7d0, 25.5d0, 1023d0) = 1552.93372863d0

implicit none
real*8 sa, t, p, gp, gtt, gtp, gpp, gsw_svel

gsw_svel = errorreturn

gtt = gsw_g(0, 2, 0, sa, t, p)
if(gtt == errorreturn .or. gtt >= 0d0) return
gtp = gsw_g(0, 1, 1, sa, t, p)
if(gtp == errorreturn) return
gpp = gsw_g(0, 0, 2, sa, t, p)
if(gpp == errorreturn) return
gp = gsw_g(0, 0, 1, sa, t, p)
if(gp == errorreturn .or. gp <= 0d0) return
if(gtp * gtp >= gtt * gpp) return

gsw_svel = gp * sqrt(gtt/(gtp**2d0 - gtt*gpp))
                                      
return
end function


!============================================================================
function gsw_ctmp_from_ptmp0(sa, tpot)
!============================================================================

!CONSERVATIVE TEMPERATURE FROM POTENTIAL TEMPERATURE OF SEAWATER  
!BASED ON GSW_CT_FROM_PT0_SI BUT USES GSW_G
!
!SA         : ABSOLUTE SALINITY                        [G/KG]
!TPOT       : POTENTIAL TEMPERATURE WITH               [DEGC]
!             REFERENCE PRESSURE OF 101325 PA
!
!RESULT     : CONSERVATIVE TEMPERATURE                 [DEGC]

!CHECKVALUE: gsw_ctmp_from_ptmp0(35.7d0, 25.5d0) = 25.4805463842d0

implicit none
real*8 sa, tpot, gsw_ctmp_from_ptmp0
real*8 cp0, pr, ct

gsw_ctmp_from_ptmp0 = errorreturn

cp0 = 3991.86795711963d0     ! IN CONCRETE 02/12/08

pr = 0d0

ct = gsw_enthalpy(sa, tpot, pr)
if (ct == errorreturn) return

gsw_ctmp_from_ptmp0 = ct / cp0

end function


!============================================================================
function gsw_ptmp0_from_ctmp(sa, ct)
!============================================================================

! POTENTIAL TEMPERATURE OF SEAWATER FROM CONSERVATIVE TEMPERATURE
! THIS CODE IS CONSISTENT WITH SEA_PTMP0_FROM_CTMP_SI IN SEA_5A
! BUT USED GSW ROUTINES AND DIFFERENT UNITS
!
! SA                  : ASOLUTE SALINITY                   [G/KG]
! CT                  : CONSERVATIVE TEMPERATURE           [DEGC]
!
! RESULT              : POTENTIAL TEMPERATURE WITH         [DEGC]
!                       REFERENCE PRESSURE OF  0 DBAR

!CHECKVALUE: gsw_ptmp0_from_ctmp(35.7d0, 25.4805463842239) = 25.5d0

implicit none
real*8 SA, ct, gsw_ptmp0_from_ctmp
integer nloops, n
real*8 s1, ct1, p0, cp0 
real*8 a0, a1, a2, a3, a4, a5, b0, b1, b2, b3
real*8 a5ct, b3ct, ct_factor, th0_num, rec_th0_den
real*8 th0, ct0, dth_dct, theta, dct, dct_dth
real*8 gtt

gsw_ptmp0_from_ctmp = errorreturn

cp0 = 3991.86795711963d0    

s1 = SA*35.d0/35.16504d0; ct1 = ct; p0 = 0.d0

a0 = -1.446013646344788d-2;     b0 =  1.000000000000000d+0
a1 = -3.305308995852924d-3;     b1 =  6.506097115635800d-4
a2 =  1.062415929128982d-4;     b2 =  3.830289486850898d-3
a3 =  9.477566673794488d-1;     b3 =  1.247811760368034d-6
a4 =  2.166591947736613d-3
a5 =  3.828842955039902d-3

a5ct = a5*ct1; b3ct = b3*ct1
ct_factor = (a3+a4*s1+a5ct)
th0_num = a0+s1*(a1+a2*s1)+ct1*ct_factor
rec_th0_den = 1.d0/(b0+b1*s1+ct1*(b2+b3ct))
th0 = th0_num*rec_th0_den 
ct0 = gsw_ctmp_from_ptmp0(sa,th0)
if(ct0 == errorreturn) return
dth_dct = (ct_factor+a5ct-(b2+b3ct+b3ct)*th0)*rec_th0_den 
theta = th0-(ct0-ct)*dth_dct

nloops = 1                  ! default

!    NOTE: nloops = 1 gives theta with a maximum error of & 

n = 0; 

do while(n.le.nloops)
    dct = gsw_ctmp_from_ptmp0(sa,theta)
    if(dct == errorreturn) return
    dct = dct - ct
    gtt = gsw_g(0, 2, 0, sa, theta, p0)
    if(gtt == errorreturn .or. gtt >= 0d0) return
    dct_dth = -(theta+273.15d0)*gtt / cp0
    if(dct_dth == 0d0) return
    theta = theta - dct/dct_dth
    n = n+1
end do

gsw_ptmp0_from_ctmp = theta

end function


!============================================================================
function gsw_ptmp(sa, t, p, pr)
!============================================================================
                               
!THIS FUNCTION COMPUTES POTENTIAL TEMPERATURE OF SEAWATER
!BASED ON SEA_POTTEMP_SI BUT USES GSW_G

!OUTPUT:
!THETA(S,T,P,PR) POTENTIAL TEMPERATURE OF SEAWATER IN DEGC

!INPUT:
!SA     ABSOLUTE SALINITY IN G/KG
!T      ABSOLUTE IN-SITU TEMPERATURE IN DEGC
!P      SEA PRESSURE IN DBAR
!PR     REFERENCE SEA PRESSURE IN DBAR

!CHECKVALUE: gsw_ptmp(35.7d0, 25.2720983155409d0, 0d0, 1023d0) = 25.5
!CHECKVALUE: gsw_ptmp(35.7d0, 25.5d0, 1023d0, 0d0) = 25.2720983155d0

implicit none
integer maxit
real*8 gsw_ptmp, sa, t, p, pr
real*8 sp, s, eps, theta

gsw_ptmp = errorreturn

if(p == pr) then
  gsw_ptmp = t
  return
end if

s = gsw_entropy(sa, t, p)  !SPECIFIC ENTROPY IN-SITU

!use David's initial estimate for theta
sp = sa * 35d0 / 35.16504d0
theta = t+(p-pr)*( 8.65483913395442d-6  - &
            sp  *  1.41636299744881d-6  - &
         (p+pr) *  7.38286467135737d-9  + &
             t  *(-8.38241357039698d-6  + &
            sp  *  2.83933368585534d-8  + &
             t  *  1.77803965218656d-8  + &
         (p+pr) *  1.71155619208233d-10))

maxit = 100
eps = 1d-12

!RUN THE ITERATION
gsw_ptmp = gsw_pottemp_iteration(sa, s, pr, theta, maxit, eps)

end function


!============================================================================
function gsw_pden(sa, t, p, pr)
!============================================================================

! POTENTIAL DENSITY OF SEAWATER
!
! SA                  : ASOLUTE SALINITY                   [G/KG]
! T                   : TEMPERATURE                        [DEGC]
! P                   : SEA (GAUGE) PRESSURE               [DBAR]
! PR                  : REFERENCE (GAUGE) PRESSURE         [DBAR]
! 
! RESULT              : POTENTIAL DENSITY                  [KG/M^3]

!CHECKVALUE: gsw_pden(35.7d0, 25.5d0, 1023d0, 0d0) = 1023.66254941d0

implicit none
real*8 sa, t, p, pr, gsw_pden, theta

gsw_pden = errorreturn

theta = gsw_ptmp(sa,t,p,pr);
if(theta == errorreturn) return

gsw_pden = gsw_dens(sa,theta,pr);

return
end function


!===========================================================================
function gsw_pottemp_iteration(sa, eta, p, theta, maxit, eps)
!============================================================================

!BASED ON POTTEMP_ITERATION FROM SEA_3B

!OUTPUT:
!THETA =   POTENTIAL TEMPERATURE OF SEAWATER IN DEGC
!          I.E. THE TEMPERATURE THAT SOLVES ETA_SI = GSW_ENTROPY_SI(SA, T, P)

!INPUT:
!SA     ABSOLUTE SALINITY IN G/KG
!ETA    ENTROPY IN J/(KG K)
!P      (REFERENCE) SEA PRESSURE IN DBAR
!T      (POTENTIAL) SEA TEMPERATURE IN DEGC, INITIAL VALUE
!MAXIT  MAX. NUMBER OF ITERATIONS
!EPS    REQUIRED TOLERANCE IN K

implicit none
integer i, maxit
real*8 gsw_pottemp_iteration, sa, eta, p, eps
real*8 s, theta, cp, dt

gsw_pottemp_iteration = errorreturn

if(eps <= 0d0) return
if(maxit <= 0d0) return

do i = 1, maxit

  !get entropy and its first derivative for newton iteration
  s = gsw_entropy(sa, theta, p)
  if(s == errorreturn) return
  cp = gsw_cp(sa, theta, p)
  if(cp == errorreturn .or. cp == 0d0) return
  
  !next temperature improvement step
  dt = (273.15d0 + theta) * (eta - s) / cp
  theta = theta + dt
  if(theta <= -273.15d0) return
  
  if(abs(dt) < eps) then
    gsw_pottemp_iteration = theta
    return
  end if
  
enddo

end function


!==========================================================================
function gsw_alpha_ct(sa, t, p)
!==========================================================================

! THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE

! SA_GKG                : ABSOLUTE SALINITY                  [G/KG]
! T_C                   : IN SITU TEMPERATURE                [DEG C]
! P_DBAR                : SEA PRESSURE                       [DBAR]
!
! GSW_ALPHA_CT        : THERMAL EXPANSION COEFFICIENT      [1/DEG C]
!                         WRT CONSERVATIVE TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!GSW_ALPHA_CT(35, 20, 1000) = 2.69418609861E-04

implicit none
real*8 sa, t, p, gsw_alpha_ct
real*8 sa_si, p_si, t_si

!CONVERT PRACTICAL INPUT UNITS TO BASIC SI UNITS
sa_si = sa * 0.001d0
p_si = SO_pressure_si + p * 10000d0
t_si = SO_temperature_si + t

!AND CALL THE SI VERSION
gsw_alpha_ct = sea_alpha_ct_si(sa_si, t_si, p_si)

end function


!==========================================================================
function gsw_alpha_pt0(sa,t,p)
!==========================================================================

!THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE

!SA_GKG                 : ABSOLUTE SALINITY                  [G/KG]
!T_C                    : IN SITU TEMPERATURE                [DEG C]
!P_DBAR                 : SEA PRESSURE                       [DBAR]
!
!GSW_ALPHA_PT0_SI       : THERMAL EXPANSION COEFFICIENT      [1/DEG C]
!                         WRT POTENTIAL TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!GSW_ALPHA_PT0(35, 20, 1000) = 2.6975373331|7E-04

implicit none
real*8 gsw_alpha_pt0, sa, t, p
real*8 sa_si, t_si, p_si

!CONVERT PRACTICAL INPUT UNITS TO BASIC SI UNITS
sa_si = sa * 0.001d0
p_si = SO_pressure_si + p * 10000d0
t_si = SO_temperature_si + t

!and call the SI version
gsw_alpha_pt0 = sea_alpha_pt0_si(sa_si, t_si, p_si)

end function


!==========================================================================
function gsw_beta_ct(sa, t, p)
!==========================================================================

! HALINE CONTRACTION COEFFICIENT OF SEAWATER WRT CONSTANT CONSERVATIVE TEMPERATURE

! SA_GKG                 : ABSOLUTE SALINITY                  [G/KG]
! T_C                    : IN SITU TEMPERATURE                [DEG C]
! P_DBAR                 : SEA PRESSURE                       [DBAR]

! RESULT                 : HALINE CONTRACTION COEFFICIENT     [KG/G]
!                          WRT CONSTANT CONSERVATIVE TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!GSW_BETA_CT(35, 20, 1000) = 7.23213672954E-04

real*8 gsw_beta_ct, sa, t, p
real*8 sa_si, p_si, t_si, bt

gsw_beta_ct = ErrorReturn

!CONVERT PRACTICAL INPUT UNITS TO BASIC SI UNITS
sa_si = sa * 0.001d0
p_si = SO_pressure_si + p * 10000d0
t_si = SO_temperature_si + t

!CALL THE SI VERSION
bt = sea_beta_ct_si(sa_si, t_si, p_si)
if (bt == errorreturn) return 

!AND CONVERT THE RESULT TO KG/G, THE UNIT OF 1/SA_GKG
gsw_beta_ct = bt * 0.001d0 !1 KG/KG = 0.001 KG/G

end function


!==========================================================================
function gsw_beta_pt0(sa, t, p)
!==========================================================================

! HALINE CONTRACTION COEFFICIENT OF SEAWATER WRT CONSTANT POTENTIAL TEMPERATURE

! SA_GKG                 : ABSOLUTE SALINITY                  [G/KG]
! T_C                    : IN SITU TEMPERATURE                [\260C]
! P_DBAR                 : SEA PRESSURE                       [DBAR]

! RESULT                 : HALINE CONTRACTION COEFFICIENT     [KG/G]
!                          WRT CONSTANT POTENTIAL TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!GSW_BETA_PT0(35, 20, 1000) = 7.31582583383E-04

real*8 gsw_beta_pt0, sa, t, p
real*8 sa_si, p_si, t_si, bt

gsw_beta_pt0 = ErrorReturn

!CONVERT PRACTICAL INPUT UNITS TO BASIC SI UNITS
sa_si = sa * 0.001d0
p_si = SO_pressure_si + p * 10000d0
t_si = SO_temperature_si + t

!CALL THE SI VERSION
bt = sea_beta_pt0_si(sa_si, t_si, p_si)
if (bt == errorreturn) return

!AND CONVERT THE RESULT TO KG/G, THE UNIT OF 1/SA_GKG
gsw_beta_pt0 = bt * 0.001d0 !1 KG/KG = 0.001 KG/G

end function







!==========================================================================
function gsw_cabb_ct(sa,t,p)
!==========================================================================

!CABBELING COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE 

!SA                  : ABSOLUTE SALINITY                  [G/KG]
!T                   : IN SITU TEMPERATURE                [DEGC]
!P                   : SEA (GAUGE) PRESSURE               [DBAR]
!
!RESULT              : CABBELING COEFFICIENT              [1/(DEGC**2)] 
!                      WRT CONSERVATIVE TEMPERATURE

!check value:
!gsw_cabb_ct(35, 20, 1d3) = 0.8969073|83083E-05  

implicit none
real*8 sa, t, p, sa_si, t_si, p_si, gsw_cabb_ct

sa_si = 1d-3 * sa
t_si = t + 273.15d0
p_si = 101325d0 + 1d4 * p

gsw_cabb_ct = sea_cabb_ct_si(sa_si, t_si, p_si)

end function


!==========================================================================
function gsw_cabb_pt0(sa,t,p)
!==========================================================================

!CABBELING COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE 

!SA                  : ABSOLUTE SALINITY                  [G/KG]
!T                   : IN SITU TEMPERATURE                [DEGC]
!P                   : SEA (GAUGE) PRESSURE               [DBAR]
!
!RESULT              : CABBELING COEFFICIENT              [1/(DEGC**2)] 
!                      WRT POTENTIAL TEMPERATURE

!check value:
!gsw_cabb_pt0(35, 20, 1d3) = 0.8759631|54048E-05  

implicit none
real*8 sa, t, p, sa_si, t_si, p_si, gsw_cabb_pt0

sa_si = 1d-3 * sa
t_si = t + 273.15d0
p_si = 101325d0 + 1d4 * p

gsw_cabb_pt0 = sea_cabb_pt0_si(sa_si, t_si, p_si)

end function


!==========================================================================
function gsw_thrmb_ct(sa,t,p)
!==========================================================================

!THERMOBARIC COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE 
!
!SA                  : ABSOLUTE SALINITY                  [G/KG]
!T                   : IN SITU TEMPERATURE                [DEGC]
!P                   : SEA (GAUGE) PRESSURE               [DBAR]
!
!RESULT              : THERMOBARIC COEFFICIENT            [1/(DEGC DBAR)]
!                      WRT CONSERVATIVE TEMPERATURE

!check value:
!gsw_thrmb_ct(35, 20, 1d3) = 0.1727083|65652E-07  

implicit none
real*8 sa, t, p, sa_si, t_si, p_si, gsw_thrmb_ct

sa_si = 1d-3 * sa
t_si = t + 273.15d0
p_si = 101325d0 + 1d4 * p

gsw_thrmb_ct = 1d4 * sea_thrmb_ct_si(sa_si, t_si, p_si)

end function

!==========================================================================
function gsw_thrmb_pt0(sa,t,p)
!==========================================================================

!THERMOBARIC COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE 
!
!SA                  : ABSOLUTE SALINITY                  [G/KG]
!T                   : IN SITU TEMPERATURE                [DEGC]
!P                   : SEA (GAUGE) PRESSURE               [DBAR]
!
!RESULT              : THERMOBARIC COEFFICIENT            [1/(DEGC DBAR)]
!                      WRT POTENTIAL TEMPERATURE

!check value:
!gsw_thrmb_pt0(35, 20, 1d3) = 0.1709450|45984E-07  

implicit none
real*8 sa, t, p, sa_si, t_si, p_si, gsw_thrmb_pt0

sa_si = 1d-3 * sa
t_si = t + 273.15d0
p_si = 101325d0 + 1d4 * p

gsw_thrmb_pt0 = 1d4 * sea_thrmb_pt0_si(sa_si, t_si, p_si)

end function

!=====================================================
function gsw_asal_from_psal(sp,lon0,lat0,p0)
!=====================================================

! CONVERT PRACTICAL SALINITY TO ABSOLUTE SALINITY 
!
! SP                  : PRACTICAL SALINITY                 [PSU, PSS]
! P0                  : SEA PRESSURE                       [DBAR]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! RESULT              : ABSOLUTE SALINITY                  [G/KG]

!CHECK VALUE:
!GSW_ASAL_FROM_PSAL(35.527515065427778, 201, -21, 1023) = 35.7
!GSW_ASAL_FROM_PSAL(35, 180, 40, 2d3) = 35.1890932889958041
!GSW_ASAL_FROM_PSAL(8, 20, 57, 0) = 8.10483771428571406

implicit none
real*8 sp, lon0, lat0, p_si, p0, gsw_asal_from_psal

p_si = 101325d0 + 1d4 * p0
gsw_asal_from_psal = sa_from_sp(sp, lon0, lat0, p_si)

return
end function

!=====================================================
function gsw_psal_from_asal(sa,lon0,lat0,p0)
!=====================================================

! CONVERT ABSOLUTE SALINITY TO PRACTICAL SALINITY 
!
! SA                  : ABSOLUTE SALINITY                  [G/KG]
! P0                  : SEA PRESSURE                       [DBAR]
! LONGS0              : LONGITUDE                          [DEG E]     
! LATS0               : LATITUDE                           [DEG N]
!
! GSW_PSAL_FROM_ASAL  : PRACTICAL SALINITY                 [PSU, PSS]

!CHECK VALUE
!GSW_PSAL_FROM_ASAL(35.7, 201, -21, 1023) = 35.527515065427778


implicit none
real*8 sa_si, sa, lon0, lat0, p_si, p0, gsw_psal_from_asal

p_si = 101325d0 + 1d4 * p0
sa_si = 1d-3 * sa
gsw_psal_from_asal = sp_from_sa(sa_si, lon0, lat0, p_si)

return
end function


end module gsw_library_5
