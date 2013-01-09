module sea_5a

!#########################################################################

!THIS MODULE IMPLEMENTS SEAWATER PROPERTIES RELATED TO POTENTIAL TEMPERATURE 
!AND CONSERVATIVE TEMPERATURE. 

!#########################################################################

!IMPLEMENTATION IN FORTRAN BY D.G. WRIGHT AND D.R. JACKETT
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
!     SEA_3C,      FILE SEA_3C.F90

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
use sea_3c

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_alpha_ct_si, sea_alpha_pt0_si, sea_alpha_t_si, &
          sea_beta_ct_si,  sea_beta_pt0_si,  sea_beta_t_si,  &
          sea_cabb_ct_si, sea_cabb_pt0_si, &
          sea_ctmp_from_ptmp0_si, sea_ptmp0_from_ctmp_si, &
          sea_thrmb_ct_si, sea_thrmb_pt0_si

contains

!==========================================================================
function sea_ctmp_from_ptmp0_si(sa_si,tpot_si)
!==========================================================================

!CONSERVATIVE TEMPERATURE FROM POTENTIAL TEMPERATURE OF SEAWATER  
!
!SA_SI         : ABSOLUTE SALINITY                        [KG/KG]
!TPOT_SI       : POTENTIAL TEMPERATURE WITH               [K]
!                REFERENCE PRESSURE OF 101325 PA
!
!RESULT        : CONSERVATIVE TEMPERATURE                 [K]

!CHECK VALUE:
!SEA_CTMP_FROM_PTMP0_SI(0.035, 300) =   0.300010069445E+03

implicit none
real*8 sa_si, tpot_si, sea_ctmp_from_ptmp0_si
real*8 cp0, pr_si, ct

sea_ctmp_from_ptmp0_si = errorreturn

cp0 = 3991.86795711963d0     ! IN CONCRETE 02/12/08

pr_si = 101325d0

ct = sea_enthalpy_si(sa_si, tpot_si, pr_si)
if (ct == errorreturn) return

sea_ctmp_from_ptmp0_si = 273.15d0 + ct / cp0

end function


!==========================================================================
function sea_ptmp0_from_ctmp_si(sa_si,ct_si)
!==========================================================================

!POTENTIAL TEMPERATURE OF SEAWATER FROM CONSERVATIVE TEMPERATURE
!
!SA_SI               : ABSOLUTE SALINITY                  [KG/KG]
!CT_SI               : CONSERVATIVE TEMPERATURE           [K]
!
!RESULT              : POTENTIAL TEMPERATURE WITH         [K]
!                      REFERENCE PRESSURE OF  101325 PA

!CHECK VALUE:
!SEA_PTMP0_FROM_CTMP_SI(0.035, 300.010069445) =   300

implicit none
integer nloops, n
real*8 sa_si, ct_si, s1, ct1, p0, sea_ptmp0_from_ctmp_si, cp0 
real*8 a0, a1, a2, a3, a4, a5, b0, b1, b2, b3
real*8 a5ct, b3ct, ct_factor, th0_num, rec_th0_den
real*8 th0, ct0, dth_dct, theta, dct, dct_dth, factor

sea_ptmp0_from_ctmp_si = errorreturn

cp0 = 3991.86795711963d0    
s1 = sa_si*35d0/35.16504d0; ct1 = ct_si-273.15d0; p0 = 101325d0

a0 = -1.446013646344788d-2;     b0 =  1.000000000000000d+0
a1 = -3.305308995852924d-3;     b1 =  6.506097115635800d-4
a2 =  1.062415929128982d-4;     b2 =  3.830289486850898d-3
a3 =  9.477566673794488d-1;     b3 =  1.247811760368034d-6
a4 =  2.166591947736613d-3
a5 =  3.828842955039902d-3

a5ct = a5*ct1; b3ct = b3*ct1
ct_factor = (a3+a4*s1+a5ct)
th0_num = a0+s1*(a1+a2*s1)+ct1*ct_factor
rec_th0_den = 1d0/(b0+b1*s1+ct1*(b2+b3ct))
th0 = th0_num*rec_th0_den

ct0 = sea_ctmp_from_ptmp0_si(sa_si,273.15d0+th0) 
if(ct0 == errorreturn) return

dth_dct = (ct_factor+a5ct-(b2+b3ct+b3ct)*th0)*rec_th0_den 
theta = th0-(ct0-ct_si)*dth_dct

nloops = 1   !DEFAULT
             !NOTE: NLOOPS = 1 GIVES THETA WITH A MAXIMUM ERROR OF & !dbg truncated sentence
n = 0; 
do while(n <= nloops)
  factor = sea_ctmp_from_ptmp0_si(sa_si,273.15d0+theta)
  if(factor == errorreturn) return
  dct = factor-ct_si

  factor = sea_g_si(0,2,0,sa_si,273.15d0+theta,p0)
  if(factor == errorreturn) return
  dct_dth = -(theta+273.15d0)*factor/cp0   

  theta = theta - dct/dct_dth
  n = n+1
end do

sea_ptmp0_from_ctmp_si = theta + 273.15d0

end function


!==========================================================================
function sea_alpha_ct_si(sa_si,t_si,p_si)
!==========================================================================

!THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!SEA_ALPHA_CT_SI        : THERMAL EXPANSION COEFFICIENT      [1/K]
!                         WRT CONSERVATIVE TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_ALPHA_CT_SI(0.035, 300, 1E8) = 3.91772847589E-04

implicit none
character *13 key
real*8 sa_si, t_si, p_si, sea_alpha_ct_si
real*8 pr_si, ah, cp0

cp0 = 3991.86795711963d0     ! IN CONCRETE 02/12/08

sea_alpha_ct_si = errorreturn

pr_si = 101325d0

key = "t"
ah = sea_eta_expansion_h_si(sa_si, t_si, p_si, pr_si, key)
if (ah == errorreturn) return

sea_alpha_ct_si = cp0 * ah

end function


!==========================================================================
function sea_alpha_pt0_si(sa_si,t_si,p_si)
!==========================================================================

!THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE

!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!SEA_ALPHA_PT0_SI       : THERMAL EXPANSION COEFFICIENT      [1/K]
!                         WRT POTENTIAL TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_ALPHA_PT0_SI(0.035, 300, 1E8) = 3.92515634064E-04

implicit none
character*13 key
real*8 sa_si, t_si, p_si, sea_alpha_pt0_si
real*8 pr_si

pr_si = 101325d0
key = "t"
sea_alpha_pt0_si = sea_eta_expansion_theta_si(sa_si, t_si, p_si, pr_si, key)

end function


!==========================================================================
function sea_alpha_t_si(sa_si,t_si,p_si)
!==========================================================================

!THERMAL EXPANSION COEFFICIENT OF SEAWATER WRT IN SITU TEMPERATURE
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!SEA_ALPHA_T_SI         : THERMAL EXPANSION COEFFICIENT      [1/K]
!                         WRT IN SITU TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_ALPHA_T_SI(0.035, 300, 1E8) = 3.73608885178E-04

implicit none
real*8 sa_si, t_si, p_si, sea_alpha_t_si

sea_alpha_t_si = sea_g_expansion_t_si(sa_si, t_si, p_si)

end function


!==========================================================================
function sea_beta_ct_si(sa_si,t_si,p_si)
!==========================================================================

!HALINE CONTRACTION COEFFICIENT OF SEAWATER WRT CONSTANT CONSERVATIVE TEMPERATURE
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!RESULT                 : HALINE CONTRACTION COEFFICIENT     [KG/KG]
!                         WRT CONSTANT CONSERVATIVE TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_BETA_CT_SI(0.035, 300, 1E8) = 0.649596383654

implicit none
character*13 key
real*8 sa_si, t_si, p_si, sea_beta_ct_si
real*8 pr_si

pr_si = 101325d0
key = "t"
sea_beta_ct_si = sea_eta_contraction_h_si(sa_si, t_si, p_si, pr_si, key)

end function


!==========================================================================
function sea_beta_pt0_si(sa_si,t_si,p_si)
!==========================================================================

!HALINE CONTRACTION COEFFICIENT OF SEAWATER WRT CONSTANT POTENTIAL TEMPERATURE 
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!RESULT                 : HALINE CONTRACTION COEFFICIENT     [KG/KG]
!                         WITH CONSTANT POTENTIAL TEMPERATURE

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_BETA_PT0_SI(0.035, 300, 1E8) = 0.663973579411

implicit none
character*13 key
real*8 sa_si, t_si, p_si, sea_beta_pt0_si
real*8 pr_si

pr_si = 101325d0
key = "t"
sea_beta_pt0_si = sea_eta_contraction_theta_si(sa_si, t_si, p_si, pr_si, key)

end function


!==========================================================================
function sea_beta_t_si(sa_si,t_si,p_si)
!==========================================================================

!HALINE CONTRACTION COEFFICIENT OF SEAWATER WRT IN SITU TEMPERATURE 
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : ABSOLUTE PRESSURE                  [PA]
!
!RESULT                 : HALINE CONTRACTION COEFFICIENT     [KG/KG]

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_BETA_T_SI(0.035, 300, 1E8) = 0.666238827368

implicit none
real*8 sa_si, t_si, p_si, sea_beta_t_si

sea_beta_t_si = sea_g_contraction_t_si(sa_si, t_si, p_si)

end function


!==========================================================================
function sea_cabb_ct_si(sa_si,t_si,p_si)
!==========================================================================

!CABBELING COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE 

!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : SEA (GAUGE) PRESSURE               [PA]
!
!RESULT                 : CABBELING COEFFICIENT              [1/(K**2)] 
!                         WRT CONSERVATIVE TEMPERATURE

!CHECK VALUE:
!SEA_CABB_CT_SI(0.035, 300, 1E5) = 8.612525|67438267E-06  !DEFAULT TOLERANCE
!SEA_CABB_CT_SI(0.035, 300, 1E5) = 8.6125256|4750277E-06  !TOLERANCE RESET TO -1D-12


implicit none
!character*13 key0 
!real*8 eps0
real*8 sa_si, t_si, p_si, sea_cabb_ct_si
real*8 pr0, th, th_l, th_u, t_l, t_u, alpha, beta, alpha_ct, ratio
real*8 ct, dct, ct_l, ct_u, dsa, sa_l, sa_u, alpha_sa, beta_sa, factor1, factor2

!EXPERIMENTATION WITH DIFFERENT ERROR TOLERANCES DID NOT LEAD TO SYSTEMATIC IMPROVEMENT
!key0="tol_liq_dens"
!eps0 = get_it_ctrl_density(key0)
!call set_it_ctrl_density(key0, -1d-12)

sea_cabb_ct_si = errorreturn

pr0 = 101325d0; 

th = sea_pottemp_si(sa_si,t_si,p_si,pr0);    ! NOTE: CT = CONST <==> THETA = CONST
if (th == errorreturn) return
ct = sea_ctmp_from_ptmp0_si(sa_si,th)  
if (ct == errorreturn) return

!CT DERIVATIVE: The modification made here is to ensure centered differences
dct = 1d-3;     
               
ct_l = ct - dct; ct_u = ct + dct; 

th_l = sea_ptmp0_from_ctmp_si(sa_si,ct_l)
if (th_l == errorreturn) return
th_u =  sea_ptmp0_from_ctmp_si(sa_si,ct_u)
if (th_u == errorreturn) return

t_l = sea_pottemp_si(sa_si,th_l,pr0,p_si)
if (t_l == errorreturn) return
t_u = sea_pottemp_si(sa_si,th_u,pr0,p_si)
if (t_u == errorreturn) return

alpha = sea_alpha_ct_si(sa_si,t_si,p_si)
if(alpha == errorreturn) return
beta = sea_beta_ct_si(sa_si,t_si,p_si)
if(beta == errorreturn .or. beta == 0d0) return
ratio = alpha/beta

factor1 = sea_alpha_ct_si(sa_si,t_u,p_si)
if(factor1 == errorreturn) return
factor2 = sea_alpha_ct_si(sa_si,t_l,p_si)
if(factor2 == errorreturn) return

!calculate d(alpha_ct)/dct
if (ct_u == ct_l) return
alpha_ct = (factor1-factor2)/(ct_u-ct_l)

!SA DERIVATIVES
dsa = min(sa_si,1d-5);
if(sa_si >= dsa)then
  sa_l = sa_si - dsa; sa_u = sa_si + dsa; 
elseif(sa_si >= 0d0)then
  sa_l = 0d0; sa_u = dsa; 
else
  print *, '****    sa_si<0 in sea_cabb_ct    ****' 
  sea_cabb_ct_si = errorreturn
  return
end if

if (sa_u == sa_l) return

factor1 = sea_alpha_ct_si(sa_u,t_si,p_si)
if(factor1 == errorreturn) return
factor2 = sea_alpha_ct_si(sa_l,t_si,p_si)
if(factor2 == errorreturn) return
alpha_sa = (factor1-factor2)/(sa_u-sa_l)

factor1 = sea_beta_ct_si(sa_u,t_si,p_si)
if(factor1 == errorreturn) return
factor2 = sea_beta_ct_si(sa_l,t_si,p_si)
if(factor2 == errorreturn) return
beta_sa =  (factor1-factor2)/(sa_u-sa_l)


!CABBELING COEFFICIENT
sea_cabb_ct_si = alpha_ct + ratio*(2d0*alpha_sa - ratio*beta_sa)

!call set_it_ctrl_density(key0, eps0)

end function


!==========================================================================
function sea_cabb_pt0_si(sa_si,t_si,p_si)
!==========================================================================

!CABBELING COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE 

!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : SEA (GAUGE) PRESSURE               [PA]
!
!RESULT                 : CABBELING COEFFICIENT              [1/(K**2)] 
!                         WRT POTENTIAL TEMPERATURE

!CHECK VALUE:
!SEA_CABB_PT0_SI(0.035, 300, 1E5) = 8.3387453|7690444E-06  !DEFAULT TOLERANCES
!SEA_CABB_PT0_SI(0.035, 300, 1E5) = 8.3387453|0126243E-06  !TOLERANCE RESET TO -1D-12



implicit none
!character*13 key0
!real*8 eps0
real*8 sa_si, t_si, p_si, sea_cabb_pt0_si
real*8 pr0, th, dth, th_l, th_u, t_l, t_u, alpha, beta, alpha_pt0, ratio
real*8 dsa, sa_l, sa_u, alpha_sa, beta_sa, factor1, factor2

!EXPERIMENTATION WITH DIFFERENT ERROR TOLERANCES DID NOT LEAD TO SYSTEMATIC IMPROVEMENT
!key0="tol_liq_dens"
!eps0 = get_it_ctrl_density(key0)
!call set_it_ctrl_density(key0,-1d-12)
 
sea_cabb_pt0_si = errorreturn

pr0 = 101325d0; 

sea_cabb_pt0_si = errorreturn
th = sea_pottemp_si(sa_si,t_si,p_si,pr0);
if (th == errorreturn) return

!TH DERIVATIVE
dth = 1d-3
th_l = th - dth; th_u = th + dth; 

t_l = sea_pottemp_si(sa_si,th_l,pr0,p_si)
if (t_l == errorreturn) return
t_u = sea_pottemp_si(sa_si,th_u,pr0,p_si)
if (t_u == errorreturn .or. t_u == t_l) return

alpha = sea_alpha_pt0_si(sa_si,t_si,p_si) 
if(alpha == errorreturn) return
beta = sea_beta_pt0_si(sa_si,t_si,p_si)
if(beta == errorreturn .or. beta == 0d0) return
ratio = alpha/beta

factor1 = sea_alpha_pt0_si(sa_si,t_u,p_si)
if(factor1 == errorreturn) return
factor2 = sea_alpha_pt0_si(sa_si,t_l,p_si)
if(factor2 == errorreturn) return
alpha_pt0 = (factor1-factor2)/(th_u-th_l)

!SA DERIVATIVES
dsa = min(sa_si,1d-5)
if(sa_si >= dsa) then
  sa_l = sa_si - dsa; sa_u = sa_si + dsa; 
elseif(sa_si >= 0d0)then
  sa_l = 0d0; sa_u = dsa; 
else
    print *, '****    sa_si<0 in sea_cabb_pt0    ****'; stop
end if

if (sa_u == sa_l) return

factor1 = sea_alpha_pt0_si(sa_u,t_si,p_si)
if(factor1 == errorreturn) return
factor2 = sea_alpha_pt0_si(sa_l,t_si,p_si)
if(factor2 == errorreturn) return
alpha_sa = (factor1-factor2)/(sa_u-sa_l)

beta_sa =  (sea_beta_pt0_si(sa_u,t_si,p_si)-sea_beta_pt0_si(sa_l,t_si,p_si))/(sa_u-sa_l)

!CABBELING COEFFICIENT
sea_cabb_pt0_si = alpha_pt0 + ratio*(2d0*alpha_sa - ratio*beta_sa)

!call set_it_ctrl_density(key0, eps0)

end function


!==========================================================================
function sea_thrmb_ct_si(sa_si,t_si,p_si)
!==========================================================================

!THERMOBARIC COEFFICIENT OF SEAWATER WRT CONSERVATIVE TEMPERATURE 
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : SEA (GAUGE) PRESSURE               [PA]
!
!RESULT                 : THERMOBARIC COEFFICIENT            [1/(K PA)]
!                         WRT CONSERVATIVE TEMPERATURE

!CHECK VALUE:
!SEA_THRMB_CT_SI(0.035, 300, 1E5) = 1.4810927|1668362E-12  !DEFAULT TOLERANCES
!SEA_THRMB_CT_SI(0.035, 300, 1E5) = 1.4810927|5172403E-12  !TOLERANCE RESET TO -1D-12


implicit none
!character*13 key0
!real*8 eps0
real*8 sa_si, t_si, p_si, sea_thrmb_ct_si
real*8 pr0, theta, dp, p_l, p_u, t_l, t_u
real*8 beta, alpha_u, beta_u, alpha_l, beta_l

!EXPERIMENTATION WITH DIFFERENT ERROR TOLERANCES DID NOT LEAD TO SYSTEMATIC IMPROVEMENT
!key0="tol_liq_dens"
!eps0 = get_it_ctrl_density(key0)
!call set_it_ctrl_density(key0, -1d-12)

sea_thrmb_ct_si = errorreturn

dp = 1d3
pr0 = 101325d0; 
theta = sea_pottemp_si(sa_si,t_si,p_si,pr0)
if (theta == errorreturn) return

if(p_si >= dp) then
  p_l = p_si - dp; p_u = p_si + dp;
else
  p_l = 0d0; p_u = dp;
end if
if (p_u == p_l) return

t_l = sea_pottemp_si(sa_si,theta,pr0,p_l);    ! CT = CONST ==> THETA = CONST
if (t_l == errorreturn) return

t_u = sea_pottemp_si(sa_si,theta,pr0,p_u)
if (t_u == errorreturn) return

beta = sea_beta_ct_si(sa_si,t_si,p_si)
if (beta == errorreturn) return

alpha_u = sea_alpha_ct_si(sa_si,t_u,p_u)
if (alpha_u == errorreturn) return
beta_u = sea_beta_ct_si(sa_si,t_u,p_u)
if(beta_u == errorreturn .or. beta_u == 0d0) return

alpha_l = sea_alpha_ct_si(sa_si,t_l,p_l)
if (alpha_l == errorreturn) return
beta_l = sea_beta_ct_si(sa_si,t_l,p_l)
if(beta_l == errorreturn .or. beta_l == 0d0) return

sea_thrmb_ct_si = beta * (alpha_u/beta_u - alpha_l/beta_l)/(p_u - p_l)

!call set_it_ctrl_density(key0, eps0)

end function


!==========================================================================
function sea_thrmb_pt0_si(sa_si,t_si,p_si)
!==========================================================================

!THERMOBARIC COEFFICIENT OF SEAWATER WRT POTENTIAL TEMPERATURE 
!
!SA_SI                  : ABSOLUTE SALINITY                  [KG/KG]
!T_SI                   : IN SITU TEMPERATURE                [K]
!P_SI                   : SEA (GAUGE) PRESSURE               [PA]
!
!RESULT                 : THERMOBARIC COEFFICIENT            [1/(K PA)]
!                         WRT POTENTIAL TEMPERATURE

!CHECK VALUE:
!SEA_THRMB_PT0_SI(0.035, 300, 1E5) = 1.4594101|0702991E-12   !DEFAULT TOLERANCES
!SEA_THRMB_PT0_SI(0.035, 300, 1E5) = 1.4594101|3482853E-12   !TOLERANCE RESET TO -1D-12


implicit none
!character*13 key0
!real*8 eps0
real*8 sa_si, t_si, p_si, sea_thrmb_pt0_si
real*8 pr0, theta, dp, p_l, p_u, t_l, t_u
real*8 beta, alpha_u, beta_u, alpha_l, beta_l

!EXPERIMENTATION WITH DIFFERENT ERROR TOLERANCES DID NOT LEAD TO SYSTEMATIC IMPROVEMENT
!key0="tol_liq_dens"
!eps0 = get_it_ctrl_density(key0)
!call set_it_ctrl_density(key0, -1d-12)

sea_thrmb_pt0_si = errorreturn

dp = 1d3
pr0 = 101325d0

theta = sea_pottemp_si(sa_si,t_si,p_si,pr0)
if (theta == errorreturn) return

if(p_si >= dp) then
  p_l = p_si - dp; p_u = p_si + dp;
else
  p_l = 0d0; p_u = dp;
end if
if (p_u == p_l) return

t_l = sea_pottemp_si(sa_si,theta,pr0,p_l)
if (t_l == errorreturn) return

t_u = sea_pottemp_si(sa_si,theta,pr0,p_u)
if (t_u == errorreturn) return

beta = sea_beta_pt0_si(sa_si,t_si,p_si)
if (beta == errorreturn) return

alpha_u = sea_alpha_pt0_si(sa_si,t_u,p_u)
if (alpha_u == errorreturn) return
beta_u = sea_beta_pt0_si(sa_si,t_u,p_u)
if(beta_u == errorreturn) return
if(beta_u == errorreturn .or. beta_u == 0d0) return

alpha_l = sea_alpha_pt0_si(sa_si,t_l,p_l)
if (alpha_l == errorreturn) return
beta_l = sea_beta_pt0_si(sa_si,t_l,p_l)
if(beta_l == errorreturn .or. beta_l == 0d0) return

sea_thrmb_pt0_si = beta * (alpha_u/beta_u - alpha_l/beta_l)/(p_u - p_l)

!call set_it_ctrl_density(key0, eps0)

end function

end module sea_5a
