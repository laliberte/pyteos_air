module flu_if97_5

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS POTENTIALS OF LIQUID WATER AND VAPOUR AND ITS
!FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE AND
!DENSITY AS DEFINED IN IAPWS-IF97 FOR THE REGIONS 1 AND 2:

!REVISED RELEASE ON THE IAPWS INDUSTRIAL FORMULATION 1997
!FOR THE THERMODYNAMIC PROPERTIES OF WATER AND STEAM
!(THE REVISION ONLY RELATES TO THE EXTENSION OF REGION 5 TO 50 MPA)
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!LUCERNE , SWITZERLAND, AUGUST 2007

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

!THIS MODULE REQUIRES THE LIBRARY MODULE
!     CONSTANTS_0, FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: chk_iapws97_table, fit_liq_density_if97_si,  fit_liq_g_if97_si, &
          fit_vap_density_if97_si, fit_vap_g_if97_si

!COEFFICIENTS OF REGION 1
integer, private :: i1i(34), j1i(34)
!COEFFICIENTS OF REGION 2, IDEAL PART
integer, private :: j0i(9), iri(43), jri(43)
!COEFFICIENTS OF REGION 2, RESIDUAL PART
real*8,  private :: n0i(9), n1i(34), nri(43)

contains

!=========================================================================
function fit_liq_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 1 (LIQUID) IN IAPWS-IF97

!OUTPUT: 
!FIT_LIQ_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!DRV_P: ORDER OF THE PRESSURE DERIVATIVE (0-2)
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 fit_liq_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 1386d0
pu = 16530000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_SI = 461.51805

fit_liq_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_1

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0); g = gamma_1(0, 0, tau, psi)
         if(g == errorreturn) return
         fit_liq_g_if97_si = rt * g

      case (1); g_p = gamma_1(0, 1, tau, psi)
         if(g_p == errorreturn) return
         fit_liq_g_if97_si = rt * g_p / pu

      case (2); g_pp = gamma_1(0, 2, tau, psi)
         if(g_pp == errorreturn) return
         fit_liq_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0)
         g = gamma_1(0, 0, tau, psi)
         if(g == errorreturn) return
         g_t = gamma_1(1, 0, tau, psi)
         if(g_t == errorreturn) return
         fit_liq_g_if97_si = r * (g - tau * g_t)

      case (1) 
         g_p = gamma_1(0, 1, tau, psi)
         if(g_p == errorreturn) return
         g_tp = gamma_1(1, 1, tau, psi)
         if(g_tp == errorreturn) return
         fit_liq_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
         g_tt = gamma_1(2, 0, tau, psi)
         if(g_tt == errorreturn) return
         fit_liq_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function

!=========================================================================
function fit_vap_g_if97_si(drv_t, drv_p, t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE GIBBS FUNCTION G(T,P) AND ITS 1ST AND 2ND DERIVATIVES
!WITH RESPECT TO TEMPERATURE AND PRESSURE, AS DEFINED FOR THE REGION 2 (VAPOUR) IN IAPWS-IF97

!OUTPUT: 
!FIT_VAP_G_IF97_SI: SPECIFIC GIBBS ENERGY IN J/KG OR ITS T-P DERIVATIVE

!INPUT:  
!DRV_T: ORDER OF THE TEMPERATURE DERIVATIVE (0-2)
!ORDER OF THE PRESSURE DERIVATIVE (0-2)
!T_SI:  ABSOLUTE TEMPERATURE IN K
!P_SI:  ABSOLUTE PRESSURE IN PA

implicit none
                          
integer drv_t, drv_p
real*8 fit_vap_g_if97_si
real*8 t_si, p_si
real*8 g, g_t, g_tt
real*8 g_p, g_tp, g_pp
real*8 psi, tau, rt
real*8 tu, pu, r

tu = 540d0
pu = 1000000d0
r = 461.526d0       !J KG-1 K-1  SPECIFIC GAS CONSTANT
                    !NOTE THIS DEVIATES FROM GAS_CONSTANT_H2O_Si = 461.51805

fit_vap_g_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_if97_2

psi = p_si / pu
tau = tu / t_si
rt = r * t_si

select case (drv_t)

  case (0)
    select case (drv_p)

      case (0) 
         g = gamma_2(0, 0, tau, psi)
         if(g == errorreturn) return
         fit_vap_g_if97_si = rt * g

      case (1) 
         g_p = gamma_2(0, 1, tau, psi)
         if(g_p == errorreturn) return
         fit_vap_g_if97_si = rt * g_p / pu

      case (2) 
         g_pp = gamma_2(0, 2, tau, psi)
         if(g_pp == errorreturn) return
         fit_vap_g_if97_si = rt * g_pp / pu ** 2

      case default; return
    end select

  case (1)
    select case (drv_p)

      case (0) 
         g = gamma_2(0, 0, tau, psi)
         if(g == errorreturn) return
         g_t = gamma_2(1, 0, tau, psi)
         if(g_t == errorreturn) return
         fit_vap_g_if97_si = r * (g - tau * g_t)

      case (1) 
         g_p = gamma_2(0, 1, tau, psi)
         if(g_p == errorreturn) return
         g_tp = gamma_2(1, 1, tau, psi)
         if(g_tp == errorreturn) return
         fit_vap_g_if97_si = r * (g_p - tau * g_tp) / pu

      case default; return
    end select

  case (2)
    select case (drv_p)

      case (0) 
         g_tt = gamma_2(2, 0, tau, psi)
         if(g_tt == errorreturn) return
         fit_vap_g_if97_si = r * tau ** 2 * g_tt / t_si

      case default; return
    end select

  case default; return
end select

end function

!=========================================================================
function fit_liq_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 1 (LIQUID) IN IAPWS-IF97

!OUTPUT:  
!FIT_LIQ_DENSITY_IF97_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 fit_liq_density_if97_si
real*8 t_si, p_si
real*8 d

fit_liq_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = fit_liq_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

fit_liq_density_if97_si = 1d0 / d

end function

!=========================================================================
function fit_vap_density_if97_si(t_si, p_si)
!=========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, IN THE REGION 2 (VAPOUR) IN IAPWS-IF97

!OUTPUT:  
!FIT_VAP_DENSITY_IF97_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

implicit none

real*8 fit_vap_density_if97_si
real*8 t_si, p_si
real*8 d

fit_vap_density_if97_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = fit_vap_g_if97_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

fit_vap_density_if97_si = 1d0 / d

end function

!=========================================================================
function gamma_0(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_0 AS GIVEN IN TABLE 13 OF IF-97

implicit none

integer i, k, drv_t, drv_p
real*8 tau, psi
real*8 gamma_0
real*8 g, pwrt

gamma_0 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

g = 0d0
if(drv_t == 0) then
  if(drv_p == 0) then
    g = log(psi)
  else
    g = 1d0 / psi
    do k = 2, drv_p
      g = dfloat(1 - k) * g / psi
    enddo
  end if
end if

if(drv_p == 0) then
  do i = 1, 9
    pwrt = tau ** (j0i(i) - drv_t)
    do k = 0, drv_t - 1
      pwrt = pwrt * dfloat(j0i(i) - k)
    enddo
    g = g + n0i(i) * pwrt
  enddo
end if

gamma_0 = g

end function

!=========================================================================
function gamma_1(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN TABLE 4 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 tau, psi
real*8 gamma_1
real*8 g, pp, tt
real*8 pwrt, pwrp

gamma_1 = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

pp = 7.1d0 - psi
tt = tau - 1.222d0

g = 0d0
do i = 1, 34

  if(tt == 0d0) then
    if(j1i(i) == drv_t) then
       pwrt = 1d0
    elseif(j1i(i) > drv_t) then
       pwrt = 0d0
    else
       return
    endif
  else
    pwrt = tt ** (j1i(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(j1i(i) - k)
  enddo

  if(pp == 0d0) then
    if(i1i(i) == drv_p) then
       pwrp = 1d0
    elseif(i1i(i) > drv_p) then
       pwrp = 0d0
    else
       return
    endif
  else
    pwrp = pp ** (i1i(i) - drv_p)
  end if
  do k = 0, drv_p - 1
    pwrp = -pwrp * dfloat(i1i(i) - k)
  enddo

  g = g + n1i(i) * pwrp * pwrt

enddo

gamma_1 = g

end function

!=========================================================================
function gamma_2(drv_t, drv_p, tau, psi)
!=========================================================================

!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA AS GIVEN IN EQ. 15 OF IF-97

implicit none

integer drv_t, drv_p
real*8 tau, psi
real*8 gamma_2
real*8 g0, gr

call init_if97_2

gamma_2 = errorreturn

g0 = gamma_0(drv_t, drv_p, tau, psi)
if(g0 == errorreturn) return

gr = gamma_r(drv_t, drv_p, tau, psi)
if(gr == errorreturn) return

gamma_2 = g0 + gr

end function

!=========================================================================
function gamma_r(drv_t, drv_p, tau, psi)
!=========================================================================
                         
!THIS FUNCTION IMPLEMENTS THE DERIVATIVES OF GAMMA_R AS GIVEN IN TABLE 14 OF IF-97

implicit none

integer drv_t, drv_p
integer i, k
real*8 gamma_r
real*8 tau, psi
real*8 g, tt
real*8 pwrt, pwrp

gamma_r = errorreturn

if(psi <= 0d0) return
if(tau <= 0d0) return
if(drv_t < 0) return
if(drv_p < 0) return

tt = tau - 0.5d0

g = 0d0
do i = 1, 43

  if(tt == 0d0) then
    if(jri(i) == drv_t) then
       pwrt = 1d0
    elseif(jri(i) > drv_t) then
       pwrt = 0d0
    else
       return
    endif
  else
    pwrt = tt ** (jri(i) - drv_t)
  end if
  do k = 0, drv_t - 1
    pwrt = pwrt * dfloat(jri(i) - k)
  enddo

  pwrp = psi ** (iri(i) - drv_p)
  do k = 0, drv_p - 1
    pwrp = pwrp * dfloat(iri(i) - k)
  enddo

  g = g + nri(i) * pwrp * pwrt

enddo

gamma_r = g

end function

!=========================================================================
subroutine init_if97_1
!=========================================================================

implicit none
integer i

if(n1i(1) == 0.14632971213167d0) return

!TABLE-2. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE DIMENSIONLESS GIBBS FREE ENERGY
!FOR REGION 1, EQ. (7)
!I II JI NI I II JI NI
i = 1; i1i(i)  = 0;  j1i(i) = -2;  n1i(i) = 0.14632971213167d0
i = 2; i1i(i)  = 0;  j1i(i) = -1;  n1i(i) = -0.84548187169114d0
i = 3; i1i(i)  = 0;  j1i(i) = 0;   n1i(i) = -3.756360367204d0
i = 4; i1i(i)  = 0;  j1i(i) = 1;   n1i(i) = 3.3855169168385d0
i = 5; i1i(i)  = 0;  j1i(i) = 2;   n1i(i) = -0.95791963387872d0
i = 6; i1i(i)  = 0;  j1i(i) = 3;   n1i(i) = 0.15772038513228d0
i = 7; i1i(i)  = 0;  j1i(i) = 4;   n1i(i) = -0.016616417199501d0
i = 8; i1i(i)  = 0;  j1i(i) = 5;   n1i(i) = 8.1214629983568d-04
i = 9; i1i(i)  = 1;  j1i(i) = -9;  n1i(i) = 2.8319080123804d-04
i = 10; i1i(i) = 1;  j1i(i) = -7;  n1i(i) = -6.0706301565874d-04
i = 11; i1i(i) = 1;  j1i(i) = -1;  n1i(i) = -0.018990068218419d0
i = 12; i1i(i) = 1;  j1i(i) = 0;   n1i(i) = -0.032529748770505d0
i = 13; i1i(i) = 1;  j1i(i) = 1;   n1i(i) = -0.021841717175414d0
i = 14; i1i(i) = 1;  j1i(i) = 3;   n1i(i) = -5.283835796993d-05  
i = 15; i1i(i) = 2;  j1i(i) = -3;  n1i(i) = -4.7184321073267d-04
i = 16; i1i(i) = 2;  j1i(i) = 0;   n1i(i) = -3.0001780793026d-04
i = 17; i1i(i) = 2;  j1i(i) = 1;   n1i(i) = 4.7661393906987d-05
i = 18; i1i(i) = 2;  j1i(i) = 3;   n1i(i) = -4.4141845330846d-06
i = 19; i1i(i) = 2;  j1i(i) = 17;  n1i(i) = -7.2694996297594d-16
i = 20; i1i(i) = 3;  j1i(i) = -4;  n1i(i) = -3.1679644845054d-05
i = 21; i1i(i) = 3;  j1i(i) = 0;   n1i(i) = -2.8270797985312d-06
i = 22; i1i(i) = 3;  j1i(i) = 6;   n1i(i) = -8.5205128120103d-10
i = 23; i1i(i) = 4;  j1i(i) = -5;  n1i(i) = -2.2425281908d-06
i = 24; i1i(i) = 4;  j1i(i) = -2;  n1i(i) = -6.5171222895601d-07
i = 25; i1i(i) = 4;  j1i(i) = 10;  n1i(i) = -1.4341729937924d-13
i = 26; i1i(i) = 5;  j1i(i) = -8;  n1i(i) = -4.0516996860117d-07
i = 27; i1i(i) = 8;  j1i(i) = -11; n1i(i) = -1.2734301741641d-09
i = 28; i1i(i) = 8;  j1i(i) = -6;  n1i(i) = -1.7424871230634d-10
i = 29; i1i(i) = 21; j1i(i) = -29; n1i(i) = -6.8762131295531d-19
i = 30; i1i(i) = 23; j1i(i) = -31; n1i(i) = 1.4478307828521d-20
i = 31; i1i(i) = 29; j1i(i) = -38; n1i(i) = 2.6335781662795d-23
i = 32; i1i(i) = 30; j1i(i) = -39; n1i(i) = -1.1947622640071d-23
i = 33; i1i(i) = 31; j1i(i) = -40; n1i(i) = 1.8228094581404d-24
i = 34; i1i(i) = 32; j1i(i) = -41; n1i(i) = -9.3537087292458d-26

end subroutine

!=========================================================================
subroutine init_if97_2
!=========================================================================

implicit none
integer i

if(n0i(1) == -9.6927686500217d0) return

!TABLE 10. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE IDEAL-GAS PART GAMMA_0 OF THE
!DIMENSIONLESS GIBBS FREE ENERGY FOR REGION 2, EQ. (16)
!I JI NI I JI NI
i = 1; j0i(i) = 0;  n0i(i) = -9.6927686500217d0
i = 2; j0i(i) = 1;  n0i(i) = 10.086655968018d0
i = 3; j0i(i) = -5; n0i(i) = -0.005608791128302d0
i = 4; j0i(i) = -4; n0i(i) = 0.071452738081455d0
i = 5; j0i(i) = -3; n0i(i) = -0.40710498223928d0
i = 6; j0i(i) = -2; n0i(i) = 1.4240819171444d0
i = 7; j0i(i) = -1; n0i(i) = -4.383951131945d0
i = 8; j0i(i) = 2;  n0i(i) = -0.28408632460772d0
i = 9; j0i(i) = 3;  n0i(i) = 0.021268463753307d0
!IF(EQ. (16) IS INCORPORATED INTO EQ. (18), INSTEAD OF THE-VALUES FOR N0I(1) AND N0I(1)
!GIVEN ABOVE, THE FOLLOWING VALUES
!N0I(1) = -0.96937268393049D1
!N0I(2) = 0.10087275970006D2
!SHOULD BE USED


!TABLE 11. NUMERICAL VALUES OF THE COEFFICIENTS AND EXPONENTS OF THE
!RESIDUAL PART GAMMA_R OF THE DIMENSIONLESS GIBBS FREE ENERGY FOR
!REGION 2, EQ. (17)
!I II JI NI
i = 1;  iri(i) = 1;  jri(i) = 0;  nri(i) = -1.7731742473213d-03
i = 2;  iri(i) = 1;  jri(i) = 1;  nri(i) = -0.017834862292358d0
i = 3;  iri(i) = 1;  jri(i) = 2;  nri(i) = -0.045996013696365d0
i = 4;  iri(i) = 1;  jri(i) = 3;  nri(i) = -0.057581259083432d0
i = 5;  iri(i) = 1;  jri(i) = 6;  nri(i) = -0.05032527872793d0
i = 6;  iri(i) = 2;  jri(i) = 1;  nri(i) = -3.3032641670203d-05
i = 7;  iri(i) = 2;  jri(i) = 2;  nri(i) = -1.8948987516315d-04
i = 8;  iri(i) = 2;  jri(i) = 4;  nri(i) = -3.9392777243355d-03
i = 9;  iri(i) = 2;  jri(i) = 7;  nri(i) = -0.043797295650573d0
i = 10; iri(i) = 2;  jri(i) = 36; nri(i) = -2.6674547914087d-05
i = 11; iri(i) = 3;  jri(i) = 0;  nri(i) = 2.0481737692309d-08
i = 12; iri(i) = 3;  jri(i) = 1;  nri(i) = 4.3870667284435d-07
i = 13; iri(i) = 3;  jri(i) = 3;  nri(i) = -3.227767723857d-05
i = 14; iri(i) = 3;  jri(i) = 6;  nri(i) = -1.5033924542148d-03
i = 15; iri(i) = 3;  jri(i) = 35; nri(i) = -0.040668253562649d0
i = 16; iri(i) = 4;  jri(i) = 1;  nri(i) = -7.8847309559367d-10
i = 17; iri(i) = 4;  jri(i) = 2;  nri(i) = 1.2790717852285d-08
i = 18; iri(i) = 4;  jri(i) = 3;  nri(i) = 4.8225372718507d-07
i = 19; iri(i) = 5;  jri(i) = 7;  nri(i) = 2.2922076337661d-06
i = 20; iri(i) = 6;  jri(i) = 3;  nri(i) = -1.6714766451061d-11
i = 21; iri(i) = 6;  jri(i) = 16; nri(i) = -2.1171472321355d-03
i = 22; iri(i) = 6;  jri(i) = 35; nri(i) = -23.895741934104d0
i = 23; iri(i) = 7;  jri(i) = 0;  nri(i) = -5.905956432427d-18
i = 24; iri(i) = 7;  jri(i) = 11; nri(i) = -1.2621808899101d-06
i = 25; iri(i) = 7;  jri(i) = 25; nri(i) = -0.038946842435739d0
i = 26; iri(i) = 8;  jri(i) = 8;  nri(i) = 1.1256211360459d-11
i = 27; iri(i) = 8;  jri(i) = 36; nri(i) = -8.2311340897998d0
i = 28; iri(i) = 9;  jri(i) = 13; nri(i) = 1.9809712802088d-08
i = 29; iri(i) = 10; jri(i) = 4;  nri(i) = 1.0406965210174d-19
i = 30; iri(i) = 10; jri(i) = 10; nri(i) = -1.0234747095929d-13
i = 31; iri(i) = 10; jri(i) = 14; nri(i) = -1.0018179379511d-09
i = 32; iri(i) = 16; jri(i) = 29; nri(i) = -8.0882908646985d-11
i = 33; iri(i) = 16; jri(i) = 50; nri(i) = 0.10693031879409d0
i = 34; iri(i) = 18; jri(i) = 57; nri(i) = -0.33662250574171d0
i = 35; iri(i) = 20; jri(i) = 20; nri(i) = 8.9185845355421d-25
i = 36; iri(i) = 20; jri(i) = 35; nri(i) = 3.0629316876232d-13
i = 37; iri(i) = 20; jri(i) = 48; nri(i) = -4.2002467698208d-06
i = 38; iri(i) = 21; jri(i) = 21; nri(i) = -5.9056029685639d-26
i = 39; iri(i) = 22; jri(i) = 53; nri(i) = 3.7826947613457d-06
i = 40; iri(i) = 23; jri(i) = 39; nri(i) = -1.2768608934681d-15
i = 41; iri(i) = 24; jri(i) = 26; nri(i) = 7.3087610595061d-29
i = 42; iri(i) = 24; jri(i) = 40; nri(i) = 5.5414715350778d-17
i = 43; iri(i) = 24; jri(i) = 58; nri(i) = -9.436970724121d-07

end subroutine

!=========================================================================
subroutine chk_iapws97_table(number)
!=========================================================================

!THIS FUNCTION RETURNS SELECTED CHECK TABLES OF IAPWS-IF97
!IN COMPARISON TO VALUES COMPUTED FROM THIS MODULE
!AS A FORMATTED STRING, TO BE PRINTED INTO A WINDOW OR FILE,
!FOR EASY VERIFICATION OF THIS IMPLEMENTATION ON A GIVEN PLATFORM
!
!IAPWS-IF97:
!REVISED RELEASE ON THE IAPWS INDUSTRIAL FORMULATION 1997
!FOR THE THERMODYNAMIC PROPERTIES OF WATER AND STEAM
!(THE REVISION ONLY RELATES TO THE EXTENSION OF REGION 5 TO 50 MPA)
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!LUCERNE, SWITZERLAND, AUGUST 2007
!AVAILABLE AT: HTTP://WWW.IAPWS.ORG

implicit none
integer number

if(number == 5) then
   call chk_iapws97_table5
elseif(number == 15) then
   call chk_iapws97_table15
else
   write(*,*)"ERROR: 5 AND 15 ARE THE ONLY LEGAL TABLE NUMBERS FOR IF15"
endif

return
end subroutine

!=========================================================================
subroutine chk_iapws97_table5
!=========================================================================

!IAPWS-IF97 CHECK TABLE:
!REVISED RELEASE ON THE IAPWS INDUSTRIAL FORMULATION 1997
!FOR THE THERMODYNAMIC PROPERTIES OF WATER AND STEAM
!(THE REVISION ONLY RELATES TO THE EXTENSION OF REGION 5 TO 50 MPA)
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!LUCERNE , SWITZERLAND, AUGUST 2007
!AVAILABLE AT: HTTP://WWW.IAPWS.ORG

!TABLE 5. THERMODYNAMIC PROPERTY VALUES CALCULATED FROM EQ. (7) FOR SELECTED VALUES OF T AND P
!                  T = 300 K,              T = 300 K,          T = 500 K,
!                  P = 3 MPA               P = 80 MPA          P = 3 MPA
!V  / (M3/KG)      0.100 215 168 E-2       0.971 180 894 E-3   0.120 241 800 E-2
!H  / (KJ/KG)      0.115 331 273 E3        0.184 142 828 E3    0.975 542 239 E3
!U  / (KJ/KG)      0.112 324 818 E3        0.106 448 356 E3    0.971 934 985 E3
!S  / (KJ/KG/K)    0.392 294 792           0.368 563 852       0.258 041 912 E1
!CP / (KJ/KG/K)    0.417 301 218 E1        0.401 008 987 E1    0.465 580 682 E1
!W  / (M/S)        0.150 773 921 E4        0.163 469 054 E4    0.124 071 337 E4

implicit none
integer i, j
real*8 p, t, d(3)
real*8 g(3), g_t(3), g_tt(3), g_tp(3), g_p(3), g_pp(3)
real*8 chooset(3), choosep(3)

                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "IMPLEMENTATION OF IAPWS-IF97 IN FORTRAN"
write(*,*) "FOR PUBLICATION IN OCEAN SCIENCE, 2009"
write(*,*) "D. WRIGHT, BIO, VERSION ", version
write(*,*) " "

write(*,*) "FUNCTION VALUES AS GIVEN IN TABLE 5 OF IF-97:"
write(*,*) "TABLE 5. THERMODYNAMIC PROPERTY VALUES CALCULATED FROM EQ. (7)"
write(*,*) "FOR SELECTED VALUES OF T AND P"
write(*,*) " "

write(*,*) "               T = 300 K           T = 300 K           T = 500 K"
write(*,*) "               P = 3 MPA           P = 80 MPA          P = 3 MPA"

chooset(1) = 300d0
chooset(2) = 300d0
chooset(3) = 500d0
choosep(1) = 3d0
choosep(2) = 80d0
choosep(3) = 3d0

do i = 1, 6
  write(*,*)' '
  select case (i)
    case (1);  write(*,*) "v/(M3/KG)      0.100215168e-2      0.971180894e-3      0.120241800e-2"
    case (2);  write(*,*) "h/(KJ/KG)      0.115331273e3       0.184142828e3       0.975542239e3"
    case (3);  write(*,*) "u/(KJ/KG)      0.112324818e3       0.106448356e3       0.971934985e3"
    case (4);  write(*,*) "s/(KJ/KG/K)    0.392294792         0.368563852         0.258041912e1"
    case (5);  write(*,*) "cp/(KJ/KG/K)   0.417301218e1       0.401008987e1       0.465580682e1"
    case (6);  write(*,*) "w/(M/S)        0.150773921e4       0.163469054e4       0.124071337e4"
  end select
  
  do j = 1, 3
    !GET TEMPERATURE AND PRESSURE INPUT
    t = chooset(j)
    p = choosep(j) * 1d6

    g(j)    = fit_liq_g_if97_si(0, 0, t, p)
    g_t(j)  = fit_liq_g_if97_si(1, 0, t, p)
    g_p(j)  = fit_liq_g_if97_si(0, 1, t, p)
    g_tt(j) = fit_liq_g_if97_si(2, 0, t, p)
    g_tp(j) = fit_liq_g_if97_si(1, 1, t, p)
    g_pp(j) = fit_liq_g_if97_si(0, 2, t, p)

  enddo

  select case (i)
    case (1);  do j = 1,3; d(j) = g_p(j); enddo
    case (2);  do j = 1,3; d(j) = 0.001d0 * (g(j) - chooset(j) * g_t(j)); enddo
    case (3);  do j = 1,3; d(j) = 0.001d0 * (g(j) - chooset(j) * g_t(j) - 1d6 * choosep(j) * g_p(j)); enddo
    case (4);  do j = 1,3; d(j) = -0.001d0 * g_t(j); enddo
    case (5);  do j = 1,3; d(j) = -0.001d0 * chooset(j) * g_tt(j); enddo
    case (6)  
      do j = 1,3
        d(j) = g_tp(j) * g_tp(j) - g_tt(j) * g_pp(j)
        if(d (j) /= 0d0) then
          d(j) = g_tt(j) / d(j)
          if(d(j) >= 0d0) then
            d(j) = g_p(j) * sqrt(d(j))
          else
            d(j) = errorreturn
          end if
        else
          d(j) = errorreturn
        end if
      enddo
  end select

  write(*,"('COMPUTED NOW:'2x,e16.9,x,e19.9,x,e19.9)") d(1),d(2),d(3)

enddo

end subroutine

!=========================================================================
subroutine chk_iapws97_table15
!=========================================================================

!IAPWS-IF97 CHECK TABLE:
!REVISED RELEASE ON THE IAPWS INDUSTRIAL FORMULATION 1997
!FOR THE THERMODYNAMIC PROPERTIES OF WATER AND STEAM
!(THE REVISION ONLY RELATES TO THE EXTENSION OF REGION 5 TO 50 MPA)
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!LUCERNE , SWITZERLAND, AUGUST 2007
!AVAILABLE AT: HTTP://WWW.IAPWS.ORG

!TABLE 15. THERMODYNAMIC PROPERTY VALUES CALCULATED FROM EQ. (15) FOR SELECTED VALUES OF T AND P A
!                      T = 300 K,       T = 700 K,       T = 700 K
!                      P = 0.0035 MPA   P = 0.0035 MPA   P = 30 MPA
!V  / (M3/KG)          0.394 913 866 E2 0.923 015 898 E2 0.542 946 619 E-2
!H  / (KJ/KG)          0.254 991 145 E4 0.333 568 375 E4 0.263 149 474 E4
!U  / (KJ/KG)          0.241 169 160 E4 0.301 262 819 E4 0.246 861 076 E4
!S  / (KJ/KG/K)        0.852 238 967 E1 0.101 749 996 E2 0.517 540 298 E1
!CP / (KJ/KG/K)        0.191 300 162 E1 0.208 141 274 E1 0.103 505 092 E2
!W  / (M/S)            0.427 920 172 E3 0.644 289 068 E3 0.480 386 523 E3

implicit none
integer i, j
real*8 t, p, d(3)
real*8 g(3), g_t(3), g_tt(3), g_tp(3), g_p(3), g_pp(3)
real*8 chooset(3), choosep(3)


                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "IMPLEMENTATION OF IAPWS-IF97 IN FORTRAN"
write(*,*) "FOR PUBLICATION IN OCEAN SCIENCE, 2009"
write(*,*) "D. WRIGHT, BIO, VERSION ", version
write(*,*) " "

write(*,*) "FUNCTION VALUES AS GIVEN IN TABLE 15 OF IF-97:"
write(*,*) "TABLE 15. THERMODYNAMIC PROPERTY VALUES CALCULATED FROM EQ. (15)"
write(*,*) "FOR SELECTED VALUES OF T AND P"
write(*,*) " "

write(*,*) "               T = 300 K           T = 700 K           T = 700 K"
write(*,*) "               P = 3500 PA         P = 3500 PA         P = 30 MPA"

chooset(1) = 300d0
chooset(2) = 700d0
chooset(3) = 700d0
choosep(1) = 3500d0
choosep(2) = 3500d0
choosep(3) = 30000000d0

do i = 1, 6

  write(*,*) " "
  select case (i)
    case (1);  write(*,*) "v/(M3/KG)      0.394913866e2       0.923015898e2       0.542946619e-2"
    case (2);  write(*,*) "h/(KJ/KG)      0.254991145e4       0.333568375e4       0.263149474e4"
    case (3);  write(*,*) "u/(KJ/KG)      0.241169160e4       0.301262819e4       0.246861076e4"
    case (4);  write(*,*) "s/(KJ/KG/K)    0.852238967e1       0.101749996e2       0.517540298e1"
    case (5);  write(*,*) "cp/(KJ/KG/K)   0.191300162e1       0.208141274e1       0.103505092e2"
    case (6);  write(*,*) "w/(M/S)        0.427920172e3       0.644289068e3       0.480386523e3"
  end select
  
  do j = 1, 3
    t = chooset(j)
    p = choosep(j)

    g(j)    = fit_vap_g_if97_si(0, 0, t, p)
    g_t(j)  = fit_vap_g_if97_si(1, 0, t, p)
    g_p(j)  = fit_vap_g_if97_si(0, 1, t, p)
    g_tt(j) = fit_vap_g_if97_si(2, 0, t, p)
    g_tp(j) = fit_vap_g_if97_si(1, 1, t, p)
    g_pp(j) = fit_vap_g_if97_si(0, 2, t, p)
  enddo

  select case (i)
    case (1); do j=1,3; d(j) = g_p(j); enddo
    case (2); do j=1,3; d(j) = 0.001d0 * (g(j) - chooset(j) * g_t(j)); enddo
    case (3); do j=1,3; d(j) = 0.001d0 * (g(j) - chooset(j) * g_t(j) - choosep(j) * g_p(j)); enddo
    case (4); do j=1,3; d(j) = -0.001d0 * g_t(j); enddo
    case (5); do j=1,3; d(j) = -0.001d0 * chooset(j) * g_tt(j); enddo
    case (6) 
              do j=1,3 
                 d(j) = g_tp(j) * g_tp(j) - g_tt(j) * g_pp(j)
                 if(d(j) /= 0d0) then
                   d(j) = g_tt(j) / d(j)
                   if(d(j) >= 0d0) then
                     d(j) = g_p(j) * sqrt(d(j))
                   else
                     d(j) = errorreturn
                   end if
                 else
                   d(j) = errorreturn
                 end if
              enddo
  end select

  write(*,"('COMPUTED NOW:'2x,e16.9,x,e19.9,x,e19.9)") d(1),d(2),d(3)

enddo

end subroutine

end module flu_if97_5
