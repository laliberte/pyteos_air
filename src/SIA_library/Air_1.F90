module air_1

!#########################################################################

!THIS MODULE IMPLEMENTS THE HELMHOLTZ POTENTIAL OF DRY AIR AND ITS
!FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE AND
!DENSITY AS DEFINED IN:

!LEMMON, E.W., JACOBSEN, R.T., PENONCELLO, S.G., FRIEND, D.G.:
!THERMODYNAMIC PROPERTIES OF AIR AND MIXTURES OF NITROGEN, ARGON AND OXYGEN
!FROM 60 TO 2000 K AT PRESSURES TO 2000 MPA.
!J. PHYS. CHEM. REF. DATA 29(2000)331-362.

!THIS MODULE ALSO IMPLEMENTS THE CROSS-VIRIAL COEFFICIENTS OF HUMID AIR
!AND THEIR FIRST AND SECOND DERIVATIVES WITH RESPECT TO TEMPERATURE.

!THE 2ND MIXED VIRIAL COEFFICIENT IS DEFINED IN:

!HARVEY, A.H., HUANG, P.H.:  
!FIRST-PRINCIPLES CALCULATION OF THE AIR-WATER SECOND VIRIAL COEFFICIENT. 
!INT. J. THERMOPHYS. 28(2007),556-565.

!THE 3RD MIXED VIRIAL COEFFICIENTS ARE DEFINED IN:

!HYLAND, R.W., WEXLER, A.:
!FORMULATIONS FOR THE THERMODYNAMIC PROPERTIES OF DRY AIR FROM 173.15 K TO 473.15 K,
!AND OF SATURATED MOIST AIR FROM 173.15 K TO 372.15 K, AT PRESSURES UP TO 5 MPA.
!ASHRAE TRANSACT. 89(1983),520-535.

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., !WRIGHT, D.G., HERRMANN, S.:
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

public :: air_baw_m3mol, air_caaw_m6mol2, air_caww_m6mol2, &
          dry_f_si, dry_init_clear, dry_init_Lemmon2000

real*8 ::  mw = molar_mass_h2o_si      !MOLAR MASS OF H2O IN KG/MOL
real*8 ::  ma = molar_mass_air_L2000   !MOLAR MASS OF AIR IN KG/MOL USED BY LEMMON et al. 2000

!ARRAYS DEFINED IN TABLES 12 AND 13 OF LEMON ET AL. (2000)
!COEFFICIENTS AND EXPONENTS FOR THE FUNDAMENTAL EQUATION (3) FOR AIR

real*8 :: nk0(13), nk(19), ik(19), jk(19), lk(19)

contains

!==========================================================================
function dry_f_si(drv_t, drv_d, t_si, d_si)
!==========================================================================
         
!THIS FUNCTION IMPLEMENTS THE HELMHOLTZ POTENTIAL OF DRY AIR AND
!ITS 1ST AND 2ND DERIVATIVES WITH RESPECT TO TEMPERATURE AND DENSITY
!AS PUBLISHED IN:

!LEMMON, E.W., JACOBSEN, R.T., PENONCELLO, S.G., FRIEND, D.G.
!THERMODYNAMIC PROPERTIES OF AIR AND MIXTURES OF NITROGEN, ARGON AND OXYGEN
!FROM 60 TO 2000 K AT PRESSURES TO 2000 MPA.
!J. PHYS. CHEM. REF. DATA 29(2000)331-362.

!HERE, IN CONTRAST TO THE ORIGINAL ARTICLE, THE COEFFICIENTS OF ABSOLUTE
!ENERGY AND ABSOLUTE ENTROPY ARE ADJUSTED TO THE REFERENCE STATE CONDITION
!OF ZERO ENTROPY AND ZERO ENTHALPY AT THE STANDARD OCEAN STATE, 0C AND 101325 PA.

!OUTPUT: 
!dry_f_si:  HELMHOLTZ ENERGY OR ITS DERIVATIVES IN J/(KG K^DRV_T) (M3/KG)^DRV_D)


!INPUT PARAMETERS:
!DRV_T:     ORDER OF THE TEMPERATURE DERIVATIVE, DRV_T < 3
!DRV_D:     ORDER OF THE DENSITY DERIVATIVE, DRV_D + DRV_T < 3
!T_SI:      ABSOLUTE ITS-90 TEMPERATURE IN K
!D_SI:      DENSITY IN KG/M3

!CHECK VALUES, V1.0:
!dry_f_si( 0, 0, 300,0.001) =-696384.495286
!dry_f_si( 1, 0, 300,0.001) =-2124.98673214
!dry_f_si( 0, 1, 300,0.001) = 86135114.6608
!dry_f_si( 2, 0, 300,0.001) =-2.39299064929
!dry_f_si( 1, 1, 300,0.001) = 287117.623822
!dry_f_si( 0, 2, 300,0.001) =-86135137748.2

!CHECK VALUES, V1.1:
!dry_f_si( 0, 0, 300,0.001) =-696239.965190
!dry_f_si( 1, 0, 300,0.001) =-2124.55145456
!dry_f_si( 0, 1, 300,0.001) = 86114714.9596
!dry_f_si( 2, 0, 300,0.001) =-2.39242390806
!dry_f_si( 1, 1, 300,0.001) = 287049.624545
!dry_f_si( 0, 2, 300,0.001) =-86114738036.1

integer :: drv_t, drv_d
real*8 :: dry_f_si
real*8 :: t_si, d_si
real*8 :: tr, dr, ra, f, tau, del, a, a_t, a_tt, a_d, a_td, a_dd

ra = gas_constant_air_L2000          !SPECIFIC GAS CONSTANT OF AIR IN J KG-1 K-1 USED BY LEMON ET AL. 2000

tr = 132.6312d0                      !REDUCING TEMPERATURE IN K
dr = 10447.7d0 * ma                  !REDUCING DENSITY IN KG M-3

dry_f_si = errorreturn
f = errorreturn

call init_iapws10

if(check_limits == 1) then
  if(t_si < dry_air_tmin .or. t_si > dry_air_tmax) return
  if(d_si <= dry_air_dmin .or. d_si > dry_air_dmax) return
else
  if (t_si <= 0d0) return
  if (d_si <= 0d0) return
endif

tau = tr / t_si  !REDUCED TEMPERATURE
del = d_si / dr  !REDUCED DENSITY

!HELMHOLTZ FUNCTION AND ITS DERIVATIVES
select case (drv_t)
  case (0)
    select case (drv_d)
      case (0)
        a = log(del) + alpha_ideal(tau) + alpha_res(0, 0, tau, del)
        f = ra * t_si * a
      case (1)
        a_d = 1d0 / del + alpha_res(0, 1, tau, del)
        f = ra * t_si * a_d / dr
      case (2)
        a_dd = -1d0 / del ** 2 + alpha_res(0, 2, tau, del)
        f = ra * t_si * a_dd / dr ** 2
      case default; return  !dbg
    end select
  case (1)
    select case (drv_d)
      case (0)
        a = log(del) + alpha_ideal(tau) + alpha_res(0, 0, tau, del)
        a_t = alpha_ideal_t(tau) + alpha_res(1, 0, tau, del)
        f = ra * (a - tau * a_t)
      case (1)
        a_d = 1d0 / del + alpha_res(0, 1, tau, del)
        a_td = alpha_res(1, 1, tau, del)
        f = ra * (a_d - tau * a_td) / dr
      case default; return  !dbg
    end select
  case (2)
    select case (drv_d)
      case (0)
        a_tt = alpha_ideal_tt(tau) + alpha_res(2, 0, tau, del)
        f = ra * tau ** 2 * a_tt / t_si
      case default; return  !dbg
    end select
end select

dry_f_si = f

end function

!==========================================================================
function air_baw_m3mol(drv_t, t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE 2ND VIRIAL COEFFICIENT BAW AS DEFINED IN

!HARVEY, A.H., HUANG, P.H.: FIRST-PRINCIPLES CALCULATION OF THE AIR-WATER SECOND
!VIRIAL COEFFICIENT. INT. J. THERMOPHYS. 28(2007)556-565.

!OUTPUT:  
!air_baw_m3mol = MOLAR VIRIAL COEFFICIENT BAW OR ITS DERIVATIVES IN M3/(MOL K^DRV_T)

!INPUT:   
!DRV_T ORDER OF THE TEMPERATURE DERIVATIVE
!T_SI  ABSOLUTE TEMPERATURE IN K

real*8 :: air_baw_m3mol
integer :: drv_t
real*8 :: t_si
real*8 :: c1, c2, c3, tr, d1, d2, d3
real*8 :: t, b

c1 = 66.5687d0; d1 = -0.237d0
c2 = -238.834d0; d2 = -1.048d0
c3 = -176.755d0; d3 = -3.183d0
tr = 100d0 !K

air_baw_m3mol = errorreturn

if (t_si <= 0d0) return

t = t_si / tr

select case (drv_t)
  case (0); b = c1 * t ** d1 + c2 * t ** d2 + c3 * t ** d3
  case (1)
    b = c1 * d1 * t ** (d1 - 1d0) + c2 * d2 * t ** (d2 - 1d0) &
        + c3 * d3 * t ** (d3 - 1d0)
    b = b / tr
  case (2)
    b = c1 * d1 * (d1 - 1d0) * t ** (d1 - 2d0) + &
        c2 * d2 * (d2 - 1d0) * t ** (d2 - 2d0) + &
        c3 * d3 * (d3 - 1d0) * t ** (d3 - 2d0)
    b = b / tr ** 2
  case default; return
end select

air_baw_m3mol = b * 0.000001d0 !FROM CM3/MOL TO M3/MOL

end function

!==========================================================================
function air_caaw_m6mol2(drv_t, t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE 3RD VIRIAL COEFFICIENT CAAW AS DEFINED IN

!HYLAND, R.W., WEXLER, A.:
!FORMULATIONS FOR THE THERMODYNAMIC PROPERTIES OF DRY AIR FROM 173.15 K TO 473.15 K,
!AND OF SATURATED MOIST AIR FROM 173.15 K TO 372.15 K, AT PRESSURES UP TO 5 MPA.
!ASHRAE TRANSACT. 89(1983)520-535.

!OUTPUT:  
!AIR_CAAW_M6MOL2 = MOLAR VIRIAL COEFFICIENT CAAW OR ITS DERIVATIVES IN M6/(MOL2 K^DRV_T)

!INPUT:   
!DRV_T ORDER OF THE TEMPERATURE DERIVATIVE
!T_SI  ABSOLUTE TEMPERATURE IN K

real*8 :: air_caaw_m6mol2
integer :: drv_t
real*8 :: t_si
real*8 :: c1, c2, c3, c4, c5, tr
real*8 :: t, c

c1 = 482.737d0
c2 = 105678d0
c3 = -65639400d0
c4 = 29444200000d0
c5 = -3193170000000d0
tr = 1d0 !K

air_caaw_m6mol2 = errorreturn

if (t_si <= 0d0) return

t = t_si / tr

select case (drv_t)
  case (0); c = c1 + c2 / t + c3 / t ** 2 + c4 / t ** 3 + c5 / t ** 4
  case (1); c = -c2 / t ** 2 - 2d0 * c3 / t ** 3 - 3d0 * c4 / t ** 4 &
                - 4d0 * c5 / t ** 5
  case (2); c = 2d0 * c2 / t ** 3 + 6d0 * c3 / t ** 4 &
                + 12d0 * c4 / t ** 5 + 20d0 * c5 / t ** 6
  case default; return
end select

air_caaw_m6mol2 = c * 1d-12 !FROM CM6/MOL2 TO M6/MOL2

end function

!==========================================================================
function air_caww_m6mol2(drv_t, t_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE 3RD VIRIAL COEFFICIENT CAWW AS DEFINED IN

!HYLAND, R.W., WEXLER, A.:
!FORMULATIONS FOR THE THERMODYNAMIC PROPERTIES OF DRY AIR FROM 173.15 K TO 473.15 K,
!AND OF SATURATED MOIST AIR FROM 173.15 K TO 372.15 K, AT PRESSURES UP TO 5 MPA.
!ASHRAE TRANSACT. 89(1983)520-535.

!OUTPUT:  
!air_caww_m6mol2 = MOLAR VIRIAL COEFFICIENT CAWW OR ITS DERIVATIVES IN M6/(MOL2 K^DRV_T)

!INPUT:   
!DRV_T ORDER OF THE TEMPERATURE DERIVATIVE
!T_SI  ABSOLUTE TEMPERATURE IN K

real*8 air_caww_m6mol2
integer :: drv_t
real*8 :: t_si
real*8 c1, c2, c3, c4, tr
real*8 t, c

c1 = -10.728876d0
c2 = 3478.02d0
c3 = -383383d0
c4 = 33406000d0
tr = 1d0 !K

air_caww_m6mol2 = errorreturn

if (t_si <= 0d0) return

t = t_si / tr

c = -1d-6 * exp(c1 + c2 / t + c3 / t ** 2 + c4 / t ** 3)
select case (drv_t)
  case (0)
  case (1) 
    c = c * (-c2 / t ** 2 - 2d0 * c3 / t ** 3 - 3d0 * c4 / t ** 4)
  case (2)
    c = c * ((-c2 / t ** 2 - 2d0 * c3 / t ** 3  &
        - 3d0 * c4 / t ** 4) ** 2 + 2d0 * c2 / t ** 3  &
        + 6d0 * c3 / t ** 4 + 12d0 * c4 / t ** 5)
  case default; return
end select

air_caww_m6mol2 = c  

end function

!==========================================================================
subroutine init_iapws10()
!==========================================================================

integer :: k 

if (nk0(1) == 0.00000006057194d0) return

!ARRAYS DEFINED IN TABLE 12 OF LEMMON ET AL. (2000)
!COEFFICIENTS FOR THE IDEAL GAS EXPRESSIONS
!-------------------------------------------------------------------------------
k = 1;  nk0(k) = 0.00000006057194d0
k = 2;  nk0(k) = -0.0000210274769d0
k = 3;  nk0(k) = -0.000158860716d0
k = 4;  nk0(k) = 9.7450251743948d0   !-13.841928076 'original Lemmon et al. 2000
k = 5;  nk0(k) = 10.0986147428912d0  ! 17.275266575 'original Lemmon et al. 2000
k = 6;  nk0(k) = -0.00019536342d0
k = 7;  nk0(k) = 2.490888032d0
k = 8;  nk0(k) = 0.791309509d0
k = 9;  nk0(k) = 0.212236768d0
k = 10; nk0(k) = -0.197938904d0

k = 11; nk0(k) = 25.36365d0
k = 12; nk0(k) = 16.90741d0
k = 13; nk0(k) = 87.31279d0

!ARRAYS DEFINED IN TABLE 13 OF LEMMON ET AL. (2000)
!COEFFICIENTS AND EXPONENTS FOR THE EQUATION OF STATE FOR AIR
!-------------------------------------------------------------------------------
k = 1;  ik(k) = 1d0;  jk(k) = 0d0;    lk(k) = 0d0; nk(k) = 0.118160747229d0
k = 2;  ik(k) = 1d0;  jk(k) = 0.33d0; lk(k) = 0d0; nk(k) = 0.713116392079d0
k = 3;  ik(k) = 1d0;  jk(k) = 1.01d0; lk(k) = 0d0; nk(k) = -1.61824192067d0
k = 4;  ik(k) = 2d0;  jk(k) = 0d0;    lk(k) = 0d0; nk(k) = 0.0714140178971d0
k = 5;  ik(k) = 3d0;  jk(k) = 0d0;    lk(k) = 0d0; nk(k) = -0.0865421396646d0
k = 6;  ik(k) = 3d0;  jk(k) = 0.15d0; lk(k) = 0d0; nk(k) = 0.134211176704d0
k = 7;  ik(k) = 4d0;  jk(k) = 0d0;    lk(k) = 0d0; nk(k) = 0.0112626704218d0
k = 8;  ik(k) = 4d0;  jk(k) = 0.2d0;  lk(k) = 0d0; nk(k) = -0.0420533228842d0
k = 9;  ik(k) = 4d0;  jk(k) = 0.35d0; lk(k) = 0d0; nk(k) = 0.0349008431982d0
k = 10; ik(k) = 6d0;  jk(k) = 1.35d0; lk(k) = 0d0; nk(k) = 0.000164957183186d0

k = 11; ik(k) = 1d0;  jk(k) = 1.6d0;  lk(k) = 1d0; nk(k) = -0.101365037912d0
k = 12; ik(k) = 3d0;  jk(k) = 0.8d0;  lk(k) = 1d0; nk(k) = -0.17381369097d0
k = 13; ik(k) = 5d0;  jk(k) = 0.95d0; lk(k) = 1d0; nk(k) = -0.0472103183731d0
k = 14; ik(k) = 6d0;  jk(k) = 1.25d0; lk(k) = 1d0; nk(k) = -0.0122523554253d0
k = 15; ik(k) = 1d0;  jk(k) = 3.6d0;  lk(k) = 2d0; nk(k) = -0.146629609713d0
k = 16; ik(k) = 3d0;  jk(k) = 6d0;    lk(k) = 2d0; nk(k) = -0.0316055879821d0
k = 17; ik(k) = 11d0; jk(k) = 3.25d0; lk(k) = 2d0; nk(k) = 0.000233594806142d0
k = 18; ik(k) = 1d0;  jk(k) = 3.5d0;  lk(k) = 3d0; nk(k) = 0.0148287891978d0
k = 19; ik(k) = 3d0;  jk(k) = 15d0;   lk(k) = 3d0; nk(k) = -0.00938782884667d0


!dynamic adjustment to the reference state conditions:
!vanishing entropy s = 0 and enthalpy h = f + Ts + pv = 0 at 273.15 K, 101325 Pa
!related density is 1.29275937468394 kg m-3, the reducing temperature 132.6312 K
!nk0(4) = 0
!nk0(4) = -dry_f_si(1, 0, 273.15, 1.29275937468394) * Molar_mass_air_L2000 / (Gas_constant_molar_L2000)
!Debug.Print nk0(4); -dry_f_si(1, 0, 273.15, 1.29275937468394)

!nk0(5) = 0
!nk0(5) = -(dry_f_si(0, 0, 273.15, 1.29275937468394) + 101325 / 1.29275937468394) 
!         * Molar_mass_air_L2000 / (132.6312 * Gas_constant_molar_L2000)
!Debug.Print nk0(5); dry_f_si(0, 0, 273.15, 1.29275937468394) + 101325 / 1.29275937468394

!output:
! 9.7450251743948 -1.05459846463211E-12
! 10.0986147428912  2.75164779850456E-10

end subroutine

!==========================================================================
subroutine dry_init_clear()
!==========================================================================

!clears the initialization of dry_f_si()

nk0(1) = 0d0

end subroutine

!==========================================================================
subroutine dry_init_Lemmon2000()
!==========================================================================

!sets original adjustable coefficients of Lemmon et al (2000)

call init_iapws10

nk0(4) = -13.841928076d0 !original Lemmon et al. 2000
nk0(5) = 17.275266575d0  !original Lemmon et al. 2000

end subroutine

!==========================================================================
function alpha_ideal(tau)
!==========================================================================

!TEMPERATURE PART OF THE IDEAL-GAS DIMENSIONLESS HELMHOLTZ FUNCTION OF AIR

integer :: k
real*8 :: alpha_ideal
real*8 :: tau, a

alpha_ideal = errorreturn

if (tau <= 0d0) return

a = 0d0
do k = 1, 5
  a = a + nk0(k) * tau ** (k - 4)
enddo
a = a + nk0(6) * tau ** 1.5d0 + nk0(7) * log(tau)
a = a + nk0(8) * log(1d0 - exp(-nk0(11) * tau))
a = a + nk0(9) * log(1d0 - exp(-nk0(12) * tau))
a = a + nk0(10) * log(2d0 / 3d0 + exp(nk0(13) * tau))

alpha_ideal = a

end function

!==========================================================================
function alpha_ideal_t(tau)
!==========================================================================

!TEMPERATURE DERIVATIVE OF THE IDEAL-GAS DIMENSIONLESS HELMHOLTZ FUNCTION OF AIR

integer :: k
real*8 :: alpha_ideal_t
real*8 :: tau, a

alpha_ideal_t = errorreturn

if (tau <= 0d0) return

a = 0d0
do k = 1, 5
  a = a + nk0(k) * (k - 4d0) * tau ** (k - 5)
enddo
a = a + 1.5d0 * nk0(6) * sqrt(tau) + nk0(7) / tau
a = a + nk0(8) * nk0(11) / (exp(nk0(11) * tau) - 1d0)
a = a + nk0(9) * nk0(12) / (exp(nk0(12) * tau) - 1d0)
a = a + nk0(10) * nk0(13) / (2d0 * exp(-nk0(13) * tau) / 3d0 + 1d0)

alpha_ideal_t = a

end function

!==========================================================================
function alpha_ideal_tt(tau)
!==========================================================================

!2ND TEMPERATURE DERIVATIVE OF THE IDEAL-GAS DIMENSIONLESS HELMHOLTZ FUNCTION OF AIR

integer k
real*8 alpha_ideal_tt
real*8 :: tau, a, e

alpha_ideal_tt = errorreturn

if (tau <= 0d0) return

a = 0d0
do k = 1, 3
  a = a + nk0(k) * (k - 4d0) * (k - 5d0) * tau ** (k - 6d0)
enddo
a = a + 0.75d0 * nk0(6) / sqrt(tau) - nk0(7) / tau ** 2
e = exp(nk0(11) * tau)
a = a - nk0(8) * nk0(11) ** 2 * e / (e - 1d0) ** 2
e = exp(nk0(12) * tau)
a = a - nk0(9) * nk0(12) ** 2 * e / (e - 1d0) ** 2
e = exp(-nk0(13) * tau)
a = a + nk0(10) * nk0(13) ** 2 * e / (2d0 * e / 3d0 + 1d0) ** 2 * 2d0 / 3d0

alpha_ideal_tt = a

end function

!==========================================================================
function alpha_res(drv_tau, drv_del, tau, del)
!==========================================================================

!THIS IMPLEMENTS EQUATION (28) AND ITS DERIVATIVES WITH RESPECT TO TAU AND DEL OF
!LEMMON, E.W., JACOBSEN, R.T., PENONCELLO, S.G., FRIEND, D.G.
!THERMODYNAMIC PROPERTIES OF AIR AND MIXTURES OF NITROGEN, ARGON AND OXYGEN
!FROM 60 TO 2000 K AT PRESSURES TO 2000 MPA.
!J. PHYS. CHEM. REF. DATA 29(2000)331-362.

integer :: drv_tau, drv_del, k
real*8 alpha_res
real*8 :: tau, del, a, a_tau, a_del, i, j, l

alpha_res = errorreturn

a = 0d0

do k = 1, 19
  i = ik(k)
  j = jk(k)
  l = lk(k)
  select case (drv_tau)
    case (0); a_tau = tau ** j
    case (1); a_tau = j * tau ** (j - 1)
    case (2); a_tau = j * (j - 1d0) * tau ** (j - 2)
    case default; return
  end select
  if (k < 11) then
    select case (drv_del)
      case (0); a_del = del ** i
      case (1); a_del = i * del ** (i - 1)
      case (2); a_del = i * (i - 1d0) * del ** (i - 2)
      case default; return
    end select
  else
    select case (drv_del)
      case (0); a_del = del ** i * exp(-del ** l)
      case (1); a_del = (i - l * del ** l) * del ** (i - 1d0) * exp(-del ** l)
      case (2); a_del = (i * (i - 1d0) - (l + 2d0 * i - 1d0) * l * del ** l &
                      + (l * del ** l) ** 2) * del ** (i - 2d0) * exp(-del ** l)
      case default; return
    end select
  end if
  a = a + nk(k) * a_tau * a_del
enddo

alpha_res = a

end function

end module air_1
