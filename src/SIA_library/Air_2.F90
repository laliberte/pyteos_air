module air_2

!#########################################################################

!THIS MODULE IMPLEMENTS THE AIR-WATER MIXING HELMHOLTZ FUNCTION
!AIR_F_MIX_SI(DRV_A, DRV_T, DRV_D, A_SI, A_SI, T_SI, D_SI),
!AND ITS FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO AIR FRACTION,
!TEMPERATURE AND DENSITY,

!THE HELMHOLTZ POTENTIAL OF HUMID AIR AS A FUNCTION OF THE DRY-AIR FRACTION
!IN KG/KG, ABSOLUTE TEMPERATURE IN K AND HUMID-AIR DENSITY IN KG/M3,
!AIR_F_SI(DRV_A, DRV_T, DRV_D, A_SI, A_SI, T_SI, D_SI),
!AND ITS FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO AIR FRACTION,
!TEMPERATURE AND DENSITY,

!AND 

!THERMODYNAMIC PROPERTIES OF HUMID AIR COMPUTED FROM THE hELMHOLTZ POTENTIAL.

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

!THIS MODULE REQUIRES THE LIBRARY MODULES
!     CONSTANTS_0,  FILE CONSTANTS_0.F90
!     FLU_1,        FILE FLU_1.F90
!     AIR_1,        FILE AIR_1.F90

!#########################################################################

use constants_0
use air_1
use flu_1

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: air_f_si, air_f_cp_si, air_f_cv_si, &
          air_f_enthalpy_si, air_f_entropy_si, &
          air_f_expansion_si, air_f_gibbs_energy_si, &
          air_f_internal_energy_si, air_f_kappa_s_si, &
          air_f_kappa_t_si, air_f_lapserate_si, &
          air_f_mix_si, air_f_pressure_si, &
          air_f_soundspeed_si, chk_iapws10_table

contains

!==========================================================================
function air_f_mix_si(drv_a, drv_t, drv_d, a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE CROSS-COMPONENT INTERACTION CONTRIBUTION TO 
!THE SPECIFIC  HELMHOLTZ ENERGY OF HUMID AIR AND ITS DERIVATIVES WITH RESPECT 
!TO THE AIR FRACTION, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!air_f_mix_si:  MIXING HELMHOLTZ ENERGY OR ITS DERIVATIVES IN J KG-1 K^-DRV_T (KG/M3)^-DRV_D

!INPUT:  
!DRV_A  ORDER OF DERIVATIVE WITH RESPECT TO THE AIR MASS FRACTION
!DRV_T  ORDER OF DERIVATIVE WITH RESPECT TO TEMPERATURE
!DRV_D  ORDER OF DERIVATIVE WITH RESPECT TO HUMID-AIR DENSITY
!A_SI   MASS FRACTION OF AIR IN HUMID AIR IN KG/KG
!T_SI   ABSOLUTE TEMPERATURE IN K
!D_SI   HUMID-AIR DENSITY IN KG/M3

!CHECK VALUES, V1.0 AND V1.1:
!air_f_mix_si( 0, 0, 0, 0.9,300,1) =-25.8379179182
!air_f_mix_si( 1, 0, 0, 0.9,300,1) = 233.827370434
!air_f_mix_si( 0, 1, 0, 0.9,300,1) = 0.164195952060
!air_f_mix_si( 0, 0, 1, 0.9,300,1) =-26.2357498619
!air_f_mix_si( 2, 0, 0, 0.9,300,1) = 500.273928155
!air_f_mix_si( 1, 1, 0, 0.9,300,1) =-1.53932744055
!air_f_mix_si( 1, 0, 1, 0.9,300,1) = 241.520643317
!air_f_mix_si( 0, 2, 0, 0.9,300,1) =-0.687329742959E-3
!air_f_mix_si( 0, 1, 1, 0.9,300,1) = 0.172192606103
!air_f_mix_si( 0, 0, 2, 0.9,300,1) =-0.795663887493


real*8 air_f_mix_si
integer :: drv_a, drv_t, drv_d
real*8 :: raw, a_si, t_si, d_si
real*8 :: baw, caaw, caww
real*8 :: bawt, caawt, cawwt

real*8 ::  mw = molar_mass_h2o_si   !MOLAR MASS OF H2O IN KG/MOL
real*8 ::  ma = molar_mass_air_si   !MOLAR MASS OF AIR IN KG/MOL
raw = gas_constant_molar_si / (ma * mw)

air_f_mix_si = errorreturn

if (drv_a < 0) return
if (drv_t < 0) return
if (drv_d < 0) return
if (drv_a + drv_t + drv_d > 2) return

if(a_si < 0d0 .or. a_si > 1d0) return

if(check_limits == 1) then
  if(t_si < mix_air_tmin .or. t_si > mix_air_tmax) return
  if(d_si <= mix_air_dmin .or. d_si > mix_air_dmax) return 
else
  if (t_si <= 0d0 .or. &
  d_si <= 0d0) return
endif

!TEMPERATURE DERIVATIVES OF THE VIRIAL COEFFICIENTS
select case (drv_t)
  case (0)
     baw = air_baw_m3mol(0, t_si)
     If(baw == ErrorReturn) return
     bawt = 2d0 * raw * t_si * baw 

     caaw = air_caaw_m6mol2(0, t_si)
     If(caaw == ErrorReturn) return
     caawt = 1.5d0 * Raw * t_si * caaw / ma

     caww = air_caww_m6mol2(0, t_si)
     If(caww == ErrorReturn) return
     cawwt = 1.5d0 * Raw * t_si * caww / mw

  case (1) 
     baw = air_baw_m3mol(1, t_si)
     If(baw == ErrorReturn) return
     bawt = 2d0 * Raw * t_si * baw
     baw = air_baw_m3mol(0, t_si)
     If(baw == ErrorReturn) return
     bawt = (bawt + 2d0 * Raw * baw)

     caaw = air_caaw_m6mol2(1, t_si)
     If(caaw == ErrorReturn) return
     caawt = 1.5d0 * Raw * t_si * caaw / ma
     caaw = air_caaw_m6mol2(0, t_si)
     If(caaw == ErrorReturn) return
     caawt = (caawt + 1.5d0 * Raw * caaw / ma)

     caww = air_caww_m6mol2(1, t_si)
     If(caww == ErrorReturn) return
     cawwt = 1.5d0 * Raw * t_si * caww / mw
     caww = air_caww_m6mol2(0, t_si)
     If(caww == ErrorReturn) return
     cawwt = (cawwt + 1.5d0 * Raw * caww / mw)

  case (2)
     baw = air_baw_m3mol(2, t_si)
     If(baw == ErrorReturn) return
     bawt = 2d0 * Raw * t_si * baw
     baw = air_baw_m3mol(1, t_si)
     If(baw == ErrorReturn) return
     bawt = (bawt + 4d0 * Raw * baw)

     caaw = air_caaw_m6mol2(2, t_si)
     If(caaw == ErrorReturn) return
     caawt = 1.5d0 * Raw * t_si * caaw / ma
     caaw = air_caaw_m6mol2(1, t_si)
     If(caaw == ErrorReturn) return
     caawt = caawt + 1.5d0 * Raw * 2d0 * caaw / ma

     caww = air_caww_m6mol2(2, t_si)
     If(caww == ErrorReturn) return
     cawwt = 1.5d0 * Raw * t_si * caww / mw
     caww = air_caww_m6mol2(1, t_si)
     If(caww == ErrorReturn) return
     cawwt = cawwt + 1.5d0 * Raw * 2d0 * caww / mw
  case default; return
end select

!DENSITY DERIVATIVES OF THE VIRIAL COEFFICIENTS
select case (drv_d)
  case (0)
     bawt  = bawt * d_si
     caawt = caawt * d_si ** 2
     cawwt = cawwt * d_si ** 2
  case (1)
     bawt  = bawt
     caawt = 2d0 * caawt * d_si
     cawwt = 2d0 * cawwt * d_si
  case (2)
     bawt  = 0d0
     caawt = 2d0 * caawt
     cawwt = 2d0 * cawwt
  case default;  return
end select

!AIR-FRACTION DERIVATIVES OF THE VIRIAL COEFFICIENTS
select case (drv_a)
  case (0)
     bawt  = bawt * a_si * (1d0 - a_si)
     caawt = caawt * a_si ** 2 * (1d0 - a_si)
     cawwt = cawwt * a_si * (1d0 - a_si) ** 2
  case (1)
     bawt  = bawt * (1d0 - 2d0 * a_si)
     caawt = caawt * a_si * (2d0 - 3d0 * a_si)
     cawwt = cawwt * (1d0 - 4d0 * a_si + 3d0 * a_si ** 2)
  case (2)
     bawt  = -2d0 * bawt
     caawt = caawt * (2d0 - 6d0 * a_si)
     cawwt = cawwt * (-4d0 + 6d0 * a_si)
  case default; return
end select

air_f_mix_si = bawt + caawt + cawwt

end function

!==========================================================================
function air_f_si(drv_a, drv_t, drv_d, a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE HELMHOLTZ POTENTIAL OF HUMID AIR AND ITS FIRST AND SECOND
!DERIVATIVES WITH RESPECT TO THE AIR FRACTION, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT:  
!AIR_F_SI = HELMHOLTZ ENERGY OR ITS DERIVATIVES IN J KG-1 K^-DRV_T (KG/M3)^-DRV_D

!INPUT:   
!DRV_A  ORDER OF THE AIR FRACTION DERIVATIVE
!DRV_T  ORDER OF THE TEMPERATURE DERIVATIVE
!DRV_D  ORDER OF THE DENSITY DERIVATIVE
!A_SI MASS FRACTION OF DRY AIR IN HUMID AIR IN KG/KG
!T_SI ABSOLUTE TEMPERATURE IN K
!D_SI DENSITY OF HUMID AIR IN KG/M3

!CHECK VALUES, V1.0:
!air_f_si( 0, 0, 0, 0.9,300,1) =-95024.7900038
!air_f_si( 1, 0, 0, 0.9,300,1) =-205630.937996
!air_f_si( 0, 1, 0, 0.9,300,1) =-940.150652241
!air_f_si( 0, 0, 1, 0.9,300,1) = 91193.7358460
!air_f_si( 2, 0, 0, 0.9,300,1) = 1447791.10862
!air_f_si( 1, 1, 0, 0.9,300,1) = 7443.19342188
!air_f_si( 1, 0, 1, 0.9,300,1) =-48827.5293333
!air_f_si( 0, 2, 0, 0.9,300,1) =-2.96533251321
!air_f_si( 0, 1, 1, 0.9,300,1) = 312.124501058
!air_f_si( 0, 0, 2, 0.9,300,1) =-91439.8037142

!CHECK VALUES, V1.1:
!air_f_si( 0, 0, 0, 0.9,300,1) =-95019.5943231
!air_f_si( 1, 0, 0, 0.9,300,1) =-205645.554995
!air_f_si( 0, 1, 0, 0.9,300,1) =-940.175394023
!air_f_si( 0, 0, 1, 0.9,300,1) = 91175.3848662
!air_f_si( 2, 0, 0, 0.9,300,1) = 1447768.46379
!air_f_si( 1, 1, 0, 0.9,300,1) = 7443.09771950
!air_f_si( 1, 0, 1, 0.9,300,1) =-48847.9096826
!air_f_si( 0, 2, 0, 0.9,300,1) =-2.96482218054
!air_f_si( 0, 1, 1, 0.9,300,1) = 312.063110700
!air_f_si( 0, 0, 2, 0.9,300,1) =-91421.4440689


real*8 :: air_f_si
integer :: drv_a, drv_t, drv_d
real*8 :: a_si, t_si, d_si
real*8 :: f, da, dw
real*8 :: fw, fa, fm

air_f_si = errorreturn

if (drv_a < 0) return
if (drv_t < 0) return
if (drv_d < 0) return
if (drv_a + drv_t + drv_d > 2) return

if (a_si < 0d0 .or. a_si > 1d0) return
if ((a_si == 0d0 .or. a_si == 1d0) .and. drv_a > 0) return  
if (t_si <= 0d0) return
if (d_si <= 0d0) return

da = a_si * d_si            !PARTIAL DENSITY OF AIR
dw = (1d0 - a_si) * d_si    !PARTIAL DENSITY OF VAPOUR

select case (drv_a)

  case (0)
    f = 0d0
    if(dw > 0d0)then
      fw = flu_f_si(drv_t, drv_d, t_si, dw)
      if(fw == errorreturn) return
      f = f + (1d0 - a_si) ** (drv_d + 1) * fw
    endif
    if(da > 0d0)then
      fa = dry_f_si(drv_t, drv_d, t_si, da)
      if(fa == errorreturn) return
      f = f + a_si ** (drv_d + 1) * fa
    endif
    if(dw > 0d0 .and. da > 0d0)then
      fm = air_f_mix_si(0, drv_t, drv_d, a_si, t_si, d_si)
      if(fm == errorreturn) return
      f = f + fm
    endif

  case (1)
    select case (drv_d)

      case (0)
        f = 0d0
        if(dw > 0d0)then
          fw = flu_f_si(drv_t, 0, t_si, dw) 
          if(fw == errorreturn) return
          f = f - fw
          fw = flu_f_si(drv_t, 1, t_si, dw)
          if(fw == errorreturn) return
          f = f - dw * fw
        endif
        if(da > 0d0)then
          fa = dry_f_si(drv_t, 0, t_si, da)
          if(fa == errorreturn) return
          f = f + fa
          fa = dry_f_si(drv_t, 1, t_si, da)
          if(fa == errorreturn) return
          f = f + da * fa
        endif
        if(dw > 0d0 .and. da > 0d0)then
          fm = air_f_mix_si(1, drv_t, 0, a_si, t_si, d_si)
          if(fm == errorreturn) return
          f = f + fm
        endif

      case (1)
        f = 0d0
        if(dw > 0d0)then
          fw = flu_f_si(drv_t, 1, t_si, dw)
          if(fw == errorreturn) return
          f = f - (1d0 - a_si) * 2d0 * fw
          fw = flu_f_si(drv_t, 2, t_si, dw)
          if(fw == errorreturn) return
          f = f - (1d0 - a_si) * dw * fw
        endif
        if(da > 0d0)then
          fa =  dry_f_si(drv_t, 1, t_si, da)
          if(fa == errorreturn) return
          f = f + a_si * 2d0 * fa
          fa =  dry_f_si(drv_t, 2, t_si, da)
          if(fa == errorreturn) return
          f = f + a_si * da * fa
        endif
        if(dw > 0d0 .and. da > 0d0)then
          fm = air_f_mix_si(1, drv_t, drv_d, a_si, t_si, d_si)
          if(fm == errorreturn) return
          f = f + fm
        endif
      case default; return
    end select

  case (2)
    f = 0d0
    if(dw > 0d0)then
      fw = flu_f_si(drv_t, 1, t_si, dw)
      if(fw == errorreturn) return
      f = f + 2d0 * d_si * fw
      fw = flu_f_si(drv_t, 2, t_si, dw)
      if(fw == errorreturn) return
      f = f + dw * d_si * fw
    endif
    if(da > 0d0)then
      fa = dry_f_si(drv_t, 1, t_si, da)
      if(fa == errorreturn) return
      f = f + 2d0 * d_si * fa
      fa = dry_f_si(drv_t, 2, t_si, da)
      if(fa == errorreturn) return
      f = f + da * d_si * fa
    endif
    if(dw > 0d0 .and. da > 0d0)then
      fm = air_f_mix_si(2, drv_t, 0, a_si, t_si, d_si)
      if(fm == errorreturn) return
      f = f + fm
    endif
  case default; return
end select

air_f_si = f

end function

!==========================================================================
function air_f_cp_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISOBARIC HEAT CAPACITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!CP = T * (DETA/DT)_P  ISOBARIC HEAT CAPACITY IN J/(KG K)

!INPUT: 
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_cp_si(0.9,300,1) = 1210.95501131
!CHECK VALUE, V1.1: air_f_cp_si(0.9,300,1) = 1210.74031058

real*8 air_f_cp_si
real*8 a_si, t_si, d_si
real*8 f_d, f_tt, f_td, f_dd, x

air_f_cp_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
if (f_tt == errorreturn) return
f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
if (f_td == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

x = 2d0 * f_d + d_si * f_dd
if (x == 0d0) return

air_f_cp_si = t_si * (f_td ** 2 * d_si / x - f_tt)

end function

!==========================================================================
function air_f_cv_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISOCHORIC HEAT CAPACITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!CV = T * (DETA/DT)_D  ISOCHORIC HEAT CAPACITY IN J/(KG K)

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_cv_si(0.9,300,1) = 889.599753963
!CHECK VALUE, V1.1: air_f_cv_si(0.9,300,1) = 889.446654163

real*8 :: air_f_cv_si
real*8 :: a_si, t_si, d_si
real*8 :: f_tt

air_f_cv_si = errorreturn

f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
if (f_tt == errorreturn) return

air_f_cv_si = -t_si * f_tt

end function

!==========================================================================
function air_f_enthalpy_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ENTHALPY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY


!OUTPUT: 
!H = F - T*(DF/DT) + D*(DF/DD) ENTHALPY IN J/KG

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_enthalpy_si(0.9,300,1) = 278214.141515
!CHECK VALUE, V1.1: air_f_enthalpy_si(0.9,300,1) = 278208.408750

real*8 :: air_f_enthalpy_si
real*8 :: a_si, t_si, d_si
real*8 :: f, f_t, f_d

air_f_enthalpy_si = errorreturn

f = air_f_si(0, 0, 0, a_si, t_si, d_si)
if (f == errorreturn) return
f_t = air_f_si(0, 1, 0, a_si, t_si, d_si)
if (f_t == errorreturn) return
f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return

air_f_enthalpy_si = f - t_si * f_t + d_si * f_d

end function

!==========================================================================
function air_f_entropy_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ENTROPY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!ETA = - (DF/DT)_A_D ENTROPY in J/(KG K)

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_entropy_si(0.9,300,1) = 940.150652241
!CHECK VALUE, V1.1: air_f_entropy_si(0.9,300,1) = 940.175394023

real*8 :: air_f_entropy_si
real*8 :: a_si, t_si, d_si
real*8 :: f_t

air_f_entropy_si = errorreturn

f_t = air_f_si(0, 1, 0, a_si, t_si, d_si)
if (f_t == errorreturn) return

air_f_entropy_si = -f_t

end function

!==========================================================================
function air_f_expansion_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE THERMAL EXPANSION COEFFICIENT OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!ALPHA = -(1/D) * (DD/DT)_A_P  THERMAL EXPANSION COEFFICIENT IN 1/K

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_expansion_si(0.9,300,1) = 3.43191318698E-03
!CHECK VALUE, V1.1: air_f_expansion_si(0.9,300,1) = 3.43193033077E-03

real*8 :: air_f_expansion_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d, f_td, f_dd, x

air_f_expansion_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
if (f_td == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

x = 2d0 * f_d + d_si * f_dd
if (x == 0d0) return

air_f_expansion_si = f_td / x

end function

!==========================================================================
function air_f_gibbs_energy_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE GIBBS ENERGY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!G = F + D*(DF/DD)   GIBBS ENERGY IN J/KG

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_gibbs_energy_si(0.9,300,1) = -3831.05415778
!CHECK VALUE, V1.1: air_f_gibbs_energy_si(0.9,300,1) = -3844.20945693

real*8 :: air_f_gibbs_energy_si
real*8 :: a_si, t_si, d_si
real*8 :: f, f_d

air_f_gibbs_energy_si = errorreturn

f = air_f_si(0, 0, 0, a_si, t_si, d_si)
if (f == errorreturn) return
f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return

air_f_gibbs_energy_si = f + d_si * f_d

end function

!==========================================================================
function air_f_internal_energy_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE INTERNAL ENERGY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!U = F - T*(DF/DT) INTERNAL ENERGY IN J/KG

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_internal_energy_si(0.9,300,1) = 187020.405669
!CHECK VALUE, V1.1: air_f_internal_energy_si(0.9,300,1) = 187033.023884

real*8 :: air_f_internal_energy_si
real*8 :: a_si, t_si, d_si
real*8 :: f, f_t

air_f_internal_energy_si = errorreturn

f = air_f_si(0, 0, 0, a_si, t_si, d_si)
if (f == errorreturn) return
f_t = air_f_si(0, 1, 0, a_si, t_si, d_si)
if (f_t == errorreturn) return

air_f_internal_energy_si = f - t_si * f_t

end function

!==========================================================================
function air_f_kappa_s_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISENTROPIC COMPRESSIBILITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!KAPPA_S = (1/D) * (DD/DP)_A_ETA  ISENTROPIC COMPRESSIBILITY IN 1/PA

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_kappa_s_si(0.9,300,1) = 8.07746486644E-06
!CHECK VALUE, V1.1: air_f_kappa_s_si(0.9,300,1) = 8.07913626816E-06

real*8 :: air_f_kappa_s_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d, f_tt, f_td, f_dd, x

air_f_kappa_s_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
if (f_tt == errorreturn) return
f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
if (f_td == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

x = d_si ** 2 * (f_tt * (2d0 * f_d + d_si * f_dd) - d_si * f_td ** 2)
if (x == 0d0) return

air_f_kappa_s_si = f_tt / x

end function

!==========================================================================
function air_f_kappa_t_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ISOTHERMAL COMPRESSIBILITY OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!KAPPA_T = (1/D) * (DD/DP)_A,T  ISOTHERMAL COMPRESSIBILITY IN 1/PA

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_kappa_t_si(0.9,300,1) = 1.09953341546E-05
!CHECK VALUE, V1.1: air_f_kappa_t_si(0.9,300,1) = 1.09975521396E-05

real*8 :: air_f_kappa_t_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d, f_dd, x

air_f_kappa_t_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

x = d_si ** 2 * (2d0 * f_d + d_si * f_dd)
if (x == 0d0) return

air_f_kappa_t_si = 1d0 / x

end function

!==========================================================================
function air_f_lapserate_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE "DRY" ADIABATIC LAPSE RATE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!GAMMA = (DT/DP)_A_ETA  "DRY" ADIABATIC LAPSE RATE IN K/PA

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_lapserate_si(0.9,300,1) = 8.50216520413E-04
!CHECK VALUE, V1.1: air_f_lapserate_si(0.9,300,1) = 8.50371537341E-04

real*8 :: air_f_lapserate_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d, f_tt, f_td, f_dd, x

air_f_lapserate_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
if (f_tt == errorreturn) return
f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
if (f_td == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

x = d_si * (f_tt * (2d0 * f_d + d_si * f_dd) - d_si * f_td ** 2)
if (x == 0d0) return

air_f_lapserate_si = -f_td / x

end function

!==========================================================================
function air_f_pressure_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE PRESSURE OF HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!P = D^2*(DF/DD) PRESSURE IN PA

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_pressure_si(0.9,300,1) = 91193.7358460
!CHECK VALUE, V1.1: air_f_pressure_si(0.9,300,1) = 91175.3848662

real*8 :: air_f_pressure_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d

air_f_pressure_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return

air_f_pressure_si = d_si ** 2 * f_d

end function

!==========================================================================
function air_f_soundspeed_si(a_si, t_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE SOUND SPEED IN HUMID AIR
!FROM THE MASS FRACTION OF AIR, TEMPERATURE AND HUMID-AIR DENSITY

!OUTPUT: 
!C = SQRT[(DP/DD)_A_ETA]  SOUND SPEED IN M/S

!INPUT:  
!A_SI      MASS FRACTION OF DRY AIR IN HUMID AIR
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      HUMID-AIR DENSITY IN KG/M3

!CHECK VALUE, V1.0: air_f_soundspeed_si(0.9,300,1) = 351.853974518
!CHECK VALUE, V1.1: air_f_soundspeed_si(0.9,300,1) = 351.817577078

real*8 :: air_f_soundspeed_si
real*8 :: a_si, t_si, d_si
real*8 :: f_d, f_tt, f_td, f_dd, x

air_f_soundspeed_si = errorreturn

f_d = air_f_si(0, 0, 1, a_si, t_si, d_si)
if (f_d == errorreturn) return
f_tt = air_f_si(0, 2, 0, a_si, t_si, d_si)
if (f_tt == errorreturn) return
f_td = air_f_si(0, 1, 1, a_si, t_si, d_si)
if (f_td == errorreturn) return
f_dd = air_f_si(0, 0, 2, a_si, t_si, d_si)
if (f_dd == errorreturn) return

if (f_tt == 0d0) return
x = 2d0 * d_si * f_d + d_si ** 2 * (f_dd - f_td ** 2 / f_tt)
if (x < 0d0) return

air_f_soundspeed_si = sqrt(x)

end function


!=========================================================================
subroutine chk_iapws10_table(number)
!=========================================================================

!THIS SUBROUTINE PRINTS SELECTED CHECK TABLES OF IAPWS-10
!FOR EASY VERIFICATION OF THIS IMPLEMENTATION ON A GIVEN PLATFORM

!GUIDELINE ON AN EQUATION OF STATE FOR HUMID AIR IN CONTACT WITH SEAWATER AND ICE,
!CONSISTENT WITH THE IAPWS FORMULATION 2008 FOR THE THERMODYNAMIC PROPERTIES OF SEAWATER
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!NIAGARA FALLS, CANADA, SEPTEMBER 2010 IN PREPARATION

integer number

select case (number)
  case(13); call chk_iapws10_table13()
  case(14); call chk_iapws10_table14()
  case(15); call chk_iapws10_table15()
  case default; write(*,*) "Tables 13, 14 and 15 are implemented here"; return
end select

end subroutine


!=========================================================================
subroutine chk_iapws10_table13()
!=========================================================================

!TABLE 13:  NUMERICAL CHECK VALUES FOR THE HELMHOLTZ FUNCTION F AND ITS DERIVATIVES,
!TABLE 6, FOR SATURATED AIR, WITH RESPECT TO ICE AT 200 K, AND WITH RESPECT TO LIQUID 
!WATER AT 300 AND 400 K. THE NUMERICAL FUNCTIONS EVALUATED HERE AT GIVEN POINTS
!(A, T, D) ARE DEFINED IN TABLES 5, 6 AND 12.
!
!Quantity       Value               Value              Value              Unit
!
!A              0.892247719         0.977605798        0.825565291        kg kg-1
!T              200                 300                400                K
!d              0.163479657E-4      0.114614216E+1     0.793354063E+1     kg m-3
!p              0.999999998         0.100000000E+6     0.100000000E+7     Pa
!fAV           -0.682093392E+6     -0.927718178E+5     0.240345570E+5     J kg-1
!(df/dA)_T,d   -0.572680404E+6     -0.263453864E+3     0.311096733E+6     J kg-1
!(df/dT)_A,d   -0.405317966E+4     -0.296711481E+3    -0.106891931E+4     J kg-1
!(df/dd)_A,T    0.374173101E+10     0.761242496E+5     0.158878781E+5     m3 kg-1
!(d2f/dA2)_T,d  0.920967684E+6      0.624886233E+7     0.113786423E+7     J kg-1
!(d2f/dAdT)_d   0.915653743E+4      0.822733446E+4     0.702631471E+4     J kg-1
!(d2f/dAdd)    -0.213442099E+10    -0.450004399E+5    -0.727972651E+4     m3 kg-1
!(d2f/dT2)_d   -0.394011921E+1     -0.244742952E+1    -0.222449294E+1     J kg-1
!(d2f/dTdd)     0.187087034E+8      0.254456302E+3     0.414350772E+2     m3 kg-1
!(d2f/dd2)_T   -0.228880603E+15    -0.664465525E+5    -0.201886184E+4     m3 kg-1
!h              0.189712231E+6      0.834908383E+5     0.577649408E+6     J kg-1
!g             -0.620923701E+6     -0.552260595E+4     0.150081684E+6     J kg-1
!s              0.405317966E+4      0.296711481E+3     0.106891931E+4     J kg-1
!muW           -0.109950917E+6     -0.526505193E+4    -0.106748981E+6     J kg-1
!cp             0.109387397E+4      0.102681324E+4     0.123552454E+4     J kg-1
!w              0.291394959E+3      0.349234196E+3     0.416656820E+3     m s-1
!

real*8 a, t, d, q(3)
integer i, j

                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) " Implementation of IAPWS-10 in Fortran"
write(*,*) " for Publication in Ocean Science, 2009"
write(*,*) " D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) " NUMERICAL CHECK VALUES FOR THE HELMHOLTZ FUNCTION AND ITS DERIVATIVES"
write(*,*) " FOR SATURATED AIR, WRT ICE AT 200 K AND WRT LIQUID WATER AT 300 AND 400 K"
write(*,*) " FUNCTION VALUES AS GIVEN IN TABLE 13 OF IAPWS-10."
write(*,*) " TOP ROW GIVES TABULATED VALUES, BOTTOM ROW COMPUTED NOW"
write(*,*) " "

write(*,*)               "Quantity        Value            Value            Value           Unit"
write(*,*) " "

do i = 1, 20
  select case(i)
    case(1);  write(*,*) "A               0.892247719      0.977605798      0.825565291     kg kg-1"
    case(2);  write(*,*) "T               200              300              400             K"
    case(3);  write(*,*) "d               0.163479657E-4   0.114614216E+1   0.793354063E+1  kg m-3"
    case(4);  write(*,*) "p               0.999999998      0.100000000E+6   0.100000000E+7  Pa"
    case(5);  write(*,*) "f              -0.682093392E+6  -0.927718178E+5   0.240345570E+5  J kg-1"
    case(6);  write(*,*) "(df/dA)_T,d    -0.572680404E+6  -0.263453864E+3   0.311096733E+6  J kg-1"
    case(7);  write(*,*) "(df/dT)_A,d    -0.405317966E+4  -0.296711481E+3  -0.106891931E+4  J kg-1"
    case(8);  write(*,*) "(df/dd)_A,T     0.374173101E+10  0.761242496E+5   0.158878781E+5  m3 kg-1"
    case(9);  write(*,*) "(d2f/dA2)_T,d   0.920967684E+6   0.624886233E+7   0.113786423E+7  J kg-1"
    case(10); write(*,*) "(d2f/dAdT)_d    0.915653743E+4   0.822733446E+4   0.702631471E+4  J kg-1"
    case(11); write(*,*) "(d2f/dAdd)_T   -0.213442099E+10 -0.450004399E+5  -0.727972651E+4  m3 kg-1"
    case(12); write(*,*) "(d2f/dT2)_d    -0.394011921E+1  -0.244742952E+1  -0.222449294E+1  J kg-1"
    case(13); write(*,*) "(d2f/dTdd)      0.187087034E+8   0.254456302E+3   0.414350772E+2  m3 kg-1"
    case(14); write(*,*) "(d2f/dd2)_T    -0.228880603E+15 -0.664465525E+5  -0.201886184E+4  m3 kg-1"
    case(15); write(*,*) "h               0.189712231E+6   0.834908383E+5   0.577649408E+6  J kg-1"
    case(16); write(*,*) "g              -0.620923701E+6  -0.552260595E+4   0.150081684E+6  J kg-1"
    case(17); write(*,*) "s               0.405317966E+4   0.296711481E+3   0.106891931E+4  J kg-1"
    case(18); write(*,*) "muW            -0.109950917E+6  -0.526505193E+4  -0.106748981E+6  J kg-1"
    case(19); write(*,*) "cp              0.109387397E+4   0.102681324E+4   0.123552454E+4  J kg-1"
    case(20); write(*,*) "w               0.291394959E+3   0.349234196E+3   0.416656820E+3  m s-1"
  end select

  do j = 1, 3
    if(j == 1)then
       a = 0.892247719d0
       t = 200d0
       d = 0.163479657d-4
    elseif (j == 2) then
       a = 0.977605798d0
       t = 300d0
       d = 0.114614216d+1
    elseif(j == 3) then
       a = 0.825565291d0
       t = 400d0
       d = 0.793354063d+1
    endif
    q(j) = 0d0
    select case(i)
      case(1); q(j) = a
      case(2); q(j) = t
      case(3); q(j) = d
      case(4); q(j) = air_f_pressure_si(a, t, d)
      case(5); q(j) = air_f_si(0, 0, 0, a, t, d)
      case(6); q(j) = air_f_si(1, 0, 0, a, t, d)
      case(7); q(j) = air_f_si(0, 1, 0, a, t, d)
      case(8); q(j) = air_f_si(0, 0, 1, a, t, d)
      case(9); q(j) = air_f_si(2, 0, 0, a, t, d)
      case(10); q(j) = air_f_si(1, 1, 0, a, t, d)
      case(11); q(j) = air_f_si(1, 0, 1, a, t, d)
      case(12); q(j) = air_f_si(0, 2, 0, a, t, d)
      case(13); q(j) = air_f_si(0, 1, 1, a, t, d)
      case(14); q(j) = air_f_si(0, 0, 2, a, t, d)
      case(15); q(j) = air_f_enthalpy_si(a, t, d)
      case(16); q(j) = air_f_gibbs_energy_si(a, t, d)
      case(17); q(j) = air_f_entropy_si(a, t, d)
      case(18); q(j) = air_f_gibbs_energy_si(a, t, d)
                if(a > 0d0) q(j) = q(j) - a * air_f_si(1, 0, 0, a, t, d)
      case(19); q(j) = air_f_cp_si(a, t, d)
      case(20); q(j) = air_f_soundspeed_si(a, t, d)
    end select
  enddo
  write(*,'( 16x, E16.9,x,E16.9,x,E16.9)')q(1), q(2), q(3)
  write(*,*)' '
enddo
write(*,*)' '
write(*,*)' '
write(*,*)' '

end subroutine

!=========================================================================
subroutine chk_iapws10_table14()
!=========================================================================

!TABLE 14  Numerical check values for the dry-air and water-vapour parts of
!the Helmholtz function f and its derivatives, Table 6, for saturated air, 
!A = Asat, with respect to ice at 200 K, and with respect to liquid water 
!at 300 and 400 K. The numerical functions evaluated here at given
!points (A, T, d) are defined in Tables 7 and 8.
 
!Quantity         Value               Value              Value              Unit
!A                0.892224944         0.977600624        0.825531379        kg kg-1
!T                200                 300                400                K
!d                0.163445112E-04     0.114587678E+1    0.793198757E+1    kg m-3
!Ad               0.145829805901E-4   0.112020985516E+1  0.654810463687E+1  kg m-3
!fA              -0.740041144E+6     -0.916103453E+5     0.895561286E+5     J kg-1
!(dfA/dT)_dA     -0.304774177E+4     -0.108476220E+3     0.193271394E+3     J kg-1 K-1
!(dfA/ddA)_T      0.393583654E+10     0.768326795E+5     0.175560114E+5     m3 kg-1
!(d2fA/dT2)_dA   -0.357677878E+1     -0.239319940E+1    -0.181809877E+1     J kg-1 K-2
!(d2fA/dTddA)     0.196791837E+8      0.256683306E+3     0.442769673E+2     m3 kg-1 K-1
!(d2fA/(ddA)2)_T -0.269828549E+15    -0.685917373E+5    -0.267635928E+4     m3 kg-1 Pa-1
!(1-A)d           0.176153059E-5      0.256669391E-1     0.138388485E+1     kg m-3`
!fV              -0.202254351E+6     -0.143157426E+6    -0.285137534E+6     J kg-1
!(dfV/dT)_Ad     -0.123787544E+5     -0.851598213E+4    -0.705288048E+4     J kg-1 K-1
!(dfV/ddV)_T      0.523995674E+11     0.538480619E+7     0.129645039E+6     m3 kg-1
!(d2fV/dT2)_Ad   -0.694877601E+1     -0.480817011E+1    -0.411710659E+1     J kg-1 K-2
!(d2fV/dTddV)     0.262001885E+9      0.181489502E+5     0.361784086E+3     m3 kg-1 K-1
!(d2fV/(ddV)2)_T -0.297466671E+17    -0.210184992E+9    -0.965539462E+5     m3 kg-1 Pa-1

real*8  a, t, d, q(3)
integer i, j

a = errorreturn
t = errorreturn
d = errorreturn

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) " Implementation of IAPWS-10 in Fortran"
write(*,*) " for Publication in Ocean Science, 2009"
write(*,*) " D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) " NUMERICAL CHECK VALUES FOR THE DRY-AIR AND WATER VAPOUR PARTS OF THE"
write(*,*) " HELMHOLTZ FUNCTION OF SATURATED AIR, AND ITS DERIVATIVES"
write(*,*) " WRT ICE at 200 K AND WRT LIQUID WATER at 300 AND 400 K"
write(*,*) " FUNCTION VALUES AS GIVEN IN TABLE 14 OF IAPWS-10:"
write(*,*) " TOP ROW GIVES TABULATED VALUES, BOTTOM ROW COMPUTED NOW"
write(*,*) " "

write(*,*)               "Quantity         Value            Value            Value           Unit"
write(*,*) " "

do i = 1, 17
  select case(i)
    case(1);  write(*,*) "A                0.892247719      0.977605798      0.825565291     kg kg-1"
    case(2);  write(*,*) "T                200              300              400             K"
    case(3);  write(*,*) "d                0.163479657E-4   0.114614216E+1   0.793354063E+1  kg m-3"
    case(4);  write(*,*) "dA = d*A         0.145864351E-4   0.112047522E+1   0.654965578E+1  kg m-3"
    case(5);  write(*,*) "fA              -0.740041144E+6  -0.916103453E+5   0.895561286E+5  J kg-1"
    case(6);  write(*,*) "(dfA/dT)_dA     -0.304774177E+4  -0.108476220E+3   0.193271394E+3  J kg-1 K-1"
    case(7);  write(*,*) "(dfA/ddA)_T      0.393583654E+10  0.768326795E+5   0.175560114E+5  m3 kg-1"
    case(8);  write(*,*) "(d2fA/dT2)_dA   -0.357677878E+1  -0.239319940E+1  -0.181809877E+1  J kg-1 K-2"
    case(9);  write(*,*) "(d2fA/dTddA)     0.196791837E+8   0.256683306E+3   0.442769673E+2  m3 kg-1 K-1"
    case(10); write(*,*) "(d2fA/(ddA)2)_T -0.269828549E+15 -0.685917373E+5  -0.267635928E+4  m3 kg-1 Pa-1"
    case(11); write(*,*) "(1-A)d           0.176153059E-5   0.256669391E-1   0.138388485E+1  kg m-3"
    case(12); write(*,*) "fV              -0.202254351E+6  -0.143157426E+6  -0.285137534E+6  J kg-1"
    case(13); write(*,*) "(dfV/dT)_Ad     -0.123787544E+5  -0.851598213E+4  -0.705288048E+4  J kg-1 K-1"
    case(14); write(*,*) "(dfV/ddV)_T      0.523995674E+11  0.538480619E+7   0.129645039E+6  m3 kg-1"
    case(15); write(*,*) "(d2fV/dT2)_Ad   -0.694877601E+1  -0.480817011E+1  -0.411710659E+1  J kg-1 K-2"
    case(16); write(*,*) "(d2fV/dTddV)     0.262001885E+9   0.181489502E+5   0.361784086E+3  m3 kg-1 K-1"
    case(17); write(*,*) "(d2fV/(ddV)2)_T -0.297466671E+17 -0.210184992E+9  -0.965539462E+5  m3 kg-1 Pa-1"
  end select

  do j = 1, 3
    if(j == 1) then
       a = 0.892247719d0
       t = 200d0
       d = 0.163479657d-4
    elseif(j == 2) then
       a = 0.977605798d0
       t = 300d0
       d = 0.114614216d1
    elseif(j == 3) then
       a = 0.825565291d0
       t = 400d0
       d = 0.793354063d1
    endif
    q(j) = 0d0
    select case(i)
      case(1); q(j) = a
      case(2); q(j) = t
      case(3); q(j) = d
      case(4); q(j) = a * d
      case(5); q(j) = dry_f_si(0, 0, t, a*d)
      case(6); q(j) = dry_f_si(1, 0, t, a*d)
      case(7); q(j) = dry_f_si(0, 1, t, a*d)
      case(8); q(j) = dry_f_si(2, 0, t, a*d)
      case(9); q(j) = dry_f_si(1, 1, t, a*d)
      case(10); q(j) = dry_f_si(0, 2, t, a*d)
      case(11); q(j) = (1d0 - a)*d
      case(12); q(j) = flu_f_si(0, 0, t, (1d0 - a) * d)
      case(13); q(j) = flu_f_si(1, 0, t, (1d0 - a) * d)
      case(14); q(j) = flu_f_si(0, 1, t, (1d0 - a) * d)
      case(15); q(j) = flu_f_si(2, 0, t, (1d0 - a) * d)
      case(16); q(j) = flu_f_si(1, 1, t, (1d0 - a) * d)
      case(17); q(j) = flu_f_si(0, 2, t, (1d0 - a) * d)
    end select
  enddo
  write(*,'(17x, E16.9,x,E16.9,x,E16.9)') q(1), q(2), q(3)
  write(*,*)' '
  
enddo
write(*,*)' '
write(*,*)' '
write(*,*)' '

end subroutine


!=========================================================================
subroutine chk_iapws10_table15()
!=========================================================================

!TABLE 15  Numerical check values for the cross-virial part, fmix, and the cross-virial
!coefficients, BAW, CAAW, CAWW, of the Helmholtz function, fav, and its derivatives,
!Table 6, for saturated air, A = Asat, with respect to ice at 200 K, and with respect
!to liquid water at 300 and 400 K. The numerical functions evaluated here at given
!points (A, T, d) are defined in Tables 10 and 11.
 
!Quantity         Value               Value               Value              Unit
!A                0.892224944         0.977600624         0.825531379        kg kg-1
!T                200                 300                 400                K
!D                0.163445112E-04     0.114587678E+1      0.793198757E+1     kg m-3"
!fM              -0.786211837111E-3  -0.711673565972E+1	 -0.161985033872E+3  J kg-1
!(dfM/dA)_T,d     0.641377589024E-2   0.311768500757E+3	  0.830802876130E+3  J kg-1
!(dfM/dT)_A,d     0.456427011454E-5   0.441245367098E-1	  0.178961265299E+1  J kg-1 K-1
!(dfM/dd)_A,T    -0.481026562160E+2  -0.623171267739E+1	 -0.223365431713E+2  m3 kg-1
!(d2fM/dA2)_T,d   0.163518396765E-1   0.534139178420E+3   0.135815609516E+4  J kg-1
!(d2fM/dAdT)_d   -0.372355250910E-4  -0.195026097643E+1	 -0.916586082354E+1  J kg-1 K-1
!(d2fM/dAdd)_T    0.392414346187E+3   0.274152648677E+3	  0.125823777783E+3  m3 kg-1
!(d2fM/dT2)_A,d  -0.378866038500E-7  -0.148782305257E-3	 -0.536718535916E-2  J kg-1 K-2
!(d2fM/dTdd)_A    0.279261664337      0.390100461800E-1	  0.249618216264     m3 kg-1 K-1
!(d2fM/dd2)_A,T  -0.192118914304E+2  -0.366162709444E-1	 -0.482803925450     m3 kg-1 Pa-1
!BAW             -0.784874277752E-4  -0.295672747428E-4	 -0.100804610474E-4  m3 mol-1
!BAW_T            0.848076624222E-6   0.280097360438E-6	  0.135021228495E-6  m3 mol-1 K-1
!BAW_TT          -0.122622146106E-7  -0.242599241306E-8	 -0.839901728946E-9  m3 mol-1 K-2
!CAAW             0.105493575000E-8   0.801977740741E-9	  0.672018171875E-9  m6 mol-2
!CAAW_T          -0.152535000000E-11 -0.196103456790E-11 -0.812416406250E-12 m6 mol-2 K-1
!CAAW_TT         -0.113436375000E-12  0.170055637860E-13  0.683147460938E-14 m6 mol-2 K-2
!CAWW            -0.349872634207E-5  -0.115552783680E-6	 -0.200806020909E-7  m6 mol-2
!CAWW_T           0.188025052349E-6   0.261363277754E-8	  0.274535402840E-9  m6 mol-2 K-1
!CAWW_TT         -0.124996855887E-7  -0.751334581804E-10 -0.491763909891E-11 m6 mol-2 K-2

real*8  a, t, d, q(3)
integer i, j

                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) " Implementation of IAPWS-10 in Fortran"
write(*,*) " for Publication in Ocean Science, 2009"

write(*,*) " Function values as given in Table 15 of IAPWS-10:"
write(*,*) " TOP ROW GIVES TABULATED VALUES, BOTTOM ROW COMPUTED NOW"
write(*,*) " "

write(*,*)              "Quantity        Value            Value            Value            Unit"

do i = 1, 20
  select case(i)
    case(1);  write(*,*)"A               0.892224944      0.977600624      0.825531379      kg kg-1"
    case(2);  write(*,*)"T               200              300              400              K"
    case(3);  write(*,*)"d               0.163445112E-4   0.114587678E+1   0.793198757E+1   kg m-3"
    case(4);  write(*,*)"fM             -0.786211837E-3  -0.711673566E+1  -0.161985034E+3   J kg-1"
    case(5);  write(*,*)"(dfm/dA)_T,d    0.641377589E-2   0.311768501E+3   0.830802876E+3   J kg-1"
    case(6);  write(*,*)"(dfm/dT)_A,d    0.456427011E-5   0.441245368E-1   0.178961265E+1   J kg-1 K-1"
    case(7);  write(*,*)"(dfm/dd)_A,T   -0.481026562E+2  -0.623171268E+1  -0.223365432E+2   m3 kg-1"
    case(8);  write(*,*)"(d2fm/dA2)_T,d  0.163518397E-1   0.534139178E+3   0.135815610E+4   J kg-1"
    case(9);  write(*,*)"(d2fm/dAdT)_d  -0.372355251E-4  -0.195026098E+1  -0.916586082E+1   J kg-1 K-1"
    case(10); write(*,*)"(d2fm/dAdd)_T   0.392414346E+3   0.274152649E+3   0.125823778E+3   m3 kg-1"
    case(11); write(*,*)"(d2fm/dT2)_A,d -0.378866039E-7  -0.148782305E-3  -0.536718536E-2   J kg-1 K-2"
    case(12); write(*,*)"(d2fm/dTdd)_A   0.279261664      0.390100462E-1   0.249618216      m3 kg-1 K-1"
    case(13); write(*,*)"(d2fm/dd2)_A,T -0.192118914E+2  -0.366162709E-1  -0.482803925      m3 kg-1 Pa-1"
    case(14); write(*,*)"BAW            -0.784874278E-4  -0.295672747E-4  -0.100804610E-4   m3 mol-1"
    case(15); write(*,*)"BAW_T           0.848076624E-6   0.280097360E-6   0.135021228E-6   m3 mol-1 K-1"
    case(16); write(*,*)"BAW_TT         -0.122622146E-7  -0.242599241E-8  -0.839901729E-9   m3 mol-1 K-2"
    case(17); write(*,*)"CAAW            0.105493575E-8   0.801977741E-9   0.672018172E-9   m6 mol-2"
    case(18); write(*,*)"CAAW_T         -0.152535000E-11 -0.196103457E-11 -0.812416406E-12  m6 mol-2 K-1"
    case(19); write(*,*)"CAAW_TT        -0.113436375E-12  0.170055638E-13  0.683147461E-14  m6 mol-2 K-2"
    case(20); write(*,*)"CAWW           -0.349872634E-5  -0.115552784E-6  -0.200806021E-7   m6 mol-2"
    case(21); write(*,*)"CAWW_T          0.188025052E-6   0.261363278E-8   0.274535403E-9   m6 mol-2 K-1"
    case(22); write(*,*)"CAWW_TT        -0.124996856E-7  -0.751334582E-10 -0.491763910E-11  m6 mol-2 K-2"
  end select

  do j = 1, 3
    if(j == 1) then
       a = 0.892224944d0
       t = 200d0
       d = 0.0000163445112d0
    elseif(j == 2) then
       a = 0.977600624d0
       t = 300d0
       d = 1.14587678d0
    elseif(j == 3) then
       a = 0.825531379d0
       t = 400d0
       d = 7.93198757d0
    endif
    q(j) = 0d0
    select case(i)
      case(1); q(j) = a
      case(2); q(j) = t
      case(3); q(j) = d
      case(4); q(j) = air_f_mix_si(0, 0, 0, a, t, d)
      case(5); q(j) = air_f_mix_si(1, 0, 0, a, t, d)
      case(6); q(j) = air_f_mix_si(0, 1, 0, a, t, d)
      case(7); q(j) = air_f_mix_si(0, 0, 1, a, t, d)
      case(8); q(j) = air_f_mix_si(2, 0, 0, a, t, d)
      case(9); q(j) = air_f_mix_si(1, 1, 0, a, t, d)
      case(10); q(j) = air_f_mix_si(1, 0, 1, a, t, d)
      case(11); q(j) = air_f_mix_si(0, 2, 0, a, t, d)
      case(12); q(j) = air_f_mix_si(0, 1, 1, a, t, d)
      case(13); q(j) = air_f_mix_si(0, 0, 2, a, t, d)
      case(14); q(j) = air_baw_m3mol(0, t)
      case(15); q(j) = air_baw_m3mol(1, t)
      case(16); q(j) = air_baw_m3mol(2, t)
      case(17); q(j) = air_caaw_m6mol2(0, t)
      case(18); q(j) = air_caaw_m6mol2(1, t)
      case(19); q(j) = air_caaw_m6mol2(2, t)
      case(20); q(j) = air_caww_m6mol2(0, t)
      case(21); q(j) = air_caww_m6mol2(1, t)
      case(22); q(j) = air_caww_m6mol2(2, t)
    end select
  enddo
  write(*,'(16x, E16.9,x,E16.9,x,E16.9)') q(1), q(2), q(3)
  write(*,*)' '
  
enddo
write(*,*)' '
write(*,*)' '
write(*,*)' '

end subroutine

end module air_2
