module flu_1

!#########################################################################

!THIS MODULE IMPLEMENTS THE HELMHOLTZ POTENTIAL OF FLUID WATER AND ITS
!FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE AND
!DENSITY AS DEFINED IN IAPWS-95:

!RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC PROPERTIES OF
!ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!FREDERICIA, DENMARK, SEPTEMBER 1996,
!REVISED DOORWERTH, THE NETHERLANDS, SEPTEMBER 2009

!WITH AN EXTENSION TO 50 - 130 K OF THE IDEAL-GAS PART, AS DESCRIBED IN

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

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
!     CONSTANTS_0,    FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: chk_iapws95_table6, chk_iapws95_table7, flu_f_si

!ARRAYS DEFINED IN TABLE 1 OF IAPWS-95
!NUMERICAL VALUES OF THE COEFFICIENTS AND PARAMETERS OF THE IDEAL-GAS
!PART OF THE DIMENSIONLESS HELMHOLTZ FUNCTION AND ENERGY, EQ.(5)

real*8 n0(8)
real*8 gam0(8)

!ARRAYS DEFINED IN TABLE 2
!NUMERICAL VALUES OF THE COEFFICIENTS AND PARAMETERS OF THE
!RESIDUAL PART OF THE DIMENSIONLESS HELMHOLTZ FUNCTION AND ENERGY, EQ.(6)

real*8 ai(56), bi(56), ci(56), di(56), ni(56), ti(56)
real*8 aa(56), bb(56), cc(56), dd(56)
real*8 alf(56), bet(56), gam(56), eps(56)

!CP EXTENSION BELOW 130 K
real*8 ee

!CRITICAL POINT
real*8 :: rhoc = cp_density_si    !kg/m3
real*8 :: tc = cp_temperature_si  !K

Contains

!==========================================================================
function flu_f_si(drv_t, drv_d, t_si, d_si)
!==========================================================================
                         
!RETURN VALUE:
!FLU_F_SI:          DERIVATIVE OF THE HELMHOLTZ FUNCTION IN THE BASIC SI UNIT

!INPUT PARAMETERS:
!DRV_T:             ORDER OF THE PARTIAL TEMPERATURE DERIVATIVE
!DRV_D:             ORDER OF THE PARTIAL DENSITY DERIVATIVE
!T_SI:              ABSOLUTE TEMPERATURE ITS-90 IN K
!D_SI:              DENSITY IN KG/M3

!CHECK VALUES:
!flu_f_si(0,0,300,1000) = -5351.74115204
!flu_f_si(1,0,300,1000) = -390.904170767
!flu_f_si(0,1,300,1000) = 7.8330013559|7
!flu_f_si(2,0,300,1000) = -13.6840204925
!flu_f_si(1,1,300,1000) = 0.639359046588
!flu_f_si(0,2,300,1000) = 2.24824656167

implicit none

integer drv_t, drv_d
real*8 flu_f_si, t_si, d_si
real*8 del, tau, RT, f, ft, fd, fdt, ftt, r

call init_iapws95

!CONST R = GAS_CONSTANT_H2O_SI  !SPECIFIC GAS CONSTANT OF H2O IN J/(KG K)
r = gas_constant_h2o_iapws95    !VALUE USED IN THE IAPWS-95 RELEASE

flu_f_si = errorreturn

!exclude illegal calls:
if(check_limits == 1) then
  if(t_si < flu_tmin .or. t_si > flu_tmax) return
  if(d_si <= flu_dmin .or. d_si > flu_dmax) return
else
  if(t_si <= 0d0 .or. &
     d_si <= 0d0) return
endif

!reduce input values:
tau = tc / t_si
del = d_si / rhoc

RT = R * t_si

Call Init_IAPWS95

!EXPRESSIONS FROM TABLE 4 OF THE IAPWS ADVISORY NOTE #3, REVISED, BERLIN 2008:
select case(drv_t)

  case (0)
    select case (drv_d)
      case (0); flu_f_si = RT * (phi0_(del, tau) + phir_(del, tau))
      case (1); flu_f_si = RT / rhoc * (phi0_d(del) + phir_d(del, tau))
      case (2); flu_f_si = RT / rhoc ** 2 * (phi0_dd(del) + phir_dd(del, tau))
      case default; return
    end select
    
  case (1)
    select case (drv_d)
      case (0); f = phi0_(del, tau) + phir_(del, tau)
                ft = phi0_t(tau) + phir_t(del, tau)
                flu_f_si = R * (f - tau * ft)
      case (1); fd = phi0_d(del) + phir_d(del, tau)
                fdt = phi0_dt() + phir_dt(del, tau)
                flu_f_si = R / rhoc * (fd - tau * fdt)
      case default; return
    end select

  case (2)
    select case (drv_d)
      case (0); ftt = phi0_tt(tau) + phir_tt(del, tau)
                flu_f_si = R / t_si * tau ** 2 * ftt
      case default; return
    end select

  case default; return

end select

end function

!==========================================================================
function delta_(i, del, tau)
!==========================================================================

!FUNCTION DELTA AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 delta_, del, tau 

!AWAY FROM AND AT THE CRITICAL DENSITY
delta_ = theta_(i, del, tau) ** 2 + bb(i) * abs(del - 1d0) ** (2d0 * ai(i))

end function

!==========================================================================
function delta_d(i, del, tau)
!==========================================================================

!FUNCTION DELTA_D AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 delta_d, del, tau

if(del == 1d0) then  !at the critical density
  delta_d = 0d0
else
  delta_d = (del - 1d0) * (aa(i) * theta_(i, del, tau) * 2d0 / bet(i) &
            * abs(del - 1d0) ** (1d0 / bet(i) - 2d0) &
            + 2d0 * bb(i) * ai(i) * abs(del - 1d0) ** (2d0 * ai(i) - 2d0))
end if

end function

!==========================================================================
function delta_dd(i, del, tau)
!==========================================================================

!FUNCTION DELTA_DD AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 delta_dd, del, tau

if(del == 1d0) then  !at the critical density
  delta_dd = 0
else
  delta_dd = delta_d(i, del, tau) / (del - 1d0) + &
             4d0 * bb(i) * ai(i) * (ai(i) - 1d0) * abs(del - 1d0) ** (2d0 * ai(i) - 2d0) + &
             2d0 * (aa(i) / bet(i)) ** 2 * abs(del - 1d0) ** (2d0 / bet(i) - 2d0) + &
             2d0 * theta_(i, del, tau) * aa(i) / bet(i) * (1d0 / bet(i) - 2d0) * &
             abs(del - 1d0) ** (1d0 / bet(i) - 2d0)
end if

end function

!==========================================================================
function deltab_d(i, del, tau)
!==========================================================================

!FUNCTION DELTA^BI_D AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 deltab_d, del, tau, delta

delta = delta_(i, del, tau)

if(delta == 0d0) then
  deltab_d = 0d0
else
  deltab_d = bi(i) * delta ** (bi(i) - 1d0) * delta_d(i, del, tau)
end if

end function

!==========================================================================
function deltab_dt(i, del, tau)
!==========================================================================

!FUNCTION DELTA^BI_DT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 deltab_dt, del, tau, delta

delta = delta_(i, del, tau)

if(delta == 0d0) then
  deltab_dt = 0d0
else
  deltab_dt = -2d0 * bi(i) * delta ** (bi(i) - 2d0) * &
              (aa(i) / bet(i) * delta * (del - 1d0) * &
              abs(del - 1d0) ** (1d0 / bet(i) - 2d0) + &
              theta_(i, del, tau) * (bi(i) - 1d0) * delta_d(i, del, tau))
end if

end function

!==========================================================================
function deltab_t(i, del, tau)
!==========================================================================

!FUNCTION DELTA^BI_T AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 deltab_t, del, tau, delta

delta = delta_(i, del, tau)
if(delta == 0d0) then
  deltab_t = 0d0
else
  deltab_t = -2d0 * theta_(i, del, tau) * bi(i) * delta ** (bi(i) - 1d0)
end if

end function

!==========================================================================
function deltab_tt(i, del, tau)
!==========================================================================

!FUNCTION DELTA^BI_TT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 deltab_tt, del, tau, delta

delta = delta_(i, del, tau)
if(delta == 0d0) then
  deltab_tt = 0d0
else
  deltab_tt = 2d0 * bi(i) * delta ** (bi(i) - 2d0) * &
              (delta + 2d0 * theta_(i, del, tau) ** 2 * (bi(i) - 1d0))
end if
end function

!==========================================================================
function deltab_dd(i, del, tau)
!==========================================================================

!FUNCTION DELTA^BI_DD AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 deltab_dd, del, tau, delta

delta = delta_(i, del, tau)
if(delta == 0d0) then
  deltab_dd = 0d0
else
  deltab_dd = bi(i) * (delta ** (bi(i) - 1d0) * delta_dd(i, del, tau) + &
            (bi(i) - 1d0) * delta ** (bi(i) - 2d0) * delta_d(i, del, tau) ** 2)
end if

end function

!==========================================================================
function theta_(i, del, tau)
!==========================================================================

!FUNCTION THETA AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 theta_, del, tau

if(del == 1d0) then  !at the critical density
  theta_ = 1d0 - tau
else
  theta_ = 1d0 - tau + aa(i) * abs(del - 1d0) ** (1d0 / bet(i))
end if

end function

!==========================================================================
function psi_(i, del, tau)
!==========================================================================

!FUNCTION PSI AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_, del, tau

psi_ = exp(-cc(i) * (del - 1d0) ** 2 - dd(i) * (tau - 1d0) ** 2)

end function

!==========================================================================
function psi_t(i, del, tau)
!==========================================================================

!FUNCTION PSI_T AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_t, del, tau

psi_t = -2d0 * dd(i) * (tau - 1d0) * psi_(i, del, tau)

end function

!==========================================================================
function psi_tt(i, del, tau)
!==========================================================================

!FUNCTION PSI_TT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_tt, del, tau

psi_tt = 2d0 * dd(i) * (2d0 * dd(i) * (tau - 1d0) ** 2 - 1d0) * psi_(i, del, tau)

end function

!==========================================================================
function psi_d(i, del, tau)
!==========================================================================

!FUNCTION PSI_D AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_d, del, tau

psi_d = -2d0 * cc(i) * (del - 1d0) * psi_(i, del, tau)

end function

!==========================================================================
function psi_dt(i, del, tau)
!==========================================================================

!FUNCTION PSI_DT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_dt, del, tau

psi_dt = 4d0 * cc(i) * dd(i) * (del - 1d0) * (tau - 1d0) * psi_(i, del, tau)

end function

!==========================================================================
function psi_dd(i, del, tau)
!==========================================================================

!FUNCTION PSI_DD AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 psi_dd, del, tau

psi_dd = 2d0 * cc(i) * (2d0 * cc(i) * (del - 1d0) ** 2 - 1d0) * psi_(i, del, tau)

end function

!==========================================================================
function phi0_(del, tau)
!==========================================================================

!FUNCTION PHI0 AS DEFINED IN ROW #1 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
integer i
real*8 phi0_, del, tau, e, phi

e = tc / 130d0

phi = log(del) + n0(1) + n0(2) * tau + n0(3) * log(tau)
Do i = 4, 8
  phi = phi + n0(i) * log(1d0 - exp(-gam0(i) * tau))
enddo

if(tau > e) then  !extension below 130 K
  phi = phi + ee * (-0.5d0 / tau - 3d0 / e ** 2 * (tau + e) * log(tau / e) &
        + 4.5d0 * tau / e ** 2 + 0.5d0 * tau ** 2 / e ** 3 - 4.5d0 / e)
end if

phi0_ = phi
 
end function

!==========================================================================
function phi0_t(tau)
!==========================================================================

!FUNCTION PHI0_T AS DEFINED IN ROW #4 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
integer i
real*8 phi0_t, tau, e, phi

e = tc / 130d0

phi = n0(2) + n0(3) / tau
Do i = 4, 8
  phi = phi + n0(i) * gam0(i) * (1d0 / (1d0 - exp(-gam0(i) * tau)) - 1d0)
enddo

if(tau > e) then  !extension below 130 K
  phi = phi + ee * (0.5d0 / tau ** 2 - 3d0 / e ** 2 * log(tau / e) &
        - 3d0 / (e * tau) + tau / e ** 3 + 1.5d0 / e ** 2)
end if

phi0_t = phi

end function

!==========================================================================
function phi0_tt(tau)
!==========================================================================

!FUNCTION PHI0_TT AS DEFINED IN ROW #5 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
integer i
real*8 phi0_tt, tau, e, phi

e = tc / 130d0

phi = -n0(3) / tau ** 2
Do i = 4, 8
  phi = phi - n0(i) * gam0(i) ** 2 * exp(-gam0(i) * tau) / (1d0 - exp(-gam0(i) * tau)) ** 2
enddo

if(tau > e) then  !extension below 130 K
  phi = phi + ee * (-1d0 / tau + 1d0 / e) ** 3
end if

phi0_tt = phi

end function

!==========================================================================
function phi0_d(del)
!==========================================================================

!FUNCTION PHI0_D AS DEFINED IN ROW #2 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
real*8 phi0_d, del

phi0_d = 1d0 / del

end function

!==========================================================================
function phi0_dd(del)
!==========================================================================

!FUNCTION PHI0_DD AS DEFINED IN ROW #3 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
real*8 phi0_dd, del

phi0_dd = -1d0 / del ** 2

end function

!==========================================================================
function phi0_dt()
!==========================================================================

!FUNCTION PHI0_DT AS DEFINED IN ROW #6 OF TABLE 4 OF IAPWS-95:
!THE IDEAL-GAS PART PHI0 OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
real*8 phi0_dt

phi0_dt = 0d0

end function

!==========================================================================
function phir_(del, tau)
!==========================================================================

!FUNCTION PHIR AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 phir_, del, tau, phi

phi = 0d0

Do i = 1, 7
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i)
enddo

Do i = 8, 51
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * exp(-del ** ci(i))
enddo

Do i = 52, 54
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2)
enddo

Do i = 55, 56
  phi = phi + ni(i) * delta_(i, del, tau) ** bi(i) * del * psi_(i, del, tau)
enddo

phir_ = phi

end function

!==========================================================================
function phir_d(del, tau)
!==========================================================================

!FUNCTION PHIR_D AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 phir_d, del, tau, phi, psi

phi = 0d0

Do i = 1, 7
  phi = phi + ni(i) * di(i) * del ** (di(i) - 1d0) * tau ** ti(i)
enddo

Do i = 8, 51
  phi = phi + ni(i) * exp(-del ** ci(i)) * (del ** (di(i) - 1d0) * tau ** ti(i) &
        * (di(i) - ci(i) * del ** ci(i)))
enddo

Do i = 52, 54
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2) * &
        (di(i) / del - 2d0 * alf(i) * (del - eps(i)))
enddo

Do i = 55, 56
  psi = psi_(i, del, tau)
  phi = phi + ni(i) * (delta_(i, del, tau) ** bi(i) * &
        (psi + del * psi_d(i, del, tau)) + deltab_d(i, del, tau) * del * psi)
enddo

phir_d = phi

end function

!==========================================================================
function phir_dd(del, tau)
!==========================================================================

!FUNCTION PHIR_DD AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 phir_dd, del, tau, phi, psi, psid, tmp

phi = 0d0

Do i = 1, 7
  phi = phi + ni(i) * di(i) * (di(i) - 1d0) * del ** (di(i) - 2d0) * tau ** ti(i)
enddo

Do i = 8, 51
  tmp = (del ** (di(i) - 2d0) * tau ** ti(i) * &
        ((di(i) - ci(i) * del ** ci(i)) * (di(i) - 1 - ci(i) * del ** ci(i)) &
        - ci(i) ** 2 * del ** ci(i)))
  phi = phi + ni(i) * exp(-del ** ci(i)) * tmp
enddo

Do i = 52, 54
  phi = phi + ni(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2) * &
        (-2d0 * alf(i) * del ** di(i) + 4d0 * alf(i) ** 2 * del ** di(i) &
        * (del - eps(i)) ** 2 - 4d0 * di(i) * alf(i) * del ** (di(i) - 1d0) &
        * (del - eps(i)) + di(i) * (di(i) - 1d0) * del ** (di(i) - 2d0))
enddo

Do i = 55, 56
  psi = psi_(i, del, tau)
  psid = psi_d(i, del, tau)
  phi = phi + ni(i) * (delta_(i, del, tau) ** bi(i) * (2d0 * psid + del * &
        psi_dd(i, del, tau)) + 2d0 * deltab_d(i, del, tau) * (psi + del * psid) + &
        deltab_dd(i, del, tau) * del * psi)
enddo

phir_dd = phi

end function

!==========================================================================
function phir_t(del, tau)
!==========================================================================

!FUNCTION PHIR_T AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES

integer i
real*8 phir_t, del, tau, phi
logical fastmode

!Jackett code variables
real*8 delta, delta2, delta3, delta4, delta6, delta13, delta14
real*8 delta1_2,  delta_tau, delta_tau2, delta_tau3, delta_tau4, delta2_tau
real*8 tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11, tau13, tau12 
real*8 tau15, tau16, tau21, tau22, tau23, tau43, tau44, tau49, tau50
real*8 taup375, taup5, taup75, taup875
real*8 tau1_2, tau2_2, tau3_2, tau_1
real*8 cd1, cd2, cd3, cd4, cd6, cdt1, cdt2, cdt3, cdt4
real*8 edt1, edt2, edt3, edt4, the_rest

fastmode = .False.
!flu_entropy_si(647, 358) = 4320.923066755d0  correct

!fastmode = .True.
!flu_entropy_si(647, 358) = 4320.92306292204d0  less correct

if(fastmode) then 
   GoTo 100  !JackettCode
endif

phi = 0d0

Do i = 1, 7
  phi = phi + ni(i) * ti(i) * del ** di(i) * tau ** (ti(i) - 1d0)
enddo

Do i = 8, 51
  phi = phi + ni(i) * ti(i) * del ** di(i) * tau ** (ti(i) - 1d0) * exp(-del ** ci(i))
enddo

Do i = 52, 54
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2) * &
        (ti(i) / tau - 2d0 * bet(i) * (tau - gam(i)))
enddo

Do i = 55, 56
  phi = phi + ni(i) * del * (deltab_t(i, del, tau) * psi_(i, del, tau) + &
        delta_(i, del, tau) ** bi(i) * psi_t(i, del, tau))
enddo

phir_t = phi

return


!======================================================================================
!fast version adopted from David Jackett!s Fortran code
!======================================================================================

100 Continue  !JackettCode:

!!    FIRST DERIVATIVE (WITH RESPECT TO TAU) OF THE RESIDUAL PART OF THE
!!    DIMENSIONLESS HELMHOLTZ FREE ENERGY OF FLUID WATER, AS IN
!!
!!    THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!!    FREDERICIA, DENMARK, SEPTEMBER 1996
!!
!!    DELTA           : DIMENSIONLESS DENSITY
!!    TAU             : DIMENSIONLESS TEMPERATURE
!
!!    DPHIRDTAU       : FIRST DERIVATIVE OF PHIR (WRT TAU)
!!
!!
!!    DAVID JACKETT   : 10/04/08

delta = del
delta2 = delta * delta; delta3 = delta * delta2; delta4 = delta * delta3; delta6 = delta3 * delta3
delta13 = delta3 * delta4 * delta6; delta14 = delta * delta13

delta1_2 = (-1d0 + delta) * (-1d0 + delta)

tau2 = tau * tau;   tau3 = tau * tau2; tau4 = tau * tau3; tau5 = tau * tau4; tau6 = tau * tau5
tau7 = tau * tau6;   tau8 = tau * tau7; tau9 = tau * tau8; tau10 = tau * tau9; tau11 = tau2 * tau9
tau12 = tau * tau11; tau13 = tau4 * tau9; tau15 = tau4 * tau11
tau16 = tau8 * tau8; tau21 = tau5 * tau16; tau22 = tau9 * tau13; tau23 = tau * tau22
tau43 = tau21 * tau22; tau44 = tau22 * tau22; tau49 = tau5 * tau44; tau50 = tau44 * tau6

taup375 = tau ** 0.375d0; taup5 = Sqrt(tau); taup75 = taup375 * taup375; taup875 = taup375 * taup5

tau1_2 = (-1.21d0 + tau) * (-1.21d0 + tau); tau2_2 = (-1.25d0 + tau) * (-1.25d0 + tau)
tau3_2 = (-1d0 + tau) * (-1d0 + tau)

cd1 = delta * (1.26917960652461d-08 * tau10 + &
      delta * (1.46328660576816d-06 * tau3 + &
      delta2 * (-6.2639586912454d-10 * delta2 - 1.72265340970684d-11 * tau12)))
cd1 = delta * (-6.83472393793296d-05 * tau8 + &
      delta2 * (1.68752938055664d-03 * tau2 + &
      delta2 * (-6.243460902854d-05 * tau3 + cd1)))
cd1 = delta * (0.64299473945004d0 * tau3 + delta * &
      (tau * (-0.080185657851614d0 + 5.11464493842302d-06 * tau11) + cd1))
cd1 = delta * (tau3 * (-2.6742628923186d0 + tau2 * (1.2260286570579d0 - 7.94551260476244d-04 * tau6)) + &
      delta * (-0.19232721156002d0 - 1.2854521501719d0 * tau4 + cd1))

cd2 = delta * (-1.3110854674024d-04 * delta2 * tau7 + &
      tau5 * (1.19733431877246d-02 + 1.42832774917413d-03 * tau3))
cd2 = cd2 + tau * (-0.05810467201917d0 + &
      tau * (tau * (-0.081573946054816d0 - 1.32432400509872d-02 * tau4) + 0.115845256722618d0))
cd2 = delta * (tau5 * (0.123167645375688d0 - 1.3636435110343d0 * tau4) + &
      delta * (0.14180634400617d0 * tau9 + &
      delta2 * (8.3326504880713d-03 + cd2)))
cd2 = (tau2 * (1.49907440972418d-02 + tau4 * (-0.219510904987843d0 - 7.4315929710341d0 * tau3)) + &
      delta * (4.780732991548d0 * tau9 + cd2))
cd2 = delta * (-0.75555206362524d0 * tau6 + delta * (0.017611491008752d0 + &
      (1.99190656507914d0 - 4.0247669763528d0 * tau) * tau8 + &
      delta * (5.8083399985759d0 * tau9 + delta * cd2)))

cd3 = delta3 * (0.697817851580976d0 * tau15 + delta * &
      ((0.76986812020283d0 - 1.76612855042628d0 * tau) * tau21 + 0.516264378636138d0 * delta * tau22))
cd4 = -6.2689710414685d-04 * delta14 * tau9
cd6 = delta3 * (-2.78555592828225d-08 * tau49 + &
      delta3 * tau43 * (-8.75851607593952d0 + tau2 * (14.6176487721395d0 - 5.9205912129905d0 * tau4)))

cdt1 = delta3 * (-11332.6263571691d0 + (20843.127003345d0 - 9463.8420713343d0 * tau) * tau)
edt1 = exp(-20d0 * delta1_2 - 150d0 * tau1_2)
cdt2 = delta3 * tau3 * (-10085.261736678d0 + tau * (-1575822.14635594d0 + 1260657.71708475d0 * tau))
edt2 = exp(-20d0 * delta1_2 - 250d0 * tau2_2)

delta_tau = (1d0 + 0.32d0 * delta1_2 ** 1.66666666666667d0 - tau); delta_tau2 = delta_tau * delta_tau 
delta2_tau = delta * delta_tau
delta_tau3 = (0.2d0 * delta1_2 ** 3.5d0 + delta_tau2); delta_tau4 = delta * tau_1

if(delta_tau3 == 0d0) then
  cdt3 = 0
else
  cdt3 = delta_tau3 ** 0.85d0 * ((0.252868894564308d0 * delta2_tau) / delta_tau3 + 208.244971994136d0 * delta_tau4)
end if

edt3 = exp(-28d0 * delta1_2 - 700d0 * tau3_2)

if(delta_tau3 == 0d0) then
  cdt4 = 0
else
  cdt4 = delta_tau3 ** 0.95d0 * ((-0.604316106690436d0 * delta2_tau) / delta_tau3 - 508.897774055104d0 * delta_tau4)
end if

edt4 = exp(-32d0 * delta1_2 - 800d0 * tau3_2)
        
the_rest = -8.7803203303561d0 * delta + 8.8089493102134d-03 * delta4 &
         - (6.2667739677615d-03 * delta * taup5) / tau2 + &
         (-2.93249068829929d-03 * delta3 * taup375 + delta2 * (0.15901254672709d0 * taup5 - &
         0.196091503945185d0 * taup75) + 6.90879303824745d0 * delta * taup875) / tau
      
phir_t = cd1 * exp(-delta) + cd2 * exp(-delta2) + cd3 * exp(-delta3) + cd4 *  &
         exp(-delta4) + cd6 * exp(-delta6) + cdt1 * edt1 + cdt2 * edt2 + cdt3 *  &
         edt3 + cdt4 * edt4 + the_rest

end function

!==========================================================================
function phir_dt(del, tau)
!==========================================================================

!FUNCTION PHIR_DT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
integer i
real*8 phir_dt, del, tau, phi, psi, psit

phi = 0

Do i = 1, 7
  phi = phi + ni(i) * di(i) * ti(i) * del ** (di(i) - 1) * tau ** (ti(i) - 1)
enddo

Do i = 8, 51
  phi = phi + ni(i) * ti(i) * exp(-del ** ci(i)) * del ** (di(i) - 1) * tau ** (ti(i) - 1) * &
        (di(i) - ci(i) * del ** ci(i))
enddo

Do i = 52, 54
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2) * &
        (di(i) / del - 2d0 * alf(i) * (del - eps(i))) * &
        (ti(i) / tau - 2 * bet(i) * (tau - gam(i)))
enddo

Do i = 55, 56
  psi = psi_(i, del, tau)
  psit = psi_t(i, del, tau)
  phi = phi + ni(i) * (delta_(i, del, tau) ** bi(i) * (psit + del * psi_dt(i, del, tau)) + &
        del * deltab_d(i, del, tau) * psit + &
        deltab_t(i, del, tau) * (psi + del * psi_d(i, del, tau)) + &
        deltab_dt(i, del, tau) * del * psi)
enddo

phir_dt = phi

end function

!==========================================================================
function phir_tt(del, tau)
!==========================================================================

!FUNCTION PHIR_TT AS DEFINED IN TABLE 5 OF IAPWS-95:
!THE RESIDUAL PART PHIR OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY AND ITS DERIVATIVES
 
integer i
real*8 phir_tt, del, tau, phi

phi = 0

Do i = 1, 7
  phi = phi + ni(i) * ti(i) * (ti(i) - 1) * del ** di(i) * tau ** (ti(i) - 2)
enddo

Do i = 8, 51
  phi = phi + ni(i) * ti(i) * (ti(i) - 1) * del ** di(i) * tau ** (ti(i) - 2) * exp(-del ** ci(i))
enddo

Do i = 52, 54
  phi = phi + ni(i) * del ** di(i) * tau ** ti(i) * &
        exp(-alf(i) * (del - eps(i)) ** 2 - bet(i) * (tau - gam(i)) ** 2) * &
        ((ti(i) / tau - 2 * bet(i) * (tau - gam(i))) ** 2 - ti(i) / tau ** 2 - 2 * bet(i))
enddo

Do i = 55, 56
  phi = phi + ni(i) * del * (deltab_tt(i, del, tau) * psi_(i, del, tau) + &
        2d0 * deltab_t(i, del, tau) * psi_t(i, del, tau) + &
        delta_(i, del, tau) ** bi(i) * psi_tt(i, del, tau))
enddo

phir_tt = phi

end function

!==========================================================================
function chk_flu_phi(drv_t, drv_d, t_si, d_si, key)
!==========================================================================

integer drv_t, drv_d 
real*8 chk_flu_phi, t_si, d_si
character (*) key

call init_iapws95

if(key == "0") then
   chk_flu_phi = chk_flu_phi0(drv_t, drv_d, t_si, d_si)
elseif(key == "r") then
   chk_flu_phi = chk_flu_phir(drv_t, drv_d, t_si, d_si)
else
   write(*,*)'invalid key specified for chk_flu'
   chk_flu_phi = errorreturn
endif

end function

!==========================================================================
function chk_flu_phi0(drv_t, drv_d, t_si, d_si)
!==========================================================================

integer drv_t, drv_d 
real*8 chk_flu_phi0, t_si, d_si, del, tau

tau = tc / t_si
del = d_si / rhoc

chk_flu_phi0 = errorreturn

select case (drv_t)

  case (0)
    select case (drv_d)
      case (0); chk_flu_phi0 = phi0_(del, tau)
      case (1); chk_flu_phi0 = phi0_d(del)
      case (2); chk_flu_phi0 = phi0_dd(del)
      case default; return
    end select
    
  case (1)
    select case (drv_d)
      case (0); chk_flu_phi0 = phi0_t(tau)
      case (1); chk_flu_phi0 = phi0_dt()
      case default; return
    end select

  case (2)
    select case (drv_d)
      case (0); chk_flu_phi0 = phi0_tt(tau)
      case default; return
    end select

  case default; return
end select

end function

!==========================================================================
function chk_flu_phir(drv_t, drv_d, t_si, d_si)
!==========================================================================

integer drv_t, drv_d
real*8 chk_flu_phir, t_si, d_si, del, tau

tau = tc / t_si
del = d_si / rhoc

chk_flu_phir = errorreturn

select case (drv_t)

  case (0)
    select case (drv_d)
      case (0); chk_flu_phir = phir_(del, tau)
      case (1); chk_flu_phir = phir_d(del, tau)
      case (2); chk_flu_phir = phir_dd(del, tau)
      case default; return
    end select
    
  case (1)
    select case (drv_d)
      case (0); chk_flu_phir = phir_t(del, tau)
      case (1); chk_flu_phir = phir_dt(del, tau)
      case default; return
    end select

  case (2)
    select case (drv_d)
      case (0); chk_flu_phir = phir_tt(del, tau)
      case default; return
    end select

  case default; return

end select

end function

!==========================================================================
subroutine Init_IAPWS95()
!==========================================================================

integer i

if(di(1) == 1d0) return

!COEFFICIENTS DEFINED IN TABLE 1 OF IAPWS-95
!NUMERICAL VALUES OF THE COEFFICIENTS AND PARAMETERS OF THE IDEAL-GAS
!PART OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY, EQ.(5)
!-------------------------------------------------------------------------------

!ADJUSTMENT TO THE IAPWS-95 REFERENCE STATE CONDITIONS AT THE TRIPLE POINT
!ORIGINAL VALUES IN THE IAPWS-95 RELEASE
!I = 1; N0(I) = -8.32044648201D0
!I = 2; N0(I) = 6.6832105268D0

!IMPROVED VALUES OF THE ADJUSTABLE COEFFICIENTS OBTAINED FROM A QUADRUPLE PRECISION IMPLEMENTATION
i = 1; n0(I) = -8.32044648374969D0
i = 2; n0(I) = 6.68321052759323D0

i = 3; n0(I) = 3.00632D0
i = 4; n0(I) = 0.012436D0; gam0(I) = 1.28728967D0
i = 5; n0(I) = 0.97315D0;  gam0(I) = 3.53734222D0
i = 6; n0(I) = 1.2795D0;   gam0(I) = 7.74073708D0
i = 7; n0(I) = 0.96956D0;  gam0(I) = 9.24437796D0
i = 8; n0(I) = 0.24873D0;  gam0(I) = 27.5075105D0
!-------------------------------------------------------------------------------

!COEFFICIENTS DEFINED IN TABLE 2 OF IAPWS-95
!NUMERICAL VALUES OF THE COEFFICIENTS AND PARAMETERS OF THE
!RESIDUAL PART OF THE DIMENSIONLESS HELMHOLTZ FREE ENERGY, EQ.(6)
!-------------------------------------------------------------------------------
i = 1;             di(i) = 1;  ti(i) = -0.5d0;  ni(i) = 0.012533547935523d0
i = 2;             di(i) = 1;  ti(i) = 0.875d0; ni(i) = 7.8957634722828d0
i = 3;             di(i) = 1;  ti(i) = 1d0;     ni(i) = -8.7803203303561d0
i = 4;             di(i) = 2;  ti(i) = 0.5d0;   ni(i) = 0.31802509345418d0
i = 5;             di(i) = 2;  ti(i) = 0.75d0;  ni(i) = -0.26145533859358d0
i = 6;             di(i) = 3;  ti(i) = 0.375d0; ni(i) = -7.8199751687981d-03
i = 7;             di(i) = 4;  ti(i) = 1d0;     ni(i) = 8.8089493102134d-03
i = 8;  ci(i) = 1; di(i) = 1;  ti(i) = 4d0;     ni(i) = -0.66856572307965d0
i = 9;  ci(i) = 1; di(i) = 1;  ti(i) = 6d0;     ni(i) = 0.20433810950965d0
i = 10; ci(i) = 1; di(i) = 1;  ti(i) = 12d0;    ni(i) = -6.6212605039687d-05
i = 11; ci(i) = 1; di(i) = 2;  ti(i) = 1d0;     ni(i) = -0.19232721156002d0
i = 12; ci(i) = 1; di(i) = 2;  ti(i) = 5d0;     ni(i) = -0.25709043003438d0
i = 13; ci(i) = 1; di(i) = 3;  ti(i) = 4d0;     ni(i) = 0.16074868486251d0
i = 14; ci(i) = 1; di(i) = 4;  ti(i) = 2d0;     ni(i) = -0.040092828925807d0
i = 15; ci(i) = 1; di(i) = 4;  ti(i) = 13d0;    ni(i) = 3.9343422603254d-07
i = 16; ci(i) = 1; di(i) = 5;  ti(i) = 9d0;     ni(i) = -7.5941377088144d-06
i = 17; ci(i) = 1; di(i) = 7;  ti(i) = 3d0;     ni(i) = 5.6250979351888d-04
i = 18; ci(i) = 1; di(i) = 9;  ti(i) = 4d0;     ni(i) = -1.5608652257135d-05
i = 19; ci(i) = 1; di(i) = 10; ti(i) = 11d0;    ni(i) = 1.1537996422951d-09
i = 20; ci(i) = 1; di(i) = 11; ti(i) = 4d0;     ni(i) = 3.6582165144204d-07
i = 21; ci(i) = 1; di(i) = 13; ti(i) = 13d0;    ni(i) = -1.3251180074668d-12
i = 22; ci(i) = 1; di(i) = 15; ti(i) = 1d0;     ni(i) = -6.2639586912454d-10
i = 23; ci(i) = 2; di(i) = 1;  ti(i) = 7d0;     ni(i) = -0.10793600908932d0
i = 24; ci(i) = 2; di(i) = 2;  ti(i) = 1d0;     ni(i) = 0.017611491008752d0
i = 25; ci(i) = 2; di(i) = 2;  ti(i) = 9d0;     ni(i) = 0.22132295167546d0
i = 26; ci(i) = 2; di(i) = 2;  ti(i) = 10d0;    ni(i) = -0.40247669763528d0
i = 27; ci(i) = 2; di(i) = 3;  ti(i) = 10d0;    ni(i) = 0.58083399985759d0
i = 28; ci(i) = 2; di(i) = 4;  ti(i) = 3d0;     ni(i) = 4.9969146990806d-03
i = 29; ci(i) = 2; di(i) = 4;  ti(i) = 7d0;     ni(i) = -0.031358700712549d0
i = 30; ci(i) = 2; di(i) = 4;  ti(i) = 10d0;    ni(i) = -0.74315929710341d0
i = 31; ci(i) = 2; di(i) = 5;  ti(i) = 10d0;    ni(i) = 0.4780732991548d0
i = 32; ci(i) = 2; di(i) = 6;  ti(i) = 6d0;     ni(i) = 0.020527940895948d0
i = 33; ci(i) = 2; di(i) = 6;  ti(i) = 10d0;    ni(i) = -0.13636435110343d0
i = 34; ci(i) = 2; di(i) = 7;  ti(i) = 10d0;    ni(i) = 0.014180634400617d0
i = 35; ci(i) = 2; di(i) = 9;  ti(i) = 1d0;     ni(i) = 8.3326504880713d-03
i = 36; ci(i) = 2; di(i) = 9;  ti(i) = 2d0;     ni(i) = -0.029052336009585d0
i = 37; ci(i) = 2; di(i) = 9;  ti(i) = 3d0;     ni(i) = 0.038615085574206d0
i = 38; ci(i) = 2; di(i) = 9;  ti(i) = 4d0;     ni(i) = -0.020393486513704d0
i = 39; ci(i) = 2; di(i) = 9;  ti(i) = 8d0;     ni(i) = -1.6554050063734d-03
i = 40; ci(i) = 2; di(i) = 10; ti(i) = 6d0;     ni(i) = 1.9955571979541d-03
i = 41; ci(i) = 2; di(i) = 10; ti(i) = 9d0;     ni(i) = 1.5870308324157d-04
i = 42; ci(i) = 2; di(i) = 12; ti(i) = 8d0;     ni(i) = -1.638856834253d-05
i = 43; ci(i) = 3; di(i) = 3;  ti(i) = 16d0;    ni(i) = 0.043613615723811d0
i = 44; ci(i) = 3; di(i) = 4;  ti(i) = 22d0;    ni(i) = 0.034994005463765d0
i = 45; ci(i) = 3; di(i) = 4;  ti(i) = 23d0;    ni(i) = -0.076788197844621d0
i = 46; ci(i) = 3; di(i) = 5;  ti(i) = 23d0;    ni(i) = 0.022446277332006d0
i = 47; ci(i) = 4; di(i) = 14; ti(i) = 10d0;    ni(i) = -6.2689710414685d-05
i = 48; ci(i) = 6; di(i) = 3;  ti(i) = 50d0;    ni(i) = -5.5711118565645d-10
i = 49; ci(i) = 6; di(i) = 6;  ti(i) = 44d0;    ni(i) = -0.19905718354408d0
i = 50; ci(i) = 6; di(i) = 6;  ti(i) = 46d0;    ni(i) = 0.31777497330738d0
i = 51; ci(i) = 6; di(i) = 6;  ti(i) = 50d0;    ni(i) = -0.11841182425981d0
!------------------------------------------------------------------------------------------------------------------------------
i = 52; di(i) = 3;  ti(i) = 0; ni(i) = -31.306260323435d0; alf(i) = 20; bet(i) = 150; gam(i) = 1.21d0; eps(i) = 1
i = 53; di(i) = 3;  ti(i) = 1; ni(i) = 31.546140237781d0;  alf(i) = 20; bet(i) = 150; gam(i) = 1.21d0; eps(i) = 1
i = 54; di(i) = 3;  ti(i) = 4; ni(i) = -2521.3154341695d0; alf(i) = 20; bet(i) = 250; gam(i) = 1.25d0; eps(i) = 1
!------------------------------------------------------------------------------------------------------------------------------
i = 55; ai(i) = 3.5d0; bi(i) = 0.85d0; bb(i) = 0.2d0;  ni(i) = -0.14874640856724d0
        cc(i) = 28;    dd(i) = 700;    aa(i) = 0.32d0; bet(i) = 0.3d0
i = 56; ai(i) = 3.5d0; bi(i) = 0.95d0; bb(i) = 0.2d0;  ni(i) = 0.31806110878444d0  
        cc(i) = 32;    dd(i) = 800;    aa(i) = 0.32d0; bet(i) = 0.3d0
!------------------------------------------------------------------------------------------------------------------------------

!CP EXTENSION BELOW 130 K
ee = 0.278296458178592d0

end subroutine

!********************************************************************************
subroutine chk_iapws95_table6
!********************************************************************************

!THIS TABLE IS USED TO COMPARE COMPUTED RESULTS WITH VALUES GIVEN IN TABLE A1 OF OS2008.
!OR EQUIVALENTLY IAPWS95_TABLE7

!IAPWS-95:
!REVISED RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC
!PROPERTIES OF ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE.
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM.
!DOORWERTH, THE NETHERLANDS, SEPTEMBER 2009
!AVAILABLE AT: HTTP://WWW.IAPWS.ORG

implicit none

!FUNCTION VALUES AS GIVEN IN TABLE 6 OF IAPWS-95, EDITORIAL VERSION OF BERLIN, 2008:
!VALUES FOR THE IDEAL-GAS PART PHI0, EQ. (5), AND FOR THE RESIDUAL PART PHIR , EQ. (6), OF THE
!DIMENSIONLESS HELMHOLTZ FREE ENERGY TOGETHER WITH THE CORRESPONDING DERIVATIVESA
!FOR T = 500 K AND RHO = 838.025 KG/M3

!PHI0_ =     0.20479773347960E1    PHIR_ =    -0.34269320568156E1
!PHI0_D =    0.38423674711375      PHIR_D =   -0.36436665036388
!PHI0_DD =  -0.147637877832556     PHIR_DD =   0.85606370097461
!PHI0_T =    0.90461110617524E1    PHIR_T =   -0.58140343523842E1
!PHI0_TT =  -0.193249185013052E1   PHIR_TT =  -0.223440736884336E1
!PHI0_DT =   0                     PHIR_DT =  -0.11217691467031E1

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "**************************************************************"
write(*,*) "**************************************************************"
write(*,*) "**************************************************************"
write(*,*) " "
write(*,*) "Implementation of IAPWS-95 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "TABLE A1 OF OS2008"
write(*,*) "UPDATED FUNCTION VALUES FOR TABLE 6 OF IAPWS95"
write(*,*) "VALUES FOR THE IDEAL-GAS PART, PHI0, EQ.(5),"
write(*,*) "FOR T = 500 K AND RHO = 838.025 KG/M3"
write(*,*) " "

write(*,*) "                    TABLE VALUE             CALCULATED VALUE" 
write(*, '(7x, "phi0       0.20479773347960e1 ", 6x, e21.15)')  chk_flu_phi(0, 0, 500d0, 838.025d0, "0")
write(*, '(7x, "phi0_d     0.38423674711375   ", 6x, e21.15)')  chk_flu_phi(0, 1, 500d0, 838.025d0, "0")
write(*, '(7x, "phi0_dd   -0.147637877832556  ", 6x, e21.15)')  chk_flu_phi(0, 2, 500d0, 838.025d0, "0")
write(*, '(7x, "phi0_t     0.90461110617524e1 ", 6x, e21.15)')  chk_flu_phi(1, 0, 500d0, 838.025d0, "0")
write(*, '(7x, "phi0_tt   -0.193249185013052e1", 6x, e21.15)')  chk_flu_phi(2, 0, 500d0, 838.025d0, "0")
write(*, '(7x, "phi0_dt    0.0                ", 6x, e21.15)')  chk_flu_phi(1, 1, 500d0, 838.025d0, "0")
write(*,*) " "

write(*,*) "VALUES FOR THE RESIDUAL PART, PHIR, EQ.(5)"
write(*,*) "FOR T = 500 K AND RHO = 838.025 KG/M3"
write(*,*) " "

write(*,*) "                    TABLE VALUE             CALCULATED VALUE"
write(*, '(7x, "phir      -0.34269320568156e1 ", 6x, e21.15)')  chk_flu_phi(0, 0, 500d0, 838.025d0, "r")    
write(*, '(7x, "phir_d    -0.36436665036388   ", 6x, e21.15)')  chk_flu_phi(0, 1, 500d0, 838.025d0, "r")    
write(*, '(7x, "phir_dd    0.85606370097461   ", 6x, e21.15)')  chk_flu_phi(0, 2, 500d0, 838.025d0, "r")    
write(*, '(7x, "phir_t    -0.58140343523842e1 ", 6x, e21.15)')  chk_flu_phi(1, 0, 500d0, 838.025d0, "r")    
write(*, '(7x, "phir_tt   -0.223440736884336e1", 6x, e21.15)')  chk_flu_phi(2, 0, 500d0, 838.025d0, "r")    
write(*, '(7x, "phir_dt   -0.11217691467031e1 ", 6x, e21.15)')  chk_flu_phi(1, 1, 500d0, 838.025d0, "r")    
write(*,*) " "

end subroutine


!==========================================================================
subroutine chk_iapws95_table7
!==========================================================================

!THIS TABLE IS USED TO COMPARE COMPUTED RESULTS WITH VALUES GIVEN IN TABLE A2 OF OS2008.
!OR EQUIVALENTLY IAPWS95_TABLE7

!IAPWS-95:
!REVISED RELEASE ON THE IAPWS FORMULATION 1995 FOR THE THERMODYNAMIC
!PROPERTIES OF ORDINARY WATER SUBSTANCE FOR GENERAL AND SCIENTIFIC USE.
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM.
!DOORWERTH, THE NETHERLANDS, SEPTEMBER 2009
!AVAILABLE AT: HTTP://WWW.IAPWS.ORG


integer irow, ir, irmin, irmax
real*8 tr(11), rhor(11)
real*8 print1, print2, print3, print4, print5, print6
real*8 t, rho, p, cv, c, sigma
real*8 ft, fd, ftt, fdd, ftd

!FUNCTION VALUES AS GIVEN IN TABLE A2 OF OS2008 AS UPDATES FOR TABLE 7 OF IAPWS-95
!THERMODYNAMIC PROPERTY VALUES IN THE SINGLE-PHASE REGION FOR SELECTED VALUES OF T AND RHO

!T   RHO           P                CV              W               S
!K   KG /M3        MPa              KJ /KG /K       M /S            KJ /KG /K
!300 0.996 556 0E3 0.992 418 352E1  0.413 018 112E1 0.150 151 914E4 0.393 062 643
!    0.100 530 8E4 0.200 022 515E2  0.406 798 347E1 0.153 492 501E4 0.387 405 401
!    0.118 820 2E4 0.700 004 704E3  0.346 135 580E1 0.244 357 992E4 0.132 609 616
!500 0.435 000 0   0.999 679 423E.1 0.150 817 541E1 0.548 314 253E3 0.794 488 271E1
!    0.453 200 0E1 0.999 938 125    0.166 991 025E1 0.535 739 001E3 0.682 502 725E1
!    0.838 025 0E3 0.100 003 858E2  0.322 106 219E1 0.127 128 441E4 0.256 690 919E1
!    0.108 456 4E4 0.700 000 405E3  0.307 437 693E1 0.241 200 877E4 0.203 237 509E1
!647 0.358 000 0E3 0.220 384 756E2  0.618 315 728E1 0.252 145 078E3 0.432 092 307E1
!900 0.241 000 0   0.100 062 559    0.175 890 657E1 0.724 027 147E3 0.916 653 194E1
!    0.526 150 0E2 0.200 000 690E2  0.193 510 526E1 0.698 445 674E3 0.659 070 225E1
!    0.870 769 0E3 0.700 000 006E3  0.266 422 350E1 0.201 933 608E4 0.417 223 802E1

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of IAPWS-95 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "TABLE A2 OF OS2008"
write(*,*) "UPDATED FUNCTION VALUES FOR TABLE 7 OF IAPWS-95:"
write(*,*) "THERMODYNAMIC PROPERTY OF H2O IN THE SINGLE-PHASE REGION FOR SELECTED VALUES OF T AND RHO"
write(*,*) " "

write(*,*) "ROW T    RHO           P              CV             C             SIGMA    "
write(*,*) "    K    KG/M3         MPa            KJ/(KG K)      M/S           KJ/(KG K)"
write(*,*) " "

write(*,*) "1   300  0.9965560e3   0.992418352e-1 0.413018112e1  0.150151914e4 0.393062643"
write(*,*) "2        0.1005308e4   0.200022515e2  0.406798347e1  0.153492501e4 0.387405401"
write(*,*) "3        0.1188202e4   0.700004704e3  0.346135580e1  0.244357992e4 0.132609616"
write(*,*) "4   500  0.4350000     0.999679423e-1 0.150817541e1  0.548314253e3 0.794488271e1"
write(*,*) "5        0.4532000e1   0.999938125    0.166991025e1  0.535739001e3 0.682502725e1"
write(*,*) "6        0.8380250e3   0.100003858e2  0.322106219e1  0.127128441e4 0.256690919e1"
write(*,*) "7        0.1084564e4   0.700000405e3  0.307437693e1  0.241200877e4 0.203237509e1"
write(*,*) "8   647  0.3580000e3   0.220384756e2  0.618315728e1  0.252145078e3 0.432092307e1"
write(*,*) "9   900  0.2410000     0.100062559    0.175890657e1  0.724027147e3 0.916653194e1"
write(*,*) "10       0.5261500e2   0.200000690e2  0.193510526e1  0.698445674e3 0.659070225e1"
write(*,*) "11       0.8707690e3   0.700000006e3  0.266422350e1  0.201933608e4 0.417223802e1"

tr(1) = 300d0
rhor(1) = 0.9965560d3
tr(2) = 300d0
rhor(2) = 0.1005308d4
tr(3) = 300d0
rhor(3) = 0.1188202d4
tr(4) = 500d0
rhor(4) = 0.4350000d0
tr(5) = 500d0
rhor(5) = 0.4532000d1
tr(6) = 500d0
rhor(6) = 0.8380250d3
tr(7) = 500d0
rhor(7) = 0.1084564d4
tr(8) = 647d0
rhor(8) = 0.3580000d3
tr(9) = 900d0
rhor(9) = 0.2410000d0
tr(10) = 900d0
rhor(10) = 0.5261500d2
tr(11) = 900d0
rhor(11) = 0.8707690d3

100 continue
#ifdef interactive
write(*,*) ' '
write(*,*) "ENTER THE ROW TO CHECK: 0 TO EXIT; 999 FOR ALL ROWS; 99 TO ENTER T AND RHO"
read(*,*) irow
#else
irow = 999
#endif

if(irow == 0) then
  return
elseif(irow == 999) then
  irmin=1
  irmax=11
elseif(irow == 99)then
  irmin=1
  irmax=1
else
  irmin=irow
  irmax=irow
endif
    
write(*,*) " "
write(*,*) " "
write(*,*) "VALUES CALCULATED NOW"

  do ir = irmin, irmax
    t=tr(ir)
    rho=rhor(ir)
    if (irow == 99) then
       write(*,*)'ENTER T AND RHO'
       read(*,*)print1, print2
       t = print1
       rho = print2
    endif
  
!COMPUTE PRESSURE
    fd = flu_f_si(0, 1, t, rho)
    p = rho**2 * fd
  
!COMPUTE CV
    ftt = flu_f_si(2, 0, t, rho)
    cv = -t * ftt
  
!COMPUTE SOUND SPEED
    fdd = flu_f_si(0, 2, t, rho)
    ftd = flu_f_si(1, 1, t, rho)
    c = sqrt(rho**2 * (ftt * fdd - ftd**2) / ftt + 2d0 * rho * fd)
    
!COMPUTE ENTROPY
    ft = flu_f_si(1, 0, t, rho)
    sigma = -ft
  
    print1=t
    print2=rho
    print3=1d-6*p
    print4=1d-3*cv
    print5=c
    print6=1d-3*sigma
    write(*,10) ir, print1, print2, print3, print4, print5, print6
10  format(i2, 3x, f4.0, 1x, e14.9, 1x, e14.9, 1x, e14.9, 1x, e14.9, &
      1x, e14.9)
    !write(6,*)1d-6*p
    !write(6,*)1d-3*cv
    !write(6,*)c
    !write(6,*)1d-3*sigma
  enddo

  if(irow /= 999)go to 100

end subroutine


end module flu_1
