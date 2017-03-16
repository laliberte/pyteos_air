module liq_f03_5

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS POTENTIALS OF LIQUID WATER AND ITS
!FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE AND
!DENSITY, AS DEFINED CONSISTENTLY WITH IAPWS-95 FOR THE OCEANOGRAPHIC
!RANGE -6 - 40 C, 0 - 100 MPA, IN THE ARTICLE:

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58 (2003) 43-115

!AND IN THE IAPWS-09 DOCUMENT:
!SUPPLEMENTARY RELEASE ON A COMPUTATIONALLY EFFICIENT THERMODYNAMIC FORMULATION
!FOR LIQUID WATER FOR OCEANOGRAPHIC USE.
!THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!DOORWERTH, THE NETHERLANDS, SEPTEMBER 2009

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

public :: chk_iapws09_table6, &
          fit_liq_cp_f03_si, fit_liq_density_f03_si, &
          fit_liq_expansion_f03_si, fit_liq_g_f03_si, &
          fit_liq_kappa_t_f03_si, fit_liq_soundspeed_f03_si

real*8, private :: gc03(0:7, 0:6)

contains

!==========================================================================
function fit_liq_g_f03_si(drv_t, drv_p, t_si, p_si)
!==========================================================================

!THIS FUNCTION IMPLEMENTS THE GIBBS FUNCTION OF PURE WATER AS DEFINED IN
!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58 (2003) 43-115

!CHECK VALUES:
!FIT_LIQ_G_F03_SI(0,0,300,1E5) =  -0.526505056073E+04
!FIT_LIQ_G_F03_SI(1,0,300,1E5) =  -0.393062597709E+03
!FIT_LIQ_G_F03_SI(0,1,300,1E5) =   0.100345554745E-02
!FIT_LIQ_G_F03_SI(2,0,300,1E5) =  -0.139354762020E+02
!FIT_LIQ_G_F03_SI(1,1,300,1E5) =   0.275754520492E-06
!FIT_LIQ_G_F03_SI(0,2,300,1E5) =  -0.452067557155E-12


implicit none

integer drv_t, drv_p
real*8 fit_liq_g_f03_si
real*8 t_si, p_si
real*8 t0, tu, p0, pu
real*8 y, z, g

t0 = celsius_temperature_si     !IN K
tu = 40d0                       !IN K
p0 = sealevel_pressure_si       !IN PA
pu = 100000000d0 !in pa

fit_liq_g_f03_si = errorreturn

if(drv_t < 0) return
if(drv_p < 0) return
if(t_si <= 0d0) return
if(p_si <= 0d0) return

call init_f03

y = (t_si - t0) / tu
z = (p_si - p0) / pu

g = polyf03_gyz(drv_t, drv_p, y, z)
if(g == errorreturn) return

fit_liq_g_f03_si = g / (tu ** drv_t * pu ** drv_p)

end function

!==========================================================================
function fit_liq_cp_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE HEAT CAPACITY CP OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!FIT_LIQ_DENSITY_F03_SI: CP IN J/(KG K)

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

!CHECK VALUE:
!FIT_LIQ_CP_F03_SI(300, 1E5) = 0.418064286060E+04

real*8 fit_liq_cp_f03_si, t_si, p_si
real*8 g_tt 

fit_liq_cp_f03_si = errorreturn

if (t_si <= 0d0) return
if (p_si <= 0d0) return

g_tt = fit_liq_g_f03_si(2, 0, t_si, p_si)

if (g_tt == errorreturn) return
if (g_tt >= 0d0) return

fit_liq_cp_f03_si = -t_si * g_tt

end function

!==========================================================================
function fit_liq_density_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE DENSITY OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!FIT_LIQ_DENSITY_F03_SI: DENSITY IN KG/M3

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

!CHECK VALUE:
!FIT_LIQ_DENSITY_F03_SI(300, 1E5) = 0.996556352243E+03

implicit none

real*8 fit_liq_density_f03_si
real*8 t_si, p_si
real*8 d

fit_liq_density_f03_si = errorreturn

if(t_si <= 0d0) return
if(p_si <= 0d0) return

d = fit_liq_g_f03_si(0, 1, t_si, p_si)

if(d == errorreturn) return
if(d <= 0d0) return

fit_liq_density_f03_si = 1d0 / d

end function

!==========================================================================
function fit_liq_expansion_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE THERMAL EXPANSION OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!FIT_LIQ_EXPANSION_F03_SI: ALPHA IN 1/K

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

!CHECK VALUE:
!FIT_LIQ_EXPANSION_F03_SI(300, 1E5) = 0.274804919056E-03

real*8 fit_liq_expansion_f03_si
real*8 t_si, p_si
real*8 g_p, g_tp

fit_liq_expansion_f03_si = errorreturn

if (t_si <= 0d0)return
if (p_si <= 0d0)return

g_p = fit_liq_g_f03_si(0, 1, t_si, p_si)
if (g_p == errorreturn) return
if (g_p <= 0d0) return

g_tp = fit_liq_g_f03_si(1, 1, t_si, p_si)
if (g_tp == errorreturn) return

fit_liq_expansion_f03_si = g_tp / g_p

end function

!==========================================================================
function fit_liq_kappa_t_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE ISOTHERMAL COMPRESSIBILITY OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!FIT_LIQ_KAPPA_T_F03_SI: KAPPA_T IN 1/PA

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

!CHECK VALUE:
!FIT_LIQ_KAPPA_T_F03_SI(300, 1E5) = 0.450510795725E-09

real*8 fit_liq_kappa_t_f03_si, t_si, p_si
real*8 g_p, g_pp

fit_liq_kappa_t_f03_si = errorreturn

if (t_si <= 0d0) return
if (p_si <= 0d0) return

g_p = fit_liq_g_f03_si(0, 1, t_si, p_si)
if (g_p == errorreturn) return
if (g_p <= 0d0) return

g_pp = fit_liq_g_f03_si(0, 2, t_si, p_si)
if (g_pp == errorreturn) return
if (g_pp >= 0d0) return

fit_liq_kappa_t_f03_si = -g_pp / g_p

end function

!==========================================================================
function fit_liq_soundspeed_f03_si(t_si, p_si)
!==========================================================================

!THIS FUNCTION RETURNS THE SOUND SPEED OF LIQUID WATER COMPUTED FROM THE
!GIBBS FUNCTION 2003 OF SEAWATER AS A FUNCTION OF TEMPERATURE
!AND PRESSURE, PUBLISHED IN

!R. FEISTEL:
!A NEW EXTENDED GIBBS THERMODYNAMIC POTENTIAL OF SEAWATER.
!PROGRESS IN OCEANOGRAPHY, 58/1 (2003) 43-115

!OUTPUT:  
!FIT_LIQ_SOUNDSPEED_F03_SI: SPEED IN M/S

!INPUT:   
!T_SI: ABSOLUTE TEMPERATURE IN K
!P_SI: ABSOLUTE PRESSURE IN PA

!CHECK VALUE:
!FIT_LIQ_SOUNDSPEED_F03_SI(300, 1e5) = 0.150152808421E+04

real*8 fit_liq_soundspeed_f03_si, t_si, p_si
real*8  g_p, g_tt, g_tp, g_pp, c

fit_liq_soundspeed_f03_si = errorreturn

if (t_si <= 0d0) return
if (p_si <= 0d0) return

g_p = fit_liq_g_f03_si(0, 1, t_si, p_si)
if (g_p == errorreturn) return
if( g_p <= 0d0) return

g_tt = fit_liq_g_f03_si(2, 0, t_si, p_si)
if (g_tt == errorreturn) return
if (g_tt >= 0d0) return

g_tp = fit_liq_g_f03_si(1, 1, t_si, p_si)
if (g_tp == errorreturn) return

g_pp = fit_liq_g_f03_si(0, 2, t_si, p_si)
if (g_pp == errorreturn) return
if (g_pp >= 0d0) return

c = g_tp * g_tp - g_tt * g_pp
if (c == 0d0) return
c = g_tt / c
if (c < 0d0) return

fit_liq_soundspeed_f03_si = g_p * sqrt(c)

end function

!==========================================================================
function polyf03_gyz(drv_y, drv_z, y, z)
!==========================================================================

!RETURNS THE VALUE OF THE POLYNOMIAL DERIVATIVE
!(D/DY)^DRV_Y (D/DZ)^DRV_Z SUM(J,K) GC03(J,K)*Y^J*Z^K

implicit none

integer drv_y, drv_z
integer j, jmax
integer k, kmax
integer l
integer maxt, maxp
real*8 polyf03_gyz, g
real*8 y, z
real*8 yj, zk
real*8 c

maxt = 7
maxp = 6

g = 0d0
if(y == 0d0) then 
   jmax = drv_y 
else 
   jmax = maxt
endif

if(z == 0d0) then 
   kmax = drv_z 
else 
   kmax = maxp
endif

yj = 1d0
do j = drv_y, jmax   !loop over powers of y

  zk = 1d0
  do k = drv_z, kmax    !loop over powers of z

    if(gc03(j, k) /= 0d0) then
      c = gc03(j, k) * yj * zk

      do l = 1, drv_y            !factors from y-derivatives
        c = c * dfloat(j - l + 1)
      enddo

      do l = 1, drv_z            !factors from z-derivatives
        c = c * dfloat(k - l + 1)
      enddo

      g = g + c
    end if

    if(k < kmax) zk = zk * z

  enddo

  if(j < jmax) yj = yj * y
enddo

polyf03_gyz = g

end function

!==========================================================================
subroutine init_f03
!==========================================================================

! G(T,P) = SUM GC03(J,K) * (T/40 C)^J * (P/100 MPA)^K
!          J,K

implicit none

real*8 g, g_t, g_p
real*8 tt, t0, tu, pt, p0, pu, yt, zt

tt = tp_temperature_si
t0 = celsius_temperature_si
tu = 40d0

pt = tp_pressure_iapws95_si !iapws-95 triple point pressure
p0 = sealevel_pressure_si
pu = 100000000d0

!WRITE(*,*)"TT, T0, TU =", tt, t0, tu
!WRITE(*,*)"PT, P0, PU =", pt, p0, pu
yt = (tt - t0) / tu
zt = (pt - p0) / pu

if(gc03(2, 0) == -12357.785933039d0) return

gc03 = 0d0

gc03(2, 0) = -12357.785933039d0
gc03(3, 0) = 736.741204151612d0
gc03(4, 0) = -148.185936433658d0
gc03(5, 0) = 58.0259125842571d0
gc03(6, 0) = -18.9843846514172d0
gc03(7, 0) = 3.05081646487967d0
gc03(0, 1) = 100015.695367145d0
gc03(1, 1) = -270.983805184062d0
gc03(2, 1) = 1455.0364540468d0
gc03(3, 1) = -672.50778314507d0
gc03(4, 1) = 397.968445406972d0
gc03(5, 1) = -194.618310617595d0
gc03(6, 1) = 63.5113936641785d0
gc03(7, 1) = -9.63108119393062d0
gc03(0, 2) = -2544.5765420363d0
gc03(1, 2) = 776.153611613101d0
gc03(2, 2) = -756.558385769359d0
gc03(3, 2) = 499.360390819152d0
gc03(4, 2) = -301.815380621876d0
gc03(5, 2) = 120.520654902025d0
gc03(6, 2) = -22.2897317140459d0
gc03(0, 3) = 284.517778446287d0
gc03(1, 3) = -196.51255088122d0
gc03(2, 3) = 273.479662323528d0
gc03(3, 3) = -239.545330654412d0
gc03(4, 3) = 152.196371733841d0
gc03(5, 3) = -55.2723052340152d0
gc03(6, 3) = 8.17060541818112d0
gc03(0, 4) = -33.3146754253611d0
gc03(1, 4) = 28.9796526294175d0
gc03(1, 4) = 28.9796526294175d0
gc03(2, 4) = -55.5604063817218d0
gc03(2, 4) = -55.5604063817218d0
gc03(3, 4) = 48.8012518593872d0
gc03(4, 4) = -26.3748377232802d0
gc03(5, 4) = 6.48190668077221d0
gc03(0, 5) = 4.20263108803084d0
gc03(1, 5) = -2.13290083518327d0
gc03(2, 5) = 4.34420671917197d0
gc03(3, 5) = -1.66307106208905d0
gc03(0, 6) = -0.546428511471039d0

!IAPWS-95 REFERENCE STATE CONDITION
!ENERGY = 0 AND ENTROPY = 0 AT THE TRIPLE POINT:
!THESE ARE ONLY USED AS A FIRST GUESS HERE SO HIGH ACCURACY IS NOT REQUIRED.
!GC03(0, 0) = 101.342743139672D0
GC03(0, 0) = 101.342743139674D0
!GC03(1, 0) = 5.90578348518236D0
gc03(1, 0) =  5.90578347909402d0

!NEVERTHELESS, THE ACCURACY CAN BE IMPROVED FOLLOWING KIYOSHI'S SUGGESTION.
!ADJUSTMENT OF COEFFICIENTS


!WRITE(*,*)"GC00, GC10: FIRST ESTIMATES = ", gc03(0,0), gc03(1,0)

g = polyf03_gyz(0, 0, yt, zt)
g_t = polyf03_gyz(1, 0, yt, zt) / tu
g_p = polyf03_gyz(0, 1, yt, zt) / pu

!WRITE(*,*)'G,G_T,G_P=',G,G_T,G_P
gc03(0, 0) = gc03(0, 0) - g + pt * g_p
gc03(1, 0) = gc03(1, 0) - g_t * tu

!WRITE(*,*)"GC00, GC10: SECOND ESTIMATES = ", gc03(0,0), gc03(1,0)

end subroutine

!==========================================================================
subroutine chk_iapws09_table6
!==========================================================================

!TABLE 6  NUMERICAL CHECK VALUES FOR THE GIBBS FUNCTION G AND ITS DERIVATIVES, TABLE 4,
!THE NUMERICAL FUNCTIONS EVALUATED HERE AT GIVEN POINTS (T, P) ARE DEFINED IN TABLES 3 AND 4.
!
!QUANTITY       VALUE              VALUE             VALUE          UNIT
!T              273.15             273.15            313.15         K
!P              101325             1E8               101325         PA
!G             0.101342743E3     0.977303868E5    -0.116198898E5    J KG-1
!(DG/DT)_P     0.147644587       0.851506346E1    -0.572365181E3    J KG-1 K-1
!(DG/DP)_T     0.100015695E-2    0.956683354E-3    0.100784471E-2   M3 KG-1
!(D2G/DT2)_P  -0.154472324E2    -0.142970174E2    -0.133463968E2    J KG-1 K-2
! D2G/DTDP    -0.677459513E-7    0.199088060E-6    0.388499694E-6   M3 KG-1 K-1
!(D2G/DP2)_T  -0.508915308E-12  -0.371527164E-12  -0.445841077E-12  M3 KG-1 Pa-1
!H             0.610136242E2     0.954044973E5     0.167616267E6    J KG-1
!F             0.183980891E-2    0.206205140E4    -0.117220097E5    J KG-1
!U            -0.403272791E2    -0.263838183E3     0.167514147E6    J KG-1
!S            -0.147644587      -0.851506346E1     0.572365181E3    J KG-1 K-1
!RHO           0.999843071E3     0.104527793E4     0.992216354E3    KG M-3
!CP            0.421941153E4     0.390523030E4     0.417942416E4    J KG-1 K-1
!W             0.140240099E4     0.157543089E4     0.152891242E4    M S-1

real*8 t, p, q(3)
integer i, j

                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of IAPWS-09 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "

write(*,*) "Function values as given in Table 6 of IAPWS-09:"
write(*,*) "Top row is tabulated value; bottom calculated now"
write(*,*) " "

write(*,*) "Quantity    Value            Value             Value          Unit" 
write(*,*) " "

do i = 1, 15
  Select Case (i)
    Case (1);  write(*,*) "T           273.15           273.15            313.15          K"
    Case (2);  write(*,*) "p           101325           1E8               101325          Pa";
               write(*,*) " "
    Case (3);  write(*,*) "g         0.101342743E+3    0.977303868E+5   -0.116198898E+5   J kg-1"
    Case (4);  write(*,*) "dg/dT     0.147644587       0.851506346E+1   -0.572365181E+3   J kg-1 K-1"
    Case (5);  write(*,*) "dg/dp     0.100015695E-2    0.956683354E-3    0.100784471E-2   m3 kg-1"
    Case (6);  write(*,*) "d2g/dT2  -0.154472324E+2   -0.142970174E+2   -0.133463968E+2   J kg-1 K-2"
    Case (7);  write(*,*) "d2g/dTdp -0.677459513E-7    0.199088060E-6    0.388499694E-6   m3 kg-1 K-1"
    Case (8);  write(*,*) "d2g/dp2  -0.508915308E-12  -0.371527164E-12  -0.445841077E-12  m3 kg-1 Pa-1"
    Case (9);  write(*,*) "h         0.610136242E+2    0.954044973E+5    0.167616267E+6   J kg-1"
    Case (10); write(*,*) "f         0.183980891E-2    0.206205140E+4   -0.117220097E+5   J kg-1"
    Case (11); write(*,*) "u        -0.403272791E+2   -0.263838183E+3    0.167514147E+6   J kg-1"
    Case (12); write(*,*) "s        -0.147644587      -0.851506346E+1    0.572365181E+3   J kg-1 K-1"
    Case (13); write(*,*) "rho       0.999843071E+3    0.104527793E+4    0.992216354E+3   kg m-3"
    Case (14); write(*,*) "cp        0.421941153E+4    0.390523030E+4    0.417942416E+4   J kg-1 K-1"
    Case (15); write(*,*) "w         0.140240099E+4    0.157543089E+4    0.152891242E+4   m s-1"
  end select

  
  if(i > 2) then
    
    do j = 1, 3
      if(j == 1)then
         t = 273.15d0
         p = 101325d0
      elseif(j == 2)then
         t = 273.15d0
         p = 1d8
      elseif(j == 3)then
         t = 313.15d0
         p = 101325d0
      endif
      q(j) = 0d0
      Select Case (i)
        Case (1); q(j) = t
        Case (2); q(j) = p
        Case (3); q(j) = fit_liq_g_f03_si(0, 0, t, p)
        Case (4); q(j) = fit_liq_g_f03_si(1, 0, t, p)
        Case (5); q(j) = fit_liq_g_f03_si(0, 1, t, p)
        Case (6); q(j) = fit_liq_g_f03_si(2, 0, t, p)
        Case (7); q(j) = fit_liq_g_f03_si(1, 1, t, p)
        Case (8); q(j) = fit_liq_g_f03_si(0, 2, t, p)
        Case (9); q(j) = fit_liq_g_f03_si(0, 0, t, p) - t * fit_liq_g_f03_si(1, 0, t, p)
        Case (10); q(j) = fit_liq_g_f03_si(0, 0, t, p) - p * fit_liq_g_f03_si(0, 1, t, p)
        Case (11); q(j) = fit_liq_g_f03_si(0, 0, t, p) - p * fit_liq_g_f03_si(0, 1, t, p) &
                                                       - t * fit_liq_g_f03_si(1, 0, t, p)
        Case (12); q(j) = -fit_liq_g_f03_si(1, 0, t, p)
        Case (13); q(j) = fit_liq_density_f03_si(t, p)
        Case (14); q(j) = fit_liq_cp_f03_si(t, p)
        Case (15); q(j) = fit_liq_soundspeed_f03_si(t, p)
      End Select
    enddo
    write(*,'(10x, e16.9, 2x, e16.9, 2x, e16.9)')q(1), q(2), q(3)
    write(*,*)" "
  end if
enddo

end subroutine

end module liq_f03_5
