module ice_1

!#########################################################################

!THIS MODULE IMPLEMENTS THE GIBBS POTENTIAL OF HEXAGONAL ICE I AND ITS
!FIRST AND SECOND PARTIAL DERIVATIVES WITH RESPECT TO TEMPERATURE AND
!PRESSURE AS DEFINED IN IAPWS-06:

!REVISED RELEASE ON AN EQUATION OF STATE i2006 FOR H2O ICE IH
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

!THIS MODULE REQUIRES 
!     CONSTANTS_0, FILE CONSTANTS_0.F90

!#########################################################################

use constants_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: chk_iapws06_table6, ice_g_si 

complex*16, private :: tk(2), rk1, rk2(3)
real*8, private :: g0(10),s0,tt,pt,pnorm,p0

contains

!==========================================================================
function ice_g_si(drv_t, drv_p, t_si, p_si)
!==========================================================================
                         
!OUTPUT:
!ICE_G_SI:          DERIVATIVE OF THE ICE IH GIBBS FUNCTION IN BASIC SI UNITS

!INPUTS:
!DRV_T:             ORDER OF THE PARTIAL DERIVATIVE WRT TEMPERATURE
!DRV_P:             ORDER OF THE PARTIAL DERIVATIVE WRT PRESSURE
!T_SI:              ABSOLUTE TEMPERATURE ITS-90 IN K
!P_SI:              PRESSURE IN PASCALS

!CHECK VALUES:
!ice_g_si(0,0,270,1E5) = -3786.74963128
!ice_g_si(1,0,270,1E5) = 1244.97335506
!ice_g_si(0,1,270,1E5) = 1.09029713624E-03
!ice_g_si(2,0,270,1E5) = -7.67955356337
!ice_g_si(1,1,270,1E5) = 1.72604208697E-07
!ice_g_si(0,2,270,1E5) = -1.27811223643E-13

implicit none

integer drv_t, drv_p
real*8 t_si, p_si
real*8 ice_g_si

ice_g_si = errorreturn

call init_iapws06

if(check_limits == 1) then
  !ICE_LIMITS
  if(t_si <= ice_tmin .or. t_si > ice_tmax .or. &
  p_si <= ice_pmin .or. p_si > ice_pmax) return
else
  if(t_si <= 0d0 .or. p_si <= 0d0) return
endif

if(drv_p == 0) then

   if(drv_t == 0) then
      ice_g_si = ice_g(t_si, p_si)
   elseif(drv_t == 1) then
      ice_g_si = ice_dgdt(t_si, p_si)
   elseif(drv_t == 2) then
      ice_g_si = ice_d2gdt2(t_si, p_si)
   elseif(drv_t == 3) then
      ice_g_si = ice_d3gdt3(t_si, p_si)
   else
      write(*,*)'ERROR: T-DERIVATIVES HIGHER THAN THIRD ORDER NOT IMPLEMENTED'
   endif

elseif(drv_p == 1) then

   if(drv_t == 0) then
      ice_g_si = ice_dgdp(t_si, p_si)
   elseif(drv_t == 1) then
      ice_g_si = ice_d2gdtdp(t_si, p_si)
   elseif(drv_t == 2) then
      ice_g_si = ice_d3gdt2dp(t_si, p_si)
   else
      write(*,*)'ERROR: DERIVATIVES HIGHER THAN THIRD ORDER NOT IMPLEMENTED'
   endif

elseif(drv_p == 2) then

   if(drv_t == 0) then
      ice_g_si = ice_d2gdp2(t_si, p_si)
   elseif(drv_t == 1) then
      ice_g_si = ice_d3gdtdp2(t_si, p_si)
   else
      write(*,*)'ERROR: DERIVATIVES HIGHER THAN THIRD ORDER NOT IMPLEMENTED'
   endif

elseif(drv_p == 3) then

   if(drv_t == 0) then
      ice_g_si = ice_d3gdp3(t_si, p_si)
   else
      write(*,*)'ERROR: DERIVATIVES HIGHER THAN THIRD ORDER NOT IMPLEMENTED'
   endif

elseif(drv_p >= 4) then

   write(*,*)'ERROR: DERIVATIVES HIGHER THAN THIRD ORDER NOT IMPLEMENTED'

endif

end function

!==========================================================================
!                                                               
!           FORTRAN IMPLEMENTATION BY CHRISTOPH GUDER,          
!              C.GUDER@THERMO.RUHR-UNI-BOCHUM.DE                
!               (CODE VERSION: 22 JUNE 2006)                
!                                                               
!                          OF THE
!
!        IAPWS RELEASE ON AN EQUATION OF STATE FOR H2O ICE IH
! THE INTERNATIONAL ASSOCIATION FOR THE PROPERTIES OF WATER AND STEAM
!                    WITNEY, UK, SEPTEMBER 2006
!                                                               
!    MORE DETAILS:                          
!                                                               
!            RAINER FEISTEL AND WOLFGANG WAGNER
!    	   A NEW EQUATION OF STATE FOR H2O ICE IH
!     J. PHYS. CHEM. REF. DATA, 35, 2 (2006) 1021-1047
!
!==========================================================================

!==========================================================================
function ice_g(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE SPECIFIC GIBBS ENERGY FOR ICE IN J/KG

implicit none
real*8 ice_g, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_g = 0d0

tn = t / tt
pn = p / pt

do i=0,4
  ice_g = ice_g + g0(i+1) * (pn-p0)**i
enddo

ice_g = ice_g - s0 * tt * tn

sr(1) = rk1
sr(2) = (0d0,0d0)
do i=0, 2
  sr(2) = sr(2) + rk2(i+1) * (pn-p0)**i
enddo

do i=1,2
  term(i) = (tk(i)-tn)*cdlog(tk(i)-tn)+(tk(i)+tn)*cdlog(tk(i)+tn) &
            -2d0*tk(i)*cdlog(tk(i))-tn**2/tk(i)
enddo

do i=1,2
  ice_g = ice_g + tt * dble(sr(i) * term(i))
enddo

end function

!==========================================================================
function ice_dgdt(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 1ST T-DERIVATIVE OF G IN J/(KG K)

implicit none
real*8 ice_dgdt, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_dgdt = -s0

tn = t / tt
pn = p / pt

sr(1) = rk1
sr(2) = (0d0,0d0)
do i=0, 2
  sr(2) = sr(2) + rk2(i+1) * (pn-p0)**i
enddo

do i=1,2
  term(i) = -cdlog(tk(i)-tn)+cdlog(tk(i)+tn)- &
            2d0*tn/tk(i)
enddo

do i=1,2
  ice_dgdt = ice_dgdt  + dble(sr(i) * term(i))
enddo

end function

!==========================================================================
function ice_dgdp(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 1ST P-DERIVATIVE OF G IN M**3/KG

implicit none
real*8 ice_dgdp, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_dgdp = 0d0

tn = t / tt
pn = p / pt

do i=1,4
  ice_dgdp = ice_dgdp + g0(i+1)*(pn-p0)**(i-1)*dble(i)*(1d0/pt)
enddo

sr(2) = (0d0,0d0)
do i=1, 2
  sr(2) = sr(2)+rk2(i+1)*(pn-p0)**(i-1)*dble(i)*1d0/pt
enddo

term(2) = (tk(2)-tn)*cdlog(tk(2)-tn)+(tk(2)+tn)*cdlog(tk(2)+tn) &
          -2d0*tk(2)*cdlog(tk(2))-tn**2d0/tk(2)

ice_dgdp = ice_dgdp + tt * dble(sr(2) * term(2))

end function

!==========================================================================
function ice_d2gdt2(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 2ND T-DERIVATIVE OF G IN J/(KG K**2)

implicit none
real*8 ice_d2gdt2, t, p, tn, pn
integer i
complex*16 sr(2),term(2)

ice_d2gdt2 = 0d0

tn = t / tt
pn = p / pt

sr(1) = rk1
sr(2) = (0d0,0d0)
do i=0, 2
  sr(2) = sr(2) + rk2(i+1) * (pn-p0)**i
enddo

do i=1,2
  term(i) = 1d0/(tk(i)-tn)+1d0/(tk(i)+tn)-2d0/tk(i)
enddo

do i=1,2
  ice_d2gdt2 = ice_d2gdt2 + 1d0/tt * dble(sr(i) * term(i))
enddo

end function

!==========================================================================
function ice_d2gdtdp(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 2ND T-P-DERIVATIVE OF G IN M**3/(KG K)

implicit none
real*8 ice_d2gdtdp, t, p, tn, pn
integer i
complex*16 sr(2),term(2)

ice_d2gdtdp = 0d0

tn = t / tt
pn = p / pt

sr(2) = (0d0,0d0)
do i=1, 2
  sr(2) = sr(2)+rk2(i+1)*(pn-p0)**(i-1)*dble(i)*1d0/pt
enddo

term(2) = -cdlog(tk(2)-tn)+cdlog(tk(2)+tn)-2d0*tn/tk(2)

ice_d2gdtdp = dble(sr(2) * term(2))

end function

!==========================================================================
function ice_d2gdp2(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 2ND P-DERIVATIVE OF G IN M**3/(KG PA)

implicit none
real*8 ice_d2gdp2, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_d2gdp2 = 0d0

tn = t / tt
pn = p / pt

do i=2,4
  ice_d2gdp2 = ice_d2gdp2+g0(i+1)*(pn-p0)**(i-2)*dble(i)*dble(i-1) &
               *(1d0/(pt**2d0))
enddo

sr(2) = rk2(3)*2d0*(1d0/pt**2)

term(2) = (tk(2)-tn)*cdlog(tk(2)-tn)+(tk(2)+tn)*cdlog(tk(2)+tn) &
          -2d0*tk(2)*cdlog(tk(2))-tn**2d0/tk(2)

ice_d2gdp2 = ice_d2gdp2 + tt * dble(sr(2) * term(2))

end function

!==========================================================================
function ice_d3gdt3(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 3RD T-DERIVATIVE OF G IN J/(KG K**3)

implicit none
real*8 ice_d3gdt3, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_d3gdt3 = 0d0

tn = t / tt
pn = p / pt

sr(1) = rk1
sr(2) = (0d0,0d0)
do i=0, 2
  sr(2) = sr(2) + rk2(i+1) * (pn-p0)**i
enddo

do i=1,2
  term(i) = (1d0/((tk(i)-tn)**2))-(1d0/((tk(i)+tn)**2))
enddo

do i=1,2
  ice_d3gdt3 = ice_d3gdt3 + 1d0/tt**2 * dble(sr(i) * term(i))
enddo

end function

!==========================================================================
function ice_d3gdt2dp(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 3RD T-T-P-DERIVATIVE OF G IN M**3/(KG K**2)

implicit none
real*8 ice_d3gdt2dp, t, p, tn, pn
complex*16 sr(2),term(2)
integer i

ice_d3gdt2dp = 0d0

tn = t / tt
pn = p / pt

sr(2) = (0d0,0d0)
do i=1, 2
  sr(2) = sr(2)+rk2(i+1)*(pn-p0)**(i-1)*dble(i)*1d0/pt
enddo

term(2) = 1d0/(tk(2)-tn)+1d0/(tk(2)+tn)-2d0/tk(2)

ice_d3gdt2dp = 1d0/tt * dble(sr(2) * term(2))

end function

!==========================================================================
function ice_d3gdtdp2(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 3RD T-P-P-DERIVATIVE OF G IN M**3/(KG PA K)

implicit none
real*8 ice_d3gdtdp2, t, p, tn, pn
complex*16 sr(2),term(2)

ice_d3gdtdp2 = 0d0

tn = t / tt
pn = p / pt

sr(2) = rk2(3)*2d0*(1d0/pt**2)

term(2) = -cdlog(tk(2)-tn)+cdlog(tk(2)+tn)- 2d0*tn/tk(2)

ice_d3gdtdp2 = dble(sr(2) * term(2))

end function

!==========================================================================
function ice_d3gdp3(t,p)
!==========================================================================

!THIS FUNCTION CALCULATES THE 3RD P-DERIVATIVE OF G IN M**3/(KG PA**2)

implicit none
real*8 ice_d3gdp3, t, p, tn, pn
integer i

ice_d3gdp3 = 0d0

tn = t / tt
pn = p / pt

do i=3,4
  ice_d3gdp3 = ice_d3gdp3+g0(i+1)*(pn-p0)**(i-3)*dble(i)*dble(i-1) &
               *dble(i-2)/pt**3
enddo
 
end function

!==========================================================================
subroutine init_iapws06
!==========================================================================

!THIS SUBROUTINE CONTAINS THE CONSTANTS FOR THE CALCULATIONS

implicit none
real*8 g01

g01 = -6.32020233335886d5

if(g0(1) == g01)return

!THE ORIGINAL VALUE OF G0(1) FROM THE IAPWS-06 RELEASE IS
!G0(1) = -632020.233449497
!THE VALUE OBTAINED FROM A QUADRUPLE PRECISION IMPLEMENTATION IS
!G0(1) = -6.32020233335886d5 AS PROVIDED IN THE REVISED RELEASE OF 2009

g0(1) = -6.32020233335886d5
g0(2) = 0.655022213658955d0
g0(3) = -1.89369929326131d-8
g0(4) = 3.39746123271053d-15
g0(5) = -5.56464869058991d-22

tk(1) = (3.68017112855051d-2,5.10878114959572d-2)
tk(2) = (0.337315741065416d0,0.335449415919309d0)

rk1 = (44.7050716285388d0,65.6876847463481d0)

rk2(1) = (-72.597457432922d0,-78.100842711287d0)
rk2(2) = (-5.57107698030123d-5,4.64578634580806d-5)
rk2(3) = (2.34801409215913d-11,-2.85651142904972d-11)

s0 = -3327.33756492168d0

tt = 273.16d0
pt = 611.657d0
pnorm = 101325.d0
p0 = pnorm / pt

end subroutine

!==========================================================================
subroutine chk_iapws06_table6
!==========================================================================

integer icol
real*8 zt(3), zp(3)
real*8 zg(3), zdgdp_t(3), zdgdt_p(3), zd2gdp2_t(3),  &
  zd2gdpdt(3), zd2gdt2_p(3), zh(3), zf(3), zu(3), zsigma(3), zcp(3), &
  zrho(3), zalpha(3), zbeta(3), zkappa_t(3), zkappa_s(3)
real*8 g, gt, gp, gtt, gtp, gpp
real*8 t, p
          
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "TABLE A4 OF OS2008 AND TABLE 6 OF IAPWS-06"
write(*,*) "Implementation of IAPWS-06 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) 'UPDATED FUNCTION VALUES FOR TABLE 6 OF IAPWS-06:'
write(*,*) 'PROPERTIES AT THE TRIPLE POINT, THE NORMAL PRESSURE MELTING POINT, '
write(*,*) 'AND AT t = 100 K, p = 100 MPa'
write(*,*) " "

write(*,*)'QUANTITY     VALUE AT            VALUE AT            &
 &VALUE AT           UNITS'
write(*,*)'             t = 273.16 K        t = 273.152519 K    &
 &t = 100 K'
write(*,*)'             p = 611.657 Pa      p = 101325 Pa       &
 &p = 100 MPa'
write(*,*) " "

write(*,*)' g           0.611784135         0.10134274069e+3   &
 &-0.222296513088e+6  J/KG'
write(*,*)'(dg/dp)_t    0.109085812737e-2   0.109084388214e-2  &
 & 0.106193389260e-2  M3/KG'
write(*,*)'(dg/dt)_p    0.122069433940e+4   0.122076932550e+4  &
 & 0.261195122589e+4  J/KG/K'
write(*,*)'(d2g/dp2)_t -0.128495941571e-12 -0.128485364928e-12 &
 &-0.941807981761e-13 M3/KG/PA'
write(*,*)' d2g/dpdt    0.174387964700e-6   0.174362219972e-6  &
 & 0.274505162488e-7  M3/KG/K'
write(*,*)'(d2g/dt2)_p -0.767602985875e+1  -0.767598233365e+1  &
 &-0.866333195517e+1  J/KG/K^2'
write(*,*)' h          -0.333444253966e+6  -0.333354873637e+6  &
 &-0.483491635676e+6  J/KG1'
write(*,*)' f          -0.55446875e-1      -0.918701567e-1     &
 &-0.328489902347e+6  J/KG'
write(*,*)' u          -0.333444921197e+6  -0.333465403393e+6  &
 &-0.589685024936e+6  J/KG'
write(*,*)' sigma      -0.122069433940e+4  -0.122076932550e+4  &
 &-0.261195122589e+4  J/KG/K'
write(*,*)' cp          0.209678431622e+4   0.209671391024e+4  &
 & 0.866333195517e+3  J/KG/K'
write(*,*)' rho         0.916709492200e+3   0.916721463419e+3  &
 & 0.941678203297e+3  KG/M3'
write(*,*)' alpha       0.159863102566e-3   0.159841589458e-3  &
 & 0.258495528207e-4  1/K'
write(*,*)' beta        0.135714764659e+7   0.135705899321e+7  &
 & 0.291466166994e+6  Pa/K'
write(*,*)' kappa_t     0.117793449348e-9   0.117785291765e-9  &
 & 0.886880048115e-10 1/Pa'
write(*,*)' kappa_s     0.114161597779e-9   0.114154442556e-9  &
 & 0.886060982687e-10 1/Pa'

zt(1) = 273.16d0
zp(1) = 611.657d0
zg(1)=            0.611784135d0
zdgdp_t(1)=       0.109085812737d-2
zdgdt_p(1)=       0.122069433940d+4 
zd2gdp2_t(1)=    -0.128495941571d-12
zd2gdpdt(1)=      0.174387964700d-6
zd2gdt2_p(1)=    -0.767602985875d+1 
zh(1)=           -0.333444253966d+6
zf(1)=           -0.55446875d-1    
zu(1)=           -0.333444921197d+6
zsigma(1)=       -0.122069433940d+4
zcp(1)=           0.209678431622d+4
zrho(1)=          0.916709492200d+3
zalpha(1)=        0.159863102566d-3
zbeta(1)=         0.135714764659d+7
zkappa_t(1)=      0.117793449348d-9
zkappa_s(1)=      0.114161597779d-9

zt(2) = 273.152519d0
zp(2) = 101325d0
zg(2)=            0.10134274069d+3
zdgdp_t(2)=       0.109084388214d-2
zdgdt_p(2)=       0.122076932550d+4
zd2gdp2_t(2)=    -0.128485364928d-12
zd2gdpdt(2)=      0.174362219972d-6
zd2gdt2_p(2)=    -0.767598233365d+1
zh(2)=           -0.333354873637d+6
zf(2)=           -0.918701567d-1 
zu(2)=           -0.333465403393d+6
zsigma(2)=       -0.122076932550d+4
zcp(2)=           0.209671391024d+4
zrho(2)=          0.916721463419d+3
zalpha(2)=        0.159841589458d-3
zbeta(2)=         0.135705899321d+7
zkappa_t(2)=      0.117785291765d-9
zkappa_s(2)=      0.114154442556d-9

zt(3) = 100d0
zp(3) = 100000000d0
zg(3)=           -0.222296513088d+6
zdgdp_t(3)=       0.106193389260d-2
zdgdt_p(3)=       0.261195122589d+4
zd2gdp2_t(3)=    -0.941807981761d-13
zd2gdpdt(3)=      0.274505162488d-7
zd2gdt2_p(3)=    -0.866333195517d+1
zh(3)=           -0.483491635676d+6
zf(3)=           -0.328489902347d+6
zu(3)=           -0.589685024936d+6
zsigma(3)=       -0.261195122589d+4
zcp(3)=           0.866333195517d+3
zrho(3)=          0.941678203297d+3
zalpha(3)=        0.258495528207d-4
zbeta(3)=         0.291466166994d+6
zkappa_t(3)=      0.886880048115d-10
zkappa_s(3)=      0.886060982687d-10


100 continue
#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'0 TO EXIT; 999 FOR ALL COLUMNS'
read(*,*) icol
if(icol == 0) return
#else
icol = 999
#endif

if(icol == 999) then

  write(*,*) " "
  write(*,*) " "
  write(*,*) "VALUES CALCULATED NOW"
  write(*,*) " "

  write(*,*)'QUANTITY     VALUE AT            VALUE AT            VALUE AT           UNITS'
  write(*,*)'             t = 273.16 K        t = 273.152519 K    t = 100 K'
  write(*,*)'             p = 611.657 Pa      p = 101325 Pa       p = 100 MPa'
  write(*,*) " "

  write(*,*)' g          ',zg(1), '       ',     zg(2), '      ',      zg(3), '    ',       'J/KG'
  write(*,*)'(dg/dp)_t   ',zdgdp_t(1), '  ',     zdgdp_t(2), '  ',     zdgdp_t(3), '   ',   'M3/KG'
  write(*,*)'(dg/dt)_p   ',zdgdt_p(1), '      ', zdgdt_p(2), '      ', zdgdt_p(3), '     ', 'J/KG/K'
  write(*,*)'(d2g/dp2)_t ',zd2gdp2_t(1), '',     zd2gdp2_t(2), '',     zd2gdp2_t(3), '',    'M3/KG/PA'
  write(*,*)' d2g/dpdt   ',zd2gdpdt(1), '    ',  zd2gdpdt(2), '  ',    zd2gdpdt(3), '  ',   'M3/KG/K'
  write(*,*)'(d2g/dt2)_p ',zd2gdt2_p(1), '    ', zd2gdt2_p(2), '    ', zd2gdt2_p(3), '    ','J/KG/K^2'
  write(*,*)' h          ',zh(1), '    ',        zh(2), '    ',        zh(3), '    ',       'J/KG1'
  write(*,*)' f          ',zf(1), '      ',      zf(2), '     ',       zf(3), '    ',       'J/KG'
  write(*,*)' u          ',zu(1), '    ',        zu(2), '    ',        zu(3), '    ',       'J/KG'
  write(*,*)' sigma      ',zsigma(1), '     ',   zsigma(2), '     ',   zsigma(3), '    ',   'J/KG/K'
  write(*,*)' cp         ',zcp(1), '     ',      zcp(2), '     ',      zcp(3), '     ',     'J/KG/K'
  write(*,*)' rho        ',zrho(1), '       ',   zrho(2), '     ',     zrho(3), '     ',    'KG/M3'
  write(*,*)' alpha      ',zalpha(1), '  ',      zalpha(2), '  ',      zalpha(3), '  ',     '1/K'
  write(*,*)' beta       ',zbeta(1), '     ',    zbeta(2), '     ',    zbeta(3), '     ',   'Pa/K'
  write(*,*)' kappa_t    ',zkappa_t(1), ' ',     zkappa_t(2), ' ',     zkappa_t(3), ' ',    '1/Pa'
  write(*,*)' kappa_s    ',zkappa_s(1), ' ',     zkappa_s(2), ' ',     zkappa_s(3), ' ',    '1/Pa'

else

  t = zt(icol)
  p = zp(icol)
  g = ice_g_si(0, 0, t, p)
  gt = ice_g_si(1, 0, t, p)
  gp = ice_g_si(0, 1, t, p)
  gtt = ice_g_si(2, 0, t, p)
  gtp = ice_g_si(1, 1, t, p)
  gpp = ice_g_si(0, 2, t, p)

  write(*,*) " "
  write(*,*) " FUNCTION VALUES AS GIVEN IN TABLE 6, COLUMN", ICOL," OF IAPWS-06"
  write(*,*) " PROPERTIES AT T =", zt(icol)," K, P =", zp(icol)," Pa"
  write(*,*) " "

  write(*,*) " QUANTITY             TABLE VALUE              COMPUTED VALUE     &
   &   UNIT"
  write(*,*) " "

  write(*,'(" g           ", 4x,e22.15,5x,e22.15,4x,"J KG-1")') &
    zg(icol), g
  write(*,'("(dg/dp)_t    ", 4x,e22.15,5x,e22.15,4x,"M3 KG-1")') &
    zdgdp_t(icol), gp
  write(*,'("(dg/dt)_p    ", 4x,e22.15,5x,e22.15,4x,"J KG-1 K-1")') &
    zdgdt_p(icol), gt
  write(*,'("(d2g/dp2)_t  ", 4x,e22.15,5x,e22.15,4x,"M3 KG-1 Pa-1")') &
    zd2gdp2_t(icol), gpp
  write(*,'(" d2g/dpdt    ", 4x,e22.15,5x,e22.15,4x,"M3 KG-1 K-1")') &
    zd2gdpdt(icol), gtp
  write(*,'("(d2g/dt2)_p  ", 4x,e22.15,5x,e22.15,4x,"J KG-1 K-2")') &
    zd2gdt2_p(icol), gtt
  write(*,'(" h           ", 4x,e22.15,5x,e22.15,4x,"J KG-1")') &
    zh(icol), g - t * gt
  write(*,'(" f           ", 4x,e22.15,5x,e22.15,4x,"J KG-1")') &
    zf(icol), g - p * gp
  write(*,'(" u           ", 4x,e22.15,5x,e22.15,4x,"J KG-1")') &
    zu(icol), g - t * gt - p * gp
  write(*,'(" sigma       ", 4x,e22.15,5x,e22.15,4x,"J KG-1 K-1")') &
    zsigma(icol), -gt
  write(*,'(" cp          ", 4x,e22.15,5x,e22.15,4x,"J KG-1 K-1")') &
    zcp(icol), -t * gtt
  write(*,'(" rho         ", 4x,e22.15,5x,e22.15,4x,"KG M-3")') &
    zrho(icol), 1d0/gp
  write(*,'(" alpha       ", 4x,e22.15,5x,e22.15,4x,"K-1")') &
    zalpha(icol), gtp/gp
  write(*,'(" beta        ", 4x,e22.15,5x,e22.15,4x,"Pa K-1")') &
    zbeta(icol), -gtp/gpp
  write(*,'(" kappa_t     ", 4x,e22.15,5x,e22.15,4x,"Pa-1")') &
    zkappa_t(icol), -gpp/gp
  write(*,'(" kappa_s     ", 4x,e22.15,5x,e22.15,4x,"Pa-1")') &
    zkappa_s(icol), (gtp**2-gtt*gpp)/(gp*gtt)
  
endif

if(icol /= 999) go to 100

end subroutine


end module ice_1
