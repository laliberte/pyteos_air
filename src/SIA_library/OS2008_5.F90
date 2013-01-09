module os2008_5

! #define interactive

!=========================================================================

!THIS MODULE IS USED TO PRODUCE COMPARISON TABLES WITH THE CHECK VALUES PUBLISHED IN
!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCI., 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/
!WITH THE IMPLEMENTATION IN THIS LIBRARY.

!==========================================================================
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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., !WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCES, 2009, IN PREPARATION.
!==========================================================================

!THIS MODULE REQUIRES
!     FLU_1,       FILE FLU_1.F90
!     ICE_1,       FILE ICE_1.F90
!     SAL_1,       FILE SAL_1.F90
!     FLU_2,       FILE FLU_2.F90
!     SAL_2,       FILE SAL_2.F90
!     FLU_3A,      FILE FLU_3A.F90
!     LIQ_VAP_4,   FILE LIQ_VAP_4.F90

!NESTED USE ASSOCIATIONS
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     MATHS_0,     FILE MATHS_0.F90

use flu_1
use ice_1
use sal_1
use flu_2
use sal_2
use flu_3a
use liq_vap_4

implicit none
private

public chk_os2008_table
 
character*16, private :: version = '8 Jun 2010'

contains

!==========================================================================
subroutine chk_os2008_table(table_number)
!==========================================================================

implicit none
character*(*) table_number

if(table_number == "2") then
   call chk_os2008_table2
elseif(table_number == "3") then
   call chk_os2008_table3
elseif(table_number == "A1") then
   call chk_os2008_tablea1  !calls chk_iapws95_table6 which is equivalent
                            !avoids need to make flu_phir and flu_phi0 public
elseif(table_number == "A2") then
   call chk_os2008_tableA2
elseif(table_number == "A3") then
   call chk_os2008_tableA3
elseif(table_number == "A4") then
   call chk_os2008_tableA4
elseif(table_number == "A5") then
   call chk_os2008_tableA5
elseif(table_number == "A6") then
   call chk_os2008_tableA6
elseif(table_number == "A7") then
   call chk_os2008_tableA7
elseif(table_number == "A8") then
   call chk_os2008_tableA8
elseif(table_number == "0") then
   stop
else
   write(*,*) "ERROR: YOU DIDN'T ENTER A LEGAL TABLE IDENTIFIER"
endif

return
end subroutine

!==========================================================================
subroutine chk_os2008_table2
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE 2 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

implicit none
real*8 ptriple, gltriple, gvtriple, dl, dv

!ptriple = 611.6547710078944d0
ptriple = liq_vap_vapourpressure_si(273.16d0)


dl = 1d0/liq_g_si(0, 1, 273.16d0, ptriple)
dv = 1d0/vap_g_si(0, 1, 273.16d0, ptriple)
gltriple = flu_f_si(0, 0, 273.16d0, dl) + ptriple/dl
gvtriple = flu_f_si(0, 0, 273.16d0, dv) + ptriple/dv
                                            
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE 2 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) "PROPERTIES OF WATER AND VAPOUR AT T_TRIPLE = 273.16 K, P =", ptriple," Pa"
write(*,*) "COMPUTED FROM THERMODYNAMIC POTENTIALS OF THE SEAWATER LIBRARY"
write(*,*) "PTRIPLE_TABLE = 611.6547710078944d0"
write(*,*) "PTRIPLE_CALCULATED = ", ptriple
write(*,*) " "

write(*,*) "QUANTITY       TABLE VALUE          COMPUTED VALUE         &
  &   UNITS"

write(*,'(" pl             611.6548         ",2x,e24.16,4x,"Pa    ")') ptriple
write(*,'(" pv             611.6547710079   ",2x,e24.16,4x,"Pa    ")') ptriple
write(*,'(" gl             0.6117817        ",2x,e24.16,4x,"J KG-1")') gltriple
write(*,'(" gv             0.6117817        ",2x,e24.16,4x,"J KG-1")') gvtriple
write(*,'(" rhol           999.7925200316   ",2x,e24.16,4x,"KG M-3")') dl
write(*,'(" rhov           4.854575724779e-3",2x,e24.16,4x,"KG M-3")') dv
write(*,*)' '

end subroutine

!==========================================================================
subroutine chk_os2008_table3
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE 3 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/


real*8 g, gt, gp, gtt, gtp, gpp
real*8 t, p, rho, temp
integer icheck
          
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE 3 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) "TABLE 3 OF OS2008"
write(*,*) "FUNCTION VALUES FOR LIQUID, VAPOUR AND ICE AT THE TRIPLE POINT"
write(*,*) "T = 273.16 K, P = 611.654771007894 Pa"
write(*,*) " "
write(*,*) "QUANTITY  LIQUID               VAPOUR               ICE                  UNITS"
write(*,'(" g         6.117817e-1          6.1178170e-1         6.11781703e-1        J KG-1")')
write(*,'(" dg/dt     5.e-11              -9.15549340929856e+3  1.22069433939648e+3  J KG-1 K-1")')
write(*,'(" dg/dp     1.00020752302525e-3  2.05991224917108e+2  1.09085812736669e-3  M3 KG-1")')
write(*,'(" d2g/dt2  -1.5448497277676e+1  -6.89834540986158    -7.67602985875191     J KG-1 K-2")')
write(*,'(" d2g/dtdp -6.797921518e-8       7.5819292780093e-1   1.74387964700076e-7  M3 KG-1 K-1")')
write(*,'(" d2g/dp2  -5.09062317110e-13   -3.36992433127456e-1 -1.28495941571693e-13 M3 KG-1 Pa-1")')
write(*,'(" h         6.117817e-1          2.50091519146570e+6 -3.3344425396784e+5   J KG-1")')
write(*,'(" f         5.e-9               -1.2599490372461e+5  -5.5446875e-2         J KG-1")')
write(*,'(" u         2.e-9                2.37491967595939e+6 -3.3344492119642e+5   J KG-1")')
write(*,'(" sigma    -5.e-11               9.15549340929856e+3 -1.22069433939648e+3  J KG-1 K-1")')
write(*,'(" rho       9.997925200316e+2    4.8545757247786e-3   9.16709492199488e+2  KG M-3")')
write(*,'(" cp        4.21991151637e+3     1.88435203215779e+3  2.09678431621667e+3  J KG-1 K-1")')
write(*,'(" alpha    -6.7965110856e-5      3.68070498200120e-3  1.59863102565909e-4  K-1")')
write(*,'(" kappa_t   5.08956696877e-10    1.63595528529462e-3  1.17793449347882e-10 Pa-1")') 
write(*,'(" kappa_s   5.08657624753e-10    1.2314112503315e-3   1.14161597778763e-10 Pa-1")') 
write(*,*)" "


999 continue
#ifdef interactive
write(*,*)'enter the column to check; 0 to exit, 999 for all'
read(*,*) icheck
if(icheck == 0) return
#else
icheck = 999
#endif


if(icheck == 1 .or. icheck == 999)then
  t = 273.16d0
  p = 611.654771007894d0
  rho=liq_density_si(t,p)

  g =  flu_f_si(0, 0, t, rho) + p/rho
  gt =  flu_f_si(1, 0, t, rho)
  gp =  1d0/rho

  temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
  temp = flu_f_si(2,0,t,rho)-rho*flu_f_si(1,1,t,rho)**2/temp
  gtt = temp

  temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
  temp = -1d0/(rho**3*temp)
  gpp = temp

  temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
  temp = flu_f_si(1,1,t,rho)/(rho*temp)
  gtp = temp

  write(*,*) " "
  write(*,*) "FUNCTION VALUES FOR LIQUID WATER AT THE TRIPLE POINT"
  write(*,*) "T = 273.16 K, P = 611.654771007894 Pa"
  write(*,*) " "
  write(*,*) "QUANTITY  TABLE VALUE           CALCULATED VALUE         UNITS"
  write(*,'(" g         6.117817e-1        ",2x,es22.15,4x,"J KG-1")') g
  write(*,'(" dg/dt     5.e-11             ",2x,es22.15,4x,"J KG-1 K-1")') gt
  write(*,'(" dg/dp     1.00020752302525e-3",2x,es22.15,4x,"M3 KG-1")') gp
  write(*,'(" d2g/dt2  -1.5448497277676e+1 ",2x,es22.15,4x,"J KG-1 K-2")') gtt
  write(*,'(" d2g/dtdp -6.797921518e-8     ",2x,es22.15,4x,"M3 KG-1 K-1")') gtp
  write(*,'(" d2g/dp2  -5.09062317110e-13  ",2x,es22.15,4x,"M3 KG-1 Pa-1")') gpp
  write(*,'(" h         6.117817e-1        ",2x,es22.15,4x,"J KG-1")') g-t*gt
  write(*,'(" f         5.e-9              ",2x,es22.15,4x,"J KG-1")') g-p*gp
  write(*,'(" u         2.e-9              ",2x,es22.15,4x,"J KG-1")') g-t*gt-p*gp
  write(*,'(" sigma    -5.e-11             ",2x,es22.15,4x,"J KG-1 K-1")') -gt
  write(*,'(" rho       9.997925200316e+2  ",2x,es22.15,4x,"KG M-3")') rho
  write(*,'(" cp        4.21991151637e+3   ",2x,es22.15,4x,"J KG-1 K-1")') -t*gtt
  write(*,'(" alpha    -6.7965110856e-5    ",2x,es22.15,4x,"K-1")') gtp/gp
  write(*,'(" kappa_t   5.08956696877e-10  ",2x,es22.15,4x,"Pa-1")') -gpp/gp
  write(*,'(" kappa_s   5.08657624753e-10  ",2x,es22.15,4x,"Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)
 
endif
          
if(icheck == 2 .or. icheck == 999)then
t = 273.16d0
p = 611.654771007894d0
rho=vap_density_si(t,p)

g =  flu_f_si(0, 0, t, rho) + p/rho
gt =  flu_f_si(1, 0, t, rho)
gp =  1d0/rho

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = flu_f_si(2,0,t,rho)-rho*flu_f_si(1,1,t,rho)**2/temp
gtt = temp

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = -1d0/(rho**3*temp)
gpp = temp

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = flu_f_si(1,1,t,rho)/(rho*temp)
gtp = temp

write(*,*) " "
write(*,*) " "
write(*,*) "FUNCTION VALUES FOR WATER VAPOUR AT THE TRIPLE POINT"
write(*,*) "T = 273.16 K, P = 611.654771007894 Pa"
write(*,*) " "
write(*,*) "QUANTITY  TABLE VALUE            CALCULATED VALUE         UNITS"
write(*,'(" g         6.1178170e-1        ",2x,es22.15,4x,"J KG-1")') g
write(*,'(" dg/dt    -9.15549340929856e+3 ",2x,es22.15,4x,"J KG-1 K-1")') gt
write(*,'(" dg/dp     2.05991224917108e+2 ",2x,es22.15,4x,"M3 KG-1")') gp
write(*,'(" d2g/dt2  -6.89834540986158    ",2x,es22.15,4x,"J KG-1 K-2")') gtt
write(*,'(" d2g/dtdp  7.5819292780093e-1  ",2x,es22.15,4x,"M3 KG-1 K-1")') gtp
write(*,'(" d2g/dp2  -3.36992433127456e-1 ",2x,es22.15,4x,"M3 KG-1 Pa-1")') gpp
write(*,'(" h         2.50091519146570e+6 ",2x,es22.15,4x,"J KG-1")') g-t*gt
write(*,'(" f        -1.2599490372461e+5  ",2x,es22.15,4x,"J KG-1")') g-p*gp
write(*,'(" u         2.37491967595939e+6 ",2x,es22.15,4x,"J KG-1")') g-t*gt-p*gp
write(*,'(" sigma     9.15549340929856e+3 ",2x,es22.15,4x,"J KG-1 K-1")') -gt
write(*,'(" rho       4.8545757247786e-3  ",2x,es22.15,4x,"KG M-3")') rho
write(*,'(" cp        1.88435203215779e+3 ",2x,es22.15,4x,"J KG-1 K-1")') -t*gtt
write(*,'(" alpha     3.68070498200120e-3 ",2x,es22.15,4x,"K-1")') gtp/gp
write(*,'(" kappa_t   1.63595528529462e-3 ",2x,es22.15,4x,"Pa-1")') -gpp/gp
write(*,'(" kappa_s   1.2314112503315e-3  ",2x,es22.15,4x,"Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)
 
endif

if(icheck == 3 .or. icheck == 999)then
t = 273.16d0
p = 611.654771007894d0

g = ice_g_si(0, 0, t, p)
gt = ice_g_si(1, 0, t, p)
gp = ice_g_si(0, 1, t, p)
gtt = ice_g_si(2, 0, t, p)
gtp = ice_g_si(1, 1, t, p)
gpp = ice_g_si(0, 2, t, p)

write(*,*) " "
write(*,*) " "
write(*,*) "FUNCTION VALUES FOR ICE AT THE TRIPLE POINT"
write(*,*) "AS GIVEN BY IAPWS06 WITH G00 CORRECTED"
write(*,*) "T = 273.16 K, P = 611.654771007894 Pa"
write(*,*) " "
write(*,*) "QUANTITY  TABLE VALUE             CALCULATED VALUE         UNITS"
write(*,'(" g         6.11781703e-1        ",2x,es22.15,4x,"J KG-1")') g
write(*,'(" dg/dt     1.22069433939648e+3  ",2x,es22.15,4x,"J KG-1 K-1")') gt
write(*,'(" dg/dp     1.09085812736669e-3  ",2x,es22.15,4x,"M3 KG-1")') gp
write(*,'(" d2g/dt2  -7.67602985875191     ",2x,es22.15,4x,"J KG-1 K-2")') gtt
write(*,'(" d2g/dp2  -1.28495941571693e-13 ",2x,es22.15,4x,"M3 KG-1 Pa-1")') gpp
write(*,'(" d2g/dpdt  1.74387964700076e-7  ",2x,es22.15,4x,"M3 KG-1 K-1")') gtp
write(*,'(" h        -3.3344425396784e+5   ",2x,es22.15,4x,"J KG-1")') g-t*gt
write(*,'(" f        -5.5446875e-2         ",2x,es22.15,4x,"J KG-1")') g-p*gp
write(*,'(" u        -3.3344492119642e+5   ",2x,es22.15,4x,"J KG-1")') g-t*gt-p*gp
write(*,'(" sigma    -1.22069433939648e+3  ",2x,es22.15,4x,"J KG-1 K-1")') -gt
write(*,'(" rho       9.16709492199488e+2  ",2x,es22.15,4x,"KG M-3")') 1d0/gp
write(*,'(" cp        2.09678431621667e+3  ",2x,es22.15,4x,"J KG-1 K-1")') -t*gtt
write(*,'(" alpha     1.59863102565909e-4  ",2x,es22.15,4x,"K-1")') gtp/gp
write(*,'(" kappa_t   1.17793449347882e-10 ",2x,es22.15,4x,"Pa-1")') -gpp/gp
write(*,'(" kappa_s   1.14161597778763e-10 ",2x,es22.15,4x,"Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)

endif

if(icheck /= 999) go to 999

end subroutine


!********************************************************************************
subroutine chk_os2008_tablea1
!********************************************************************************

!This table is used to compare computed results with values given in Table A1 of OS2008.
!It calls chk_iapws95_table6 to create the table output in order to avoid the need to make
!flu_phi0 and flu_phir public which would otherwise conflict with uniqueness of level 1 routines.

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

call chk_iapws95_table6

end subroutine


!==========================================================================
subroutine chk_os2008_tablea2
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE A2 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

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
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A2 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) "UPDATED FUNCTION VALUES FOR TABLE 7 OF IAPWS-95:"
write(*,*) "THERMODYNAMIC PROPERTY OF H2O IN THE SINGLE-PHASE REGION FOR SELECTED VALUES OF T AND RHO"
write(*,*) " "

write(*,'(A)') " ROW T     RHO            P              CV             C              SIGMA    "
write(*,'(A)') "     K     KG/M3          MPa            KJ/(KG K)      M/S            KJ/(KG K)"
write(*,'(A)') " "

write(*,'(A)') " 1   300  .9965560e+3    .992418352e-1  .413018112e+1  .150151914e+4  .393062643"
write(*,'(A)') " 2   300  .1005308e+4    .200022515e+2  .406798347e+1  .153492501e+4  .387405401"
write(*,'(A)') " 3   300  .1188202e+4    .700004704e+3  .346135580e+1  .244357992e+4  .132609616"
write(*,'(A)') " 4   500  .4350000       .999679423e-1  .150817541e+1  .548314253e+3  .794488271e+1"
write(*,'(A)') " 5   500  .4532000e+1    .999938125     .166991025e+1  .535739001e+3  .682502725e+1"
write(*,'(A)') " 6   500  .8380250e+3    .100003858e+2  .322106219e+1  .127128441e+4  .256690919e+1"
write(*,'(A)') " 7   500  .1084564e+4    .700000405e+3  .307437693e+1  .241200877e+4  .203237509e+1"
write(*,'(A)') " 8   647  .3580000e+3    .220384756e+2  .618315728e+1  .252145078e+3  .432092307e+1"
write(*,'(A)') " 9   900  .2410000       .100062559     .175890657e+1  .724027147e+3  .916653194e+1"
write(*,'(A)') "10   900  .5261500e+2    .200000690e+2  .193510526e+1  .698445674e+3  .659070225e+1"
write(*,'(A)') "11   900  .8707690e+3    .700000006e+3  .266422350e+1  .201933608e+4  .417223802e+1"

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

!==========================================================================
subroutine chk_os2008_tablea3
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE A3 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

implicit none
integer it
real*8 pl, pv, rhol, rhov
real*8 t, zt(3)
real*8 g, gt, gp
real*8 temp

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A3 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) "PROPERTIES OF H2O IN THE TWO-PHASE REGION AT SELECTED VALUES OF TEMPERATURE" 

zt(1) = 275d0
zt(2) = 450d0
zt(3) = 625d0

do it = 1,3

  t = zt(it)

#ifdef interactive
  write(*,*)' '
  write(*,*)'ENTER 1 TO SHOW RESULTS FOR T = ',t,'; 0 TO EXIT'
  read(*,*)icheck
  if(icheck /= 1) go to 999
#endif

  !if(it == 1)rhol = 9.99887406119850d2; rhov = 5.50664918504123d-3  
  !if(it == 2)rhol = 8.90341249761673d2; rhov = 4.81200360125671d0
  !if(it == 3)rhol = 5.67090385146353d2; rhov = 1.18290280451157d2
  temp = liq_vap_vapourpressure_si(t)
  rhol = liq_density_si(t,temp)
  rhov = vap_density_si(t,temp)

  pl=flu_pressure_si(t, rhol)
  g =  flu_f_si(0, 0, t, rhol) + pl/rhol
  gt =  flu_f_si(1, 0, t, rhol)
  gp =  1d0/rhol

  write(*,*) " "
  write(*,*) "PROPERTIES AT T =", zt(it), " K, P =", pl," Pa"
  write(*,*) " "

  write(*,*) "QUANTITY       TABLE VALUE           COMPUTED VALUE          &
  &  UNITS"

  if(it == 1)then
    write(*,'(" pl             6.9845117e8       ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           999.887406120     ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             7.759722016e3     ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         2.83094669595e1   ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 2)then
    write(*,'(" pl             9.32203563628e5   ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           890.34124976167   ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             7.4916158501217e5 ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         2.10865844688447e3",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 3)then
    write(*,'(" pl             1.6908269318578e7 ",2x,e24.16,4x,"Pa    ")') pl
    write(*,'(" rhol           5.670903851464e2  ",2x,e24.16,4x,"KG M-3")') rhol
    write(*,'(" hl             1.6862697594697e6 ",2x,e24.16,4x,"J KG-1")') g- t * gt
    write(*,'(" sigmal         3.8019468301114e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  endif

  write(*,*)' '

  pv=flu_pressure_si(t, rhov)
  g = flu_f_si(0, 0, t, rhov) + pv/rhov
  gt =  flu_f_si(1, 0, t, rhov)

  if(it == 1)then
    write(*,'(" pv             6.98451166701e8   ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           5.506649185041e-3 ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.5042899500405e6 ",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         9.1066012052322e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 2)then
    write(*,'(" pv             9.322035636282e5  ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           4.812003601257    ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.77441077988962e6",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         6.6092122132788e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  elseif(it == 3)then
    write(*,'(" pv             1.6908269318578e7 ",2x,e24.16,4x,"Pa    ")') pv
    write(*,'(" rhov           1.182902804512e2  ",2x,e24.16,4x,"KG M-3")') rhov
    write(*,'(" hv             2.5507162456235e6 ",2x,e24.16,4x,"J KG-1")') g - t * gt
    write(*,'(" sigmav         5.1850612079574e3 ",2x,e24.16,4x,"J KG-1 K-1")') -gt
  endif

go to 999
999 write(*,*)' '
 
enddo

end subroutine

!==========================================================================
subroutine chk_os2008_tablea4
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE A4 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

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
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A4 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*)'UPDATED FUNCTION VALUES FOR TABLE 6 OF IAPWS-06:'
write(*,*)'PROPERTIES AT THE TRIPLE POINT, THE NORMAL PRESSURE &
 &MELTING POINT, '
write(*,*)'AND AT t = 100 K, p = 100 MPa'
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

  write(*,*)'QUANTITY     VALUE AT            VALUE AT            VALUE AT          UNITS'
  write(*,*)'             t = 273.16 K        t = 273.152519 K    t = 100 K'
  write(*,*)'             p = 611.657 Pa      p = 101325 Pa       p = 100 MPa'
  write(*,*) " "

  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' g          ',zg(1),       ' ', zg(2),       ' ', zg(3),       ' J/KG'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')'(dg/dp)_t   ',zdgdp_t(1),  ' ', zdgdp_t(2),  ' ', zdgdp_t(3),  ' M3/KG'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')'(dg/dt)_p   ',zdgdt_p(1),  ' ', zdgdt_p(2),  ' ', zdgdt_p(3),  ' J/KG/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')'(d2g/dp2)_t ',zd2gdp2_t(1),' ', zd2gdp2_t(2),' ', zd2gdp2_t(3),' M3/KG/PA'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' d2g/dpdt   ',zd2gdpdt(1), ' ', zd2gdpdt(2), ' ', zd2gdpdt(3), ' M3/KG/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')'(d2g/dt2)_p ',zd2gdt2_p(1),' ', zd2gdt2_p(2),' ', zd2gdt2_p(3),' J/KG/K^2'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' h          ',zh(1),       ' ', zh(2),       ' ', zh(3),       ' J/KG1'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' f          ',zf(1),       ' ', zf(2),       ' ', zf(3),       ' J/KG'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' u          ',zu(1),       ' ', zu(2),       ' ', zu(3),       ' J/KG'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' sigma      ',zsigma(1),   ' ', zsigma(2),   ' ', zsigma(3),   ' J/KG/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' cp         ',zcp(1),      ' ', zcp(2),      ' ', zcp(3),      ' J/KG/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' rho        ',zrho(1),     ' ', zrho(2),     ' ', zrho(3),     ' KG/M3'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' alpha      ',zalpha(1),   ' ', zalpha(2),   ' ', zalpha(3),   ' 1/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' beta       ',zbeta(1),    ' ', zbeta(2),    ' ', zbeta(3),    ' Pa/K'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' kappa_t    ',zkappa_t(1), ' ', zkappa_t(2), ' ', zkappa_t(3), ' 1/Pa'
  write(*,'(A,E19.12,A,E19.12,A,E19.12,A)')' kappa_s    ',zkappa_s(1), ' ', zkappa_s(2), ' ', zkappa_s(3), ' 1/Pa'

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
   &   UNITS"

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

!==========================================================================
subroutine chk_os2008_tablea5
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE A5 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A5 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) 'PROPERTIES AT S = 0.03516504, T = 273.15 K, P = 101325 Pa'
write(*,*) " "

write(*,*)'QUANTITY  WATER PART           SALINE PART          SEAWATER VALUE      UNITS'

write(*,*)'g         1.01342742e2        -1.0134274172939e2    4.e-9               J/KG'
write(*,*)'dg/ds     0.0                  6.39974067312299e4   6.39974067312299e4  J/KG'
write(*,*)'dg/dt     1.47643376e-1       -1.47643376346e-1    -6.e-11              J/KG/K'
write(*,*)'dg/dp     1.00015693912169e-3 -2.74957224268433e-5  9.7266121669485e-4  M3/KG'
write(*,*)'d2g/dsdp  0.0                 -7.59615411515309e-4 -7.59615411515309e-4 M3/KG'
write(*,*)'d2g/dt2  -1.54473542320e1      8.5286115117592e-1  -1.459449308080e1    J/KG/K2'
write(*,*)'d2g/dtdp -6.777003179e-8       1.1928678741396e-7   5.1516755627e-8     M3/KG/K'
write(*,*)'d2g/dp2  -5.08928894643e-13    5.8153517233288e-14 -4.507753774102e-13  M3/KG/Pa'
write(*,*)'h         6.10139535e1        -6.10139534804e1      2.e-8               J/KG'
write(*,*)'f         1.8399e-3           -9.8556737654491e1   -9.855489778e1       J/KG'
write(*,*)'u        -4.03269484e1        -5.82279494055e1     -9.85548978e1        J/KG'
write(*,*)'sigma    -1.47643376e-1        1.47643376346e-1     6.e-11              J/KG/K'
write(*,*)'rho       9.9984308550433e2   -9d99                 1.02810719995401e3  KG/M3'
write(*,*)'cp        4.21944480846e3     -2.3295902344370e2    3.98648578502e3     J/KG/K'
write(*,*)'c         1.40238253109e3     -9d99                 1.44900246362e3     M/S'
write(*,*)'muw       1.01342742e2        -2.35181410932936e3  -2.250471368e3       J/KG'
write(*,*)' '

s = 0.03516504d0
t = 273.15d0
p = 101325d0

zg(1) = 1.01342742d2  	
zdgds(1) = 0d0
zdgdt(1) = 1.47643376d-1  
zdgdp(1) = 1.00015693912169d-3
zd2gdsdp(1) = 0d0
zd2gdt2(1) = -1.54473542320d1
zd2gdtdp(1) = -6.777003179d-8
zd2gdp2(1) = -5.08928894643d-13	
zh(1) = 6.10139535d1
zf(1) = 1.8399d-3
zu(1) = -4.03269484d1	
zsigma(1) = -1.47643376d-1
zrho(1) = 9.9984308550433d2
zcp(1) = 4.21944480846d3
zc(1) = 1.40238253109d3	
zmuw(1) = 1.01342742d2

zg(2) = -1.0134274172939d2      
zdgds(2) = 6.39974067312299d4 
zdgdt(2) = -1.47643376346d-1
zdgdp(2) = -2.74957224268433d-5	
zd2gdsdp(2) = -7.59615411515309d-4
zd2gdt2(2) = 8.5286115117592d-1
zd2gdtdp(2) = 1.1928678741396d-7	
zd2gdp2(2) = 5.8153517233288d-14
zh(2) = -6.10139534804d1
zf(2) = -9.8556737654491d1	
zu(2) = -5.82279494055d1
zsigma(2) = 1.47643376346d-1
zrho(2) = -9.d99	       
zcp(2) = -2.3295902344370d2  
zc(2) = -9.d99
zmuw(2) = -2.35181410932936d3

zg(3) = 4.d-9
zdgds(3) = 6.39974067312299d4
zdgdt(3) = -6.d-11          
zdgdp(3) = 9.7266121669485d-4
zd2gdsdp(3) = -7.59615411515309d-4 
zd2gdt2(3) = -1.459449308080d1   
zd2gdtdp(3) = 5.1516755627d-8    
zd2gdp2(3) = -4.507753774102d-13
zh(3) = 2.d-8             
zf(3) = -9.855489778d1   
zu(3) = -9.85548978d1   
zsigma(3) = 6.d-11         
zrho(3) = 1.02810719995401d3
zcp(3) = 3.98648578502d3  
zc(3) = 1.44900246362d3 
zmuw(3) = -2.250471368d3 
             
icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A <1 OR >3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol ==0 .or. icol >= 4) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif

write(*,*) " "
write(*,*) " "
write(*,*) "COMPARISON WITH FUNCTION VALUES AS GIVEN IN TABLE A5, COLUMN", icol, &
           "OF OS2008"
write(*,*) "PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
 &       UNITS"

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q

go to 100

end subroutine

!==========================================================================
subroutine chk_os2008_tablea6
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE A6 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A6 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) 'PROPERTIES AT S = 0.1, T = 353 K, P = 101325 Pa'
write(*,*) " "

write(*,*)'QUANTITY  WATER PART           SALINE PART          SEAWATER VALUE     UNITS'
write(*,*)'g        -4.46114968996e4      1.5087174003705e4   -2.95243228959e4    J/KG'
write(*,*)'dg/ds     0.0                  2.51957275851413e5   2.51957275851413e5 J/KG'
write(*,*)'dg/dt    -1.0737599318875e3    1.56230907404291e2  -9.175290244832e2   J/KG/K'
write(*,*)'dg/dp     1.02892955635611e-3 -5.7922728577126e-5   9.7100682777898e-4 M3/KG'
write(*,*)'d2g/dsdp  0.0                 -3.0595780244234e-4  -3.0595780244234e-4 M3/KG'
write(*,*)'d2g/dt2  -1.1888500004755e1    1.27922649315507e0  -1.0609273511600e1  J/KG/K2'
write(*,*)'d2g/dtdp  6.59051552339e-7     8.0306159575153e-7   1.462113148091e-6  M3/KG/K'
write(*,*)'d2g/dp2  -4.746728193611e-13   2.1308615424374e-13 -2.61586665117e-13  M3/KG/Pa'
write(*,*)'h         3.344257590567e5    -4.006233631001e4     2.94363422747e5    J/KG'
write(*,*)'f        -4.47157531869e4      1.5093043024178e4   -2.96227101627e4    J/KG'
write(*,*)'u         3.343215027694e5    -4.0056467289536e4    2.942650354799e5   J/KG'
write(*,*)'sigma     1.0737599318875e3   -1.56230907404291e2   9.175290244832e2   J/KG/K'
write(*,*)'rho       9.7188383191308e2   -9d99                 1.02985887574790e3 KG/M3'
write(*,*)'cp        4.1966405016784e3   -4.51566952083741e2   3.745073549595e3   J/KG/K'
write(*,*)'c         1.5544629665347e3   -9d99                 3.9612783529e3     M/S'
write(*,*)'muw      -4.46114968996e4     -1.01085535814360e4  -5.47200504810e4    J/KG'
write(*,*)' '

s = 0.1d0
t = 353d0
p = 101325d0

zg(1) = -4.46114968996d4
zdgds(1) = 0d0
zdgdt(1) = -1.0737599318875d3
zdgdp(1) = 1.02892955635611d-3
zd2gdsdp(1) = 0d0
zd2gdt2(1) = -1.1888500004755d1
zd2gdtdp(1) = 6.59051552339d-7
zd2gdp2(1) = -4.746728193611d-13
zh(1) = 3.344257590567d5
zf(1) = -4.47157531869d4
zu(1) = 3.343215027694d5
zsigma(1) = 1.0737599318875d3
zrho(1) = 9.7188383191308d2
zcp(1) = 4.1966405016784d3
zc(1) = 1.5544629665347d3
zmuw(1) = -4.46114968996d4

zg(2) = 1.5087174003705d4
zdgds(2) = 2.51957275851413d5
zdgdt(2) = 1.56230907404291d2
zdgdp(2) = -5.7922728577126d-5
zd2gdsdp(2) = -3.0595780244234d-4
zd2gdt2(2) = 1.27922649315507d0
zd2gdtdp(2) = 8.0306159575153d-7
zd2gdp2(2) = 2.1308615424374d-13
zh(2) = -4.006233631001d4
zf(2) = 1.5093043024178d4
zu(2) = -4.0056467289536d4
zsigma(2) = -1.56230907404291d2
zrho(2) = -9d99	           
zcp(2) = -4.51566952083741d2
zc(2) = -9d99	           
zmuw(2) = -1.01085535814360d4

zg(3) = -2.95243228959d4 
zdgds(3) = 2.51957275851413d5
zdgdt(3) = -9.175290244832d2
zdgdp(3) = 9.7100682777898d-4
zd2gdsdp(3) = -3.0595780244234d-4
zd2gdt2(3) = -1.0609273511600d1
zd2gdtdp(3) = 1.462113148091d-6
zd2gdp2(3) = -2.61586665117d-13
zh(3) = 2.94363422747d5
zf(3) = -2.96227101627d4
zu(3) = 2.942650354799d5
zsigma(3) = 9.175290244832d2
zrho(3) = 1.02985887574790d3
zcp(3) = 3.745073549595d3
zc(3) = 3.9612783529d3
zmuw(3) = -5.47200504810d4
             
icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A VALUE < 1 or > 3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol < 1 .or. icol > 3) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif

write(*,*) " "
write(*,*) " "
write(*,*) "FUNCTION VALUES AS GIVEN IN TABLE A6, COLUMN", icol, &
           " OF OS2008"
write(*,*) "PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
&       UNITS"

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q
 
go to 100

end subroutine

!==========================================================================
subroutine chk_os2008_tablea7
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE 3 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

implicit none
integer icol
real*8 s, t, p, rho
real*8 zg(3), zdgds(3), zdgdt(3), zdgdp(3), zd2gdsdp(3), zd2gdt2(3)
real*8 zd2gdtdp(3), zd2gdp2(3), zh(3), zf(3), zu(3), zsigma(3)
real*8 zrho(3), zcp(3), zc(3), zmuw(3)
real*8 g, gs, gt, gp, gss, gtt, gpp, gst, gsp, gtp, q

write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A7 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*)'PROPERTIES AT S = 0.1, T = 353 K, P = 10**8 Pa'
write(*,*) " "

write(*,*)'QUANTITY  WATER PART           SALINE PART           SEAWATER VALUE      UNITS'

write(*,*)'g         9.773038621954e+4   -2.60093050730637e+3   9.512945571223e+4   J/KG'
write(*,*)'dg/ds     0.0                 -5.4586158064880e+3   -5.4586158064880e+3  J/KG'
write(*,*)'dg/dt     8.5146650206e        7.5404568488117e+0    1.60551218694e+1    J/KG/K'
write(*,*)'dg/dp     9.5668332915351e-4  -2.2912384179113e-5    9.3377094497440e-4  M3/KG'
write(*,*)'d2g/dsdp  0.0                 -6.4075761854575e-4   -6.40757618545748e-4 M3/KG'
write(*,*)'d2g/dt2  -1.429698733876e+1    4.88076973942251e-1  -1.380891036482e+1   J/KG/K2'
write(*,*)'d2g/dtdp  1.99079570803e-7     4.6628441224121e-8    2.45708012027e-7    M3/KG/K'
write(*,*)'d2g/dp2  -3.715308894234e-13   3.57345735845327e-14 -3.35796315839e-13   M3/KG/Pa'
write(*,*)'h         9.54046054692e+4    -4.6606062955593e+3    9.07439991736e+4    J/KG'
write(*,*)'f         2.06205330419e+3    -3.0969208939506e+2    1.75236121479e+3    J/KG'
write(*,*)'u        -2.637274462e+2      -2.3693678776480e+3   -2.6330953238e+3     J/KG' 
write(*,*)'sigma    -8.5146650206e+0     -7.5404568488117      -1.6055121869e+1     J/KG/K'
write(*,*)'rho       1.04527796139692e+3 -9d99                  1.0709264465574e+3  KG/M3'
write(*,*)'cp        3.905222091582e+3   -1.33318225432326e+2   3.77190386615e+3    J/KG/K'
write(*,*)'c         1.575422398486e+3   -9d99                  1.621989976499e+3   M/S'
write(*,*)'muw       9.77303862195e+4    -2.4089780641266e+3    9.532140815541e+4   J/KG'
write(*,*)' '

s = 0.03516504d0
t = 273.15d0
p = 1d8

zg(1) = 9.773038621954d+4
zdgds(1) = 0.0	     
zdgdt(1) = 8.5146650206d+0
zdgdp(1) = 9.5668332915351d-4
zd2gdsdp(1) = 0.0	      
zd2gdt2(1) = -1.429698733876d+1
zd2gdtdp(1) = 1.99079570803d-7
zd2gdp2(1) = -3.715308894234d-13
zh(1) = 9.54046054692d+4
zf(1) = 2.06205330419d+3
zu(1) = -2.637274462d+2
zsigma(1) = -8.5146650206d+0
zrho(1) = 1.04527796139692d+3
zcp(1) = 3.905222091582d+3
zc(1) = 1.575422398486d+3
zmuw(1) = 9.77303862195d+4

zg(2) = -2.60093050730637d+3
zdgds(2) = -5.4586158064880d+3
zdgdt(2) = 7.5404568488117d+0
zdgdp(2) = -2.2912384179113d-5
zd2gdsdp(2) = -6.4075761854575d-4
zd2gdt2(2) = 4.88076973942251d-1
zd2gdtdp(2) = 4.6628441224121d-8
zd2gdp2(2) = 3.57345735845327d-14
zh(2) = -4.6606062955593d+3
zf(2) = -3.0969208939506d+2
zu(2) = -2.3693678776480d+3
zsigma(2) = -7.5404568488117d+0
zrho(2) = -9d99
zcp(2) = -1.33318225432326d+2
zc(2) = -9d99
zmuw(2) = -2.4089780641266d+3

zg(3) = 9.512945571223d+4
zdgds(3) = -5.4586158064880d+3
zdgdt(3) = 1.60551218694d+1
zdgdp(3) = 9.3377094497440d-4
zd2gdsdp(3) = -6.40757618545748d-4
zd2gdt2(3) = -1.380891036482d+1
zd2gdtdp(3) = 2.45708012027d-7
zd2gdp2(3) = -3.35796315839d-13
zh(3) = 9.07439991736d+4
zf(3) = 1.75236121479d+3
zu(3) = -2.6330953238d+3
zsigma(3) = -1.6055121869d+1
zrho(3) = 1.0709264465574d+3
zcp(3) = 3.77190386615d+3
zc(3) = 1.621989976499d+3
zmuw(3) = 9.532140815541d+4

icol = 0
100 continue

#ifdef interactive
write(*,*)' '
write(*,*)'ENTER THE COLUMN NUMBER THAT YOU WOULD LIKE TO CHECK'
write(*,*)'A VALUE < 1 or > 3 WILL CAUSE AN EXIT'
read(*,*) icol
#else
icol = icol + 1
#endif

if(icol < 1 .or. icol > 3) return

if(icol == 1) then
  g =   liq_g_si(0, 0, t, p)
  gs =  0d0
  gt =  liq_g_si(1, 0, t, p)
  gp =  liq_g_si(0, 1, t, p)
  gss = 0d0
  gtt = liq_g_si(2, 0, t, p)
  gpp = liq_g_si(0, 2, t, p)
  gst = 0d0
  gsp = 0d0
  gtp = liq_g_si(1, 1, t, p)
elseif(icol == 2) then
  g =   sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  sal_g_si(0, 1, 0, s, t, p)
  gp =  sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = sal_g_si(0, 2, 0, s, t, p)
  gpp = sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = sal_g_si(0, 1, 1, s, t, p)
elseif(icol == 3) then
  g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
  gs =  sal_g_si(1, 0, 0, s, t, p)
  gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
  gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
  gss = sal_g_si(2, 0, 0, s, t, p)
  gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
  gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
  gst = sal_g_si(1, 1, 0, s, t, p)
  gsp = sal_g_si(1, 0, 1, s, t, p)
  gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)
endif
if(icol == 2) then
  rho = -9d99
else
  rho = 1d0/gp
endif


write(*,*) " "
write(*,*) " "
write(*,*) "FUNCTION VALUES AS GIVEN IN TABLE A7, COLUMN", icol, &
           " OF OS2008"
write(*,*) "PROPERTIES AT S =", s,"G/KG, T =", t," K, P =", p," Pa"
write(*,*) " "

write(*,*) "QUANTITY          TABLE VALUE               CALCULATED VALUE  &
&       UNITS"

write(*,'(" g            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zg(icol), g
write(*,'("(dg/ds)_tp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgds(icol), gs
write(*,'("(dg/dt)_sp    ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zdgdt(icol), gt
write(*,'("(dg/dp)_st    ", 4x,e22.15,4x,e22.15,4x,"M3/KG")') zdgdp(icol), gp
write(*,'("(d2g/dsdp)_t  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdsdp(icol), gsp
write(*,'("(d2g/dt2)_p   ", 4x,e22.15,4x,e22.15,4x,"J/KG/K2")') zd2gdt2(icol), gtt
write(*,'("(d2g/dp2)_st  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/Pa")') zd2gdp2(icol), gpp
write(*,'("(d2g/dtdp)_s  ", 4x,e22.15,4x,e22.15,4x,"M3/KG/K")') zd2gdtdp(icol), gtp
write(*,'(" h            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zh(icol), g-t*gt
write(*,'(" f            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zf(icol), g-p*gp
write(*,'(" u            ", 4x,e22.15,4x,e22.15,4x,"J/KG")') zu(icol), g - t * gt - p * gp
write(*,'(" sigma        ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zsigma(icol), -gt
write(*,'(" rho          ", 4x,e22.15,4x,e22.15,4x,"KG/M3")') zrho(icol), rho
write(*,'(" cp           ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zcp(icol), -t*gtt
write(*,'(" muw          ", 4x,e22.15,4x,e22.15,4x,"J/KG/K")') zmuw(icol), g-s*gs

if (gtt < 0d0) then
  q = (gtp ** 2 - gtt * gpp) / gtt
  if (q > 0d0) then
    q = gp / q**0.5d0
  else
    write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 1'
    q = -9.9999999d99
  end if
else
  write(*,*)'ERROR IN DETERMINING SOUND SPEED : TYPE 2'
  q = -9.9999999d99
end if

write(*,'(" c            ", 4x,e22.15,4x,e22.15,4x,"J KG-1 K-1")') zc(icol), q
 
go to 100

end subroutine

!==========================================================================
subroutine chk_os2008_tablea8
!==========================================================================

!COMPARISON OF THE LIBRARY FUNCTION WITH CHECK VALUES
!COMPUTED IN QUADRUPLE PRECISION (128 BIT FLOATING POINT)
!AS PUBLISHED IN TABLE 3 OF

!R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER
!MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:
!A NEW STANDARD FOR OCEANOGRAPHY.
!OCEAN SCIENCE, 4, 275-291, 2008
!WWW.OCEAN-SCI.NET/4/275/2008/

integer icol
real*8 g, gt, gp, gtt, gtp, gpp, gs, gss, gst, gsp
real*8 s, t, p, rho, temp

          
icol = 0
999 CONTINUE

#ifdef interactive
write(*,*)" "
write(*,*)"TABLE A8 OF OS2008"
write(*,*)"ENTER THE COLUM TO CHECK; 0 TO EXIT; 999 FOR ALL"
write(*,*)"1 FOR LIQUID WATER AT THE STANDARD OCEAN STATE"
write(*,*)"2 FOR ICE AT THE STANDARD OCEAN STATE WITH S=0"
write(*,*)"3 FOR SEAWATER AT THE STANDARD OCEAN STATE"
read(*,*)icol
#else
icol = icol + 1
#endif

if(icol == 0 .or. icol == 4) return

if(icol == 1 .or. icol == 999) then
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) "****************************************************************"
write(*,*) " "
write(*,*) "Implementation of OS2008 in Fortran"
write(*,*) "for Publication in Ocean Science, 2009"
write(*,*) "D. Wright, BIO, VERSION ", version
write(*,*) " "
write(*,*) "COMPARISON FOR TABLE A8 OF OS2008"
write(*,*) "R. FEISTEL, D.G. WRIGHT, K. MIYAGAWA, A.H. HARVEY, J. HRUBY, D.R. JACKETT, T.J. MCDOUGALL, W.WAGNER"
write(*,*) "MUTUALLY CONSISTENT THERMODYNAMIC POTENTIALS FOR FLUID WATER, ICE AND SEAWATER:"
write(*,*) "A NEW STANDARD FOR OCEANOGRAPHY."
write(*,*) "OCEAN SCIENCE, 4, 275-291, 2008"
write(*,*) "WWW.OCEAN-SCI.NET/4/275/2008/"
write(*,*) " "
write(*,*) "FUNCTION VALUES FOR LIQUID WATER AT THE STANDARD OCEAN STATE"
write(*,*) "AS GIVEN BY IAPWS95 WITH N01, N02 CORRECTED"
write(*,*) "T = 273.15 K, P = 101325 Pa"

t = 273.15d0
p = 101325d0
rho=liq_density_si(t,p)

g =  flu_f_si(0, 0, t, rho) + p/rho
gt =  flu_f_si(1, 0, t, rho)
gp =  1d0/rho

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = flu_f_si(2,0,t,rho)-rho*flu_f_si(1,1,t,rho)**2/temp
gtt = temp

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = -1d0/(rho**3*temp)
gpp = temp

temp = 2d0*flu_f_si(0,1,t,rho)+rho*flu_f_si(0,2,t,rho)
temp = flu_f_si(1,1,t,rho)/(rho*temp)
gtp = temp

write(*,*)" "
write(*,*)"QUANTITY      TABLE VALUE          COMPUTED VALUE            UNITS"

write(*,'(" g             1.013427417e2      ",1x,e22.15,5x,"J KG-1")') g
write(*,'("(dg/dt)_p      1.476433763e-1     ",1x,e22.15,5x,"J KG-1 K-1")') gt
write(*,'("(dg/dp)_t      1.00015693912169e-3",1x,e22.15,5x,"J M3 KG-1")') gp
write(*,'("(d2g/dt2)_p   -1.544735423198e+1  ",1x,e22.15,5x,"J KG-1 K-2")') gtt
write(*,'(" d2g/dpdt     -6.777003179e-8     ",1x,e22.15,5x,"J M3 KG-1 K-1")') gtp
write(*,'("(d2g/dp2)_t   -5.089288946435e-13 ",1x,e22.15,5x,"J M3 KG-1 Pa-1")') gpp
write(*,'(" h             6.10139535e+1      ",1x,e22.15,5x,"J KG-1")') g-t*gt
write(*,'(" f             1.8399e-3          ",1x,e22.15,5x,"J KG-1")') g-p*gp
write(*,'(" u            -4.03269484e+1      ",1x,e22.15,5x,"J KG-1")') g-t*gt-p*gp
write(*,'(" sigma        -1.47643376e-1      ",1x,e22.15,5x,"J KG-1 K-1")') -gt
write(*,'(" rho           9.9984308550433e+2 ",1x,e22.15,5x,"KG M-3")') rho
write(*,'(" cp            4.21944480846e+3   ",1x,e22.15,5x,"J KG-1 K-1")') -t*gtt
write(*,'(" alpha        -6.7759397686e-5    ",1x,e22.15,5x,"J K-1")') gtp/gp
write(*,'(" kappa_t       5.08849036323e-10  ",1x,e22.15,5x,"J Pa-1")') -gpp/gp
write(*,'(" kappa_s       5.08551764928e-10  ",1x,e22.15,5x,"J Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)
write(*,'(" muw           1.013427417e+2     ",1x,e22.15,5x,"J/KG/K")') g
 
endif
 

if(icol == 2 .or. icol == 999) then
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "TABLE A8, COLUMN 2 OF OS2008"
write(*,*) " "
write(*,*) " FUNCTION VALUES FOR ICE AT THE STANDARD OCEAN STATE WITH S=0"
write(*,*) " AS GIVEN BY IAPWS06 WITH G00 CORRECTED"
write(*,*) " T = 273.15 K, P = 101325 Pa"

t = 273.15d0
p = 101325d0

g = ice_g_si(0, 0, t, p)
gt = ice_g_si(1, 0, t, p)
gp = ice_g_si(0, 1, t, p)
gtt = ice_g_si(2, 0, t, p)
gtp = ice_g_si(1, 1, t, p)
gpp = ice_g_si(0, 2, t, p)

write(*,*)" "
write(*,*)" QUANTITY     TABLE VALUE            COMPUTED VALUE            UNITS"

write(*,'(" g             9.8267598403e1       ",1x,e22.15,5x,"J KG-1")') g
write(*,'("(dg/dt)_p      1.22078866129995e3   ",1x,e22.15,5x,"J KG-1 K-1")') gt
write(*,'("(dg/dp)_t      1.09084344292644e-3  ",1x,e22.15,5x,"J M3 KG-1")') gp
write(*,'("(d2g/dt2)_p   -7.67598511156675     ",1x,e22.15,5x,"J KG-1 K-2")') gtt
write(*,'(" d2g/dpdt      1.74360824960850e-7  ",1x,e22.15,5x,"J M3 KG-1 K-1")') gtp
write(*,'("(d2g/dp2)_t   -1.28484824639762e-13 ",1x,e22.15,5x,"J M3 KG-1 Pa-1")') gpp
write(*,'(" h            -3.33360155235679e5   ",1x,e22.15,5x,"J KG-1")') g-t*gt
write(*,'(" f            -1.226211345e1        ",1x,e22.15,5x,"J KG-1")') g-p*gp
write(*,'(" u            -3.334706849475e5     ",1x,e22.15,5x,"J KG-1")') g-t*gt-p*gp
write(*,'(" sigma        -1.22078866129995e3   ",1x,e22.15,5x,"J KG-1 K-1")') -gt
write(*,'(" rho           9.16721832527382e2   ",1x,e22.15,5x,"KG M-3")') 1d0/gp
write(*,'(" cp            2.09669533322446e3   ",1x,e22.15,5x,"J KG-1 K-1")') -t*gtt
write(*,'(" alpha         1.59840374979096e-4  ",1x,e22.15,5x,"J K-1")') gtp/gp
write(*,'(" kappa_t       1.17784843895722e-10 ",1x,e22.15,5x,"J Pa-1")') -gpp/gp
write(*,'(" kappa_s       1.14154052637228e-10 ",1x,e22.15,5x,"J Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)
write(*,'(" muw           9.8267598403e1       ",1x,e22.15,5x,"J/KG/K")') g
 
endif


if(icol == 3 .or. icol == 999) then
write(*,*) " "
write(*,*) " "
write(*,*) " "
write(*,*) "TABLE A8, COLUMN 3 OF OS2008"
write(*,*) " "
write(*,*) " FUNCTION VALUES FOR SEAWATER AT THE STANDARD OCEAN STATE"
write(*,*) " T = 273.15 K, P = 101325 Pa, S = 0.03516504"

s = 0.03516504d0
t = 273.15d0
p = 101325d0
rho=liq_density_si(t,p)

g =   liq_g_si(0, 0, t, p) + sal_g_si(0, 0, 0, s, t, p)
gs =  sal_g_si(1, 0, 0, s, t, p)
gt =  liq_g_si(1, 0, t, p) + sal_g_si(0, 1, 0, s, t, p)
gp =  liq_g_si(0, 1, t, p) + sal_g_si(0, 0, 1, s, t, p)
gss = sal_g_si(2, 0, 0, s, t, p)
gtt = liq_g_si(2, 0, t, p) + sal_g_si(0, 2, 0, s, t, p)
gpp = liq_g_si(0, 2, t, p) + sal_g_si(0, 0, 2, s, t, p)
gst = sal_g_si(1, 1, 0, s, t, p)
gsp = sal_g_si(1, 0, 1, s, t, p)
gtp = liq_g_si(1, 1, t, p) + sal_g_si(0, 1, 1, s, t, p)

write(*,*)" "
write(*,*)"QUANTITY    TABLE VALUE           COMPUTED VALUE            UNITS"

write(*,'(" g           4.e-9               ",1x,e22.15,5x,"J KG-1")') g
write(*,'(" dg/ds       6.39974067312299e4  ",1x,e22.15,5x,"J KG-1")') gs
write(*,'(" dg/dt      -6.e-11              ",1x,e22.15,5x,"J KG-1 K-1")') gt
write(*,'(" dg/dp       9.7266121669485e-4  ",1x,e22.15,5x,"J M3 KG-1")') gp
write(*,'(" d2g/dsdp   -7.59615411515309e-4 ",1x,e22.15,5x,"M3 KG-1")') gsp
write(*,'(" d2g/dt2    -1.459449308080e1    ",1x,e22.15,5x,"J KG-1 K-2")') gtt
write(*,'(" d2g/dpdt    5.15167556274e-8    ",1x,e22.15,5x,"J M3 KG-1 K-1")') gtp
write(*,'(" d2g/dp2    -4.507753774102e-13  ",1x,e22.15,5x,"J M3 KG-1 Pa-1")') gpp
write(*,'(" h           2.e-8               ",1x,e22.15,5x,"J KG-1")') g-t*gt
write(*,'(" f          -9.855489778e1       ",1x,e22.15,5x,"J KG-1")') g-p*gp
write(*,'(" u          -9.85548978e1        ",1x,e22.15,5x,"J KG-1")') g-t*gt-p*gp
write(*,'(" sigma       6.e-11              ",1x,e22.15,5x,"J KG-1 K-1")') -gt
write(*,'(" rho         1.02810719995401e3  ",1x,e22.15,5x,"KG M-3")') 1d0/gp
write(*,'(" cp          3.98648578502e3     ",1x,e22.15,5x,"J KG-1 K-1")') -t*gtt
write(*,'(" alpha       5.2964747379e-5     ",1x,e22.15,5x,"J K-1")') gtp/gp
write(*,'(" kappa_t     4.634454110774e-10  ",1x,e22.15,5x,"J Pa-1")') -gpp/gp
write(*,'(" kappa_s     4.63258452069e-10   ",1x,e22.15,5x,"J Pa-1")') (gtp**2-gtt*gpp)/(gp*gtt)
write(*,'(" muw        -2.250471368e3       ",1x,e22.15,5x,"J/KG/K")') g - s*gs

endif

if(icol /= 999) go to 999

end subroutine

end module os2008_5
