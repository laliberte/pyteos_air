module values_sal

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE SAL ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     SAL_1,       FILE SAL_1.F90
!     SAL_2,       FILE SAL_2.F90

use constants_0
use convert_0
use sal_1
use sal_2

implicit none 
character*16, private :: version = '8 Jun 2010'
character (len = *), parameter :: fmt = "(A,T38,1X,E19.12,3X,A)"

contains

!==========================================================================
subroutine sal_chk_values
!==========================================================================

implicit none 

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF SEA SALT'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SAL_1'

if(IsExtension2010)then

   write(*,fmt)'sal_g_term(1,0,0,300,1d6) =', sal_g_term_si(1,0,0,300d0,1d6),&
   ' 79427.9694846'
   write(*,fmt)'sal_g_term(1,1,0,300,1d6) =', sal_g_term_si(1,1,0,300d0,1d6),&
   ' 264.759898282'
   write(*,fmt)'sal_g_term(1,0,1,300,1d6) =', sal_g_term_si(1,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,2,0,300,1d6) =', sal_g_term_si(1,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,1,1,300,1d6) =', sal_g_term_si(1,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,0,2,300,1d6) =', sal_g_term_si(1,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(2,0,0,300,1d6) =', sal_g_term_si(2,0,0,300d0,1d6),&
   ' 301223.934546274'
   write(*,fmt)'sal_g_term(2,1,0,300,1d6) =', sal_g_term_si(2,1,0,300d0,1d6),&
   ' 1558.10730392553'
   write(*,fmt)'sal_g_term(2,0,1,300,1d6) =', sal_g_term_si(2,0,1,300d0,1d6),&
   '-7.60101106780E-04'
   write(*,fmt)'sal_g_term(2,2,0,300,1d6) =', sal_g_term_si(2,2,0,300d0,1d6),&
   ' 19.0397747694'
   write(*,fmt)'sal_g_term(2,1,1,300,1d6) =', sal_g_term_si(2,1,1,300d0,1d6),&
   ' 9.35496617545E-07'
   write(*,fmt)'sal_g_term(2,0,2,300,1d6) =', sal_g_term_si(2,0,2,300d0,1d6),&
   ' 1.28077803426E-12'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(3,0,0,300,1d6) =', sal_g_term_si(3,0,0,300d0,1d6),&
   '-345570.781498'
   write(*,fmt)'sal_g_term(3,1,0,300,1d6) =', sal_g_term_si(3,1,0,300d0,1d6),&
   '-1749.19911587'
   write(*,fmt)'sal_g_term(3,0,1,300,1d6) =', sal_g_term_si(3,0,1,300d0,1d6),&
   ' 2.01608449881E-04'
   write(*,fmt)'sal_g_term(3,2,0,300,1d6) =', sal_g_term_si(3,2,0,300d0,1d6),&
   '-9.3208882362'
   write(*,fmt)'sal_g_term(3,1,1,300,1d6) =', sal_g_term_si(3,1,1,300d0,1d6),&
   ' 2.17412289837E-08'
   write(*,fmt)'sal_g_term(3,0,2,300,1d6) =', sal_g_term_si(3,0,2,300d0,1d6),&
   '-4.78063827318E-13'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(4,0,0,300,1d6) =', sal_g_term_si(4,0,0,300d0,1d6),&
   ' 1468073.64559' 
   write(*,fmt)'sal_g_term(4,1,0,300,1d6) =', sal_g_term_si(4,1,0,300d0,1d6),&
   ' 7741.24404962' 
   write(*,fmt)'sal_g_term(4,0,1,300,1d6) =', sal_g_term_si(4,0,1,300d0,1d6),&
   '-4.33952339916E-04'
   write(*,fmt)'sal_g_term(4,2,0,300,1d6) =', sal_g_term_si(4,2,0,300d0,1d6),&
   '-6.13689642370'
   write(*,fmt)'sal_g_term(4,1,1,300,1d6) =', sal_g_term_si(4,1,1,300d0,1d6),&
   '-3.50876195080E-06'
   write(*,fmt)'sal_g_term(4,0,2,300,1d6) =', sal_g_term_si(4,0,2,300d0,1d6),&
   '-6.06204383305E-13'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(5,0,0,300,1d6) =', sal_g_term_si(5,0,0,300d0,1d6),&
   '-3776969.31546'
   write(*,fmt)'sal_g_term(5,1,0,300,1d6) =', sal_g_term_si(5,1,0,300d0,1d6),&
   '-15135.6522248'
   write(*,fmt)'sal_g_term(5,0,1,300,1d6) =', sal_g_term_si(5,0,1,300d0,1d6),&
   ' 1.11272425476E-03'
   write(*,fmt)'sal_g_term(5,2,0,300,1d6) =', sal_g_term_si(5,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(5,1,1,300,1d6) =', sal_g_term_si(5,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(5,0,2,300,1d6) =', sal_g_term_si(5,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(6,0,0,300,1d6) =', sal_g_term_si(6,0,0,300d0,1d6),&
   ' 6151235.69234'
   write(*,fmt)'sal_g_term(6,1,0,300,1d6) =', sal_g_term_si(6,1,0,300d0,1d6),&
   ' 14157.0509933'
   write(*,fmt)'sal_g_term(6,0,1,300,1d6) =', sal_g_term_si(6,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,2,0,300,1d6) =', sal_g_term_si(6,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,1,1,300,1d6) =', sal_g_term_si(6,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,0,2,300,1d6) =', sal_g_term_si(6,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(7,0,0,300,1d6) =', sal_g_term_si(7,0,0,300d0,1d6),&
   '-3734033.38866'
   write(*,fmt)'sal_g_term(7,1,0,300,1d6) =', sal_g_term_si(7,1,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,0,1,300,1d6) =', sal_g_term_si(7,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,2,0,300,1d6) =', sal_g_term_si(7,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,1,1,300,1d6) =', sal_g_term_si(7,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,0,2,300,1d6) =', sal_g_term_si(7,0,2,300d0,1d6),&
   ' 0'

else    !IsExtension2010 = .false.

   write(*,fmt)'sal_g_term(1,0,0,300,1d6) =', sal_g_term_si(1,0,0,300d0,1d6),&
   ' 79427.9694846'
   write(*,fmt)'sal_g_term(1,1,0,300,1d6) =', sal_g_term_si(1,1,0,300d0,1d6),&
   ' 264.759898282'
   write(*,fmt)'sal_g_term(1,0,1,300,1d6) =', sal_g_term_si(1,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,2,0,300,1d6) =', sal_g_term_si(1,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,1,1,300,1d6) =', sal_g_term_si(1,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(1,0,2,300,1d6) =', sal_g_term_si(1,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(2,0,0,300,1d6) =', sal_g_term_si(2,0,0,300d0,1d6),&
   ' 301223.934546'
   write(*,fmt)'sal_g_term(2,1,0,300,1d6) =', sal_g_term_si(2,1,0,300d0,1d6),&
   ' 1558.10730393'
   write(*,fmt)'sal_g_term(2,0,1,300,1d6) =', sal_g_term_si(2,0,1,300d0,1d6),&
   '-7.60101106780E-04'
   write(*,fmt)'sal_g_term(2,2,0,300,1d6) =', sal_g_term_si(2,2,0,300d0,1d6),&
   ' 19.0397747694'
   write(*,fmt)'sal_g_term(2,1,1,300,1d6) =', sal_g_term_si(2,1,1,300d0,1d6),&
   ' 9.35496617545E-07'
   write(*,fmt)'sal_g_term(2,0,2,300,1d6) =', sal_g_term_si(2,0,2,300d0,1d6),&
   ' 1.28077803426E-12'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(3,0,0,300,1d6) =', sal_g_term_si(3,0,0,300d0,1d6),&
   '-345570.781498'
   write(*,fmt)'sal_g_term(3,1,0,300,1d6) =', sal_g_term_si(3,1,0,300d0,1d6),&
   '-1749.19911587'
   write(*,fmt)'sal_g_term(3,0,1,300,1d6) =', sal_g_term_si(3,0,1,300d0,1d6),&
   ' 2.01608449881E-04'
   write(*,fmt)'sal_g_term(3,2,0,300,1d6) =', sal_g_term_si(3,2,0,300d0,1d6),&
   '-9.32088823620'
   write(*,fmt)'sal_g_term(3,1,1,300,1d6) =', sal_g_term_si(3,1,1,300d0,1d6),&
   ' 2.17412289837E-08'
   write(*,fmt)'sal_g_term(3,0,2,300,1d6) =', sal_g_term_si(3,0,2,300d0,1d6),&
   '-4.78063827318E-13'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(4,0,0,300,1d6) =', sal_g_term_si(4,0,0,300d0,1d6),&
   ' 1468073.64559'
   write(*,fmt)'sal_g_term(4,1,0,300,1d6) =', sal_g_term_si(4,1,0,300d0,1d6),&
   ' 7741.24404962'
   write(*,fmt)'sal_g_term(4,0,1,300,1d6) =', sal_g_term_si(4,0,1,300d0,1d6),&
   '-4.33952339916E-04'
   write(*,fmt)'sal_g_term(4,2,0,300,1d6) =', sal_g_term_si(4,2,0,300d0,1d6),&
   '-6.13689642370'
   write(*,fmt)'sal_g_term(4,1,1,300,1d6) =', sal_g_term_si(4,1,1,300d0,1d6),&
   '-3.50876195080E-06'
   write(*,fmt)'sal_g_term(4,0,2,300,1d6) =', sal_g_term_si(4,0,2,300d0,1d6),&
   '-6.06204383305E-13'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(5,0,0,300,1d6) =', sal_g_term_si(5,0,0,300d0,1d6),&
   '-3776969.31546'
   write(*,fmt)'sal_g_term(5,1,0,300,1d6) =', sal_g_term_si(5,1,0,300d0,1d6),&
   '-15135.6522248'
   write(*,fmt)'sal_g_term(5,0,1,300,1d6) =', sal_g_term_si(5,0,1,300d0,1d6),&
   ' 1.11272425476E-03'
   write(*,fmt)'sal_g_term(5,2,0,300,1d6) =', sal_g_term_si(5,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(5,1,1,300,1d6) =', sal_g_term_si(5,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(5,0,2,300,1d6) =', sal_g_term_si(5,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(6,0,0,300,1d6) =', sal_g_term_si(6,0,0,300d0,1d6),&
   ' 6151235.69234'
   write(*,fmt)'sal_g_term(6,1,0,300,1d6) =', sal_g_term_si(6,1,0,300d0,1d6),&
   ' 14157.0509933'
   write(*,fmt)'sal_g_term(6,0,1,300,1d6) =', sal_g_term_si(6,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,2,0,300,1d6) =', sal_g_term_si(6,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,1,1,300,1d6) =', sal_g_term_si(6,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(6,0,2,300,1d6) =', sal_g_term_si(6,0,2,300d0,1d6),&
   ' 0'

   write(*,"(A)")' '
   write(*,fmt)'sal_g_term(7,0,0,300,1d6) =', sal_g_term_si(7,0,0,300d0,1d6),&
   '-3734033.38866'
   write(*,fmt)'sal_g_term(7,1,0,300,1d6) =', sal_g_term_si(7,1,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,0,1,300,1d6) =', sal_g_term_si(7,0,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,2,0,300,1d6) =', sal_g_term_si(7,2,0,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,1,1,300,1d6) =', sal_g_term_si(7,1,1,300d0,1d6),&
   ' 0'
   write(*,fmt)'sal_g_term(7,0,2,300,1d6) =', sal_g_term_si(7,0,2,300d0,1d6),&
   ' 0'

endif

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM SAL_2'

write(*,fmt)'sal_g(0,0,0,0.035,300,1d6) =', sal_g_si(0,0,0,0.035d0,300d0,1d6),&
' 127.033640309'
write(*,fmt)'sal_g(1,0,0,0.035,300,1d6) =', sal_g_si(1,0,0,0.035d0,300d0,1d6),&
' 77949.2100395'
write(*,fmt)'sal_g(0,1,0,0.035,300,1d6) =', sal_g_si(0,1,0,0.035d0,300d0,1d6),&
' 18.6360407073'
write(*,fmt)'sal_g(0,0,1,0.035,300,1d6) =', sal_g_si(0,0,1,0.035d0,300d0,1d6),&
'-2.55600080319E-05'
write(*,fmt)'sal_g(2,0,0,0.035,300,1d6) =', sal_g_si(2,0,0,0.035d0,300d0,1d6),&
' 2248200.54660'
write(*,fmt)'sal_g(1,1,0,0.035,300,1d6) =', sal_g_si(1,1,0,0.035d0,300d0,1d6),&
' 790.563810558'
write(*,fmt)'sal_g(1,0,1,0.035,300,1d6) =', sal_g_si(1,0,1,0.035d0,300d0,1d6),&
'-7.15686520588E-04'
write(*,fmt)'sal_g(0,2,0,0.035,300,1d6) =', sal_g_si(0,2,0,0.035d0,300d0,1d6),&
' 0.597842170749'
write(*,fmt)'sal_g(0,1,1,0.035,300,1d6) =', sal_g_si(0,1,1,0.035d0,300d0,1d6),&
' 2.85865076269E-08'
write(*,fmt)'sal_g(0,0,2,0.035,300,1d6) =', sal_g_si(0,0,2,0.035d0,300d0,1d6),&
' 4.09543164905E-14'

write(*,"(A)")' '
write(*,fmt)'sal_act_coeff(0.035,300,1d6) =', sal_act_coeff_si(0.035d0,300d0,1d6),&
'-0.527003008913'
write(*,fmt)'sal_act_potential(0.035,300,1d6) =', sal_act_potential_si(0.035d0,300d0,1d6),&
'-0.429940465498'
write(*,fmt)'sal_activity_w(0.035,300,1d6) =', sal_activity_w_si(0.035d0,300d0,1d6),&
' 0.981388410188'
write(*,fmt)'sal_dilution(0.035,300,1d6) =', sal_dilution_si(0.035d0,300d0,1d6),&
' 78687.0191309'
write(*,fmt)'sal_mixenth(0.01,0.035,0.6,300,1d6) =', sal_mixenthalpy_si(0.01d0,0.035d0,0.6d0,300d0,1d6),&
' 16.1539911284'
write(*,fmt)'sal_mixentr(0.01,0.035,0.6,300,1d6) =', sal_mixentropy_si(0.01d0,0.035d0,0.6d0,300d0,1d6),&
' 0.96682942261|7'
write(*,fmt)'sal_mixvol(0.01,0.035,0.6,300,1d6) =', sal_mixvolume_si(0.01d0,0.035d0,0.6d0,300d0,1d6),&
'-5.94174956892E-08'
write(*,fmt)'sal_osm_coeff(0.035,300,1d6) =', sal_osm_coeff_si(0.035d0,300d0,1d6),&
' 0.902937456585'
write(*,fmt)'sal_chempot_h2o(0.035,300,1d6) =', sal_chempot_h2o_si(0.035d0,300d0,1d6),&
'-2601.18871107'
write(*,fmt)'sal_chempot_rel(0.035,300,1d6) =', sal_chempot_rel_si(0.035d0,300d0,1d6),&
' 77949.2100395'
write(*,fmt)'sal_chem_coeff(0.035,300,1d6) =', sal_chem_coeff_si(0.035d0,300d0,1d6),&
' 2754.04566958'
write(*,fmt)'sal_saltenthalpy(0.035,300,1d6) =', sal_saltenthalpy_si(0.035d0,300d0,1d6),&
'-156107.959196'
write(*,fmt)'sal_saltentropy(0.035,300,1d6) =', sal_saltentropy_si(0.035d0,300d0,1d6),&
'-532.458305922'
write(*,fmt)'sal_saltvolume(0.035,300,1d6) =', sal_saltvolume_si(0.035d0,300d0,1d6),&
'-7.30285943768E-04'
write(*,fmt)'sal_molality(0.035) =', sal_molality_si(0.035d0),&
' 1.15493681893'

end subroutine

end module values_sal
