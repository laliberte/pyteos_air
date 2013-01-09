module values_gsw

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM 
!THE SEAWATER ROUTINES RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN 
!THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     SAL_1,         FILE SAL_1.F90
!     SAL_2,         FILE SAL_2.F90
!     LIQ_F03_5,     FILE LIQ_F03_5
!     GSW_LIBRARY_5, FILE GSW_LIBRARY_5

use constants_0
use convert_0
use maths_0
use sal_1
use sal_2
use liq_f03_5
use convert_5
use gsw_library_5

implicit none 
character*16, private :: version = '1 Feb 2011'

contains

!==========================================================================
subroutine gsw_chk_values
!==========================================================================

implicit none 
character (len = *), parameter :: fmt = "(A,T52,1X,E19.12,3X,A)"

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF SEAWATER'
write(*,"(A)")'OBTAINED WITH THE GSW LIBRARY ROUTINES BUT USING LIQ_F03_5'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,fmt)'gsw_g(0,0,0,35,26.85,0d0) =', gsw_g(0,0,0,35d0,26.85d0,0d0), &
'-5113.70064124'
write(*,fmt)'gsw_g(1,0,0,35,26.85,0d0) =', gsw_g(1,0,0,35d0,26.85d0,0d0), &
' 78.5928261339'
write(*,fmt)'gsw_g(0,1,0,35,26.85,0d0) =', gsw_g(0,1,0,35d0,26.85d0,0d0), &
'-374.452000830'
write(*,fmt)'gsw_g(0,0,1,35,26.85,0d0) =', gsw_g(0,0,1,35d0,26.85d0,0d0), &
' 0.977858058750E-03'
write(*,fmt)'gsw_g(2,0,0,35,26.85,0d0) =', gsw_g(2,0,0,35d0,26.85d0,0d0), &
' 2.24755137017'
write(*,fmt)'gsw_g(1,1,0,35,26.85,0d0) =', gsw_g(1,1,0,35d0,26.85d0,0d0), &
' 0.789935187192'
write(*,fmt)'gsw_g(1,0,1,35,26.85,0d0) =', gsw_g(1,0,1,35d0,26.85d0,0d0), &
'-0.716680931996E-06'
write(*,fmt)'gsw_g(0,2,0,35,26.85,0d0) =', gsw_g(0,2,0,35d0,26.85d0,0d0), &
'-13.3358337534'
write(*,fmt)'gsw_g(0,1,1,35,26.85,0d0) =', gsw_g(0,1,1,35d0,26.85d0,0d0), &
' 0.304607508052E-06'
write(*,fmt)'gsw_g(0,0,2,35,26.85,0d0) =', gsw_g(0,0,2,35d0,26.85d0,0d0), &
'-0.410939723950E-12'

write(*,"(A)")' '
!write(*,fmt)'gsw_asal_from_psal(35.527515065427778,201,-21,1023) = ', &
! gsw_asal_from_psal(35.527515065427778, 201d0, -21d0, 1023d0), " 35.7"
!Rank mismatch!
!write(*,fmt)'gsw_asal_from_psal(35.527515065,201,-21,1023) = ', &
! gsw_asal_from_psal(35.527515065, 201d0, -21d0, 1023d0), " 35.7"

write(*,fmt)'gsw_asal_from_psal(35,180,40,2d3) = ', &
 gsw_asal_from_psal(35d0, 180d0, 40d0, 2d3), " 35.1888478029"

write(*,fmt)'gsw_asal_from_psal(8,20,57,0d0) = ', &
 gsw_asal_from_psal(8d0, 20d0, 57d0, 0d0), " 8.13338057143"

write(*,fmt)'gsw_psal_from_asal(35.7,201,-21,1023) = ', &
 gsw_psal_from_asal(35.7d0, 201d0, -21d0, 1023d0), " 35.527515065427778"

write(*,"(A)")' '
write(*,fmt)'gsw_alpha_t(35.7,25.5,1023) = ', &
   gsw_alpha_t(35.7d0, 25.5d0, 1023d0), ' 0.000309837839319d0'

write(*,fmt)'gsw_beta_t(35.7,25.5,1023) = ', &
   gsw_beta_t(35.7d0, 25.5d0, 1023d0), ' 0.000725729797839d0'

write(*,fmt)'gsw_cp(35.7,25.5,1023) = ', &
   gsw_cp(35.7d0, 25.5d0, 1023d0), ' 3974.42541260d0'

write(*,fmt)'gsw_ctmp_from_ptmp0(35.7,25.5) = ', &
   gsw_ctmp_from_ptmp0(35.7d0, 25.5d0), ' 25.4805463842d0'

write(*,fmt)'gsw_ptmp0_from_ctmp(35.7,25.4805463842239) = ', &
   gsw_ptmp0_from_ctmp(35.7d0, 25.4805463842239d0), ' 25.5d0'

write(*,fmt)'gsw_dens(35.7,25.5,1023) = ', &
   gsw_dens(35.7d0, 25.5d0, 1023d0), ' 1027.95249316d0'

write(*,fmt)'gsw_enthalpy(35.7,25.5,1023) = ', &
   gsw_enthalpy(35.7d0, 25.5d0, 1023d0), ' 110776.712409d0'

write(*,fmt)'gsw_entropy(35.7,25.5,1023) = ', &
   gsw_entropy(35.7d0, 25.5d0, 1023d0), ' 352.818797715d0'

write(*,fmt)'gsw_kappa(35.7,25.5,1023) = ', &
   gsw_kappa(35.7d0, 25.5d0, 1023d0), ' 4.03386268546d-6'

write(*,fmt)'gsw_kappa_t(35.7,25.5,1023) = ', &
   gsw_kappa_t(35.7d0, 25.5d0, 1023d0), ' 4.10403794615d-6'

write(*,fmt)'gsw_pden(35.7,25.5,1023, 0) = ', &
   gsw_pden(35.7d0, 25.5d0, 1023d0, 0d0), ' 1023.66254941d0'

write(*,fmt)'gsw_ptmp(35.7,25.5,1023, 0) = ', &
   gsw_ptmp(35.7d0, 25.5d0, 1023d0, 0d0), ' 25.2720983155d0'

write(*,fmt)'gsw_ptmp(35.7,25.2720983155409d0,0,1023d0) = ', &
   gsw_ptmp(35.7d0, 25.2720983155409d0, 0d0, 1023d0), ' 25.5'

write(*,fmt)'gsw_specvol(35.7,25.5,1023) = ', &
   gsw_specvol(35.7d0, 25.5d0, 1023d0), ' 0.000972807602158d0'

write(*,fmt)'gsw_svel(35.7,25.5,1023) = ', &
   gsw_svel(35.7d0, 25.5d0, 1023d0), ' 1552.93372863d0'

write(*,"(A)")' '
write(*,fmt)'gsw_alpha_ct(35,20,1000) =', &
   gsw_alpha_ct(35d0, 20d0, 1000d0), ' 2.69418609861E-04'

write(*,fmt)'gsw_alpha_pt0(35,20,1000) =', &
   gsw_alpha_pt0(35d0, 20d0, 1000d0), ' 2.6975373331|7E-04'

write(*,fmt)'gsw_beta_ct(35,20,1000) =', &
   gsw_beta_ct(35d0, 20d0, 1000d0), ' 7.23213672954E-04'

write(*,fmt)'gsw_beta_pt0(35,20,1000) =', &
   gsw_beta_pt0(35d0, 20d0, 1000d0), ' 7.31582583383E-04'

write(*,"(A)")' '
write(*,fmt)'gsw_cabb_ct(35,20,1d3) =', &
   gsw_cabb_ct(35d0, 20d0, 1d3), ' 0.8969073|83083E-05'  

write(*,fmt)'gsw_cabb_pt0(35,20,1d3) =', &
   gsw_cabb_pt0(35d0, 20d0, 1d3), ' 0.8759631|54048E-05'  

write(*,fmt)'gsw_thrmb_ct(35,20,1d3) =', &
   gsw_thrmb_ct(35d0, 20d0, 1d3), ' 0.1727083|65652E-07'  

write(*,fmt)'gsw_thrmb_pt0(35,20,1d3) =', &
   gsw_thrmb_pt0(35d0, 20d0, 1d3), ' 0.1709450|45984E-07'

end subroutine


end module values_gsw
