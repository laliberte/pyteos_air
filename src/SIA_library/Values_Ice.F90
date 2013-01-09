module values_ice

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM THE ICE ROUTINES
!RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0, FILE CONSTANTS_0.F90
!     CONVERT_0,   FILE CONVERT_0.F90
!     ICE_1,       FILE ICE_1.F90
!     ICE_2,       FILE ICE_2.F90

use constants_0
use convert_0
use ice_1
use ice_2

implicit none 
character*16, private :: version = '8 Jun 2010'
character (len = *), parameter :: fmt = "(A,T32,1X,E19.12,3X,A)"

contains

!==========================================================================
subroutine ice_chk_values
!==========================================================================

implicit none

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR PROPERTIES OF PURE ICE'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_1'
write(*,fmt)'ice_g(0,0,270,1d5) =', ice_g_si(0,0,270d0,1d5),&
'-3786.74963128'
write(*,fmt)'ice_g(1,0,270,1d5) =', ice_g_si(1,0,270d0,1d5),&
' 1244.97335506'
write(*,fmt)'ice_g(0,1,270,1d5) =', ice_g_si(0,1,270d0,1d5),&
' 1.09029713624E-03'
write(*,fmt)'ice_g(2,0,270,1d5) =', ice_g_si(2,0,270d0,1d5),&
'-7.67955356337'
write(*,fmt)'ice_g(1,1,270,1d5) =', ice_g_si(1,1,270d0,1d5),&
' 1.72604208697E-07'
write(*,fmt)'ice_g(0,2,270,1d5) =', ice_g_si(0,2,270d0,1d5),&
'-1.27811223643E-13'

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'COMPARISON OF CHECK VALUES FROM ICE_2'
write(*,fmt)'ice_chempot(270,1d5) =', ice_chempot_si(270d0,1d5),&
'-3786.74963128'
write(*,fmt)'ice_cp(270,1d5) =', ice_cp_si(270d0,1d5),&
' 2073.47946211'
write(*,fmt)'ice_density(270,1d5) =', ice_density_si(270d0,1d5),&
' 917.181167192'
write(*,fmt)'ice_enthalpy(270,1d5) =', ice_enthalpy_si(270d0,1d5),&
'-339929.555499'
write(*,fmt)'ice_entropy(270,1d5) =', ice_entropy_si(270d0,1d5),&
'-1244.97335506'
write(*,fmt)'ice_expansion(270,1d5) =', ice_expansion_si(270d0,1d5),&
' 1.58309329594E-04'
write(*,fmt)'ice_helmholtz_energy(270,1d5) =', ice_helmholtz_energy_si(270d0,1d5),&
'-3895.77934490'
write(*,fmt)'ice_internal_energy(270,1d5) =', ice_internal_energy_si(270d0,1d5),&
'-340038.585212'
write(*,fmt)'ice_kappa_s(270,1d5) =', ice_kappa_s_si(270d0,1d5),&
' 1.13667916416E-10'
write(*,fmt)'ice_kappa_t(270,1d5) =', ice_kappa_t_si(270d0,1d5),&
' 1.17226047281E-10'
write(*,fmt)'ice_lapserate(270,1d5) =', ice_lapserate_si(270d0,1d5),&
' 2.24758128545E-08'
write(*,fmt)'ice_p_coefficient(270,1d5) =', ice_p_coefficient_si(270d0,1d5),&
' 1350462.06254'
write(*,fmt)'ice_specific_volume(270,1d5) =', ice_specific_volume_si(270d0,1d5),&
' 1.09029713624E-03'

end subroutine

end module values_ice
