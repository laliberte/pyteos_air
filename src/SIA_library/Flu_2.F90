module flu_2

!#########################################################################

!THIS MODULE IMPLEMENTS THERMODYNAMIC PROPERTIES OF PURE WATER
!AS FUNCTIONS OF ABSOLUTE TEMPERATURE IN K AND DENSITY IN KG/M3,
!COMPUTED FROM THE HELMHOLTZ FUNCTION FLU_F_SI(DRV_T, DRV_D, T_SI, D_SI)

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

!FEISTEL, R., KRETZSCHMAR, H.-J., SPAN, R., HAGEN, E., !WRIGHT, D.G., JACKETT, D.R.:
!THERMODYNAMIC PROPERTIES OF SEA AIR.
!OCEAN SCIENCE DISCUSSION 6(2009)2193-2325.

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES
!     CONSTANTS_0,     FILE CONSTANTS_0.F90
!     FLU_1 (IAPWS95), FILE FLU_1.F90 

!#########################################################################

use constants_0
use flu_1

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: flu_cp_si, flu_cv_si, flu_enthalpy_si, flu_entropy_si, &
          flu_expansion_si, flu_gibbs_energy_si, flu_internal_energy_si, &
          flu_kappa_s_si, flu_kappa_t_si, flu_lapserate_si, flu_pressure_si, &
          flu_soundspeed_si

contains

!*******************************************************************
function flu_cp_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_cp_si
real*8 temp, f_d, f_dd, f_td, f_tt

!OUTPUT:
!CP(T,D) = T * (DETA/DT)_P 
!        = ISOBARIC HEAT CAPACITY IN J/(KG*K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K,
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_cp_si(300,1000) = 4159.37519963

flu_cp_si = errorreturn

f_d  = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return
f_td = flu_f_si(1, 1, t_si, d_si)
if(f_td == errorreturn) return
f_tt = flu_f_si(2, 0, t_si, d_si)
if(f_tt == errorreturn) return

temp = 2d0*f_d + d_si*f_dd
if(temp == 0d0) return

flu_cp_si = t_si * (f_td**2*d_si/temp - f_tt)

end function

!*******************************************************************
function flu_cv_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_cv_si
real*8 f_tt

!OUTPUT:
!CV(T,D) = T * (DETA/DT)_D 
!        = ISOCHORIC HEAT CAPACITY IN J/(KG*K)

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K,
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_cv_si(300,1000) = 4105.20614776

flu_cv_si = errorreturn

f_tt = flu_f_si(2, 0, t_si, d_si)
if(f_tt == errorreturn) return

flu_cv_si = -t_si * f_tt

end function

!*******************************************************************
function flu_enthalpy_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_enthalpy_si
real*8 f, f_t, f_d

!OUTPUT:
!H(T,D) = F - T*(DF/DT)_D + D*(DF/DD)_T ENTHALPY IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_enthalpy_si(300,1000) = 119752.511434

flu_enthalpy_si = errorreturn

f = flu_f_si(0, 0, t_si, d_si)
if(f == errorreturn) return
f_t = flu_f_si(1, 0, t_si, d_si)
if(f_t == errorreturn) return
f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return

flu_enthalpy_si = f - t_si * f_t + d_si * f_d

end function

!*******************************************************************
function flu_entropy_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_entropy_si
real*8 f_t

!OUTPUT:
!S(T,D) = - (DF/DT)_D ENTROPY IN J/KG K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_entropy_si(300,1000) = 390.904170767

flu_entropy_si = errorreturn

f_t = flu_f_si(1, 0, t_si, d_si)
if(f_t == errorreturn) return

flu_entropy_si = -flu_f_si(1, 0, t_si, d_si)

end function

!*******************************************************************
function flu_expansion_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_expansion_si
real*8 temp, f_d, f_dd, f_td

!OUTPUT:
!THERMAL EXPANSION COEFFICIENT IN 1/K

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_expansion_si(300,1000) = 2.82413312531E-04

flu_expansion_si = errorreturn

f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_td = flu_f_si(1, 1, t_si, d_si)
if(f_td == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return

temp = 2d0*f_d + d_si*f_dd
if(temp == 0d0) return

flu_expansion_si = f_td/temp

end function

!*******************************************************************
function flu_gibbs_energy_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_gibbs_energy_si
real*8 f, f_d

!OUTPUT:
!G(T,D) = F + D*(DF/DD)_T   GIBBS ENERGY IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K,
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_gibbs_energy_si(300,1000) = 2481.26020392

flu_gibbs_energy_si = errorreturn

f = flu_f_si(0, 0, t_si, d_si)
if(f == errorreturn) return
f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return

flu_gibbs_energy_si = f + d_si * f_d

end function

!*******************************************************************
function flu_internal_energy_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_internal_energy_si
real*8 f, f_t

!OUTPUT:
!U(T,D) = F - T*(DF/DT)_D INTERNAL ENERGY IN J/KG

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_internal_energy_si(300,1000) = 111919.510078

flu_internal_energy_si = errorreturn

f = flu_f_si(0, 0, t_si, d_si)
if(f == errorreturn) return
f_t = flu_f_si(1, 0, t_si, d_si)
if(f_t == errorreturn) return

flu_internal_energy_si = f - t_si * f_t

end function

!*******************************************************************
function flu_kappa_s_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_kappa_s_si
real*8 f_d, f_tt, f_td, f_dd, temp

!OUTPUT:
!KAPPA_S(T,D) = (1/D) * (DD/DP)_ETA 
!             = ISENTROPIC COMPRESSIBILITY IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_kappa_s_si(300,1000) = 4.35960581171E-10

flu_kappa_s_si = errorreturn

f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_tt = flu_f_si(2, 0, t_si, d_si)
if(f_tt == errorreturn) return
f_td = flu_f_si(1, 1, t_si, d_si)
if(f_td == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return

temp = f_tt*(2d0*f_d+d_si*f_dd)-d_si*f_td**2
if(temp == 0d0) return

flu_kappa_s_si = f_tt/d_si**2 / temp

end function

!*******************************************************************
function flu_kappa_t_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_kappa_t_si
real*8 temp, f_d, f_dd

!OUTPUT:
!ISOTHERMAL COMPRESSIBILITY IN 1/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_kappa_t_si(300,1000) = 4.41713172024E-10

flu_kappa_t_si = errorreturn

f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return

temp = 2d0*f_d + d_si*f_dd
if(temp == 0d0 .or. d_si == 0d0) return

flu_kappa_t_si = 1d0/(d_si**2*temp)

end function

!*******************************************************************
function flu_lapserate_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_lapserate_si
real*8 temp, f_d, f_tt, f_dd, f_td

!OUTPUT:
!GAMMA(T,D) = (DT/DP)_ETA  
!           = ADIABATIC LAPSE RATE IN K/PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_lapserate_si(300,1000) = 2.0369403983|3E-08

flu_lapserate_si = errorreturn

f_d  = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_tt = flu_f_si(2, 0, t_si, d_si)
if(f_tt == errorreturn) return
f_td = flu_f_si(1, 1, t_si, d_si)
if(f_td == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return

temp = d_si*f_td**2-f_tt*(2d0*f_d + d_si*f_dd)
if(temp == 0d0 .or. d_si == 0d0) return

flu_lapserate_si = (f_td/d_si)/temp

end function

!*******************************************************************
function flu_pressure_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_pressure_si
real*8 f_d

!OUTPUT:
!P(T,D) = D^2*(df/dD)_T PRESSURE IN PA

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_pressure_si(300,1000) = 7833001.35597

flu_pressure_si = errorreturn

f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return

flu_pressure_si = d_si**2 * f_d

end function

!*******************************************************************
function flu_soundspeed_si(t_si, d_si)
!*******************************************************************

implicit none 

real*8 t_si, d_si, flu_soundspeed_si
real*8 f_d, f_tt, f_td, f_dd, c2

!OUTPUT:
!C(T,D) = SQRT[ (DP/DD)_ETA ] SOUND SPEED IN M/S

!INPUTS:
!T_SI      ABSOLUTE TEMPERATURE IN K
!D_SI      ABSOLUTE DENSITY IN KG/M3

!CHECK VALUE: flu_soundspeed_si(300,1000) = 1514.52479780

flu_soundspeed_si = errorreturn

f_d = flu_f_si(0, 1, t_si, d_si)
if(f_d == errorreturn) return
f_tt = flu_f_si(2, 0, t_si, d_si)
if(f_tt == errorreturn) return
f_td = flu_f_si(1, 1, t_si, d_si)
if(f_td == errorreturn) return
f_dd = flu_f_si(0, 2, t_si, d_si)
if(f_dd == errorreturn) return

c2 = 2d0*d_si*f_d + d_si**2*(f_dd - f_td**2/f_tt)

if(c2 < 0d0)return

flu_soundspeed_si = sqrt(c2)

end function

end module flu_2
