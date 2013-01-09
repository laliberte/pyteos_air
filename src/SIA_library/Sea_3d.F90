module sea_3d

!#########################################################################

!THIS MODULE COMPUTES THE PROPERTIES OF SEAWATER AS FUNCTIONS OF TEMPERATURE,
!PRESSURE AND DENSITY, IN PARTICULAR, THE ABSOLUTE SALINITY

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

!#########################################################################

!THIS MODULE REQUIRES THE LIBRARY MODULES:
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90
!     MATHS_0,       FILE MATHS_0.F90
!     FLU_1,         FILE FLU_1.F90
!     SAL_1,         FILE SAL_1.F90
!     SAL_2,         FILE SAL_2.F90
!     FLU_3A,        FILE FLU_3A.F90

!#########################################################################

use constants_0
use sal_2
use flu_3a

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: sea_sa_si, set_it_ctrl_salinity

!CONTROL PARAMETERS OF THE SALINITY ITERATION
integer ctrl_initialized

integer ctrl_mode_salinity
integer ctrl_loop_maximum
real*8 ctrl_init_salinity
real*8 ctrl_eps_exit_salinity

contains

!=========================================================================
function sea_sa_si(t_si, p_si, d_si)
!==========================================================================

!THIS FUNCTION COMPUTES THE ABSOLUTE SALINITY SA_SI IN KG/KG FROM
!THE ABSOLUTE TEMPERATURE T_SI IN K, THE ABSOLUTE PRESSURE P_SI IN PA
!AND THE SEAWATER DENSITY D_SI IN KG/M3

!CHECK VALUE WITH DEFAULT SETTINGS:
!SEA_SA_SI(273.15, 101325, 1028) = 3.50315257709E-02

implicit none

integer maxit
real*8 s, eps, sea_sa_si, t_si, p_si, d_si

sea_sa_si = errorreturn

if(t_si <= 0d0 .or. &
  p_si <= 0d0 .or. &
  d_si <= 0d0) return

call init_it_ctrl_salinity

Select Case (ctrl_mode_salinity)
  Case (0);     s = SO_salinity_si
  Case Default; s = ctrl_init_salinity
End Select

if(ctrl_loop_maximum == 0)then
  maxit = 100
elseif(ctrl_loop_maximum == -1)then
  sea_sa_si = s; return
elseif(ctrl_loop_maximum > 0)then
  maxit = ctrl_loop_maximum
else
  return
endif

eps = ctrl_eps_exit_salinity
if(eps == 0d0)return

!RUN THE ITERATION
sea_sa_si = salinity_iteration(t_si, p_si, d_si, s, maxit, eps)

end function

!=========================================================================
function salinity_iteration(t_si, p_si, d_si, s, maxit, eps)
!==========================================================================

!THIS FUNCTION CALCULATES THE ABSOLUTE SALINITY S FROM GIVEN VALUES
!OF THE ABSOLUTE TEMPERATURE, T_SI, IN K, ABSOLUTE PRESSURE, P_SI,
!IN PA AND THE SEAWATER DENSITY, D_SI, IN KG/M3

!OUTPUT: SALINITY_ITERATION IS THE RESULT FOR THE ABSOLUTE SALINITY IN KG/KG

!INPUT: T_SI:  ABSOLUTE TEMPERATURE IN K
!       P_SI:  ABSOLUTE PRESSURE IN PA
!       D_SI:  SEAWATER DENSITY IN KG/M3
!       S:     INITIAL SALINITY ESTIMATE IN KG/KG
!       MAXIT: MAXIUMUM NUMBER OF ITERATION STEPS
!       EPS:   SALINITY TOLERANCE, EPS < 0 MEANS RELATIVE ERROR

implicit none

integer it, maxit
real*8 salinity_iteration, t_si, p_si, d_si, s, eps
real*8 ds, dv, gs_p, gs_sp

salinity_iteration = errorreturn

if(eps <= 0d0) return
if(maxit <= 0d0) return

if(check_limits == 1) then
  !FLU_LIMITS
  if(t_si < flu_tmin .or. t_si > flu_tmax .or. &
  d_si < flu_dmin .or. d_si > flu_dmax) return
  !SAL_LIMITS
  if(t_si < sal_tmin .or. t_si > sal_tmax .or. &
  p_si < sal_pmin .or. p_si > sal_pmax) return
else
  if(t_si <= 0d0 .or. &
  p_si <= 0d0 .or. &
  d_si <= 0d0) return
endif

check_limits = check_limits - 1

!PURE WATER DENSITY
dv = liq_density_si(t_si, p_si)
If(dv == errorreturn .or. dv == 0) go to 999

!SALINE PART OF THE GIVEN SPECIFIC VOLUME
dv = 1d0 / d_si - 1d0 / dv

do it = 1, maxit

  !DERIVATIVES FOR NEWTON ITERATION
  gs_p = sal_g_si(0, 0, 1, s, t_si, p_si)
  If(gs_p == errorreturn) go to 999
  gs_sp = sal_g_si(1, 0, 1, s, t_si, p_si)
  If(gs_sp == errorreturn) go to 999
  If(gs_sp == 0) go to 999
  

  !COMPUTE SALINITY CORRECTION
  ds = (dv - gs_p) / gs_sp
  !UPDATE SALINITY ESTIMATE
  s = s + ds
  If(s < 0 .or. s >= 1) go to 999
  
  If((eps > 0d0 .and. abs(ds) < eps) .or. &
     (eps < 0d0 .and. abs(ds) < eps * s)) Then
    salinity_iteration = s
    go to 999
  End If
  
enddo

999 continue

check_limits = check_limits + 1

if(check_limits == 1) then
  if(s < sal_smin .or. s > sal_smax) then
    salinity_iteration = errorreturn
  endif
endif

end function

!==========================================================================
subroutine init_it_ctrl_salinity()
!==========================================================================

implicit none

if(ctrl_initialized == -1) return

ctrl_initialized = -1

!SET DEFAULT VALUES AND MODES FOR SALINITY ITERATION
ctrl_loop_maximum = 100
ctrl_mode_salinity = 0                !DEFAULT: SA = SO_SALINITY_SI = 36.16504 G/KG
ctrl_init_salinity = SO_salinity_si
ctrl_eps_exit_salinity = 0.0000001d0  !DEFAULT = 0.0001 G/KG

end subroutine

!==========================================================================
subroutine set_it_ctrl_salinity(key, value)
!==========================================================================

!THIS SUB SETS CONTROL PARAMETERS FOR THE NEWTON ITERATION USED TO COMPUTE
!POTENTIAL TEMPERATURE FROM REFERENCE PRESSURE

!KEY             VALUE
!IT_STEPS        0           SET ITERATION NUMBER TO DEFAULT (100)
!IT_STEPS        N > 0       SET ITERATION NUMBER TO N
!IT_STEPS       -1           DO NOT ITERATE, USE INITIAL VALUE

!INIT_SA         0           USE DEFAULT SA = 35.16504 G/KG TO START
!INIT_SA         T > 0       USE VALUE S AS SALINITY TO START

!TOL_SA          0           USE DEFAULT EXIT ACCURACY FOR SALINITY (0.1 MK)
!TOL_SA          EPS > 0     USE EPS AS EXIT ACCURACY FOR SALINITY

implicit none
character*13 key
real*8 value

call init_it_ctrl_salinity

select case(trim(key))

  case("it_steps")                !ITERATION STEPS
    if(value == 0d0)then
      ctrl_loop_maximum = 100     !DEFAULT = 100
    elseif(value < 0d0)then
      ctrl_loop_maximum = -1
    else
      ctrl_loop_maximum = value
    endif

  case("init_sa")                 !START SALINITY
    if(value == 0)then
      ctrl_mode_salinity = 0      !DEFAULT: SA = 35.16504 G/KG
    elseif(value > 0)then
      ctrl_mode_salinity = 1
      ctrl_init_salinity = value
    endif

  case("tol_sa")                  !REQUIRED SALINITY TOLERANCE
      if(value == 0d0)then
        ctrl_eps_exit_salinity = 0.0000001d0 !DEFAULT = 0.1 MG/KG
      elseif(value > 0)then
        ctrl_eps_exit_salinity = value
      endif

  case default; write(*,*)key,' not supported'

end select

end subroutine

end module
