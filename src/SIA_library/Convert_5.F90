module convert_5

!#########################################################################

!THIS MODULE IMPLEMENTS CONVERSION FORMULAS BETWEEN DIFFERENT SCALES AND UNITS

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

!THIS MODULE REQUIRES THE LIBRARY MODULE
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90

!#########################################################################

use constants_0
use convert_0

implicit none
private

character*16, private :: version = '8 Jun 2010'

public :: cnv_pressure, cnv_salinity, cnv_temperature

contains

!==========================================================================
function cnv_pressure(unit_out, p_in , unit_in, lat0)
!==========================================================================

!NOTE: UNIT_IN AND UNIT_OUT ARE CASE SENSITIVE. mPa IS NOT THE SAME AS MPA.
!      ALL UNITS MUST BE IN CAPITAL LETTERS

! CHECK VALUES:
! CNV_PRESSURE("MPA",4000,"DBAR")               = 40.101325
! CNV_PRESSURE("DBAR",10331325,"PA")            = 0.102300000000E+04
! CNV_PRESSURE("PA",1023,"DBAR")                = 0.103313250000E+08
! CNV_PRESSURE("TORR",10331325,"PA")            = 0.774913101406E+05
! CNV_PRESSURE("PA",0.774913101406d5,"TORR")    = 0.103313250000E+08
! CNV_PRESSURE("KGF",10331325,"PA")             = 0.105350196040E+06
! CNV_PRESSURE("PA",0.105350196040d6,"KGF")     = 0.103313250000E+08
! CNV_PRESSURE("ATM",10331325,"PA")             = 0.105350196040E+06
! CNV_PRESSURE("PA",0.105350196040d6,"ATM")     = 0.103313250000E+08
! CNV_PRESSURE("LBF/IN2",10331325,"PA")         = 0.149842272437E+04
! CNV_PRESSURE("PA",0.149842272437d4,"LBF/IN2") = 0.103313250000E+08
! CNV_PRESSURE("PSI",10331325,"PA")             = 0.149842272437E+04
! CNV_PRESSURE("PA",0.149842272437d4,"PSI")     = 0.103313250000E+08
! CNV_PRESSURE("M", 1023d0, "DBAR")             = 1011.94563591
! CNV_PRESSURE("DBAR", 1011.94563591d0, "M")    = 1023

character*(*) unit_out, unit_in
real*8 cnv_pressure, p_in
real*8 p_si, lat
real*8, optional :: lat0

cnv_pressure = ErrorReturn
p_si = ErrorReturn

if(.not. present(lat0))then
   lat = 45d0
else
   lat = lat0
endif

!CONVERT INPUT TO SI
Select Case(unit_in)
  Case ("PA"); p_si = p_in
  Case ("HPA"); p_si = p_in * 100d0
  Case ("MPA"); p_si = p_in * 1000000d0
  Case ("KPA"); p_si = p_in * 1000d0
  Case ("DBAR"); p_si = p_in * 10000d0 + 101325d0  !SEA PRESSURE ASSUMED
  Case ("MBAR"); p_si = p_in * 100d0
  Case ("BAR"); p_si = p_in * 100000d0
  Case ("KBAR"); p_si = p_in * 100000000d0
  Case ("TORR"); p_si = p_in * 101325d0 / 760d0
  Case ("MMHG"); p_si = p_in * 101325d0 / 760d0
  Case ("KGF"); p_si = p_in * 98.0665d0 !ATM = KG/CM2 GAUGE
  Case ("ATM"); p_si = p_in * 98.0665d0 !ATM = KG/CM2 GAUGE
  Case ("LBF/IN2"); p_si = p_in * 6894.8d0  !POUND PER SQUARE INCH
  Case ("PSI"); p_si = p_in * 6894.8d0  !POUND PER SQUARE INCH
  Case ("M"); p_si = 101325d0 + 1d4 * sw_pres(p_in, lat)
  Case default; write(*,*)'default P input--dbg'; return
End Select

!CONVERT SI TO OUTPUT
Select Case(unit_out)
  Case ("PA"); cnv_pressure = p_si
  Case ("HPA"); cnv_pressure = p_si / 100d0
  Case ("MPA"); cnv_pressure = p_si / 1000000d0
  Case ("KPA"); cnv_pressure = p_si / 1000d0
  Case ("DBAR"); cnv_pressure = (p_si - 101325d0) / 10000d0 !SEA PRESSURE ASSUMED
  Case ("MBAR"); cnv_pressure = p_si / 100d0
  Case ("BAR"); cnv_pressure = p_si / 100000d0
  Case ("KBAR"); cnv_pressure = p_si / 100000000d0
  Case ("TORR"); cnv_pressure = p_si / 101325d0 * 760d0
  Case ("MMHG"); cnv_pressure = p_si / 101325d0 * 760d0
  Case ("KGF"); cnv_pressure = p_si / 98.0665d0 !ATM = KG/CM2 GAUGE
  Case ("ATM"); cnv_pressure = p_si / 98.0665d0 !ATM = KG/CM2 GAUGE
  Case ("LBF/IN2"); cnv_pressure = p_si / 6894.8d0  !POUND PER SQUARE INCH
  Case ("PSI"); cnv_pressure = p_si / 6894.8d0  !POUND PER SQUARE INCH
  Case ("M"); cnv_pressure =  sw_dpth((p_si - 101325d0) / 10000d0, lat) !SEA PRESSURE ASSUMED
  Case default; write(*,*)'default P output--dbg';  return
End Select
                         
End Function

!==========================================================================
function cnv_temperature(unit_out, t_in, unit_in)
!==========================================================================

!NOTE: UNIT_IN AND UNIT_OUT ARE CASE SENSITIVE. DegC IS NOT THE SAME AS DEGC.
!      ALL UNITS MUST BE IN CAPITAL LETTERS

! CHECK VALUES:

! CNV_TEMPERATURE("DEGF(T48)",300,"K(T48)")           = 0.803300000000E+02
! CNV_TEMPERATURE("DEGF(T68)",300,"K(T48)")           = 0.803139948158E+02
! CNV_TEMPERATURE("DEGF(T90)",300,"K(T48)")           = 0.803018730496E+02
! CNV_TEMPERATURE("DEGC(T48)",300,"K(T48)")           = 0.268500000000E+02
! CNV_TEMPERATURE("DEGC(T68)",300,"K(T48)")           = 0.26841108231|2E+02
! CNV_TEMPERATURE("DEGC(T90)",300,"K(T48)")           = 0.26834373916|7E+02
! CNV_TEMPERATURE("K(T68)",300,"K(T48)")              = 0.299991108231E+03
! CNV_TEMPERATURE("K(T90)",300,"K(T48)")              = 0.29998437391|7E+03
! CNV_TEMPERATURE("DEGF(T68)",299.991108231,"K(T68)") = 0.803139916338E+02
! CNV_TEMPERATURE("DEGF(T90)",299.991108231,"K(T68)") = 0.803018698685E+02
! CNV_TEMPERATURE("DEGC(T68)",299.991108231,"K(T68)") = 0.268411064632E+02
! CNV_TEMPERATURE("DEGC(T90)",299.991108231,"K(T68)") = 0.268343721491E+02
! CNV_TEMPERATURE("K(T90)",299.991108231,"K(T68)")    = 0.299984372149E+03
! CNV_TEMPERATURE("DEGF(T90)",299.984372149,"K(T90)") = 0.803018698682E+02
! CNV_TEMPERATURE("DEGC(T90)",299.984373917,"K(T90)") = 0.268343739170E+02

character*(*) unit_out, unit_in
real*8 cnv_temperature, t_in, t_c
real*8 t_si

if(unit_out == unit_in) then; cnv_temperature = t_in; return; endif

cnv_temperature = ErrorReturn

!CONVERT INPUT TO T90
Select Case(unit_in)
  !Note: ISO(1993), ANNEX A, P. 86, GIVES T_F/DEGF = 32 + 1.8 * T_C/DEGC
  Case ("DEGF(T48)"); t_c = (t_in -32d0)/1.8d0; t_si = t_si_from_t48_k(t_k_from_t_c(t_c))
  Case ("DEGF(T68)"); t_c = (t_in -32d0)/1.8d0; t_si = t_si_from_t68_k(t_k_from_t_c(t_c))
  Case ("DEGF(T90)"); t_c = (t_in -32d0)/1.8d0; t_si = t_k_from_t_c(t_c)
  Case ("DEGC(T48)"); t_si = t_si_from_t48_k(t_k_from_t_c(t_in))
  Case ("DEGC(T68)"); t_si = t_si_from_t68_k(t_k_from_t_c(t_in))
  Case ("DEGC(T90)"); t_si = t_k_from_t_c(t_in)
  Case ("K(T48)"); t_si = t_si_from_t48_k(t_in)
  Case ("K(T68)"); t_si = t_si_from_t68_k(t_in)
  Case ("K(T90)"); t_si = t_in
  Case default; write(*,*)'TEMPERATURE INPUT UNIT NOT SUPPORTED'; return
End Select

!CONVERT T90 TO OUTPUT
Select Case(unit_out)
  Case ("DEGF(T48)")
    if(unit_in == "K(T48)")then
       cnv_temperature = 32d0 + 1.8d0 * t_c_from_t_k(t_in)
    elseif(unit_in == "DEGC(T48)")then
       cnv_temperature = 32d0 + 1.8d0 * t_in
    else
       write(*,*)'t48_f is not a valid output except for t48_k or t48_c input'; return
    endif
  Case ("DEGF(T68)"); cnv_temperature = 32d0 + 1.8d0 * t_c_from_t_k(t68_k_from_t_si(t_si))
  Case ("DEGF(T90)"); cnv_temperature = 32d0 + 1.8d0 * t_c_from_t_k(t_si)
  Case ("DEGC(T48)")
    if(unit_in == "K(T48)")then
       cnv_temperature = t_c_from_t_k(t_in)
    else
       write(*,*)'t48_c is not a valid output except for t48_k input'; return
    endif
  Case ("DEGC(T68)"); cnv_temperature = t_c_from_t_k(t68_k_from_t_si(t_si))
  Case ("DEGC(T90)"); cnv_temperature = t_c_from_t_k(t_si)
  Case ("K(T48)"); write(*,*)'t48_k is not a valid output'; return
  Case ("K(T68)"); cnv_temperature = t68_k_from_t_si(t_si)
  Case ("K(T90)"); cnv_temperature = t_si
  Case default; write(*,*)'TEMPERATURE OUTPUT UNIT NOT SUPPORTED';  return
End Select
                         
End Function

!==========================================================================
function cnv_salinity(unit_out, s_in, unit_in, t0_si, p0_si, lon0, lat0)
!==========================================================================

!NOTE: UNIT_IN AND UNIT_OUT ARE CASE SENSITIVE. Cl IS NOT THE SAME AS CL.

! CHECK VALUES:
! CNV_SALINITY("CL",0.0357,"KG/KG(ABS)",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.196660177563E+02
! CNV_SALINITY("KN",0.0357,"KG/KG(ABS)",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.355271620502E+02
! CNV_SALINITY("CND",0.0357,"KG/KG(ABS)",t_si=273.15+25.5,
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.127556680822E+01
! CNV_SALINITY("PSS",0.0357,"KG/KG(ABS)",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.355276443777E+02
! CNV_SALINITY("KG/KG(REF)",0.0357d0,"KG/KG(ABS)",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.356951724471E-01
! CNV_SALINITY("KG/KG(ABS)",0.0356951724471,"KG/KG(REF)",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.0357
! CNV_SALINITY("CL",35.52764437773386,"PSS")         = 0.196660177563E+02
! cnv_salinity("KN",35.52764437773386,"PSS")         = 0.355271620502E+02
! CNV_SALINITY("CND",35.52764437773386,"PSS",
!              t_si=273.15+25.5,p_si=101325+1023d4)  = 0.127556680822E+01
! CNV_SALINITY("KG/KG(REF)",35.52764437773386,"PSS") = 0.356951724471E-01
! CNV_SALINITY("KG/KG(ABS)",35.52764437773386,"PSS",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.0357
! CNV_SALINITY("KG/KG(ABS)", 35d0, "PSS", &
!              lon0=179d0, lat0=40d0, p0_si=101325d0+2d7), "3.51888478029404E-2"
! CNV_SALINITY("KG/KG(ABS)", 8d0, "PSS", &
!              lon0=20d0, lat0=57d0, p0_si=101325d0), "8.13338057142857E-03"
! CNV_SALINITY("PSS",0.355271620502E+02,"KN")        = 0.355276443777E+02
! CNV_SALINITY("KG/KG(ABS)",0.355271620502E+02,"KN",
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.357000000000E-01
! CNV_SALINITY("PSS",0.127556680822E+01,"CND",
!              t_si=273.15+25.5,p_si=101325+1023d4)  = 0.355226965418E+02
! CNV_SALINITY("KG/KG(ABS)",0.127556680822E+01,"CND",t_si=273.15+25.5,
!              lon0=201,lat0=-21,p_si=101325+1023d4)  = 0.0357

character*(*) unit_out, unit_in
real*8 cnv_salinity, s_in
real*8 t, p
real*8, optional :: lat0, lon0, t0_si, p0_si
real*8 :: t_si, p_si, ploc_si, lon, lat
real*8 ups, ucl, c35, t68_c, p_dbar, psal, s_si, c_ratio

if (unit_in == unit_out)then
   cnv_salinity = s_in
   return
endif

ups = SO_salinity_si / 35d0  !SR = SP * ups
ucl = 1.80655d0 * ups        !SR = Cl * ucl
c35 = 4.2914d0               !C(35,15,0) = 42.914 mS/CM = 4.2914 S/M,  CULKIN & SMITH

if(unit_in == unit_out) then
   cnv_salinity = s_in
   return
endif

cnv_salinity = ErrorReturn
s_si = ErrorReturn

if(present(lon0) .and. present(lat0) .and. present(p0_si)) then
   lon = lon0
   lat = lat0
   ploc_si = p0_si
else
   lon = -20d0
   lat = 50d0
   ploc_si = 101325d0
endif
if(present(t0_si)) then
   t_si = t0_si
else
   t_si = 298.15d0
endif
if(present(p0_si)) then
   p_si = p0_si
else
   p_si = 101325d0
endif
t = t_si - 273.15d0
p = (p_si - 101325d0)/1d4

!CONVERT INPUT TO REFERENCE SALINITY IN KG/KG

Select Case(unit_in)

  !Chlorinity
  Case ("CL"); s_si = s_in * ucl

  !Knudsen Salinity
  Case ("KN"); s_si = (s_in - 0.03d0) / 1.805d0 * ucl

  !Conductivity
  Case ("MS/CM")  !CONDUCTIVITY IN mS/CM = MMHO/CM
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert conductivity to salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     s_si = sal78_from_cnd(s_in, t, p) * ups
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     c_ratio = 0.1 * s_in / c35
     s_si = sal78_from_cnd(c_ratio, t68_c, p_dbar) * ups
  Case ("MMHO/CM")
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert conductivity to salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     s_si = sal78_from_cnd(s_in, t, p) * ups
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     c_ratio = 0.1 * s_in / c35
     s_si = sal78_from_cnd(c_ratio, t68_c, p_dbar) * ups
  Case ("S/M")
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert conductivity to salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     c_ratio = s_in / c35
     s_si = sal78_from_cnd(c_ratio, t68_c, p_dbar) * ups
  Case ("CND")
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert conductivity to salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     c_ratio = s_in
     s_si = sal78_from_cnd(c_ratio, t68_c, p_dbar) * ups

  !PRACTICAL SALINITY
  Case ("PSU"); s_si = s_in * ups
  Case ("PSS"); s_si = s_in * ups
  Case ("ONE"); s_si = s_in * ups

  !REFERENCE SALINITY
  Case ("KG/KG(REF)"); s_si = s_in
  Case ("G/KG(REF)"); s_si = 0.001d0 * s_in


  !ABSOLUTE SALINITY
  Case ("KG/KG(ABS)")
     if(.not. (present(lat0) .and. present(lon0) .and. present(p0_si))) then
        write(*,*) 'Position and pressure are needed to convert from or to Absolute Salinity'
        write(*,*) 'in order to determine the Absolute Salinity anomaly'
        write(*,*) 'Proceeding value from surface North Atlantic (near zero)'
     endif
     s_si = sp_from_sa(s_in, lon, lat, ploc_si) * ups
  Case ("G/KG(ABS)")
     if(.not. (present(lat0) .and. present(lon0) .and. present(p0_si))) then
        write(*,*) 'Position and pressure are needed to convert from or to Absolute Salinity'
        write(*,*) 'in order to determine the Absolute Salinity anomaly'
        write(*,*) 'Proceeding value from surface North Atlantic (near zero)'
     endif
     s_si = sp_from_sa(0.001d0 * s_in, lon, lat, ploc_si) * ups

  Case default; return
End Select


!CONVERT SI REFERENCE SALINITY TO OUTPUT
Select Case(unit_out)

  !Chlorinity
  Case ("CL"); cnv_salinity = s_si / ucl

  !Knudsen Salinity
  Case ("KN"); cnv_salinity = 0.03d0 + 1.805d0 * s_si / ucl

  !Conductivity
  Case ("MS/CM")  !conductivity in mS/cm = mmho/cm
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     psal = s_si / ups
     cnv_salinity = 10 * c35 * cnd_from_sal78(psal, t68_c, p_dbar)
  Case ("MMHO/CM")
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     psal = s_si / ups
     cnv_salinity = 10 * c35 * cnd_from_sal78(psal, t68_c, p_dbar)
  Case ("S/M")  !conductivity in S/m
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     psal = s_si / ups
     cnv_salinity = c35 * cnd_from_sal78(psal, t68_c, p_dbar)
  Case ("CND")    !conductivity ratio
     if(.not. present(p0_si)) then
        write(*,*) 'Warning: pressure is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with atmospheric pressure'
     endif
     if(.not. present(t0_si)) then
        write(*,*) 'Warning: temperature is needed to convert between conductivity and salinity'
        write(*,*) 'Proceeding with T90 = 273.15 + 25'
     endif
     t68_c = cnv_temperature("DEGC(T68)", t_si, "K(T90)")
     p_dbar = cnv_pressure("DBAR", p_si, "PA")
     psal = s_si / ups
     cnv_salinity = cnd_from_sal78(psal, t68_c, p_dbar)

  Case ("PSU"); cnv_salinity = s_si / ups
  Case ("PSS"); cnv_salinity = s_si / ups
  Case ("ONE"); cnv_salinity = s_si / ups

  !Reference Salinity
  Case ("KG/KG(REF)"); cnv_salinity = s_si
  Case ("G/KG(REF)"); cnv_salinity = 1d3 * s_si

  !Absolute Salinity
  Case ("KG/KG(ABS)")
     if(.not. (present(lat0) .and. present(lon0) .and. present(p0_si))) then
        write(*,*) 'Position and pressure are needed to convert from or to Absolute Salinity'
        write(*,*) 'in order to determine the Absolute Salinity anomaly'
        write(*,*) 'Proceeding value from surface North Atlantic (near zero)'
     endif
     psal = s_si / ups
     cnv_salinity = sa_from_sp(psal, lon, lat, ploc_si)
  Case ("G/KG(ABS)")
     if(.not. (present(lat0) .and. present(lon0) .and. present(p0_si))) then
        write(*,*) 'Position and pressure are needed to convert from or to Absolute Salinity'
        write(*,*) 'in order to determine the Absolute Salinity anomaly'
        write(*,*) 'Proceeding value from surface North Atlantic (near zero)'
     endif
     psal = s_si / ups
     cnv_salinity = 1d3 * sa_from_sp(psal, lon, lat, ploc_si)

  Case default; return
End Select
                         
End Function

!==========================================================================
function t_k_from_t_c(t_c)
!==========================================================================

!THIS FUNCTION COMPUTES TEMPERATURE IN K FROM TEMPERATURE IN deg C (ITS-90, IPTS-68, IPTS-48)

implicit none
real*8 t_k_from_t_c, t_c

t_k_from_t_c = 273.15d0 + t_c

end function

!==========================================================================
function t_c_from_t_k(t_k)
!==========================================================================

!THIS FUNCTION COMPUTES TEMPERATURE IN C FROM TEMPERATURE IN K (ITS-90, IPTS-68, IPTS-48)

implicit none
real*8 t_c_from_t_k, t_k

t_c_from_t_k = t_k - 273.15d0

end function

!==========================================================================
function t68_k_from_t48_k(t48_k)
!==========================================================================

!THIS FUNCTION RETURNS IPTS-68 TEMPERATURE IN K COMPUTED FROM IPTS-48 TEMPERATURE IN K.

!VB CODE ADAPTED FROM:

!****************************************************************************
!
!                          RUHR-UNIVERSITAET BOCHUM
!                        LEHRSTUHL FUER THERMODYNAMIK
!                          PROF. DR.-ING. W. WAGNER
!                              D - 44780 BOCHUM
!                             TEL.: 0234/700-3033
!                             FAX : 0234/7094 -163
!
!****************************************************************************

! PROGRAMM ZUR TEMPERATURUMRECHNUNG: IPTS-48 IN IPTS-68

! DER WERT !T48_K! IN DER IPTS-48 (IN KELVIN) WIRD VOM RUFENDEN PROGRAMM UEBENOMMEN
! UND IN T68_K_FROM_T48_K (IN KELVIN) WIEDER UEBERGEBEN.

implicit none
real*8 t68_k_from_t48_k, t48_k
integer i, k
real*8 a48a4, b48, c48c4, a48a68, b48b68, a68, b68
real*8 t48, t68, t68alt, tstri
real*8 hg1, hg2, hg3, fkt
real*8 a(21)

a48a4 = 0.003984517d0
b48 = -0.000000585502d0
c48c4 = -4.35716d-12
a48a68 = -0.00000019539d0
b48b68 = 0.0000000019539d0
a68 = 0.003984712d0
b68 = -0.000000587456d0

a(1) = 273.15d0
a(2) = 250.846209678803d0
a(3) = 135.099869965d0
a(4) = 52.78567590085d0
a(5) = 27.6768548854105d0
a(6) = 39.1053205376684d0
a(7) = 65.5613230578069d0
a(8) = 80.8035868559867d0
a(9) = 70.5242118234052d0
a(10) = 44.7847589638966d0
a(11) = 21.2525653556078d0
a(12) = 7.67976358170846d0
a(13) = 2.1368945938285d0
a(14) = 0.459843348928069d0
a(15) = 7.63614629231648d-02
a(16) = 9.69328620373121d-03
a(17) = 9.23069154007008d-04
a(18) = 6.38116590952654d-05
a(19) = 3.02293237874619d-06
a(20) = 8.7755139130376d-08
a(21) = 1.17702613125477d-09

if (t48_k < 90d0) then
  !MSGBOX "ACHTUNG: KEINE UMRECHNUNG T48T68 FUER T48_K < 90 K !"
  t68_k_from_t48_k = t48_k
  return
end if

if (t48_k > 273.15d0) goto 10

t48 = t48_k - 273.15d0
hg1 = 1d0 + a48a4*t48 + b48*t48*t48 + c48c4*t48*t48*t48*(t48 - 100d0)
hg1 = log(hg1)
t68 = 0d0
do k = 1, 21
  i = 22 - k
  t68 = t68 * hg1 + a(i)
enddo

t68_k_from_t48_k = t68
return

10    continue
t48 = t48_k - 273.15d0
hg1 = a68 + 2d0 * b68 * t48
hg2 = a48a68 * t48 + b48b68 * t48 * t48
t68 = t48
t68alt = t48
tstri = t48

11    continue
fkt = 0.045d0*tstri/100d0*(tstri/100d0 - 1d0)  &
      *(tstri/419.58d0 - 1d0)*(tstri/630.74d0 - 1d0)
hg3 = b68 * (t68 - t48 - fkt) * (t68 - t48 - fkt)
t68 = t48 + fkt + hg2 / hg1 - hg3 / hg1
if (abs(t68 - t68alt) < 1d-12) goto 12  !dbg NOTE 1D-5 REPLACED BY 1D-12
t68alt = t68
tstri = t68
goto 11

12    continue
t68_k_from_t48_k = t68 + 273.15d0

end function

!==========================================================================
function t_si_from_t48_k(t48_k)
!==========================================================================

implicit none 
real*8 t_si_from_t48_k, t48_k

t_si_from_t48_k = t_si_from_t68_k(t68_k_from_t48_k(t48_k))

end function

!==========================================================================
function t_si_from_t68_k(t68_k)
!==========================================================================

!THIS FUNCTION RETURNS ITS-90 TEMPERATURE IN K COMPUTED FROM IPTS-68 TEMPERATURE IN K.

!VB CODE ADAPTED FROM:

!****************************************************************************
!
!                          RUHR-UNIVERSITAET BOCHUM
!                        LEHRSTUHL FUER THERMODYNAMIK
!                          PROF. DR.-ING. W. WAGNER
!                              D - 44780 BOCHUM
!                             TEL.: 0234/700-3033
!                             FAX : 0234/7094 -163
!
!****************************************************************************
!
!   REGION 1:  13.8 K  < T68_K <   83.8 K
!
!   REGION 2:  83.8 K  < T68_K <  903.75 K
!
!   REGION 3: 903.89 K < T68 < 1337.58 K
!
!   REGION 4:     1337.58 < T68
!
!   INPUT:  T68_K      TEMPERATURE IN K (IPTS-68)
!
!   OUTPUT: T90_K_FROM_T68_K   TEMPERATURE IN K (ITS-90)
!
! FROM:
!
!  SUPPLEMENTARY INFORMATION FOR THE ITS-90, BUREAU INTERNATIONAL DES POIDS ET
!  MESURES, PAV BRETEUIL, 92312 SEVRES CEDEX, FRANCE, 1990.
!
!  R.L.RUSBY: THE CONVERSION OF THERMAL REFERENCE VALUES TO THE
!             ITS-90. J. CHEM. THERMODYNAMICS, 1991, 23, 1153 - 1161.
!
!  RUSBY ET AL.: INTERNATIONAL REPORT. REVISED VALUES FOR (T90-T68) FROM
!                630C TO 1064C. METROLOGIA, 1994, 31, 149 - 153.
!
!****************************************************************************


implicit none
real*8 t_si_from_t68_k, t68_k
integer j
real*8 deltat, tn
real*8 a(0:12), b(1:8), c(0:5)


a(0) = -0.005903d0
a(1) = 0.008174d0
a(2) = -0.061924d0
a(3) = -0.193388d0
a(4) = 1.490793d0
a(5) = 1.252347d0
a(6) = -9.835868d0
a(7) = 1.411912d0
a(8) = 25.277595d0
a(9) = -19.183815d0
a(10) = -18.437089d0
a(11) = 27.000895d0
a(12) = -8.716324d0

b(1) = -0.148759d0
b(2) = -0.267408d0
b(3) = 1.08076d0
b(4) = 1.269056d0
b(5) = -4.089591d0
b(6) = -1.871251d0
b(7) = 7.438081d0
b(8) = -3.536296d0

c(0) = 78.687209d0
c(1) = -0.47135991d0
c(2) = 0.0010954715d0
c(3) = -0.0000012357884d0
c(4) = 6.7736583d-10
c(5) = -1.4458081d-13

deltat = 0d0

if (t68_k < 13.8d0) goto 999

!*** REGION 4 ***
if (t68_k >= 1337.58d0) then
  deltat = -0.25d0 * (t68_k / 1337.33d0) ** 2
  goto 1000

!*** REGION 3 ***
!DIE GLEICHUNG IST IN ABHAENGIGKEIT VON T90 FORMULIERT, FUER DIE
!BERECHNUNG VON DELTAT WIRD HIER T68=T90 GESETZT /AK 16.7.96
elseif (t68_k >= 903.89d0) then
  tn = t68_k - 273.15d0
  do j = 1, 5
    deltat = deltat + c(j) * tn ** j
  enddo

  deltat = deltat + c(0)
  goto 1000

! *** REGION 2 ***
elseif (t68_k >= 83.8d0) then
  tn = (t68_k - 273.15d0) / 630d0
  do j = 1, 8
    deltat = deltat + b(j) * tn ** j
  enddo
  goto 1000

!*** REGION 1 ***
elseif (t68_k >= 13.8d0) then
  tn = (t68_k - 40d0) / 40d0
  do j = 1, 12
    deltat = deltat + a(j) * tn ** j
  enddo

  deltat = deltat + a(0)
  goto 1000
end if

999   write(*,*) "temperature out of range of validity, no conversion."
      t_si_from_t68_k = t68_k
      return

1000  t_si_from_t68_k = t68_k + deltat

end function

!==========================================================================
function t68_k_from_t_si(t_si)
!==========================================================================

!THIS FUNCTION RETURNS IPTS-68 TEMPERATURE IN K COMPUTED FROM ITS-90 TEMPERATURE IN K.

!VB CODE ADAPTED FROM:

!****************************************************************************
!
!                          RUHR-UNIVERSITAET BOCHUM
!                        LEHRSTUHL FUER THERMODYNAMIK
!                          PROF. DR.-ING. W. WAGNER
!                              D - 44780 BOCHUM
!                             TEL.: 0234/700-3033
!                             FAX : 0234/7094 -163
!
!****************************************************************************
!
!   REGION 1:  13.8 K  < T_SI <   83.8 K
!
!   REGION 2:  83.8 K  < T_SI <  903.75 K
!
!   REGION 3: 903.765 K < T_SI < 1337.33 K
!
!   REGION 4:     1337.33 < T_SI
!
!   INPUT:  T_SI      TEMPERATURE IN K (ITS-90)
!
!   OUTPUT: T68_K_FROM_T_SI   TEMPERATURE IN K (IPTS-68)
!
! FROM:
!
!  SUPPLEMENTARY INFORMATION FOR THE ITS-90, BUREAU INTERNATIONAL DES POIDS ET
!  MESURES, PAV BRETEUIL, 92312 SEVRES CEDEX, FRANCE, 1990.
!
!  R.L.RUSBY: THE CONVERSION OF THERMAL REFERENCE VALUES TO THE
!             ITS-90. J. CHEM. THERMODYNAMICS, 1991, 23, 1153 - 1161.
!
!  RUSBY ET AL.: INTERNATIONAL REPORT. REVISED VALUES FOR (T90-T68) FROM
!                630C TO 1064C. METROLOGIA, 1994, 31, 149 - 153.
!
!****************************************************************************

implicit none 
integer j
real*8 t68_k_from_t_si, t_si
real*8 deltat, tn
real*8 a(0:12), b(1:8), c(0:5)

a(0) = -0.005903d0
a(1) = 0.008174d0
a(2) = -0.061924d0
a(3) = -0.193388d0
a(4) = 1.490793d0
a(5) = 1.252347d0
a(6) = -9.835868d0
a(7) = 1.411912d0
a(8) = 25.277595d0
a(9) = -19.183815d0
a(10) = -18.437089d0
a(11) = 27.000895d0
a(12) = -8.716324d0

b(1) = -0.148759d0
b(2) = -0.267408d0
b(3) = 1.08076d0
b(4) = 1.269056d0
b(5) = -4.089591d0
b(6) = -1.871251d0
b(7) = 7.438081d0
b(8) = -3.536296d0

c(0) = 78.687209d0
c(1) = -0.47135991d0
c(2) = 0.0010954715d0
c(3) = -0.0000012357884d0
c(4) = 6.7736583d-10
c(5) = -1.4458081d-13

deltat = 0d0
if (t_si < 13.8d0) goto 999

!*** REGION 4 ***
if (t_si >= 1337.33d0) then
  deltat = -0.25d0 * (t_si / 1337.33d0) ** 2
  goto 1000

!*** REGION 3 ***
elseif (t_si >= 903.765d0) then
  tn = t_si - 273.15d0
  do j = 1, 5
    deltat = deltat + c(j) * tn ** j
  enddo
  deltat = deltat + c(0)
  goto 1000

!*** REGION 2 ***
elseif (t_si >= 83.8d0) then
  tn = (t_si - 273.15d0) / 630d0
  do j = 1, 8
    deltat = deltat + b(j) * tn ** j
  enddo
  goto 1000

!*** REGION 1 ***
elseif (t_si >= 13.8d0) then
  tn = (t_si - 40d0) / 40d0
  do j = 1, 12
    deltat = deltat + a(j) * tn ** j
  enddo

  deltat = deltat + a(0)
  goto 1000
end if

999   write(*,*) "temperature out of range of validity, no conversion."
      t68_k_from_t_si = t_si
      return

1000  t68_k_from_t_si = t_si - deltat

end function

!==========================================================================
function p_si_from_psea_pa(psea_pa)
!==========================================================================

!THIS FUNCTION COMPUTES ABSOLUTE PRESSURE IN PA FROM SEA PRESSURE IN PA

implicit none
real*8 p_si_from_psea_pa, psea_pa

p_si_from_psea_pa = psea_pa + 101325d0

end function

!==========================================================================
function psea_pa_from_p_si(p_si)
!==========================================================================

!THIS FUNCTION COMPUTES SEA PRESSURE IN PA FROM ABSOLUTE PRESSURE IN PA

implicit none 
real*8 psea_pa_from_p_si, p_si

psea_pa_from_p_si = p_si - 101325d0

end function

!==========================================================================
function psea_pa_from_psea_bar(psea_bar)
!==========================================================================

!THIS FUNCTION COMPUTES SEA PRESSURE IN PA FROM ABSOLUTE PRESSURE IN PA

implicit none
real*8 psea_pa_from_psea_bar, psea_bar

psea_pa_from_psea_bar = psea_bar * 100000d0

end function

!==========================================================================
function psea_bar_from_psea_pa(psea_pa)
!==========================================================================

!THIS FUNCTION COMPUTES SEA PRESSURE IN PA FROM ABSOLUTE PRESSURE IN PA

implicit none
real*8 psea_bar_from_psea_pa, psea_pa

psea_bar_from_psea_pa = psea_pa * 0.00001d0

end function

!=====================================================
function sal78_from_cnd(cnd, t, p)
!=====================================================

!THE CONDUCTIVITY RATIO (CND)=1.0000000 FOR SALINITY=35 PSS-78
!TEMPERATURE=15.0 DEG. CELSIUS AND ATMOSPHERIC PRESSURE.
!FUNCTION TO CONVERT CONDUCTIVITY RATIO TO SALINITY

!REFERENCES: ALSO LOCATED IN UNESCO REPORT NO. 37 1981
!PRACTICAL SALINITY SCALE 1978: E.L. LEWIS IEEE OCEAN ENG. JAN. 1980

!----------------------------------------------------------
!UNITS:
!  PRESSURE      P      DECIBARS
!  TEMPERATURE   T      DEG CELSIUS IPTS-68
!  CONDUCTIVITY  CND    RATIO
!  SALINITY      SAL78  PSS-78
!----------------------------------------------------------
!CHECKVALUES:
!  SAL78=40.00000 for CND=1.8880911555787682, T=40 DEG C, P=10000 DECIBARS
!----------------------------------------------------------
!SAL78 RATIO: RETURNS ZERO FOR CONDUCTIVITY RATIO: < 0.0005
!----------------------------------------------------------
!ORIGINAL FORTRAN CODE IS FOUND IN:
!  UNESCO TECHNICAL PAPERS IN MARINE SCIENCE 44 (1983) -
!  'ALGORITHMS FOR COMPUTATION OF FUNDAMENTAL PROPERTIES OF SEAWATER'
!----------------------------------------------------------
!TRANSLATED TO OBJECT PASCAL BY:
!  DR. JAN SCHULZ, 19 MAY 2008, WWW.CODE10.INFO
!----------------------------------------------------------
!TRANSLATED BACK TO FORTRAN BY:
!  D. WRIGHT, BIO, 23 SEPT., 2009
!----------------------------------------------------------

real*8 sal78_from_cnd, cnd, t, p
real*8 dt, res, rt

sal78_from_cnd = 0d0

! ZERO SALINITY/CONDUCTIVITY TRAP
if (cnd <= 5d-4) return

dt = t - 15d0

! CONVERT CONDUCTIVITY TO SALINITY
res   = cnd
rt    = res / (rt35(t) * (1d0 + c(p) / (b(t) + a(t) * res)))
rt    = sqrt (abs(rt))
sal78_from_cnd = sal (rt, dt)

end function


!=====================================================
function cnd_from_sal78 (sal78, t, p)
!=====================================================

!THE CONDUCTIVITY RATIO (CND)=1.0000000 FOR SALINITY=35 PSS-78
!TEMPERATURE=15.0 DEG. CELSIUS AND ATMOSPHERIC PRESSURE.
!Function to convert SALINITY TO CONDUCTIVITY RATIO

!REFERENCES: ALSO LOCATED IN UNESCO REPORT NO. 37 1981
!PRACTICAL SALINITY SCALE 1978: E.L. LEWIS IEEE OCEAN ENG. JAN. 1980

!----------------------------------------------------------
!UNITS:
!  PRESSURE      P      DECIBARS
!  TEMPERATURE   T      DEG CELSIUS IPTS-68
!  CONDUCTIVITY  CND    RATIO
!  SALINITY      SAL78  PSS-78
!----------------------------------------------------------
!CHECKVALUES:
!  CND=1.8880911555787682 for SAL78=40.0000, T=40 DEG C, P=10000 DECIBARS
!----------------------------------------------------------
!CND: RETURNS ZERO FOR SALINITY: < 0.02
!----------------------------------------------------------
!ORIGINAL FORTRAN CODE IS FOUND IN:
!  UNESCO TECHNICAL PAPERS IN MARINE SCIENCE 44 (1983) -
!  'ALGORITHMS FOR COMPUTATION OF FUNDAMENTAL PROPERTIES OF SEAWATER'
!----------------------------------------------------------
!TRANSLATED TO OBJECT PASCAL BY:
!  DR. JAN SCHULZ, 19 MAY 2008, WWW.CODE10.INFO
!----------------------------------------------------------
!TRANSLATED BACK TO FORTRAN BY:
!  D. WRIGHT, BIO, 23 SEPT., 2009
!----------------------------------------------------------

integer n
real*8 cnd_from_sal78, cnd, sal78, t, p
real*8 dt, rt, si, dels, rtt, at, bt, cp, res

cnd_from_sal78 = errorreturn
cnd = 0d0

!ZERO SALINITY/CONDUCTIVITY TRAP
if (sal78 <= 0.2d0) return

dt = t - 15d0

!INVERT SALINITY TO CONDUCTIVITY BY THE
!NEWTON-RAPHSON ITERATIVE METHOD

!FIRST APPROXIMATION
rt = sqrt(sal78 / 35d0)
si = sal(rt, dt)

!ITERATION LOOP BEGINS HERE WITH A MAXIMUM OF 20 CYCLES
do n = 1, 20
   rt   = rt + (sal78 - si) / dsal(rt, dt)
   si   = sal(rt, dt)
   dels = abs(si - sal78)
   if (n.gt.1 .and. dels < 1d-20) go to 99
enddo

99 continue

!COMPUTE CONDUCTIVITY RATIO
rtt = rt35(t) * rt * rt
at  = a(t)
bt  = b(t)
cp  = c(p)
cp  = rtt * (cp + bt)
bt  = bt - rtt * at

!SOLVE QUADRATIC EQUATION FOR R: R=RT35*RT*(1+C/AR+B)
res = sqrt(abs (bt * bt + 4d0 * at * cp)) - bt

!RETURN CONDUCTIVITY
cnd_from_sal78 = 0.5d0 * res / at

end function


!=====================================================
function sal(xr, xt)
!=====================================================

!PRACTICAL SALINITY SCALE 1978 DEFINITION WITH TEMPERATURE
!CORRECTION;XT :=T-15.0; XR:=SQRT(RT)

real*8 sal, xr, xt

sal = ((((2.7081d0*xr-7.0261d0)*xr+14.0941d0)*xr+25.3851d0)*xr &
    - 0.1692d0)*xr+0.0080d0 &
    + (xt/(1d0+0.0162d0*xt))*(((((-0.0144d0*xr &
    + 0.0636d0)*xr-0.0375d0)*xr-0.0066d0)*xr &
    - 0.0056d0)*xr+0.0005d0)

end function


!=====================================================
function dsal(xr, xt)
!=====================================================

!FUNCTION FOR DERIVATIVE OF SAL(XR,XT) WITH XR

real*8 dsal, xr, xt

dsal = ((((13.5405d0 * xr - 28.1044d0) * xr + 42.2823d0) * xr + 50.7702d0) &
     * xr - 0.1692d0) + (xt / (1d0 + 0.0162d0 * xt)) * ((((-0.0720d0 &
     * xr + 0.2544d0) * xr - 0.1125d0) * xr - 0.0132d0) * xr - 0.0056d0)

end function


!=====================================================
function rt35(xt)
!=====================================================

!FUNCTION RT35: C(35,T,0)/C(35,15,0) VARIATION WITH TEMPERATURE

real*8 rt35, xt

rt35 = (((1.0031d-9 * xt - 6.9698d-7) * xt + 1.104259d-4) * xt &
     + 2.00564d-2) * xt + 0.6766097d0

end function


!=====================================================
function c(xp)
!=====================================================

!C(XP) POLYNOMIAL CORRESPONDS TO A1-A3 CONSTANTS: LEWIS 1980

real*8 c, xp

c = ((3.989d-15 * xp - 6.370d-10) * xp + 2.070d-5) * xp

end function


!=====================================================
function b(xt)
!=====================================================

real*8 b, xt

b = (4.464d-4 * xt + 3.426d-2) * xt + 1d0

end function


!=====================================================
function a(xt)
!=====================================================

!A(XT) POLYNOMIAL CORRESPONDS TO B3 AND B4 CONSTANTS: LEWIS 1980

real*8 a, xt

a = -3.107d-3 * xt + 0.4215d0

end function

!=========================================================================
function sw_dpth(P,LAT)
!=========================================================================

! SW_DPTH    DEPTH FROM PRESSURE
!=========================================================================
! SW_DPTH   IS A MODIFICATION OF THE MATLAB ROUTINE WRITTEN BY PHIL MORGAN
!           OF CSIRO FOR THE SEAWATER LIBRARY DEVELOPED BY HIM AROUND 1992.
!
! USAGE:  DPTH = SW_DPTH(P,LAT)
!
! DESCRIPTION:
!    CALCULATES DEPTH IN METRES FROM PRESSURE IN DBARS. 
!    BASED ON EQN 25, P26.  UNESCO 1983.
!
! INPUT:
!   P   = PRESSURE    [DB]
!   LAT = LATITUDE IN DECIMAL DEGRESS NORTH [-90..+90]
!
! OUTPUT:
!  DPTH = DEPTH [METRES]
!
! DISCLAIMER:
!   THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND.
!
! REFERENCES: 
!    SAUNDERS, P.M. 1981
!    "PRACTICAL CONVERSION OF PRESSURE TO DEPTH"
!    JOURNAL OF PHYSICAL OCEANOGRAPHY, 11, 573-574
!=========================================================================

real*8 p, sw_dpth
real*8 deg2rad, x, c1
real*8, optional :: lat

if(.not. present(lat))lat=45d0

sw_dpth = errorreturn

! CHECK INPUTS
if(p .lt. 0d0 .or. abs(lat) .gt. 90d0)return

! The following commented code is the original
! It is not an exact inverse of sw_press
!deg2rad = 3.1415926535d0/180d0
!c1 = +9.72659d0
!c2 = -2.2512d-5
!c3 = +2.279d-10
!c4 = -1.82d-15
!gam_dash = 2.184d-6
!x   = sin(lat*deg2rad)**2  ! CONVERT TO RADIANS
!sw_dpth = (((c4*p+c3)*p+c2)*p+c1)*p
!sw_dpth = sw_dpth / (9.780318d0*(1d0+(5.2788d-3+2.36d-5*x)*x) + gam_dash*0.5d0*p)

! The following code is simply the analyticalinverse of swpres
deg2rad = 3.1415926535d0/180d0
x       = sin(abs(lat)*deg2rad)  ! convert to radians
c1      = 5.92d-3 + x**2 * 5.25d-3
sw_dpth = ((1d0-c1)**2 - (1d0 - c1 - 4.42d-6*p)**2)/8.84d-6

return
end function

!===========================================================================
function sw_pres(depth,lat)
!===========================================================================

! SW_PRES    PRESSURE FROM DEPTH
!===========================================================================
! SW_PRES   IS A MODIFICATION OF THE MATLAB ROUTINE WRITTEN BY PHIL MORGAN
!           OF CSIRO FOR THE SEAWATER LIBRARY DEVELOPED BY HIM AROUND 1992.
!
! USAGE:  PRES = SW_PRES(DEPTH,LAT)
!
! DESCRIPTION:
!    CALCULATES PRESSURE IN DBARS FROM DEPTH IN METERS.
!
! INPUT: 
!   DEPTH = DEPTH [METRES]
!   LAT   = LATITUDE IN DECIMAL DEGRESS NORTH [-90..+90]
!
! OUTPUT:
!  PRES   = PRESSURE    [DB]
!
! DISCLAIMER:
!   THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND.
!
! REFERENCES:
!    SAUNDERS, P.M. 1981
!    "PRACTICAL CONVERSION OF PRESSURE TO DEPTH"
!    JOURNAL OF PHYSICAL OCEANOGRAPHY, 11, 573-574
!
! CHECK VALUE:
!    P=7500.00 DB for LAT=30 DEG, DEPTH=7321.45 METERS
!=========================================================================

real*8 depth, sw_pres
real*8 deg2rad, c1, x
real*8, optional :: lat

if(.not. present(lat))lat=45d0

sw_pres = errorreturn

! CHECK INPUTS
if(depth .lt. 0d0 .or. abs(lat) .gt. 90d0)return

deg2rad = 3.1415926535d0/180d0
x       = sin(abs(lat)*deg2rad)  ! convert to radians
c1      = 5.92d-3 + x**2 * 5.25d-3
sw_pres = ((1d0-c1) - sqrt((1d0-c1)**2 - 8.84d-6*depth))/4.42d-6
return
end function

end module convert_5
