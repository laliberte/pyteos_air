module values_cnv

!THIS PROGRAM PROVIDES AN EASY WAY FOR THE USER TO COMPARE RESULTS FROM 
!THE CONVERT ROUTINES RUN LOCALLY WITH THE CHECK VALUES PROVIDED IN 
!THE VARIOUS ROUTINES

!THIS PROGRAM REQUIRES
!     CONSTANTS_0,   FILE CONSTANTS_0.F90
!     CONVERT_0,     FILE CONVERT_0.F90
!     CONVERT_5,     FILE CONVERT_5.F90

use constants_0
use convert_0
use convert_5

implicit none 
character*16, private :: version = '8 Jun 2010'

contains

!==========================================================================
subroutine cnv_chk_values
!==========================================================================

implicit none 
character*13 key
character (len = *), parameter :: fmt = "(A,T74,1X,E19.12,3X,A)"

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR THE PRESSURE CONVERSION ROUTINE'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,fmt)'cnv_pressure("DBAR",10331325,"PA") = ', &
 cnv_pressure("DBAR", 10331325d0, "PA"),  &
 "1023"

write(*,fmt)'cnv_pressure("PA",1023,"DBAR") = ', &
 cnv_pressure("PA", 1023d0, "DBAR"),  &
 "10331325"

write(*,fmt)'cnv_pressure("TORR",10331325,"PA") = ', &
 cnv_pressure("TORR", 10331325d0, "PA"),  &
 "77491.3101406"

write(*,fmt)'cnv_pressure("PA",0.774913101406d5,"TORR") = ', &
 cnv_pressure("PA", 0.774913101406d5, "TORR"),  &
 "10331325"

write(*,fmt)'cnv_pressure("KGF",10331325,"PA") = ', &
 cnv_pressure("KGF", 10331325d0, "PA"),  &
 "105350.196040"

write(*,fmt)'cnv_pressure("PA",0.105350196040d6,"KGF") = ', &
 cnv_pressure("PA", 0.105350196040d6, "KGF"),  &
 "10331325"

write(*,fmt)'cnv_pressure("ATM",10331325,"PA") = ', &
 cnv_pressure("ATM", 10331325d0, "PA"),  &
 "105350.196040"

write(*,fmt)'cnv_pressure("PA",0.105350196040d6,"ATM") = ', &
 cnv_pressure("PA", 0.105350196040d6, "ATM"),  &
 "10331325"

write(*,fmt)'cnv_pressure("LBF/IN2",10331325,"PA") = ', &
 cnv_pressure("LBF/IN2", 10331325d0, "PA"),  &
 "1498.42272437"

write(*,fmt)'cnv_pressure("PA",0.149842272437d4,"LBF/IN2") = ', &
 cnv_pressure("PA", 0.149842272437d4, "LBF/IN2"),  &
 "10331325"

write(*,fmt)'cnv_pressure("PSI",10331325,"PA") = ', &
 cnv_pressure("PSI", 10331325d0, "PA"),  &
 "1498.42272437"

write(*,fmt)'cnv_pressure("PA",0.149842272437d4,"PSI") = ', &
 cnv_pressure("PA", 0.149842272437d4, "PSI"),  &
 "10331325"

write(*,fmt)'cnv_pressure("M",1023d0,"DBAR") = ', &
 cnv_pressure("M", 1023d0, "DBAR"),  &
 "1011.94563591"

write(*,fmt)'cnv_pressure("DBAR",1023d0,"M") = ', &
 cnv_pressure("DBAR", 1011.94563591d0, "M"),  &
 "1023"



write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR THE TEMPERATURE CONVERSION ROUTINE'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

write(*,fmt)'cnv_temperature("DEGF(T48)",300,"K(T48)") = ', &
 cnv_temperature("DEGF(T48)", 300d0, "K(T48)"),  &
 "80.33"

write(*,fmt)'cnv_temperature("DEGF(T68)",300,"K(T48)") = ', &
 cnv_temperature("DEGF(T68)", 300d0, "K(T48)"),  &
 "80.3139948158"

write(*,fmt)'cnv_temperature("DEGF(T90)",300,"K(T48)") = ', &
 cnv_temperature("DEGF(T90)", 300d0, "K(T48)"),  &
 "80.3018730496"

write(*,fmt)'cnv_temperature("DEGC(T48)",300,"K(T48)") = ', &
 cnv_temperature("DEGC(T48)", 300d0, "K(T48)"),  &
 "26.85"

write(*,fmt)'cnv_temperature("DEGC(T68)",300,"K(T48)") = ', &
 cnv_temperature("DEGC(T68)", 300d0, "K(T48)"),  &
 "26.841108231|2"

write(*,fmt)'cnv_temperature("DEGC(T90)",300,"K(T48)") = ', &
 cnv_temperature("DEGC(T90)", 300d0, "K(T48)"),  &
 "26.834373916|7"

write(*,fmt)'cnv_temperature("K(T68)",300,"K(T48)") = ', &
 cnv_temperature("K(T68)", 300d0, "K(T48)"),  &
 "299.991108231"

write(*,fmt)'cnv_temperature("K(T90)",300,"K(T48)") = ', &
 cnv_temperature("K(T90)", 300d0, "K(T48)"),  &
 "299.98437391|7"

write(*,fmt)'cnv_temperature("DEGF(T68)",299.991108231,"K(T68)") = ', &
 cnv_temperature("DEGF(T68)", 299.991108231d0, "K(T68)"),  &
 "80.3139916338"

write(*,fmt)'cnv_temperature("DEGF(T90)",299.991108231,"K(T68)") = ', &
 cnv_temperature("DEGF(T90)", 299.991108231d0, "K(T68)"),  &
 "80.3018698685"

write(*,fmt)'cnv_temperature("DEGC(T68)",299.991108231,"K(T68)") = ', &
 cnv_temperature("DEGC(T68)", 299.991108231d0, "K(T68)"),  &
 "26.8411064632"

write(*,fmt)'cnv_temperature("DEGC(T90)",299.991108231,"K(T68)") = ', &
 cnv_temperature("DEGC(T90)", 299.991108231d0, "K(T68)"),  &
 "26.8343721491"

write(*,fmt)'cnv_temperature("K(T90)",299.991108231,"K(T68)") = ', &
 cnv_temperature("K(T90)", 299.991108231d0, "K(T68)"),  &
 "299.984372149"

write(*,fmt)'cnv_temperature("DEGF(T90)",299.984372149,"K(T90)") = ', &
 cnv_temperature("DEGF(T90)",299.984372149d0, "K(T90)"),  &
 "80.3018698682"

write(*,fmt)'cnv_temperature("DEGC(T90)",299.984373917,"K(T90)") = ', &
 cnv_temperature("DEGC(T90)",299.984373917d0, "K(T90)"),  &
 "26.8343739170"

write(*,"(A)")' '
write(*,"(A)")' '
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'COMPARISON OF CHECK VALUES FOR THE SALINITY CONVERSION ROUTINE'
write(*,"(A)")'**************************************************************'
write(*,"(A)")'**************************************************************'
write(*,"(A)")' '
write(*,"(A)")' '

if(.true.) then
write(*,fmt)'asal_from_psal(35.52764437773386,201,-21,P0+1023d4) = ', &
 cnv_salinity("KG/KG(ABS)", 35.52764437773386d0, "PSS", &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.0357"

write(*,fmt)'asal_from_psal(35,180,40,P0+2d7) = ', &
 cnv_salinity("KG/KG(ABS)", 35d0, "PSS", &
 p0_si=101325d0+2d7, lon0=180d0, lat0=40d0), "3.51888478029E-2"

write(*,fmt)'asal_from_psal(8,20,57,101325d0) = ', &
 cnv_salinity("KG/KG(ABS)", 8d0, "PSS", &
 p0_si=101325d0, lon0=20d0, lat0=57d0), "8.13338057143E-03"

write(*,fmt)'psal_from_asal(0.0357,201,-21,P0+1023d4) = ', &
 cnv_salinity("PSS", 0.0357d0, "KG/KG(ABS)", &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "35.5276443777"
write(*,*)" "
endif


if(.false.) then
write(*,*)"***********************************************************************"
write(*,*)"The following shows the comparison with results not including silicate"
write(*,*)"***********************************************************************"
write(*,fmt)'asal_from_psal(35.52764437773386, 201, -21, P0+1023d4) = ', &
 cnv_salinity("KG/KG(ABS)", 35.52764437773386d0, "PSS", &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "3.56951724471E-02"

write(*,fmt)'asal_from_psal(35,180,40,P0+2d7) = ', &
 cnv_salinity("KG/KG(ABS)", 35d0, "PSS", &
 p0_si=101325d0+2d7, lon0=180d0, lat0=40d0), "0.03516504"       !north Pacific (without silicate)

write(*,fmt)'asal_from_psal(8,20,57,101325d0) = ', &
 cnv_salinity("KG/KG(ABS)", 8d0, "PSS", &
 p0_si=101325d0, lon0=20d0, lat0=57d0), "8.13338057143E-03"  !central Baltic Sea

write(*,fmt)'psal_from_asal(0.0357,201,-21,P0+1023d4) = ', &
 cnv_salinity("PSS", 0.0357d0, "KG/KG(ABS)", &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "35.5324492735"
write(*,*)" "
endif

write(*,fmt)'cnv_salinity("CL",0.0357,"KG/KG(ABS)",P0+1023d4,201,-21) = ', &
 cnv_salinity("CL", 0.0357d0, "KG/KG(ABS)",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.196660177563E+02"

write(*,fmt)'cnv_salinity("KN",0.0357,"KG/KG(ABS)",P0+1023d4,201,-21) = ', &
 cnv_salinity("KN", 0.0357d0, "KG/KG(ABS)",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.355271620502E+02"

write(*,fmt)'cnv_salinity("CND",0.0357,"KG/KG(ABS)",P0+1023d4,201,-21) = ', &
 cnv_salinity("CND", 0.0357d0, "KG/KG(ABS)",  &
 t0_si = 273.15d0 + 25.5d0, p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.127556680822E+01"

write(*,fmt)'cnv_salinity("PSS",0.0357,"KG/KG(ABS)",P0+1023d4,201,-21) = ', &
 cnv_salinity("PSS", 0.0357d0, "KG/KG(ABS)",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.355276443777E+02"

write(*,fmt)'cnv_salinity("KG/KG(REF)",0.0357d0,"KG/KG(ABS)",P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(REF)",0.0357d0, "KG/KG(ABS)",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.356951724471E-01"

write(*,fmt)'cnv_salinity("KG/KG(ABS)",.0356951724471,"KG/KG(REF)",P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(ABS)", 0.0356951724471d0, "KG/KG(REF)",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.0357"

write(*,fmt)'cnv_salinity("CL",35.52764437773386,"PSS",P0+1023d4,201,-21) = ', &
 cnv_salinity("CL",35.52764437773386d0, "PSS"), "0.196660177563E+02"

write(*,fmt)'cnv_salinity("KN",35.52764437773386,"PSS",101325+1023d4,201,-21) = ', &
 cnv_salinity("KN",35.52764437773386d0, "PSS"), "0.355271620502E+02"

write(*,fmt)'cnv_salinity("CND",35.52764437773386,"PSS",P0+1023d4,201,-21) = ', &
 cnv_salinity("CND",35.52764437773386d0, "PSS",  &
 t0_si = 273.15d0 + 25.5d0, p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.127556680822E+01"

write(*,fmt)'cnv_salinity("KG/KG(REF)",35.52764437773386,"PSS",P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(REF)",35.52764437773386d0, "PSS",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.356951724471E-01"

write(*,fmt)'cnv_salinity("KG/KG(ABS)",35.52764437773386,"PSS",P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(ABS)",35.52764437773386d0, "PSS",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.0357"

write(*,fmt)'cnv_salinity("PSS",0.355271620502E+02,"KN") = ', &
 cnv_salinity("PSS", 0.355271620502d2, "KN",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.355276443777E+02"

write(*,fmt)'cnv_salinity("KG/KG(ABS)",0.355271620502E+02,"KN",P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(ABS)", 0.355271620502d2, "KN",  &
 p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.0357"

write(*,fmt)'cnv_salinity("PSS",1.27540928136,"CND",25.5,P0+1023d4,201,-21) = ', &
 cnv_salinity("PSS", 0.127540928136d1, "CND",  &
 t0_si=273.15d0+25.5d0, p0_si=101325d0+1023d4), "0.355226965418E+02"

write(*,fmt)'cnv_salinity("KG/KG(ABS)",0.127556680822d1,"CND",25.5,P0+1023d4,201,-21) = ', &
 cnv_salinity("KG/KG(ABS)", 0.127556680822d1, "CND",  &
 t0_si = 273.15d0 + 25.5d0, p0_si=101325d0+1023d4, lon0=201d0, lat0=-21d0), "0.357000000001E-01"

end subroutine


end module values_cnv
