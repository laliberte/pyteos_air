!FBL October 2012, University of Toronto

!This f2py interface to the python code is structured the following way:
! 1. a subroutine for each type of inputs (e.g. a_si, t_si, p_si) is created
!    and can be called with the function name and the numerical inputs.
! 2. for each subroutine, there is a companion subroutine with the suffix _desc
!    that serves to communicate with python. It lists the available functions
!    and give a simple description of the outputs with units.
! 3. the numerical inputs are three-dimensional, to minimize calls to the fortran
!    interface when applying a function point-wise on an array. The python bindings
!    make sure that the numerical inputs are broadcast to the right size, making most
!    these functions act like typical numpy functions for arrays with 3 or less 
!    dimensions.

subroutine liq_ice_air_sat(in_name,in_t,in_p, k, l, m, out_val)
    !This subroutine is a wrapper for the Gibbs function formulation for wet air saturation.
    !Available routines:
    !
    ! Liq_Air_4a:
    ! massfraction_air

    use liq_ice_air_5b
    implicit none
    integer :: k, l, m, idk, idl, idm
    real(8) :: in_t(k,l,m), in_p(k,l,m)
    real(8) :: out_val(k,l,m)
    character (len = *) :: in_name

!f2py real(8), intent(out), dimension(k,l,m) :: out_val
!f2py real(8), intent(in), dimension(k,l,m) :: in_t, in_p
!f2py character(len = *), intent(in) :: in_name
!f2py integer, intent(in) :: k,l,m

    do idk=1,k
    do idl=1,l
    do idm=1,m
            select case(in_name)
                case( 'massfraction_air' )
                    out_val(idk,idl,idm)=liq_ice_air_massfraction_air_si(&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case default
                    write(*,*) 'Unknown function'
            end select
    enddo
    enddo
    enddo
 
end subroutine liq_ice_air_sat

subroutine liq_ice_air_g_ref(in_name,in_q,in_t,in_p, in_ref, k, l, m, out_val)
    !This subroutine is a wrapper for the Gibbs function formulation of humid air thermodynamics.
    !All of these functions require a reference value
    !Available routines:
    !
    ! Liq_Ice_Air_5
    ! pottemp
    ! pottempequi
    ! pottempequipseudo

    use liq_ice_air_5d
    implicit none
    integer :: k, l, m, idk, idl, idm
    real(8) :: in_q(k,l,m), in_t(k,l,m), in_p(k,l,m), in_ref(k,l,m)
    real(8) :: out_val(k,l,m)
    character (len = *) :: in_name

!f2py real(8), intent(out), dimension(k,l,m) :: out_val
!f2py real(8), intent(in), dimension(k,l,m) :: in_q, in_t, in_p, in_ref
!f2py integer, intent(in) :: k,l,m
!f2py character(len = *), intent(in) :: in_name

    do idk=1,k
    do idl=1,l
    do idm=1,m
            select case(in_name)
                case( 'pottemp' )
                    out_val(idk,idl,idm)=liq_ice_air_pottemp_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm), in_ref(idk,idl,idm))
                case( 'pottempequisat' )
                    out_val(idk,idl,idm)=liq_ice_air_pottempequisat_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm), in_ref(idk,idl,idm))
                case( 'pottempequiapprox' )
                    out_val(idk,idl,idm)=liq_ice_air_pottempequiapprox_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm), in_ref(idk,idl,idm))
                case( 'pottempequi' )
                    out_val(idk,idl,idm)=liq_ice_air_pottempequi_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm), in_ref(idk,idl,idm))
                case( 'pottempequipseudo' )
                    out_val(idk,idl,idm)=liq_ice_air_pottempequipseudo_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm), in_ref(idk,idl,idm))
                case default
                    write(*,*) 'Unknown function'
            end select
    enddo
    enddo
    enddo
 
end subroutine liq_ice_air_g_ref

subroutine liq_ice_air_g(in_name,in_q,in_t,in_p, k, l, m, out_val)
    !This subroutine is a wrapper for the Gibbs function formulation of humid air thermodynamics.
    !All of these functions require a reference value
    !Available routines:
    !
    ! Liq_Ice_Air_5
    ! entropy

    use liq_ice_air_5b
    use liq_ice_air_5c
    use liq_ice_air_5d
    !use liq_ice_air_5e
    implicit none
    integer :: k, l, m, idk, idl, idm
    real(8) :: in_q(k,l,m), in_t(k,l,m), in_p(k,l,m)
    real(8) :: out_val(k,l,m)
    character (len = *) :: in_name

!f2py real(8), intent(out), dimension(k,l,m) :: out_val
!f2py real(8), intent(in), dimension(k,l,m) :: in_q, in_t, in_p
!f2py integer, intent(in) :: k,l,m
!f2py character(len = *), intent(in) :: in_name

    do idk=1,k
    do idl=1,l
    do idm=1,m
            select case(in_name)
                case( 'rh_wmo' )
                    out_val(idk,idl,idm)=liq_ice_air_g_rh_wmo_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case( 'entropy' )
                    out_val(idk,idl,idm)=liq_ice_air_g_entropy_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case( 'density' )
                    out_val(idk,idl,idm)=liq_ice_air_g_density_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case( 'enthalpy' )
                    out_val(idk,idl,idm)=liq_ice_air_g_enthalpy_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case( 'cond_entropy' )
                    out_val(idk,idl,idm)=liq_ice_air_g_cond_entropy_si(in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case( 'temperatureequi' )
                    out_val(idk,idl,idm)=liq_ice_air_g_temperatureequi_si(&
                    in_q(idk,idl,idm),&
                    in_t(idk,idl,idm), in_p(idk,idl,idm))
                !case( 'temperaturewetbulb' )
                !    out_val(idk,idl,idm)=liq_ice_air_g_temperaturewetbulb_si(&
                !    in_q(idk,idl,idm),&
                !    in_t(idk,idl,idm), in_p(idk,idl,idm))
                case default
                    write(*,*) 'Unknown function'
            end select
    enddo
    enddo
    enddo
 
end subroutine liq_ice_air_g

subroutine liq_ice_air_h(in_name,in_q,in_s,in_p, k, l, m, out_val)
    !This subroutine is a wrapper for the Gibbs function formulation of humid air thermodynamics.
    !All of these functions require a reference value
    !Available routines:
    !
    ! Liq_Ice_Air_5
    ! temperature

    use liq_ice_air_5a
    implicit none
    integer :: k, l, m, idk, idl, idm
    real(8) :: in_q(k,l,m), in_s(k,l,m), in_p(k,l,m)
    real(8) :: out_val(k,l,m)
    character (len = *) :: in_name

!f2py real(8), intent(out), dimension(k,l,m) :: out_val
!f2py real(8), intent(in), dimension(k,l,m) :: in_q, in_s, in_p
!f2py integer, intent(in) :: k,l,m
!f2py character(len = *), intent(in) :: in_name

    do idk=1,k
    do idl=1,l
    do idm=1,m
            select case(in_name)
                case( 'temperature' )
                    out_val(idk,idl,idm)=liq_ice_air_h_temperature_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'density' )
                    out_val(idk,idl,idm)=liq_ice_air_h_density_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'enthalpy' )
                    out_val(idk,idl,idm)=liq_ice_air_h_enthalpy_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'gc_gv' )
                    out_val(idk,idl,idm)=liq_ice_air_h_gc_gv_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'gd_gc' )
                    out_val(idk,idl,idm)=liq_ice_air_h_gd_gc_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'gd' )
                    out_val(idk,idl,idm)=liq_ice_air_h_gd_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'gc' )
                    out_val(idk,idl,idm)=liq_ice_air_h_gc_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case( 'gh_gc' )
                    out_val(idk,idl,idm)=liq_ice_air_h_gh_gc_si(in_q(idk,idl,idm),&
                    in_s(idk,idl,idm), in_p(idk,idl,idm))
                case default
                    write(*,*) 'Unknown function'
            end select
    enddo
    enddo
    enddo
 
end subroutine liq_ice_air_h

subroutine liq_ice_air_il(in_name,in_q,in_s, k, l, m, out_val)
    !This subroutine is a wrapper for the Gibbs function formulation of humid air thermodynamics.
    !All of these functions require a reference value
    !Available routines:
    !
    ! Liq_Ice_Air_4
    ! iml
    ! ifl

    use liq_ice_air_4
    implicit none
    integer :: k, l, m, idk, idl, idm
    real(8) :: in_q(k,l,m), in_s(k,l,m)
    real(8) :: out_val(k,l,m)
    character (len = *) :: in_name

!f2py real(8), intent(out), dimension(k,l,m) :: out_val
!f2py real(8), intent(in), dimension(k,l,m) :: in_q, in_s
!f2py integer, intent(in) :: k,l,m
!f2py character(len = *), intent(in) :: in_name

    do idk=1,k
    do idl=1,l
    do idm=1,m
            select case(in_name)
                case( 'iml' )
                    out_val(idk,idl,idm)=liq_ice_air_iml_si(in_q(idk,idl,idm),&
                                                                      in_s(idk,idl,idm))
                case( 'ifl' )
                    out_val(idk,idl,idm)=liq_ice_air_ifl_si(in_q(idk,idl,idm),&
                                                                      in_s(idk,idl,idm))
                case default
                    write(*,*) 'Unknown function'
            end select
    enddo
    enddo
    enddo
 
end subroutine liq_ice_air_il
