!=======================================================================
!   This module contains the F90 kind definitions for integers and
!   reals to be used for all variables and constants
!=======================================================================

module wave_kind_mod

    implicit none

    integer, parameter :: iknd=selected_int_kind(8)
    integer, parameter :: ilongknd=selected_int_kind(18)
    integer, parameter :: rknd=selected_real_kind(14,100)

end module wave_kind_mod