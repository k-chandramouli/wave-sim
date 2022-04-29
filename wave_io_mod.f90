!=======================================================================
!     This file contains the wave_io_mod module which defines 
!     I/O unit numbers with mnemonics to avoid confusion.
!=======================================================================

module wave_io_mod

    use wave_kind_mod
    implicit none

    integer(iknd), parameter :: rmesh_unit=10
    integer(iknd), parameter :: thetamesh_unit=11
    integer(iknd), parameter :: eigenval_unit=12
    integer(iknd), parameter :: eigenvec_unit=13

end module wave_io_mod
