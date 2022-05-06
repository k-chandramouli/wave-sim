!=======================================================================
!     This file contains the wave_io_mod module which defines 
!     I/O unit numbers with mnemonics to avoid confusion.
!=======================================================================

module wave_io_mod

    use wave_kind_mod
    implicit none

    integer(iknd), parameter :: rad_mesh_unit=10
    integer(iknd), parameter :: ang_mesh_unit=11
    integer(iknd), parameter :: eigval_unit=12
    integer(iknd), parameter :: r_eigvec_unit=13
    integer(iknd), parameter :: l_eigvec_unit=14
    integer(iknd), parameter :: nml_unit=15

end module wave_io_mod
