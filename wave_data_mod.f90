!=======================================================================
!   This module holds the independent and dependent variable arrays,
!   Laplacian matrix, and other variables/constants related to the
!   computation.
!=======================================================================

module wave_data_mod

    use wave_kind_mod
    implicit none

    !------------------------------------------------------------------
    !   General information.
    !------------------------------------------------------------------

    real(rknd) :: pi,twopi

    !------------------------------------------------------------------
    !   Radial and angular step values.
    !------------------------------------------------------------------

    real(rknd) :: dr,dtheta

    !------------------------------------------------------------------
    !   Displacement array, mode frequency array and 
    !   Laplacian matrix declared as allocatables.
    !------------------------------------------------------------------

    integer(iknd) :: n_nodes
    real(rknd), dimension(:), allocatable :: init_disp ! Initial displacement array

    real(rknd), dimension(:), allocatable :: w_r ! real part of eigenvalues
    real(rknd), dimension(:), allocatable :: w_im ! imaginary part of eigenvalues
    real(rknd), dimension(:,:), allocatable :: v_left ! left eigenvectors
    real(rknd), dimension(:,:), allocatable :: v_right ! right eigenvectors
    real(rknd), dimension(:,:), allocatable :: lap_mat ! Laplacian matrix



    contains

    subroutine wave_data_init

        use wave_inp_mod
        implicit none

        integer(iknd) :: i_rad,j_ang

        pi=acos(-1._rknd)
        twopi=2._rknd*pi

        dr = r_max/(n_rad+1)
        dtheta = twopi/n_theta

        n_nodes = n_rad * n_theta + 1
        allocate(w_r(n_nodes))
        allocate(w_im(n_nodes))
        allocate(v_right(n_nodes,n_nodes))
        allocate(v_left(n_nodes,n_nodes))
        allocate(lap_mat(n_nodes,n_nodes))
        lap_mat = 0

        !--------------------------------------------------------------
        !   Allocate initial displacement array
        !--------------------------------------------------------------
        
        allocate(init_disp(n_nodes))
        

    end subroutine wave_data_init

end module wave_data_mod