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
    !   Radial and angular mesh arrays and step values.
    !------------------------------------------------------------------

    real(rknd) :: dr,dtheta
    real(rknd), dimension(:), allocatable :: rad_mesh
    real(rknd), dimension(:), allocatable :: ang_mesh

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

        dr = r_max/(n_rad+2)
        dtheta = twopi/n_theta

        !--------------------------------------------------------------
        !   Allocate eigenvalue, eigenvectors and Laplacian arrays
        !--------------------------------------------------------------

        n_nodes = n_rad * n_theta + 1
        allocate(w_r(n_nodes))
        allocate(w_im(n_nodes))
        allocate(v_right(n_nodes,n_nodes))
        allocate(v_left(n_nodes,n_nodes))
        allocate(lap_mat(n_nodes,n_nodes))
        lap_mat = 0

        !--------------------------------------------------------------
        !   Allocate and fill mesh arrays
        !--------------------------------------------------------------

        allocate(rad_mesh(n_rad+2))
        allocate(ang_mesh(n_theta))

        do i_rad=0,n_rad+1
            rad_mesh(i_rad+1) = (i_rad)*dr
        end do

        do j_ang=0,n_theta-1
            ang_mesh(j_ang+1) = j_ang*dtheta
        end do

        !--------------------------------------------------------------
        !   Allocate initial displacement array
        !--------------------------------------------------------------
        
        allocate(init_disp(n_nodes))
        

    end subroutine wave_data_init

end module wave_data_mod