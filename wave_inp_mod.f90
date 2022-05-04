!=======================================================================
!   This module declares all input parameters and sets their default
!   values. These items are declared in the namelist that is read in
!   the driver program.
!=======================================================================

module wave_inp_mod

    use wave_kind_mod
    implicit none

    !-----------------------------------------------------------------------
    !   Density is set based on the density_type parameter, which
    !   can be either constant or variable.
    !-----------------------------------------------------------------------

    character(16) :: density_type="constant"
    real(rknd) :: rho0=1._rknd

    !-----------------------------------------------------------------------
    !   Set the maximum radius of the membrane (in m) and the initial and
    !   final times (in s).
    !-----------------------------------------------------------------------

    real(rknd) :: r_max=1._rknd
    real(rknd) :: t_initial=0._rknd
    real(rknd) :: t_final=1._rknd

    !-----------------------------------------------------------------------
    !   Set 0 displacement boundary conditions for a clamped membrane.
    !-----------------------------------------------------------------------

    real(rknd) :: disp_bdry=0._rknd

    !-----------------------------------------------------------------------
    !   Radial and angular resolution parameters.
    !   Angular resolution must be in multiples of 4 (to resolve Laplacian
    !   at r = 0). Also set the number of modes for which eigenvalues and
    !   eigenvectors are to be computed.
    !-----------------------------------------------------------------------

    integer(iknd) :: n_rad=20_iknd
    integer(iknd) :: n_theta=16_iknd
    integer(iknd) :: num_modes=10_iknd

    !-----------------------------------------------------------------------
    !   Initial condition parameters.
    !-----------------------------------------------------------------------

    ! How can the initial conditions be set? What's the cleanest way
    ! of doing it? Initial conditions can be in terms of velocity
    ! (partial_u/partial_t = phi(r,theta)) or it can be directly assigning
    ! u(0,r,theta) = f(r,theta), which makes more sense?

end module wave_inp_mod