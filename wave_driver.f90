!=======================================================================
!   This program is the main driver program for the membrane wave
!   equation solver. 
!
!   It finds the eigenvalues and eigenvectors of
!   the 2D wave equation du/dt = c^2 * del^2(u) in polar geometry
!   using a finite difference method to compute a Laplacian matrix.  
!=======================================================================

program wave_driver

    use wave_kind_mod
    use wave_inp_mod
    use wave_io_mod
    use wave_data_mod
    implicit none

    real(rknd) :: time_st,time_en
    character :: JOBVL,JOBVR
    integer(iknd) :: LDA,LDVL,LDVR,LWORK,INFO
    real(rknd), dimension(:), allocatable :: WORK

    integer(iknd) :: i_rad,j_ang


    !-----------------------------------------------------------------------
    !   Namelist declarations for input parameters.
    !-----------------------------------------------------------------------

    namelist /wave_inp/ density_type,rho0,r_max,t_initial,t_final,&
                        n_rad,n_theta,num_modes

    !-----------------------------------------------------------------------
    !   Read namelist from file wave_input.dat.
    !-----------------------------------------------------------------------

    OPEN(UNIT=nml_unit,FILE="wave_inp.dat",STATUS="OLD",&
        FORM="FORMATTED")
    READ(UNIT=nml_unit,NML=wave_inp)
    CLOSE(nml_unit)

    !-----------------------------------------------------------------------
    !   Start the timer, call the data initialisation routine and write
    !   the radial and angular mesh values to file for plotting
    !-----------------------------------------------------------------------

    call timer(time_st)
    call wave_data_init

    OPEN(UNIT=rad_mesh_unit,FILE='rad_mesh.dat',STATUS='UNKNOWN')
    OPEN(UNIT=ang_mesh_unit,FILE='ang_mesh.dat',STATUS='UNKNOWN')
    do i_rad=1,n_rad+2
        WRITE(rad_mesh_unit,'(es12.5)') rad_mesh(i_rad)
    end do
    do j_ang=1,n_theta
        WRITE(ang_mesh_unit,'(es12.5)') ang_mesh(i_rad)
    end do

    !-----------------------------------------------------------------------
    !   Open data files to write eigenvalues and eigenvectors to.
    !-----------------------------------------------------------------------

    open(unit=eigval_unit,file='eigenvalues.dat',status='unknown')
    open(unit=r_eigvec_unit,file='right_eigenvectors.dat',status='unknown')
    open(unit=l_eigvec_unit,file='left_eigenvectors.dat',status='unknown')

    !-----------------------------------------------------------------------
    !   Call subroutine to populate FD Laplacian matrix
    !-----------------------------------------------------------------------

    call wave_fdlap_mat

    !-----------------------------------------------------------------------
    !   Set up and call LAPACK routine to find
    !   eigenvalues and eigenvectors
    !-----------------------------------------------------------------------

    JOBVL='V'
    JOBVR='V'
    
    LWORK = 5*n_nodes
    allocate(WORK(LWORK))
    LDA = n_nodes
    LDVL = n_nodes
    LDVR = n_nodes

    call dgeev(JOBVL,JOBVR,n_nodes,lap_mat,LDA,w_r,w_im,v_left,LDVL,&
               v_right,LDVR,WORK,LWORK,INFO)

    call wave_write

    write(*,*) "DGEEV result: ", INFO

    call timer(time_en)
    write(*,*) "Time taken: ", (time_en - time_st), " seconds."

end program

!=======================================================================
!   Subroutine wave_write writes the dependent variable for analysis
!   and plotting
!=======================================================================

subroutine wave_write

    use wave_kind_mod
    use wave_inp_mod
    use wave_data_mod
    use wave_io_mod
    implicit none

    integer(iknd) :: ix,jx

    do ix=1,n_nodes
        write(eigval_unit,'(es14.6)') w_r(ix)
        do jx=1,n_nodes-1
            write(r_eigvec_unit,'(es14.6)',ADVANCE="NO") v_right(jx,ix)
            write(l_eigvec_unit,'(es14.6)',ADVANCE="NO") v_left(jx,ix)
        end do
        write(r_eigvec_unit,'(es14.6)') v_right(n_nodes,ix)
        write(l_eigvec_unit,'(es14.6)') v_left(n_nodes,ix)
    end do

    return
end subroutine wave_write