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


    !-----------------------------------------------------------------------
    !   Namelist declarations for input parameters.
    !-----------------------------------------------------------------------

    !   FILL HERE

    !-----------------------------------------------------------------------
    !   Start the timer and call the data initialisation routine.
    !-----------------------------------------------------------------------

    call timer(time_st)
    call wave_data_init

    !-----------------------------------------------------------------------
    !   Read namelist from file wave_input.dat.
    !-----------------------------------------------------------------------

    !   FILL HERE

    !-----------------------------------------------------------------------
    !   Call subroutine to populate FD Laplacian matrix
    !-----------------------------------------------------------------------

    call wave_fdlap_mat

    !-----------------------------------------------------------------------
    !   Set up and call LAPACK routine to find
    !   eigenvalues and eigenvectors
    !-----------------------------------------------------------------------

    JOBVL='N'
    JOBVR='V'
    
    LWORK = 5*n_nodes
    allocate(WORK(LWORK))
    LDA = n_nodes
    LDVL = n_nodes
    LDVR = n_nodes

    call dgeev(JOBVL,JOBVR,n_nodes,lap_mat,LDA,w_r,w_im,v_left,LDVL,&
               v_right,LDVR,WORK,LWORK,INFO)

    write(*,*) "DGEEV result: ", INFO

    write(*,*) "WI: ", w_im

    write(*,*) "WR: ", w_r

    call timer(time_en)
    write(*,*) "Time taken: ", (time_en - time_st), " seconds."

end program