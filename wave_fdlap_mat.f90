!=======================================================================
!   This subroutine populates the finite difference Laplacian matrix
!   using mesh details. 
!=======================================================================

subroutine wave_fdlap_mat

    use wave_kind_mod
    use wave_inp_mod
    use wave_data_mod
    implicit none

    integer(iknd) :: k,i
    real(rknd) :: a,b,c


    a = (1._rknd/(2*dr**2))
    b = (1._rknd/dr**2)
    c = (1._rknd/(dr*dtheta)**2)
    lap_mat(1,1) = n_theta*b
    do k=2,n_theta+1
        lap_mat(1,k) = b
        lap_mat(k,k) = 2*(b + c)
        if(mod(k,n_theta) == 2) then ! angle index = 1
            lap_mat(k,k+1) = -c
            lap_mat(k,k+n_theta-1) = -c
        else if(mod(k,n_theta) == 1) then ! angle index = n_theta
            lap_mat(k,k+1-n_theta) = -c
            lap_mat(k,k-1) = -c
        else
            lap_mat(k,k+1) = -c
            lap_mat(k,k-1) = -c
        end if
        lap_mat(k,k+n_theta) = -1*(b + a)
        lap_mat(k,1) = -1*(b - a)
    end do

    do k=n_theta+2,n_nodes
        i = (k-2)/n_theta + 1 ! retrive radial index
        a = (1._rknd/(2*i*dr**2))
        b = (1._rknd/dr**2)
        c = (1._rknd/(i*dr*dtheta)**2)

        lap_mat(k,k) = 2*(b + c)
        if(mod(k,n_theta) == 2) then ! angle index = 1
            lap_mat(k,k+1) = -c
            lap_mat(k,k+n_theta-1) = -c
        else if(mod(k,n_theta) == 1) then ! angle index = n_theta
            lap_mat(k,k+1-n_theta) = -c
            lap_mat(k,k-1) = -c
        else
            lap_mat(k,k+1) = -c
            lap_mat(k,k-1) = -c
        end if
        lap_mat(k,k+n_theta) = -1*(b + a)
        lap_mat(k,k-n_theta) = -1*(b - a)
    end do

    write(*,*) "Laplacian matrix created"

end subroutine wave_fdlap_mat