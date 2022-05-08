!=======================================================================
!   This subroutine populates the finite difference Laplacian matrix
!   using mesh details. 
!=======================================================================

subroutine wave_fdlap_mat2

    use wave_kind_mod
    use wave_inp_mod
    use wave_data_mod
    implicit none

    integer(iknd) :: k,i
    real(rknd) :: a,b,c,rho
    real(rknd), external :: wave_rho

    if (density_type == "constant") then
    
        !-----------------------------------------------------------------------
        !   Fill matrix for points at r = 0 and first circumference
        !-----------------------------------------------------------------------
        a = (1._rknd/(dr**2))*c0**2
        b = (1._rknd/dr**2)*c0**2
        c = (1._rknd/(dr*dtheta)**2)*c0**2
        lap_mat(1,1) = n_theta*b
        do k=2,n_theta/2+1
            lap_mat(1,k) = b
            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,1) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k+n_theta/2) = (1._rknd/12._rknd)*(b-a)
        end do

        do k=n_theta/2+2,n_theta+1
            lap_mat(1,k) = b
            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,1) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k-n_theta/2) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the second circumference
        !-----------------------------------------------------------------------
        
        do k=n_theta+2,n_nodes-2*n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            a = (1._rknd/(i*dr**2))*c0**2
            b = (1._rknd/dr**2)*c0**2
            c = (1._rknd/(dr*dtheta)**2)*c0**2

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,1) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points from third circumference to
        !   second from last circumference
        !-----------------------------------------------------------------------
        do k=n_theta+2,n_nodes-2*n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            a = (1._rknd/(i*dr**2))*c0**2
            b = (1._rknd/dr**2)*c0**2
            c = (1._rknd/(dr*dtheta)**2)*c0**2

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k-2*n_theta) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the second from last circumference
        !-----------------------------------------------------------------------
        do k=n_nodes-2*n_theta+1,n_nodes-n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            a = (1._rknd/(i*dr**2))*c0**2
            b = (1._rknd/dr**2)*c0**2
            c = (1._rknd/(dr*dtheta)**2)*c0**2

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k-2*n_theta) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the last circumference
        !-----------------------------------------------------------------------
        do k=n_nodes-n_theta+1,n_nodes
            i = (k-2)/n_theta + 1 ! retrive radial index
            a = (1._rknd/(i*dr**2))*c0**2
            b = (1._rknd/dr**2)*c0**2
            c = (1._rknd/(dr*dtheta)**2)*c0**2

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
            lap_mat(k,k-n_theta) = (a-b)
        end do

    else if (density_type == "variable") then

        !-----------------------------------------------------------------------
        !   Fill matrix for points at r = 0 and first circumference
        !-----------------------------------------------------------------------
        rho = wave_rho(dr)
        a = (1._rknd/(dr**2))*d0**2/rho
        b = (1._rknd/dr**2)*d0**2/rho
        c = (1._rknd/(dr*dtheta)**2)*d0**2/rho
        lap_mat(1,1) = n_theta*b
        do k=2,n_theta/2+1
            lap_mat(1,k) = b
            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,1) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k+n_theta/2) = (1._rknd/12._rknd)*(b-a)
        end do

        do k=n_theta/2+2,n_theta+1
            lap_mat(1,k) = b
            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,1) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k-n_theta/2) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the second circumference
        !-----------------------------------------------------------------------

        do k=n_theta+2,n_nodes-2*n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            rho = wave_rho(i*dr)
            a = (1._rknd/(i*dr**2))*d0**2/rho
            b = (1._rknd/dr**2)*d0**2/rho
            c = (1._rknd/(dr*dtheta)**2)*d0**2/rho

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,1) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points from third circumference to
        !   second from last circumference
        !-----------------------------------------------------------------------
        do k=n_theta+2,n_nodes-2*n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            rho = wave_rho(i*dr)
            a = (1._rknd/(i*dr**2))*d0**2/rho
            b = (1._rknd/dr**2)*d0**2/rho
            c = (1._rknd/(dr*dtheta)**2)*d0**2/rho

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k+2*n_theta) = (1._rknd/12._rknd)*(b+a)
            lap_mat(k,k-2*n_theta) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the second from last circumference
        !-----------------------------------------------------------------------
        do k=n_nodes-2*n_theta+1,n_nodes-n_theta
            i = (k-2)/n_theta + 1 ! retrive radial index
            rho = wave_rho(i*dr)
            a = (1._rknd/(i*dr**2))*d0**2/rho
            b = (1._rknd/dr**2)*d0**2/rho
            c = (1._rknd/(dr*dtheta)**2)*d0**2/rho

            lap_mat(k,k) = (15._rknd/6._rknd)*b + 2*c
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
            lap_mat(k,k+n_theta) = -1*((4._rknd/3._rknd)*b + (2._rknd/3._rknd)*a)
            lap_mat(k,k-n_theta) = -1*((4._rknd/3._rknd)*b - (2._rknd/3._rknd)*a)
            lap_mat(k,k-2*n_theta) = (1._rknd/12._rknd)*(b-a)
        end do

        !-----------------------------------------------------------------------
        !   Fill matrix for points on the last circumference
        !-----------------------------------------------------------------------
        do k=n_nodes-n_theta+1,n_nodes
            i = (k-2)/n_theta + 1 ! retrive radial index
            rho = wave_rho(i*dr)
            a = (1._rknd/(i*dr**2))*d0**2/rho
            b = (1._rknd/dr**2)*d0**2/rho
            c = (1._rknd/(dr*dtheta)**2)*d0**2/rho

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
            lap_mat(k,k-n_theta) = (a-b)
        end do

    else
        write(*,*) "Density type not recognised"
    end if

    write(*,*) "Laplacian matrix created"

end subroutine wave_fdlap_mat2