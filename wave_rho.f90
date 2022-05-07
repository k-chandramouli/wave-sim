!=======================================================================
!   This subroutine populates the finite difference Laplacian matrix
!   using mesh details. 
!=======================================================================

function wave_rho(radius) result(rho)
    use wave_kind_mod
    use wave_inp_mod
    implicit none

    real(rknd), intent(in) :: radius
    real(rknd) :: rho

    if(radius <= rho_bdry) then
        rho = rho1
        return
    else
        rho = rho0
        return
    end if

end function wave_rho