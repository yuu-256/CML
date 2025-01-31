!! Two dimentional Cloud Dynamics Model implementation as Coupled Map Lattice (CML) model
!! paper: Modeling and Characterization of Cloud Dynamics, 1997, Tatsuo Yanagita and Kunihiko Kaneko
program main
    implicit none
    integer :: width, height

    !!! Make lattices 

contains
    subroutine laplacian(A, delta_A, nx, ny)
        ! Declaration
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: A(nx, ny)
        real, intent(out) :: delta_A(nx, ny)
        integer :: i, j
        
        ! Execution
        do i = 2, nx-1
            do j = 2, ny-1
                delta_A(i, j) = 0.25 * (A(i-1, j) + A(i+1,j) + A(i,j-1) + A(i,j+1) - 4 * A(i,j))
            end do
        end do

        ! Boundary condition 
        delta_A(1,:) = 0.0
        delta_A(nx,:) = 0.0
        delta_A(:,1) = 0.0
        delta_A(:,ny) = 0.0
    end subroutine laplacian
    
    !!! 境界処理してください。useもお願いします。
    subroutine buoyancy_dragging(v_y, v_y_new, E, w_l, V, nx, ny)
        ! Declaration
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: v_y(nx, ny), E(nx, ny), w_l(nx, ny), V
        real, intent(out) :: v_y_new(nx, ny)

        ! Execution
        do i = 1, nx
            do j = 1, ny
                v_y_new(i, j) = v_y(i, j) + c * (E(i+1,j) + E(i-1,j) - 2 * E(i,j)) / 2 - r * w_l(i,j) * (v_y(i,j) - V)
            end do
        end do
    end subroutine buoyancy_dragging 
       
    subroutine viscosity_pressure(v, v_new, nx, ny)
        ! Declaration
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: v(nx, ny)
        real, intent(out) :: v_new(nx, ny)

end program main
