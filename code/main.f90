module constant
    implicit none
    real, parameter :: c = 0.2      ! 
    real, parameter :: r = 0.2      ! dragging coefficient
    real, parameter :: V = 0.2      ! velocity of liquid droplet
    real, parameter :: nu = 0.2     ! 
    real, parameter :: eta = 0.2    !
    real, parameter :: lambda = 0.2
    real, parameter :: beta = 0.2   ! adiabatic expansion coefficient
    real, parameter :: alpha = 0.2  ! phase transition coefficient
    real, parameter :: Q = 0.2
    real, parameter :: A = 0.1
    real, parameter :: q = 0.1
    real, parameter :: const = 0.0
    real, parameter :: dt = 1     ! time step

end module constant

program main
    implicit none
    integer :: nx, ny
    real :: v_0(nx, ny), v_1(nx, ny), v_2(nx, ny), v_3(nx, ny)
    real :: u_0(nx, ny), u_2(nx, ny), u_3(nx, ny)
    real :: E_0(nx, ny), E_1(nx, ny), E_2(nx, ny), E_3(nx, ny)
    real :: w_l_0(nx, ny), w_l_2(nx, ny), w_l_3(nx, ny)
    real :: w_v_0(nx, ny), w_v_1(nx, ny), w_v_2(nx, ny), w_v_3(nx, ny)
    integer :: step, nsteps

    nx = 10; ny = 10
    nsteps = 200
    
    ! Read input
    call read_input(v_0, u_0, E_0, w_l_0, w_v_0, nx, ny)

    ! Update variables
    do step = 1, nsteps
        call update_variables(v_0, v_1, v_2, v_3, u_0, u_2, u_3, E_0, E_1, E_2, E_3, w_l_0, w_l_2, w_l_3, w_v_0, w_v_1, w_v_2, w_v_3, nx, ny)
        v_0 = v_3; u_0 = u_3; E_0 = E_3; w_l_0 = w_l_3; w_v_0 = w_v_3
    end do

    ! Write output
    call write_output(v_0, u_0, E_0, w_l_0, w_v_0, nx, ny)
    
contains
    !!! Process boundary conditions !!! 
    !!! Read input
    subroutine read_input()
        implicit none
        print *, 'Read input'
    end subroutine read_input

    !!! Write output
    subroutine write_output()
        implicit none
        print *, 'Write output'
    end subroutine write_output

    !!! Update variables
    subroutine update_variables(v_0, v_1, v_2, v_3, u_0, u_2, u_3, E_0, E_1, E_2, E_3, w_l_0, w_l_2, w_l_3, w_v_0, w_v_1, w_v_2, w_v_3, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(inout) :: v_0(nx, ny), v_1(nx, ny), v_2(nx, ny), v_3(nx, ny)
        real, intent(inout) :: u_0(nx, ny), u_2(nx, ny), u_3(nx, ny)
        real, intent(inout) :: E_0(nx, ny), E_1(nx, ny), E_2(nx, ny), E_3(nx, ny)
        real, intent(inout) :: w_l_0(nx, ny), w_l_2(nx, ny), w_l_3(nx, ny)
        real, intent(inout) :: w_v_0(nx, ny), w_v_1(nx, ny), w_v_2(nx, ny), w_v_3(nx, ny)

        call buoyancy_dragging(v_0, v_1, E_0, w_l_0, nx, ny)
        call viscosity_pressure(u_0, v_1, u_2, v_2, nx, ny)
        call diffusion(E_0, E_1, w_v_0, w_v_1, v_0, nx, ny)
        call phase_transition(w_v_1, w_v_2, w_l_0, w_l_2, E_1, E_2, nx, ny)
        call lagrangian_procedure(u_2, v_2, u_3, v_3, E_2, E_3, w_l_2, w_l_3, w_v_2, w_v_3, nx, ny)

    end subroutine update_variables
    
    !!! Laplacian
    subroutine laplacian(A, lap_A, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: A(nx, ny)
        real, intent(out) :: lap_A(nx, ny)
        integer :: i, j

        do i = 2, nx-1
            do j = 2, ny-1
                lap_A(i, j) = 0.25 * (A(i-1, j) + A(i+1, j) + A(i, j-1) + A(i, j+1) - 4.0 * A(i, j))
            end do
        end do

        lap_A(1, :) = 0.0; lap_A(nx, :) = 0.0
        lap_A(:, 1) = 0.0; lap_A(:, ny) = 0.0

    end subroutine laplacian

    !!! Divergence
    subroutine divergence(u, v, div_v, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), v(nx, ny)
        real, intent(out) :: div_v(nx, ny)
        integer :: i, j

        do i = 2, nx-1
            do j = 2, ny-1
                div_v(i, j) = 0.5 * ((u(i+1, j) - u(i-1, j)) + (v(i, j+1) - v(i, j-1)))
            end do
        end do

        div_v(1, :) = 0.0; div_v(nx, :) = 0.0
        div_v(:, 1) = 0.0; div_v(:, ny) = 0.0

    end subroutine divergence
    
    !!! Gradient
    subroutine gradient(div_v, grad_div_u, grad_div_v, nx, ny)
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: div_v(nx, ny)
        real, intent(out) :: grad_div_u(nx, ny), grad_div_v(nx, ny)
        integer :: i, j

        do i = 2, nx-1
            do j = 2, ny-1
                grad_div_u(i, j) = 0.5 * (div_v(i+1, j) - div_v(i-1, j))
                grad_div_v(i, j) = 0.5 * (div_v(i, j+1) - div_v(i, j-1))
            end do
        end do

        grad_div_u(1, :) = 0.0; grad_div_u(nx, :) = 0.0
        grad_div_u(:, 1) = 0.0; grad_div_u(:, ny) = 0.0
        grad_div_v(1, :) = 0.0; grad_div_v(nx, :) = 0.0
        grad_div_v(:, 1) = 0.0; grad_div_v(:, ny) = 0.0

    end subroutine gradient
    
    !!! Representation of Buoyancy and Dragging forces
    subroutine buoyauncy_dragging(v, v_new, E, w_l, nx, ny)
        use constant, only: c, r, V
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: v(nx, ny), E(nx, ny), w_l(nx, ny)
        real, intent(out) :: v_new(nx, ny)
        integer :: i, j
        
        !update v
        do i = 1, nx
            do j = 1, ny
                v_new(i, j) = v(i, j) + c * (E(i + 1, j) + E(i - 1, j) - 2 * E(i, j)) / 2 + r * w_l(i, j) * (v(i, j)  - V)
            end do
        end do

    end subroutine buoyancy_dragging
    
    !!! Representation of Viscosity and Pressure effects
    subroutine viscosity_pressure(u, v, u_new, v_new, nx, ny)
        use constant, only: nu, eta
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), v(nx, ny)
        real, intent(out) :: u_new(nx, ny), v_new(nx, ny)
        real :: lap_u(nx, ny), lap_v(nx, ny)
        real :: div_u(nx, ny), div_v(nx, ny)
        real :: grad_div_u(nx, ny), grad_div_v(nx, ny)
        integer :: i, j

        call laplacian(u, lap_u, nx, ny)
        call laplacian(v, lap_v, nx, ny)
        call divergence(u, v, div_u, nx, ny)
        call gradient(div_u, grad_div_u, grad_div_v, nx, ny)

        do i = 1, nx
            do j = 1, ny
                u_new(i, j) = u(i, j) + nu * lap_u(i, j) + eta * grad_div_u(i, j)
                v_new(i, j) = v(i, j) + nu * lap_v(i, j) + eta * grad_div_v(i, j)
            end do
        end do

    end subroutine viscosity_pressure 
    
    !!! Represntation of Thermal and vapor diffusion
    subroutine diffusion(E, E_new, w_v, w_v_new, v, nx, ny)
        use constant, only: lambda, beta
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: E(nx, ny), v(nx, ny), w_v(nx, ny)
        real, intent(out) :: E_new(nx, ny), w_v_new(nx, ny)
        real :: lap_E(nx, ny), lap_w_v(nx, ny)
        integer :: i, j

        call laplacian(E, lap_E, nx, ny)
        call laplacian(w_v, lap_w_v, nx, ny)
        
        do i = 1, nx
            do j = 1, ny
                E_new(i, j) = E(i, j) + lambda * lap_E(i, j) - beta * u(i, j)
                w_v_new(i, j) = w_v(i, j) + lambda * lap_w_v(i, j)
            end do
        end do

    end subroutine diffusion
    
    !!! Representation of phase transition
    subroutine phase_transition(w_v, w_v_new, w_l, w_l_new, E, E_new, nx, ny)
        use constant, only: alpha, Q, A, q, const, dt
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: E(nx, ny), w_l(nx, ny), w_v(nx, ny)
        real, intent(out) :: E_new(nx, ny), w_l_new(nx, ny), w_v_new(nx, ny)
        real :: W, w_eq(nx, ny)
        integer :: i, j

        do i = 1, nx
            do j = 1, ny
                W = w_l(i, j) + w_v(i, j) 
                if (A * exp(q / (E(i, j) + const)) > W) then
                    w_eq(i, j) = A * exp(q / (E(i, j) + const))
                else
                    w_eq(i, j) = W
                end if
                w_v_new(i, j) = w_v(i, j) + dt * (alpha * (w_v(i, j) - w_eq(i, j)))
                w_l_new(i, j) = w_l(i, j) - dt * (alpha * (w_v(i, j) - w_eq(i, j)))
                E_new(i, j)   = E(i, j) - dt * (Q * (2 * alpha * (w_v(i, j) - w_eq(i, j))))
            end do
        end do

    end subroutine phase_transition

    ! Allocate field variables to neighbors after update position(velocity, internal energy, liquid water, water vapor)
    subroutine lagrangian_procedure(u, v, u_new, v_new, E, E_new, w_l, w_l_new, w_v, w_v_new, nx, ny)
        use constant, only: dt, V
        implicit none
        integer, intent(in) :: nx, ny
        real, intent(in) :: u(nx, ny), v(nx, ny), E(nx, ny), w_l(nx, ny), w_v(nx, ny)
        real, intent(out) :: u_new(nx, ny), v_new(nx, ny), E_new(nx, ny), w_l_new(nx, ny), w_v_new(nx, ny)
        integer :: i, j, ni, nj
        real :: x, y, dx, dy
        real :: w11, w12, w21, w22

        ! initialize field variables
        do i = 1, nx
            do j = 1, ny
                u_new(i, j) = 0.0
                v_new(i, j) = 0.0
                E_new(i, j) = 0.0
                w_l_new(i, j) = 0.0
                w_v_new(i, j) = 0.0
            end do
        end do

        ! update field variables
        do i = 1, nx
            do j = 1, ny
                ! updated position
                x = i + u(i, j) * dt
                y = j + v(i, j) * dt
                y_l = j + (w_l(i, j) - V) * dt
                ! nearest grid point
                ni = floor(x)
                nj = floor(y)
                nj_l = floor(y_l)
                
                ! boundary condition
                if (ni < 1 .or. ni >= nx .or. nj < 1 .or. nj >= ny .or. nj_l < 1 .or. nj_l >= ny) then
                    ! do something
                end if
                
                ! relative position
                dx = x - ni
                dy = y - nj
                dy_l = y_l - nj_l
                
                ! weight
                w11 = (1 - dx) * (1 - dy)
                w12 = (1 - dx) * dy
                w21 = dx * (1 - dy)
                w22 = dx * dy
                w11_l = (1 - dx) * (1 - dy_l)
                w12_l = (1 - dx) * dy_l
                w21_l = dx * (1 - dy_l)
                w22_l = dx * dy_l
                
                ! update field variables
                ! velocity
                u_new(ni, nj) += w11 * u(i, j)
                u_new(ni+1, nj) += w21 * u(i, j)
                u_new(ni, nj+1) += w12 * u(i, j)
                u_new(ni+1, nj+1) += w22 * u(i, j) 
                v_new(ni, nj) += w11 * v(i, j)
                v_new(ni+1, nj) += w21 * v(i, j)
                v_new(ni, nj+1) += w12 * v(i, j)
                v_new(ni+1, nj+1) += w22 * v(i, j)
                ! internal energy
                E_new(ni, nj) += w11 * E(i, j)
                E_new(ni+1, nj) += w21 * E(i, j)
                E_new(ni, nj+1) += w12 * E(i, j)
                E_new(ni+1, nj+1) += w22 * E(i, j)

                ! liquid water
                w_l_new(ni, nj_l) += w11_l * w_l(i, j)
                w_l_new(ni+1, nj_l) += w21_l * w_l(i, j)
                w_l_new(ni, nj_l+1) += w12_l * w_l(i, j)
                w_l_new(ni+1, nj_l+1) += w22_l * w_l(i, j)
                ! water vapor
                w_v_new(ni, nj) += w11 * w_v(i, j)
                w_v_new(ni+1, nj) += w21 * w_v(i, j)
                w_v_new(ni, nj+1) += w12 * w_v(i, j)
                w_v_new(ni+1, nj+1) += w22 * w_v(i, j)
            end do
        end do 
        
    end subroutine lagrangian_procedure
end program main
