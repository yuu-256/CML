program velocity_evolution
  implicit none
  integer, parameter :: nx = 10, ny = 10
  real, parameter :: nu = 0.1, eta = 0.1, dt = 0.01
  real :: u(nx, ny), v(nx, ny), u_new(nx, ny), v_new(nx, ny)
  real :: lap_u(nx, ny), lap_v(nx, ny)
  real :: div_v(nx, ny), grad_div_u(nx, ny), grad_div_v(nx, ny)
  integer :: i, j

  ! 初期化（例としてランダムな値）
  call random_seed()
  call random_number(u)
  call random_number(v)

  ! 速度場の更新
  call update_velocity(u, v, u_new, v_new, nx, ny, nu, eta, dt)

  ! 結果出力
  do i = 1, nx
    do j = 1, ny
      write(*,*) "x=", i, "y=", j, " u_new=", u_new(i,j), " v_new=", v_new(i,j)
    end do
  end do

contains

  subroutine update_velocity(u, v, u_new, v_new, nx, ny, nu, eta, dt)
    implicit none
    integer, intent(in) :: nx, ny
    real, intent(in) :: u(nx, ny), v(nx, ny)
    real, intent(out) :: u_new(nx, ny), v_new(nx, ny)
    real, intent(in) :: nu, eta, dt
    real :: lap_u(nx, ny), lap_v(nx, ny)
    real :: div_v(nx, ny), grad_div_u(nx, ny), grad_div_v(nx, ny)
    integer :: i, j

    ! ラプラス作用素 Δu, Δv の計算
    call compute_laplacian(u, lap_u, nx, ny)
    call compute_laplacian(v, lap_v, nx, ny)

    ! 発散 ∇⋅v の計算
    call compute_divergence(u, v, div_v, nx, ny)

    ! 勾配 ∇(∇⋅v) の計算
    call compute_gradient(div_v, grad_div_u, grad_div_v, nx, ny)

    ! 時間発展
    do i = 2, nx-1
      do j = 2, ny-1
        u_new(i, j) = u(i, j) + dt * (nu * lap_u(i, j) + eta * grad_div_u(i, j))
        v_new(i, j) = v(i, j) + dt * (nu * lap_v(i, j) + eta * grad_div_v(i, j))
      end do
    end do

    ! 境界の処理（ゼロにする）
    u_new(1, :) = 0.0;  u_new(nx, :) = 0.0
    u_new(:, 1) = 0.0;  u_new(:, ny) = 0.0
    v_new(1, :) = 0.0;  v_new(nx, :) = 0.0
    v_new(:, 1) = 0.0;  v_new(:, ny) = 0.0

  end subroutine update_velocity

  ! ラプラス作用素 ΔA の計算
  subroutine compute_laplacian(A, lap_A, nx, ny)
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

  end subroutine compute_laplacian

  ! 発散 ∇⋅v の計算
  subroutine compute_divergence(u, v, div_v, nx, ny)
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

  end subroutine compute_divergence

  ! 勾配 ∇(∇⋅v) の計算
  subroutine compute_gradient(div_v, grad_div_u, grad_div_v, nx, ny)
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

  end subroutine compute_gradient

end program velocity_evolution

