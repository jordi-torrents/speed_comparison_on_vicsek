program levy_program
    use subroutines
    use mt95
    implicit none
    real(8) :: pi2, rho, L, eta, v0, polar
    real(8), allocatable :: pos(:,:), unitary_vel(:,:), theta(:)
    integer :: seed, N, N_reset, N_steps, N_cells, int_L, iteration, int_eta
    integer, allocatable :: list_nearest_nbr_cells(:,:)
    parameter (pi2 = 6.283185307d0)

    open(unit=645,file='input.dat')
    read(645,*) int_L
    read(645,*) v0
    read(645,*) rho
    read(645,*) N_reset
    read(645,*) N_steps
    read(645,*) seed
    close(645)

    L=dble(int_L)
    N=nint(L*L*rho)
    N_cells=int_L*int_L


    allocate(pos(N,2))
    allocate(unitary_vel(N,2))
    allocate(theta(N))
    allocate(list_nearest_nbr_cells(N_cells,9))


    call set_geometry(list_nearest_nbr_cells, int_L)
    call reset_system(seed, pos, unitary_vel, N, L)

    eta = 0.d0
    ! print*, pos(40,1)

    do iteration = 1, N_reset
        call integrate(unitary_vel, pos, N, L, N_cells, list_nearest_nbr_cells, eta, v0, int_L)
    end do

    do int_eta = 0, 100, 10
        eta = dble(int_eta)/100.d0

        do iteration = 1, N_reset
            call integrate(unitary_vel, pos, N, L, N_cells, list_nearest_nbr_cells, eta, v0, int_L)
        end do

        polar = 0.0d0

        do iteration = 1, N_steps
            call integrate(unitary_vel, pos, N, L, N_cells, list_nearest_nbr_cells, eta, v0, int_L)
            polar = polar + sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)/dble(N)
        end do

        print*, eta, polar/dble(N_steps)
    end do

  end program levy_program
