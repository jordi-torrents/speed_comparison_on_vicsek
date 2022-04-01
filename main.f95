program levy_program
    use subroutines

    real(8) :: polar, polar_i, polar2

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
    inv_L = 1.d0/L
    inv_2L = 2.d0/L


    call set_geometry()
    call reset_system()

    eta = 0.d0

    do iteration = 1, N_reset
        call integrate()
    end do

    do int_eta = 0, 100, 5
        eta = dble(int_eta)/100.d0

        do iteration = 1, N_reset
            call integrate()
        end do

        polar = 0.0d0
        polar2 = 0.0d0

        do iteration = 1, N_steps
            call integrate()
            polar_i = sqrt(sum(unitary_vel(:,1))**2+sum(unitary_vel(:,2))**2)

            polar = polar + polar_i
            polar2 = polar2 + polar_i*polar_i
        end do

        ! polar = polar/dble(N_steps*N)
        ! polar2 = polar2/dble(N_steps*N*N)

        ! print*, eta,',', polar/(dble(N)*dble(N_steps)), ',',(polar2 - (polar * polar)/dble(N_steps))/polar
        print*, &
            eta,',', &
            polar/(dble(N)*dble(N_steps)), ',', &
            sqrt(polar2/dble(N_steps) - (polar * polar)/dble(N_steps*N_steps))/dble(N), ',', &
            (polar2 - (polar * polar)/dble(N_steps))/polar

    end do

  end program levy_program
