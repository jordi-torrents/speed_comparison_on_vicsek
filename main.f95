program levy_program
    use subroutines

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


    call integrate(N_reset, .false.)

    do int_eta = 0, 100, 5
        eta = dble(int_eta)/100.d0

        call integrate(N_reset, .false.)

        polar = 0.0d0
        polar2 = 0.0d0


        call integrate(N_steps, .true.)




        print*, &
            eta,',', &
            polar/(dble(N)*dble(N_steps)), ',', &
            sqrt(polar2/dble(N_steps) - (polar * polar)/dble(N_steps*N_steps))/dble(N), ',', &
            (polar2 - (polar * polar)/dble(N_steps))/polar

    end do

  end program levy_program
