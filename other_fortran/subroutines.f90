module subroutines
  use mt95
  implicit none
contains

  subroutine integrate(unitary_vel, pos, N, L, N_cells, list_nearest_nbr_cells, eta, v0, int_L)
    integer :: N, N_cells, n_turns, int_L
    real(8) :: unitary_vel(N, 2), pos(N, 2)
    real(8) :: rand_val(N), theta(N), nbr_direction(N)
    real(8) :: noise_cnst, L, eta, v0
    integer :: list_nearest_nbr_cells(N_cells,9)
    noise_cnst = eta*6.283185307d0
    n_turns = 0

    call neighbours_direction(pos, unitary_vel, nbr_direction, N, N_cells, list_nearest_nbr_cells, int_L)

    call genrand_real2(rand_val)
    theta = nbr_direction+noise_cnst*(rand_val-0.5)
    unitary_vel(:,1) = cos(theta)
    unitary_vel(:,2) = sin(theta)
    pos(:,1) = modulo(pos(:,1) + v0*unitary_vel(:,1),L)
    pos(:,2) = modulo(pos(:,2) + v0*unitary_vel(:,2),L)
  end subroutine



  subroutine neighbours_direction(pos, unitary_vel, nbr_direction, N, N_cells, list_nearest_nbr_cells, int_L)
    integer :: i, j, k, where_is, N, cell, N_cells, int_L
    real(8) :: nbr_vel(N,2), nbr_direction(N)
    integer :: cell_list(N), header(N_cells), list_nearest_nbr_cells(N_cells,9)
    real(8) :: unitary_vel(N, 2), pos(N, 2)

    header = 0
    do i=1,N
      cell = int(pos(i,1)) + int_L*int(pos(i,2)) + 1
      cell_list(i) = header(cell)
      header(cell) = i
    end do


    do i = 1, N
      nbr_vel(i,1) = unitary_vel(i,1)
      nbr_vel(i,2) = unitary_vel(i,2)
    end do


    do where_is=1,N_cells
      i = header(where_is)
      do while (i>0)
        j = cell_list(i)
        do while (j>0)
            nbr_vel(i,:) = nbr_vel(i,:) + unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + unitary_vel(i,:)
            j = cell_list(j)
        end do
        i = cell_list(i)
      end do

    end do

    do i=1,N
      where_is = int(pos(i,1)) + int_L*int(pos(i,2)) + 1
      do k=1,4
        j = header(list_nearest_nbr_cells(where_is,k))
        do while (j>0)
            nbr_vel(i,:) = nbr_vel(i,:) + unitary_vel(j,:)
            nbr_vel(j,:) = nbr_vel(j,:) + unitary_vel(i,:)
          j = cell_list(j)
        end do
      end do
    end do


    nbr_direction=atan2(nbr_vel(:,2),nbr_vel(:,1))
  end subroutine




subroutine set_geometry(list_nearest_nbr_cells, int_L)
    integer :: cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y, nbr_cell, int_L
    integer :: nbr_X(9), nbr_Y(9), k, list_nearest_nbr_cells(int_L*int_L,9)
    nbr_X=(/-1,0,1,-1,0,1,-1, 0, 1/)
    nbr_Y=(/ 1,1,1, 0,0,0,-1,-1,-1/)

    do cell_X=0,int_L-1
      do cell_Y=0,int_L-1
        cell = cell_X + int_L*cell_Y + 1
        do k=1,9
          nbr_cell_X=cell_X+nbr_X(k)
          nbr_cell_Y=cell_Y+nbr_Y(k)

          if (nbr_cell_X < 0) then
            nbr_cell = nbr_cell_X + int_L
          else if (nbr_cell_X >= int_L) then
            nbr_cell = nbr_cell_X - int_L
          else
            nbr_cell = nbr_cell_X
          end if

          if (nbr_cell_Y < 0) then
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y + int_L) + 1
          else if (nbr_cell_Y >= int_L) then
            nbr_cell = nbr_cell + int_L*(nbr_cell_Y - int_L) + 1
          else
            nbr_cell = nbr_cell + int_L*nbr_cell_Y + 1
          end if

          list_nearest_nbr_cells(cell,k) = nbr_cell
        end do
      end do
    end do
  end subroutine

  subroutine reset_system(seed, pos, unitary_vel, N, L)
    integer :: seed, N, i
    real(8) :: rand_val(3), L
    real(8) :: theta_i, pos(N, 2), unitary_vel(N, 2)
    call genrand_init( put=seed )

    do i=1,N
        call genrand_real2(rand_val)
        pos(i,:) = rand_val(1:2)*L
        theta_i = rand_val(3)*6.283185307d0
        unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
    end do
  end subroutine

end module
