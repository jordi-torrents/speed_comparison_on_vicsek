module subroutines
  use mt19937_64

  implicit none
  real(8) :: rho, L, eta, v0, inv_L, inv_2L
  real(8), allocatable :: pos(:,:), unitary_vel(:,:)
  integer :: N, N_reset, N_steps, N_cells, int_L, i
  integer, allocatable :: list_nearest_nbr_cells(:,:)
  integer(8)::seed

contains

  subroutine integrate()
    real(8), dimension(N) :: nbr_direction(N)

    call neighbours_direction(nbr_direction)


    do i=1, N
        nbr_direction(i) = nbr_direction(i)+eta*6.283185307d0*(genrand64_real2()-0.5d0)
    end do
    ! nbr_direction = nbr_direction+eta*6.283185307d0*(rand_val-0.5)
    unitary_vel(:,1) = cos(nbr_direction)
    unitary_vel(:,2) = sin(nbr_direction)
    pos = modulo(pos + v0*unitary_vel,L)
  end subroutine



!   subroutine neighbours_direction(nbr_direction)
!     real(8) :: nbr_direction(N), integrated_vel(N,2)
!     integer :: cell, header(N_cells), cell_list(N), nbr_cell, j, nbr_cell_index, where_is(N)

!     header = 0

!     do i = 1,N
!         cell = int(pos(i,1)) + int_L*int(pos(i,2)) + 1
!         cell_list(i) = header(cell)
!         header(cell) = i
!         where_is(i) = cell;
!         integrated_vel(i,:) = unitary_vel(i,:)
!     end do

!     do i=1,N
!         cell = where_is(i)

!         do nbr_cell_index=1,9
!             nbr_cell = list_nearest_nbr_cells(cell, nbr_cell_index)
!             j = header(nbr_cell)
!             do while (j > 0)
!                 if (i<j) then
!                     if (dist_PBC(pos(i,:), pos(j,:))>1.d0) then
!                         integrated_vel(i,:) = integrated_vel(i,:) + unitary_vel(j,:)
!                         integrated_vel(j,:) = integrated_vel(j,:) + unitary_vel(i,:)
!                     end if
!                 end if
!                 j = cell_list(j)
!             end do
!         end do
!         nbr_direction(i) = atan2(integrated_vel(i,2), integrated_vel(i,1))
!     end do

!   end subroutine


  subroutine neighbours_direction(nbr_direction)
    real(8) :: nbr_direction(N), integrated_vel(N,2)
    integer :: cell, header(N_cells), cell_list(N), nbr_cell, j, nbr_cell_index

    header = 0

    do i = 1,N
        cell = int(pos(i,1)) + int_L*int(pos(i,2)) + 1
        cell_list(i) = header(cell)
        header(cell) = i
        integrated_vel(i,:) = unitary_vel(i,:)
    end do

    do cell=1, N_cells

        i = header(cell)
        do while (i>0)
            j = cell_list(i)
            do while (j>0)
                if (dist_simple(pos(i,:), pos(j,:))<1.d0) then
                    integrated_vel(i,:) = integrated_vel(i,:) + unitary_vel(j,:)
                    integrated_vel(j,:) = integrated_vel(j,:) + unitary_vel(i,:)
                end if
                j = cell_list(j)
            end do

            do nbr_cell_index=1,4
                nbr_cell = list_nearest_nbr_cells(cell, nbr_cell_index)
                j = header(nbr_cell)
                do while (j > 0)
                    if (dist_PBC(pos(i,:), pos(j,:))<1.d0) then
                        integrated_vel(i,:) = integrated_vel(i,:) + unitary_vel(j,:)
                        integrated_vel(j,:) = integrated_vel(j,:) + unitary_vel(i,:)
                    end if
                    j = cell_list(j)
                end do
            end do
            i = cell_list(i)
        end do
    end do

    do i=1, N
        nbr_direction(i) = atan2(integrated_vel(i,2), integrated_vel(i,1))
    end do

  end subroutine

subroutine set_geometry()
    integer :: cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y
    integer :: nbr_X(4), nbr_Y(4)
    allocate(pos(N,2))
    allocate(unitary_vel(N,2))
    allocate(list_nearest_nbr_cells(N_cells,4))
    nbr_X=(/-1,0,1,-1/)
    nbr_Y=(/ 1,1,1, 0/)

    do cell_X=0,int_L-1
      do cell_Y=0,int_L-1
        cell = cell_X + int_L*cell_Y + 1
        do i=1,4
          nbr_cell_X=cell_X+nbr_X(i)
          nbr_cell_Y=cell_Y+nbr_Y(i)

          if (nbr_cell_X <  0    ) nbr_cell_X = nbr_cell_X + int_L
          if (nbr_cell_X >= int_L) nbr_cell_X = nbr_cell_X - int_L

          if (nbr_cell_Y <  0    ) nbr_cell_Y = nbr_cell_Y + int_L
          if (nbr_cell_Y >= int_L) nbr_cell_Y = nbr_cell_Y - int_L

          list_nearest_nbr_cells(cell,i) = nbr_cell_X + int_L*nbr_cell_Y +1
        end do
      end do
    end do
  end subroutine

  subroutine reset_system()
    real(8) :: theta_i
    call init_genrand64(seed)

    do i=1,N
        pos(i,:) = (/genrand64_real2(), genrand64_real2()/)*L
        theta_i = genrand64_real2()*6.283185307d0
        unitary_vel(i,:) = (/cos(theta_i),sin(theta_i)/)
    end do
  end subroutine

  real(8) function dist_PBC(pos1,pos2)
    real(8) :: dx, dy, pos1(2), pos2(2)

    dx = pbc_dist(pos1(1) - pos2(1))
    dy = pbc_dist(pos1(2) - pos2(2))

    dist_PBC = dx*dx + dy*dy

  end function dist_PBC

  function dist_simple(pos1,pos2)
    real(8) :: dist_simple, dx, dy, pos1(2), pos2(2)

    dx = pos1(1) - pos2(1)
    dy = pos1(2) - pos2(2)

    dist_simple = dx*dx + dy*dy

  end function dist_simple

  function pbc_dist(x)
   real(8) :: pbc_dist, x
     pbc_dist = x - int(inv_2L*x)*L
  end function pbc_dist

end module
