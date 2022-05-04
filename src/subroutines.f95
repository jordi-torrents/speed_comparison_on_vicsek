module subroutines
  use mt19937_64

  implicit none
  real(8) :: rho, L, eta, v0, inv_L, inv_2L, polar, polar2
  real(8), allocatable :: pos(:,:), unitary_vel(:,:)
  integer :: N, N_reset, N_steps, N_cells, int_L, i
  integer, allocatable :: list_nearest_nbr_cells(:,:)
  integer(8) :: seed

contains


  subroutine integrate(steps, measure)
    real(8) :: nbr_direction(N), polar_i
    real(8) :: integrated_vel(2,N), eff_noise!, randoms(N)
    integer :: cell, header(N_cells), cell_list(N), j, steps, step!, nbr_cell_index
    logical :: measure
    eff_noise = eta*6.283185307d0

    do step=1, steps
        header = 0
        ! integrated_vel = unitary_vel

        do i = 1,N
            integrated_vel(:,i) = unitary_vel(:,i)
            cell = int(pos(1,i)) + int_L*int(pos(2,i)) + 1
            cell_list(i) = header(cell)
            header(cell) = i
        end do

        do cell=1, N_cells

            i = header(cell)
            do while (i>0)
                j = cell_list(i)
                do while (j>0)
                    if (dist_simple(pos(:,i), pos(:,j))<1.d0) then
                        integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                        integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                    end if
                    j = cell_list(j)
                end do

                ! do nbr_cell_index=1,4
                !     j = header(list_nearest_nbr_cells(cell, nbr_cell_index))
                !     do while (j > 0)
                !         if (dist_PBC(pos(:,i), pos(:,j))<1.d0) then
                !             integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                !             integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                !         end if
                !         j = cell_list(j)
                !     end do
                ! end do

                j = header(list_nearest_nbr_cells(1, cell))
                do while (j > 0)
                    if (dist_PBC(pos(:,i), pos(:,j))<1.d0) then
                        integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                        integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                    end if
                    j = cell_list(j)
                end do

                j = header(list_nearest_nbr_cells(2, cell))
                do while (j > 0)
                    if (dist_PBC(pos(:,i), pos(:,j))<1.d0) then
                        integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                        integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                    end if
                    j = cell_list(j)
                end do

                j = header(list_nearest_nbr_cells(3, cell))
                do while (j > 0)
                    if (dist_PBC(pos(:,i), pos(:,j))<1.d0) then
                        integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                        integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                    end if
                    j = cell_list(j)
                end do

                j = header(list_nearest_nbr_cells(4, cell))
                do while (j > 0)
                    if (dist_PBC(pos(:,i), pos(:,j))<1.d0) then
                        integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(:,j)
                        integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(:,i)
                    end if
                    j = cell_list(j)
                end do

                i = cell_list(i)
            end do
        end do


        do i=1, N
            nbr_direction(i) = atan2(integrated_vel(2,i), integrated_vel(1,i))+eff_noise*(genrand64_real2()-0.5)
        end do


        unitary_vel(1,:) = cos(nbr_direction)
        unitary_vel(2,:) = sin(nbr_direction)
        pos = modulo(pos + v0*unitary_vel,L)

        if (measure) then
            polar_i = sqrt(sum(unitary_vel(1,:))**2+sum(unitary_vel(2,:))**2)

            polar = polar + polar_i
            polar2 = polar2 + polar_i*polar_i
        end if
    end do

  end subroutine



subroutine set_geometry()
    integer :: cell, cell_X, cell_Y, nbr_cell_X, nbr_cell_Y
    integer :: nbr_X(4), nbr_Y(4)
    allocate(pos(2,N))
    allocate(unitary_vel(2,N))
    allocate(list_nearest_nbr_cells(4,N_cells))
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

          list_nearest_nbr_cells(i, cell) = nbr_cell_X + int_L*nbr_cell_Y +1
        end do
      end do
    end do
  end subroutine

  subroutine reset_system()
    real(8) :: theta_i
    call init_genrand64(seed)

    do i=1,N
        pos(:,i) = (/genrand64_real2(), genrand64_real2()/)*L
        theta_i = genrand64_real2()*6.283185307d0
        unitary_vel(:,i) = (/cos(theta_i),sin(theta_i)/)
    end do
  end subroutine

  pure function dist_PBC(pos1,pos2)
  real(8), INTENT(IN) :: pos1(2), pos2(2)
  real(8) :: dist_PBC, dx, dy

  dx = wrap_dist(pos1(1) - pos2(1))
!   dx = dx - int(inv_2L*dx)*L

  dy = wrap_dist(pos1(2) - pos2(2))
!   dy = dy - int(inv_2L*dy)*L

  dist_PBC = dx*dx + dy*dy

  end function dist_PBC

  pure real(8) function  wrap_dist(val)
  real(8), INTENT(IN) :: val
  wrap_dist = val - int(inv_2L*val)*L
    end function wrap_dist

  pure function dist_simple(pos1,pos2)
  real(8), INTENT(IN) :: pos1(2), pos2(2)
  real(8) :: dist_simple, dx, dy

  dx = pos1(1) - pos2(1)
  dy = pos1(2) - pos2(2)

  dist_simple = dx*dx + dy*dy

  end function dist_simple





end module




!   subroutine neighbours_direction(nbr_direction)
!     real(8) :: nbr_direction(N), integrated_vel(N,2)
!     integer :: cell, header(N_cells), cell_list(N), nbr_cell, j, nbr_cell_index, where_is(N)

!     header = 0

!     do i = 1,N
!         cell = int(pos(i,1)) + int_L*int(pos(i,2)) + 1
!         cell_list(i) = header(cell)
!         header(cell) = i
!         where_is(i) = cell;
!         integrated_vel(:,i) = unitary_vel(i,:)
!     end do

!     do i=1,N
!         cell = where_is(i)

!         do nbr_cell_index=1,9
!             nbr_cell = list_nearest_nbr_cells(cell, nbr_cell_index)
!             j = header(nbr_cell)
!             do while (j > 0)
!                 if (i<j) then
!                     if (dist_PBC(pos(i,:), pos(j,:))>1.d0) then
!                         integrated_vel(:,i) = integrated_vel(:,i) + unitary_vel(j,:)
!                         integrated_vel(:,j) = integrated_vel(:,j) + unitary_vel(i,:)
!                     end if
!                 end if
!                 j = cell_list(j)
!             end do
!         end do
!         nbr_direction(i) = atan2(integrated_vel(i,2), integrated_vel(i,1))
!     end do

!   end subroutine
