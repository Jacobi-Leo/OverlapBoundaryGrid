! ******************** START ********************
! Created by Liu Zeyu (liuzeyu271828@gmail.com)o
! Distributed under GPL v3.0
! ********************  END  ********************
module repo
  !**********************************************************
  ! repo is the repositary of the program which provide basic
  ! parameters and basic subroutines.
  !**********************************************************
  implicit none
  integer, parameter :: n=256 ! number of nodes of grid
  real, parameter :: pi = atan(1.0)*4.0, U0=1., L=1.
  real, parameter :: dt = 1.e-5

contains

  subroutine ChebGridGenarator1dRightHalf (n, nodes, lengths)
    implicit none
    integer, intent(in) :: n
    real, intent(out), dimension(0:n) :: nodes
    real, intent(out), dimension(n) :: lengths
    integer :: i

    nodes(0) = 0.
    do i=1,n
       lengths(i) = 2*sin(pi/4.0/real(n))*cos(pi/real(n)/4.0*(2*real(i)-1))
       nodes(i) = nodes(i-1) + lengths(i)
    enddo
  end subroutine ChebGridGenarator1dRightHalf

  subroutine Show1dGrid (n, nodes, lengths)
    implicit none
    integer, intent(in) :: n
    real, intent(in), dimension(n) :: lengths
    real, intent(in), dimension(0:n) :: nodes
    integer :: i
    
    100 format(' ', 'grid_node = ', ES16.9, '     grid_size = ', ES16.9)
    do i=1,n
       write(*,100) nodes(i), lengths(i)
    enddo
  end subroutine Show1dGrid

end module repo


program main
  use repo
  use fortran_gnuplot
  implicit none

  real, allocatable, dimension(:) :: grid_node, grid_size, u
  real, allocatable, dimension(:,:) :: xydata
  integer :: i

  allocate(grid_node(0:n), grid_size(n), u(n), xydata(n,2))
  call ChebGridGenarator1dRightHalf  (n, grid_node, grid_size)
  grid_node = grid_node * L
  grid_size = grid_size * L

  do i=1,n
     u(i) = u0*(L-grid_node(i)+0.5*grid_size(i))
  enddo

  

  xydata(:,1) = grid_node(1:n)
  xydata(:,2) = u(:)
  ! call f2gp (n, 0, xydata, xydata, 1, 'x', 'u', '1', '1')
  deallocate(grid_node, grid_size)
end program main

