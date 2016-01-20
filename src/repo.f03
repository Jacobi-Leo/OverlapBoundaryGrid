module repo
  !**********************************************************
  ! repo is the repositary of the program which provide basic
  ! parameters and basic subroutines.
  !**********************************************************
  implicit none
  integer, parameter :: n=32 ! number of nodes of grid
  integer, parameter :: m=5  ! number of overlapping level
  integer, parameter :: iteration=300 ! number of iteration
  real, parameter :: pi = atan(1.0)*4.0 ! some mathematical constants
  real, parameter :: U0 = 1, L = 1, nu = 1.e-5 ! some physical constant
  real, parameter :: dt = 1.9e-3 ! some computation constant
contains

  subroutine ChebGridGenarator1dRightHalf (n, nodes, lengths)
    implicit none
    !! This subroutine generate Chebshev grid on (0,1) with
    !! x=0 the most coarse and x=1 the most refined.
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

  subroutine BoundaryOverlappingLeft (n, m, grid_size)
    implicit none
    integer, intent(in) :: n, m
    real, intent(inout), dimension(0:n+1) :: grid_size
    integer :: i
    do i = 2, m
       grid_size(n-m+i) = grid_size(n-m+i) + grid_size(n-m+i-1)
    end do
  end subroutine BoundaryOverlappingLeft

  subroutine Show1dGrid (n, nodes, lengths)
    implicit none
    !! The subroutine to print out grid node coordinates
    !! and the size of each element in 1D case.
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