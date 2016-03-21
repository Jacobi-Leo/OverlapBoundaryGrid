module repo
  !**********************************************************
  ! repo is the repositary of the program which provide basic
  ! parameters and basic subroutines.
  !**********************************************************
  implicit none
  integer, parameter :: n=32 ! number of nodes of grid
  integer, parameter :: m=1 ! Number of overlapping level
  integer, parameter :: iteration=200000 ! number of iteration
  real, parameter :: pi = atan(1.0)*4.0 ! some mathematical constants
  real, parameter :: U0 = 1., L = 1., nu = 0.1 ! some physical constant
  real, parameter :: dt = 1.e-7 ! some computation constant
  character(len=*), parameter:: rawoutfile = 'raw.dat'

  
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

  subroutine UniformGridGenerator (n, nodes, lengths)
    implicit none
    integer, intent(in) :: n
    real, intent(out), dimension(0:n) :: nodes
    real, intent(out), dimension(n) :: lengths
    integer :: i
    real :: tmp

    nodes(0) = 0.
    tmp = 1. / real(n)
    do i = 1, n
       lengths(i) = tmp
       nodes(i) = tmp * real(i)
    end do
  end subroutine UniformGridGenerator
  

  subroutine BoundaryOverlappingLeft (n, m, grid_size, grid_size_new)
    implicit none
    integer, intent(in) :: n, m
    real, intent(in), dimension(0:n+1) :: grid_size
    real, intent(out), dimension(0:n+1) :: grid_size_new
    integer :: i
    grid_size_new = grid_size !!! This can be optimized
    do i = 1, m-1
       grid_size_new(m-i) = grid_size(m-i) + grid_size_new(m-i+1)
    end do
  end subroutine BoundaryOverlappingLeft
  

  subroutine BoundaryOverlappingRight (n, m, grid_size, grid_size_new)
    implicit none
    integer, intent(in) :: n, m
    real, intent(in), dimension(0:n+1) :: grid_size
    real, intent(out), dimension(0:n+1) :: grid_size_new
    integer :: i
    grid_size_new = grid_size !!! This can be optimized
    do i = n-m+2, n
       grid_size_new(i) = grid_size_new(i-1) + grid_size(i)
    end do
  end subroutine BoundaryOverlappingRight
  

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
  

  subroutine Converter (n, m, grid_size_new, u, u_new)
    ! From normal to overlapped grid
    implicit none
    integer, intent(in) :: n, m
    real, intent(in), dimension(0:n+1) :: u, grid_size_new
    real, intent(out), dimension(0:n+1) :: u_new
    integer :: i
    real :: sum

    sum = u(n-m+1)*grid_size_new(n-m+1)

    do i=0, n-m+1
       u_new(i) = u(i)
    end do

    do i=n-m+2, n
       sum = sum + u(i)*(grid_size_new(i)-grid_size_new(i-1)) ! This can be optimized
       u_new(i) = sum / grid_size_new(i)
    end do
    u_new(n+1) = u(n+1)
  end subroutine Converter
  

  subroutine ConverterReverse (n, m, grid_size_new, u, u_new)
    ! From overlapped to normal
    implicit none
    integer, intent(in) :: n, m
    real, intent(in), dimension(0:n+1) :: u, grid_size_new
    real, intent(out), dimension(0:n+1) :: u_new
    integer :: i

    do i=0, n-m+1
       u_new(i) = u(i)
    end do

    do i=n-m+2, n
       u_new(i) = (u(i)*grid_size_new(i) - u(i-1)*grid_size_new(i-1)) /&
            & (grid_size_new(i) - grid_size_new(i-1))
    end do
    u_new(n+1) = u(n+1)
  end subroutine ConverterReverse

  pure function f0 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: f0, a, b, c
    a = (3.*x0 - 2.*(x1+x2+x3))*x0
    b = x1*x2 + x2*x3 + x1*x3
    c = (x2-x0)*(x3-x0)
    f0 = (a+b)/c
  end function f0

  pure function f1 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: f1, a, b
    a = (x2-x1)*(x3-x1)
    b = (x2-x0)*(x3-x0)
    f1 = a/b
  end function f1

  pure function f2 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: f2, a, b
    a = (x1-x2)*(x3-x2)
    b = (x2-x0)*(x3-x0)
    f2 = a/b
  end function f2

  pure function f3 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: f3, a, b
    a = (x1-x3)*(x2-x3)
    b = (x2-x0)*(x3-x0)
    f3 = a/b
  end function f3

  pure function g0 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: g0, a, b
    a = 2.*(3.*x0-(x1+x2+x3))
    b = (x2-x0)*(x3-x0)
    g0 = a/b
  end function g0

  pure function g1 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: g1, a, b
    a = 2.*(2.*x1-(x2+x3))
    b = (x2-x0)*(x3-x0)
    g1 = a/b
  end function g1

  pure function g2 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: g2, a, b
    a = 2.*(2.*x2-(x1+x3))
    b = (x2-x0)*(x3-x0)
    g2 = a/b
  end function g2

  pure function g3 (x0, x1, x2, x3)
    implicit none
    real, intent(in) :: x0, x1, x2, x3
    real :: g3, a, b
    a = 2.*(2.*x3-(x2+x1))
    b = (x2-x0)*(x3-x0)
    g3 = a/b
  end function g3



end module repo
