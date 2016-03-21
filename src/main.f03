! ******************** START ********************
! Created by Liu Zeyu (liuzeyu271828@gmail.com)
! Distributed under GPL v3.0
! ********************  END  ********************

program main
  use repo
  use fortran_gnuplot
  implicit none

  real, allocatable, dimension(:) :: grid_node, grid_size, grid_size_new
  real, allocatable, dimension(:) :: u, u_new, u_tmp, u_output, u_exact
  real, allocatable, dimension(:,:) :: xydata, xydata2
  ! real, allocatable, dimension(:) :: De, Dw, Fe, Fw, Pe, Pw, aE, aW, aP
  ! real :: De, Dw, Fe, Fw, Pe, Pw, aE, aW
  real :: uplus, uminus, uxplus, uxminus
  integer :: i, j

  allocate(grid_node(0:n), grid_size(0:n+1), grid_size_new(0:n+1), &
       &   u_tmp(n), u(0:n+1), u_new(0:n+1), u_output(0:n+1), u_exact(0:n+1), &
       &   xydata(n,2), xydata2(n,2))
  
  ! allocate(De(n), Dw(n), Fe(n), Fw(n), Pe(n), &
  !      &   Pw(n), aE(n), aW(n), aP(n))

  !! Generate grid node
  ! call UniformGridGenerator (n, grid_node, grid_size(1:n))
  call ChebGridGenarator1dRightHalf (n, grid_node, grid_size(1:n))
  grid_node = grid_node * L
  grid_size = grid_size * L

  !! The 0th and (n+1)th node are auxiliary node to satisfy the 1st kind of BC
  grid_size(0) = grid_size(1)
  grid_size(n+1) = grid_size(n)
  
  open(10, status='replace', file=rawoutfile, action='write')
  open(11, status='replace', file='grid.dat', action='write')
  write(11,*) grid_node(1:n) - 0.5*grid_size(1:n)

  do i=1,n
     u_exact(i) = (1. - exp((grid_node(i)-0.5*grid_size(i) - 1.)/nu))/(1. - exp(-1./nu))
  end do
  u_exact(0) = 2. - u_exact(1)
  u_exact(n+1) = 0. - u_exact(n)

  do i=1,n
     u(i) = u0*(L-grid_node(i)+0.5*grid_size(i))
  enddo
  !! The following two equations are simplified version
  u(0) = 2. - u(1)
  u(n+1) = 0. - u(n)
  
  write(10, *) u(1:n)
  
  call BoundaryOverlappingRight (n, m, grid_size, grid_size_new)
  call Converter(n, m, grid_size_new, u, u_new)
    
  ! !! Forward Euler method
  iterations: do j = 1, iteration
     do i = 1, n-m
        uplus = findu(grid_node, u, n, i, f1, f2)
        
        if ( i == 1 ) then
           uminus = 1.
        else
           uminus = findu(grid_node, u, n, i-1, f1, f2)
        end if
        
        uxplus = findux(grid_node, u, n, i, g1, g2)
        
        if ( i == 1 ) then
           uxminus = (u(1)-u(0)) / grid_size(0)
        else
           uminus = findux(grid_node, u, n, i-1, g1, g2)
        end if
        
        u_tmp(i) = u(i) + dt/grid_size(i)*(U0*(uminus-uplus)+nu*(uxplus-uxminus))
     end do

     do i = n-m+1, n
        if ( i == n ) then
           uplus = 0.
        else
           uplus = findu(grid_node, u, n, i, f1, f2)
        end if
        
        uminus = findu(grid_node, u, n, n-m, f1, f2)
        uxplus = findux(grid_node, u, n, i, g1, g2)
        uxminus = findux(grid_node, u, n, n-m, g1, g2)
        u_tmp(i) = u_new(i) + dt/grid_size_new(i)*(U0*(uminus-uplus)+nu*(uxplus-uxminus))
     end do
     
     u_new(1:n) = u_tmp(1:n)
     u_new(0) = 2. - u_new(1)
     u_new(n+1) = 0. - u_new(n)

     call ConverterReverse (n, m, grid_size_new, u_new, u)
     write(10, *) u(1:n)
  end do iterations

  xydata(:,1) = grid_node(1:n) - 0.5*grid_size(1:n)
  xydata(:,2) = u(1:n)
  xydata2(:,1) = xydata(:,1)
  xydata2(:,2) = u_exact(1:n)
  call f2gp (n, n, xydata, xydata2, 1, 'x', 'u', 'u', 'exact')
  deallocate(grid_node, grid_size, u, u_new, u_tmp, grid_size_new)

  !! The end of the program.
  !! Output this for debug purpose, though not in the constrains of Unix philosophy.
  write(*,*) 'Program finished.'


contains
  
  pure function findu (xx, u, n, j)
    use repo, only: f1, f2
    !! Interpolation of u, according to Mr. Cai,
    !! This algorithm has been verified.
    implicit none
    integer, intent(in) :: n, j
    real, dimension(0:n+1), intent(in) :: u
    real, dimension(0:n), intent(in) :: xx
    real :: h1, h2, h3, findu
    real, dimension(-1:n+1) :: x

    x(0:n) = xx
    x(-1) = 2.*x(0) - x(1)
    x(n) = 2.*x(n) - x(n-1)

    !! the old version which should be equivalent
    ! h1 = ((x(j)-x(j-1))*(x(j)-x(j+1))) / ((x(j-2)-x(j))*(x(j-2)-x(j+1)))
    ! h3 = ((x(j-2)-x(j))*(x(j)-x(j-1))) / ((x(j-2)-x(j+1))*(-x(j-1)+x(j+1)))
    ! h2 = 1. - h1 - h3

    h1 = f2(x(j-2), x(j-1), x(j), x(j+1))
    h3 = f1(x(j+1), x(j), x(j-1), x(j-1))
    h2 = 1. - h1 - h3

    findu = u(j-1)*h1 + u(j)*h2 + u(j+1)*h3
  end function findu


  pure function findux (xx, u, n, j)
    use repo, only: g1, g2
    !! Interpolation of u_x, according to Mr. Cai.
    !! This algorithm has been verified.
    implicit none
    integer, intent(in) :: n, j    
    real, dimension(0:n+1), intent(in) :: u
    real, dimension(0:n), intent(in) :: xx
    real :: h1, h2, h3, findux
    real, dimension(-1:n+1) :: x

    x(0:n) = xx
    x(-1) = 2.*x(0) - x(1)
    x(n) = 2.*x(n) - x(n-1)

    !! One method that seems not work
    ! h1 = -2.*(x(j-1)-2.*x(j)+x(j+1)) / ((x(j-2)-x(j))*(x(j-2)-x(j+1)))
    ! h3 = 2.*(x(j-2)+x(j-1)-2.*x(j)) / ((x(j-2)-x(j+1))*(x(j+1)-x(j-1)))
    ! h2 = 0. - h1 - h3
    ! findux = u(j-1)*h1+u(j)*h2+u(j+1)*h3

    h1 = g2(x(j-2), x(j-1), x(j), x(j+1))
    h3 = g1(x(j+1), x(j), x(j-1), x(j-1))
    h2 = 0. - h1 - h3

    findux = u(j-1)*h1 + u(j)*h2 + u(j+1)*h3

  end function findux

 
end program

