! ******************** START ********************
! Created by Liu Zeyu (liuzeyu271828@gmail.com)
! Distributed under GPL v3.0
! ********************  END  ********************

program main
  use repo
  use fortran_gnuplot
  implicit none

  real, allocatable, dimension(:) :: grid_node, grid_size, u, u_tmp
  real, allocatable, dimension(:,:) :: xydata
  real, allocatable, dimension(:) :: De, Dw, Fe, Fw, Pe, Pw, aE, aW, aP
  integer :: i, j

  allocate(grid_node(0:n), grid_size(0:n+1), u_tmp(n), xydata(n,2), u(0:n+1))
  allocate(De(n), Dw(n), Fe(n), Fw(n), Pe(n), Pw(n), aE(n), aW(n), aP(n))
  call ChebGridGenarator1dRightHalf (n, grid_node, grid_size(1:n))
  grid_node = grid_node * L
  grid_size = grid_size * L
  grid_size(0) = grid_size(1)
  grid_size(n+1) = grid_size(n)

  write(*,*) grid_size(n)
  write(*,*) grid_size(n-m+1)

  do i=1,n
     u(i) = u0*(L-grid_node(i)+0.5*grid_size(i))
  enddo
  u(0) = u(1)
  u(n+1) = u(n)

  do i = 1, n
     De(i) = nu/(0.5*(grid_size(i)+grid_size(i+1)))
     Dw(i) = nu/(0.5*(grid_size(i)+grid_size(i-1)))
     !! optimization can be done here
     Fe(i) = 1.0
     Fw(i) = -1.0
     Pe(i) = Fe(i)/De(i)
     Pw(i) = Fw(i)/Dw(i)
     aE(i) = De(i)*(MAX(0., (1.-0.1*abs(Pe(i)))**5) + MAX(0., -Pe(i)))
     aW(i) = Dw(i)*(MAX(0., (1.-0.1*abs(Pw(i)))**5) + MAX(0., -Pw(i)))
     aP(i) = aE(i) + aW(i)
  enddo
  
  call BoundaryOverlappingLeft (n, m, grid_size)

  do i=n-m+2, n
     u(i) = u0*(L-grid_node(n-m+1)-0.5*grid_size(i))
  enddo
  u(0) = u(1)
  u(n+1) = u(n)

  do j = 1, iteration
     do i = 1, n-m+1
        u_tmp(i) = aE(i)*(u(i+1)-u(i))+aW(i)*(u(i-1)-u(i))*dt/grid_size(i) + u(i)
     enddo

     do i = n-m+2, n
        u_tmp(i) = aE(i)*(u(i+1)-u(i))+aW(n-m-1)*(u(i-1)-u(i))*dt/grid_size(i) + u(i)
     end do
     
     u(1:n) = u_tmp(1:n)
  enddo


  xydata(:,1) = grid_node(1:n)
  xydata(:,2) = u(1:n)
  call f2gp (n, 0, xydata, xydata, 1, 'x', 'u', '1', '1')
  deallocate(grid_node, grid_size)
end program main

