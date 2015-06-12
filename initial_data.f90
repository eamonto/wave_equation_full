
! ===========================================================================
! initial_data.f90
! ===========================================================================
! Here are initialized all the functions defined on the grid.

!     Copyright (C) 2015  Edison Montoya, eamonto@gmail.com

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

! Up to date: 12 Jun 2015					


  subroutine initial_data(phi,pi,psi,x,dx,alpha,beta,grid_points)

    use param
    use mylibrary

    implicit none

    type(dynamical_func) phi,pi,psi
    type(extra_func)     x,dx   
    type(extra_func)     alpha,beta

    integer :: grid_points

    real(double) :: amp   = 1.0D0
    real(double) :: x0    = 0.0D0
    real(double) :: sigma = 1.0D0
    integer :: i

    dx%f = dx_aux

    do i=0,grid_points
       x%f(i) = -dx%f(i)*grid_points/2.0D0 + i*dx%f(i)
    enddo

    phi%f = amp*exp(-(x%f-x0)**2/sigma**2)
    
    pi%f  = 0.0D0
    
    psi%f = -2.0D0*(x%f-x0)/sigma**2*phi%f

    alpha%f = 1.0D0

    beta%f = 0.0D0

  end subroutine initial_data
