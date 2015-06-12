
! ===========================================================================
! sources.f90
! ===========================================================================
! Here are implemented the right hand side (derivatives) of the functions
! to be evolve in time.

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


  subroutine sources(phi,pi,psi,x,dx,alpha,beta,grid_points)
    
    use mylibrary

    implicit none

    type(dynamical_func) phi,pi,psi
    type(extra_func) x,dx   
    type(extra_func) alpha,beta

    integer :: grid_points

    integer :: i

    do i=1,grid_points-1

       phi%s(i) = alpha%f(i)*pi%f(i) + beta%f(i)*psi%f(i)

       pi%s(i) = (alpha%f(i+1)*psi%f(i+1)-alpha%f(i-1)*psi%f(i-1))/(2.0D0*dx%f(i)) &
                 +(beta%f(i+1)*pi%f(i+1)  -beta%f(i-1)*pi%f(i-1))/(2.0D0*dx%f(i))

       psi%s(i) = (alpha%f(i+1)*pi%f(i+1)-alpha%f(i-1)*pi%f(i-1))/(2.0D0*dx%f(i)) &
                  +(beta%f(i+1)*psi%f(i+1)-beta%f(i-1)*psi%f(i-1))/(2.0D0*dx%f(i))
    enddo

  end subroutine sources
