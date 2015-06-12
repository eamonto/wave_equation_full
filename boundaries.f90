
! ===========================================================================
! boundaries.f90
! ===========================================================================
! Implementation of the boundaries conditions

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


  subroutine boundaries(phi,pi,psi,x,dx,alpha,beta,Np)

    use mylibrary
    use param

    implicit none

    type(dynamical_func) phi,pi,psi
    type(extra_func)     x,dx   
    type(extra_func)     alpha,beta

    integer :: Np

    real(double) :: der_pi,der_psi

    !Periodic
    if(boundary.eq.1) then

       phi%s(Np) = alpha%f(Np)*pi%f(Np) + beta%f(Np)*psi%f(Np)

       pi%s(Np) = (alpha%f(1)*psi%f(1)-alpha%f(Np-1)*psi%f(Np-1))/(2.0D0*dx%f(Np)) &
                  +(beta%f(1)*pi%f(1)  -beta%f(Np-1)*pi%f(Np-1))/(2.0D0*dx%f(Np))

       psi%s(Np) = (alpha%f(1)*pi%f(1)-alpha%f(Np-1)*pi%f(Np-1))/(2.0D0*dx%f(Np)) &
                  +(beta%f(1)*psi%f(1)-beta%f(Np-1)*psi%f(Np-1))/(2.0D0*dx%f(Np))

       phi%s(0) = phi%s(Np)
       pi%s(0)  = pi%s(Np)
       psi%s(0) = psi%s(Np) 

    !Free
    else if(boundary.eq.2) then

       der_pi  = (- pi%f(2)+4.0D0* pi%f(1)-3.0D0* pi%f(0))/(2.0D0*dx%f(0))
       der_psi = (-psi%f(2)+4.0D0*psi%f(1)-3.0D0*psi%f(0))/(2.0D0*dx%f(0))

       phi%s(0) = alpha%f(0)*pi%f(0) + beta%f(0)*psi%f(0)
       pi%s(0)  = (der_pi+der_psi)/2.0D0
       psi%s(0) = pi%s(0)

       der_pi  = ( pi%f(Np-2)-4.0D0* pi%f(Np-1)+3.0D0 *pi%f(Np))/(2.0D0*dx%f(Np))
       der_psi = (psi%f(Np-2)-4.0D0*psi%f(Np-1)+3.0D0*psi%f(Np))/(2.0D0*dx%f(Np))

       phi%s(Np) = alpha%f(Np)*pi%f(Np) + beta%f(Np)*psi%f(Np)
       pi%s(Np) = -(der_pi-der_psi)/2.0D0
       psi%s(Np) = -pi%s(Np)

    !Fix
    else if(boundary.eq.3) then

       der_pi  = (- pi%f(2)+4.0D0* pi%f(1)-3.0D0* pi%f(0))/(2.0D0*dx%f(0))
       der_psi = (-psi%f(2)+4.0D0*psi%f(1)-3.0D0*psi%f(0))/(2.0D0*dx%f(0))

       phi%s(0) = 0.0D0
       pi%s(0)  = 0.0D0
       psi%s(0) = (der_pi+der_psi)/2.0D0

       der_pi  = ( pi%f(Np-2)-4.0D0* pi%f(Np-1)+3.0D0 *pi%f(Np))/(2.0D0*dx%f(Np))
       der_psi = (psi%f(Np-2)-4.0D0*psi%f(Np-1)+3.0D0*psi%f(Np))/(2.0D0*dx%f(Np))

       phi%s(Np) = 0.0D0
       pi%s(Np)  = 0.0D0
       psi%s(Np) = (der_pi-der_psi)/2.0D0

    else

       print*,'Boundary condition not implemented!'
       stop
       
    endif

  end subroutine boundaries
