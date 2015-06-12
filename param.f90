
! ===========================================================================
! param.f90
! ===========================================================================
! Global parameters for the physical system.

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


  module param

    use mylibrary

    implicit none

    !evolution
    integer,parameter :: Ntime = 10000         !Number of time steps

    integer,parameter :: Nx = 1000             !Number of points in the grid
    real(double),parameter :: length = 10.0D0  !length

    integer :: boundary = 1                    !1=periodic 2=eigenvalue 3=fix

    character(100) :: output_dir = "test"      !output directory

    real(double),parameter :: dx_aux = length/(1.0D0*Nx)   !dx homogeneous grid
    real(double),parameter :: courant = 0.2D0  !Courant factor

    real(double) :: dt = courant*dx_aux        !time step

    real(double) :: time =0.0D0                !initial time

    integer :: every_0D = 100                  !time output
    integer :: every_1D = 100                  !spatial output

  end module param
