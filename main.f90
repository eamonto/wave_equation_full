
! ===========================================================================
! main.f90
! ===========================================================================
! Principal program, this program solves the wave equation in 1+1 dimensions,
! using the 1+1 decomposition
!
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


  program main

    use param
    use mylibrary

    implicit none

    integer :: l,k

    type(dynamical_func) phi,pi,psi
    type(extra_func)     x,dx   
    type(extra_func)     alpha,beta

    !Allocation of function that are evolve
    call allocate_dyn(phi,Nx)  !NO USER
    call allocate_dyn(pi ,Nx)  !NO USER
    call allocate_dyn(psi,Nx)  !NO USER

    !Allocation of function that are not evolve
    call allocate_extra(x ,Nx)  !NO USER
    call allocate_extra(dx,Nx)  !NO USER
    call allocate_extra(alpha,Nx)  !NO USER
    call allocate_extra(beta ,Nx)  !NO USER

    !Initialization of the functions, including the grid
    call initial_data(phi,pi,psi,x,dx,alpha,beta,Nx)  !USER

    !Created the output file for the dynamical functions "func" in  
    !the "output_dir" with name "output_file" and id "file_number"
    call create_output(phi,output_dir,'phi.x',200)     !NO USER
    call create_output(pi ,output_dir,"pi.x" ,201)     !NO USER
    call create_output(psi,output_dir,"psi.x",202)     !NO USER

    !Created the output file for the extra function "func" in  
    !the "output_dir" with name "output_file" and id "file_number"
    call create_output_extra(alpha,output_dir,'alpha.x',300)     !NO USER
    call create_output_extra(beta ,output_dir,'beta.x' ,301)     !NO USER

    !Print the output of two scalar functions on the grid
    call output_obs_obs(x%f,phi%f,phi%name,phi%id,Nx) !NO USER
    call output_obs_obs(x%f, pi%f, pi%name, pi%id,Nx) !NO USER
    call output_obs_obs(x%f,psi%f,psi%name,psi%id,Nx) !NO USER

    call output_obs_obs(x%f,alpha%f,alpha%name,alpha%id,Nx) !NO USER
    call output_obs_obs(x%f, beta%f, beta%name, beta%id,Nx) !NO USER
 
    print *,''
    print *,'    Time     '
    print *,'-------------'
    write(*,"(F12.5)") time

    do l=1,Ntime

       !Store the initial values before the integration
       call store_levels_rk4(phi)  !NO USER
       call store_levels_rk4(pi)   !NO USER
       call store_levels_rk4(psi)  !NO USER

       do k=1,4

          !Derivatives of the functions to be evolve
          call sources(phi,pi,psi,x,dx,alpha,beta,Nx)    !USER

          !Boundaries
          call boundaries(phi,pi,psi,x,dx,alpha,beta,Nx)   !USER

          !Implemetation of the Runge-Kutta 4 method 
          call evolution_rk4(k,phi,dt)  !NO USER
          call evolution_rk4(k,pi ,dt)  !NO USER
          call evolution_rk4(k,psi,dt)  !NO USER

       enddo

       time = time + dt

       if(mod(l,every_0D).eq.0) then
          write(*,"(F12.5)") time
       endif

       !Print the output of two scalar functions on the grid
       if(mod(l,every_1D).eq.0) then
          call output_obs_obs(x%f, phi%f, phi%name, phi%id,Nx)  !NO USER
          call output_obs_obs(x%f,  pi%f,  pi%name,  pi%id,Nx)  !NO USER
          call output_obs_obs(x%f, psi%f, psi%name, psi%id,Nx)  !NO USER
       endif

    enddo

    !Deallocation of function that are evolve
    call deallocate_dyn(phi)  !NO USER
    call deallocate_dyn(pi )  !NO USER
    call deallocate_dyn(psi)  !NO USER

    !Deallocation of function that are not evolve
    call deallocate_extra(x )  !NO USER
    call deallocate_extra(dx)  !NO USER
    call deallocate_extra(alpha)  !NO USER
    call deallocate_extra(beta )  !NO USER

    print *,''
    print *,'   Finish'
    print *,''

  end program main
