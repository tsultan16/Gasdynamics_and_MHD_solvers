!Driver for WENO5 MHD Solver

program WENO5MHD
use globalVariables_mod
use init_mod
use WENOflux_mod
implicit none

real::q(-3:nx+4,8)
integer::tSteps

open(unit=10,file='output_fluid.txt')
open(unit=11,file='output_mag.txt')


!Initialize spatial grid and MHD state variables
call initialize(q) 
!file output
call fileOut(q)
print*,'Done initializing.'

!Simulation loop (w/ fixed total # of time steps)
do tSteps=1,nt

  print*,'Time Step:',tSteps
  !compute dt
  print*,'Computing dt...'
  call timeStep(q,dt)
  print*,'dt=',dt

  !Runge-Kutta time integration: q^n --> q^n+1
  call RK3TVD_integrator(q,dt)

  !file output
  call fileOut(q)

end do

print*,'Simulation Completed.'

contains

subroutine RK3TVD_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4,8),dt
!Subroutine local variables
real::q1(-3:nx+4,8),q2(-3:nx+4,8)
real::flux0(0:nx,8),flux1(0:nx,8),flux2(0:nx,8)
integer::i,j

  !Step 1:
  call compute_flux(q,flux0)
  do i=1,nx
    do j=1,8
      q1(i,j)=q(i,j)-(dt/dx)*(flux0(i,j)-flux0(i-1,j))
    end do
  end do  

  call bound(q1)
  
  !Step2:
  call compute_flux(q1,flux1)
  do i=1,nx
    do j=1,8
      q2(i,j)=q1(i,j)-(1./4.)*(dt/dx)*(-3.*(flux0(i,j)-flux0(i-1,j)) &
              +(flux1(i,j)-flux1(i-1,j)) )
    end do
  end do

  call bound(q2)

  !Step3:
  call compute_flux(q2,flux2)
  do i=1,nx
    do j=1,8 
      q(i,j)=q2(i,j)-(1./12.)*(dt/dx)*(-(flux0(i,j)-flux0(i-1,j))-(flux1(i,j)-flux1(i-1,j))&
        +8.*(flux2(i,j)-flux2(i-1,j)))
    end do
  end do

  call bound(q)

end subroutine RK3TVD_integrator


subroutine bound(q)

!Subroutine I/O variables
real::q(-3:nx+4,8)
!Subroutine local variables
integer::j

!Outflow boundary conditions
do j=1,8
  q(0,j)=q(1,j)	
  q(-1,j)=q(1,j)
  q(-2,j)=q(1,j)
  q(-3,j)=q(1,j)
  q(nx+1,j)=q(nx,j)
  q(nx+2,j)=q(nx,j)
  q(nx+3,j)=q(nx,j)
  q(nx+4,j)=q(nx,j)
end do

end subroutine bound

subroutine fileOut(q)
!Subroutine I/O variables
real::q(-3:nx+4,8)
!Subroutine local variables
integer::i
real::x
real::dens,velx,vely,velz,pres,magx,magy,magz

do i=1,nx
  x=xmin+(i-1)*dx
  dens=q(i,1)
  velx=q(i,2)/q(i,1)
  vely=q(i,3)/q(i,1)
  velz=q(i,4)/q(i,1)
  magx=q(i,5)
  magy=q(i,6)
  magz=q(i,7)
  pres=(gam-1.)*(q(i,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  write(10,*) x,dens,velx,vely,velz,pres
  write(11,*) x,magy,magz
end do


end subroutine fileOut

end program WENO5MHD
