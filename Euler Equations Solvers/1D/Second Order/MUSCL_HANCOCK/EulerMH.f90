!Godunovs method (1st order) fopr Eulers Equation in 1D
!(Improved version)

program EulerMH
use RoeSolver_entropyfixed_mod
implicit none

real::start,finish
call cpu_time(start)
call init()

!Compute Solution
do i=1,nt
  print*,'TIME STEP=',i
  !compute numerical fluxes
  call computeFlux()
  !set time-step size
  dt=dx*cour/smax
  print*,'dt=',dt
  write(12,*) dt
  do j=1,nx 
    
    !update cell averages of conerved variables
    do k=1,3    
     u2(j,k)=u1(j,k)+(dt/dx)*(flux(j-1,k)-flux(j,k))
    end do
    x=xmin+(j-0.5)*dx 
    dens=u2(j,1)
    vel=u2(j,2)/u2(j,1)
    pres=(gam-1.)*(u2(j,3)-0.5*u2(j,2)*u2(j,2)/u2(j,1))
    write(11,*) x,dens,vel,pres

  end do
  !enforce outflow boundary conditions
  !do k=1,3
    !u2(0,k)=u2(1,k)	
    !u2(-1,k)=u2(1,k)
    !u2(nx+1,k)=u2(nx,k)
    !u2(nx+2,k)=u2(nx,k)
    !outflow boundary on the right and reflecting on the left
    u2(0,1)=u2(1,1)	
    u2(-1,1)=u2(1,1)
    u2(0,2)=-u2(1,2)	
    u2(-1,2)=-u2(1,2)
    u2(0,3)=u2(1,3)	
    u2(-1,3)=u2(1,3)
    do k=1,3
     u2(nx+1,k)=u2(nx,k)
     u2(nx+2,k)=u2(nx,k)
    end do	
  u1=u2

end do 


print*,'Done.'
close(unit=10)
close(unit=11)
call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

end program EulerMH
