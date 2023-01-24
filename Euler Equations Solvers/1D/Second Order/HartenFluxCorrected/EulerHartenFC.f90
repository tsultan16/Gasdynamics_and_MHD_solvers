!Harten's First Order Upwind Method for Euler's eqautions

program EulerHartenFC
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
  
  print*,'dt=',dt
  write(12,*) dt
 
  do j=1,nx    
    !update cell averages of conerved variables
    do k=1,3    
     u2(j,k)=u1(j,k)+(dt/dx)*(g(j-1,k)-g(j,k))
    end do
    x=xmin+(j-0.5)*dx 
    dens=u2(j,1)
    vel=u2(j,2)/u2(j,1)
    pres=(gam-1.)*(u2(j,3)-0.5*u2(j,2)*u2(j,2)/u2(j,1))
    write(11,*) x,dens,vel,pres
  
  end do

  !Enforce boundary conditions
  call bound()
  u1=u2

end do 


print*,'Done.'
close(unit=10)
close(unit=11)
call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start

end program EulerHartenFC
