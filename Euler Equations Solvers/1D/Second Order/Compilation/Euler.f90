!Second Order Godunov-type Methods for Euler's eqautions

program Euler
use RoeSolver_mod
implicit none


real::start,finish
call cpu_time(start)
call init()

!Compute Solution
do i=1,nt
  print*,'TIME STEP=',i
  !compute numerical fluxes
  if(flux_option==1)then
    call computeFluxMH()
  else if(flux_option==1)then
    call computeFluxHarten()
  else if(flux_option==3)then
    call computeFluxENO()
  end if

  print*,'dt=',dt
  write(12,*) dt
 

  do j=1,nx    
    !update cell averages of conerved variables
    do k=1,5    
     u2(j,k)=u1(j,k)+(dt/dx)*(g(j-1,k)-g(j,k))
    end do
    x=xmin+(j-0.5)*dx 
    dens=u2(j,1)
    velx=u2(j,2)/u2(j,1)
    vely=u2(j,3)/u2(j,1)
    velz=u2(j,4)/u2(j,1)
    pres=(gam-1.)*( u1(j,5)-0.5*dens &
         *(velx**2.+vely**2.+velz**2.) )
    write(11,*) x,dens,velx,pres
  
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

end program Euler
