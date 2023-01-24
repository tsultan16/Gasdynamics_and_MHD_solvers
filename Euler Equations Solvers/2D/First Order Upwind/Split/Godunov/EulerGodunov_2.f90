!Godunovs method (1st order) fopr Eulers Equation in 2D
!with dimensional splitting

program EulerGodunov2D
use RoeSolver2D_entropyfix_mod
implicit none

integer::i,j,k,l
character(len=3)::uniti
real::start,finish

call cpu_time(start)
call init()

!Compute Solution 
do i=1,nt
  print*,'% Complete=',(real(i)/real(nt))*100.

  !-----------------------------------------------
  !x-sweeps
  !-----------------------------------------------   
  do k=1,ny
    !compute numerical fluxes
    call computeFlux(nx,k,1)
    !set time-step size
    dt=dx*cour/smax
 
    do j=1,nx   
      !update cell averages of conserved variables
      do l=1,4    
        u2(j,k,l)=u1(j,k,l)+(dt/dx)*(xflux(j-1,l)-xflux(j,l)) 
      end do
    end do
  end do

  call bound()	
  u1=u2
  !-----------------------------------------------
  !y-sweeps
  !-----------------------------------------------  
  do j=1,nx
    !compute numerical fluxes
    call computeFlux(ny,j,2)
    !set time-step size
    dt=dx*cour/smax
    
    do k=1,ny   
      !update cell averages of conserved variables  
      u2(j,k,1)=u1(j,k,1)+(dt/dx)*(yflux(k-1,1)-yflux(k,1)) 
      u2(j,k,2)=u1(j,k,2)+(dt/dx)*(yflux(k-1,3)-yflux(k,3)) 
      u2(j,k,3)=u1(j,k,3)+(dt/dx)*(yflux(k-1,2)-yflux(k,2)) 
      u2(j,k,4)=u1(j,k,4)+(dt/dx)*(yflux(k-1,4)-yflux(k,4))       
    end do
  end do

  call bound()	
  u1=u2
  
  if(mod(i,10)==0)then

  if(i<10)then
    write(uniti,'(I1.1)') i
  else if(i>=10 .and. i<100)then
    write(uniti,'(I2.2)') i
  else if(i>=100 .and. i<1000)then
    write(uniti,'(I3.3)') i
  end if
  filename=trim('Output/t=')//trim(uniti)//trim('.txt')
  open(unit=i+1000,file=filename) 
  call fileOutput(i+1000)
  close(unit=i+1000)

  end if
  call horcut()
  call vercut()
  

end do 

close(unit=10)
close(unit=12)
close(unit=13)

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start

end program EulerGodunov2D
