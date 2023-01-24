!Godunovs method (1st order) fopr Eulers Equation in 2D
!with dimensional splitting

program EulerMH2D
use RoeSolver2D_mod
implicit none

integer::i,j,k,l
character(len=6)::uniti
real::start,finish


call cpu_time(start)
call init()

!Compute Solution 
do i=1,nt
  print*,'% Complete=',(real(i)/real(nt))*100.
  
  !compute numerical fluxes

  call computeFlux()
  print*,'dt=',dt
  !update cell averages of conserved variables
  do j=1,nx 
   do k=1,ny   
     u2(j,k,1)=u1(j,k,1)+(dt/dx)*(f(j-1,k,1)-f(j,k,1)) &
               +(dt/dy)*(g(j,k-1,1)-g(j,k,1))         
     u2(j,k,2)=u1(j,k,2)+(dt/dx)*(f(j-1,k,2)-f(j,k,2)) &
               +(dt/dy)*(g(j,k-1,3)-g(j,k,3))        
     u2(j,k,3)=u1(j,k,3)+(dt/dx)*(f(j-1,k,3)-f(j,k,3)) &
               +(dt/dy)*(g(j,k-1,2)-g(j,k,2))        
     u2(j,k,4)=u1(j,k,4)+(dt/dx)*(f(j-1,k,4)-f(j,k,4)) &
               +(dt/dy)*(g(j,k-1,4)-g(j,k,4))        
     u2(j,k,5)=u1(j,k,5)+(dt/dx)*(f(j-1,k,5)-f(j,k,5)) &
               +(dt/dy)*(g(j,k-1,5)-g(j,k,5))        
   end do
  end do

  call bound()	
  u1=u2



  !print*,'dt',dt
  
  if(mod(i,tSkip)==0)then

  if(i<10)then
    write(uniti,'(I1.1)') i
  else if(i>=10 .and. i<100)then
    write(uniti,'(I2.2)') i
  else if(i>=100 .and. i<1000)then
    write(uniti,'(I3.3)') i
  else if(i>=1000 .and. i<10000)then
    write(uniti,'(I4.3)') i
  else if(i>=10000 .and. i<100000)then
    write(uniti,'(I5.3)') i
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
close(unit=22)

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start

end program EulerMH2D
