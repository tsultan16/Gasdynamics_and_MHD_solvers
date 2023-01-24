!Godunovs method (1st order) fopr Eulers Equation in 2D
!with dimensional splitting

program EulerMH2D
use RoeSolver2D_mod
implicit none

integer::i,j,k,l,alt_flag
character(len=6)::uniti
real::start,finish

integer,parameter::tSkip=4

alt_flag=0

call cpu_time(start)
call init()

!Compute Solution 
do i=1,nt
  print*,'% Complete=',(real(i)/real(nt))*100.
  
  if(alt_flag==0)then

  !-----------------------------------------------
  !x-sweeps
  !-----------------------------------------------   
  do k=1,ny
    !compute numerical fluxes
    call computeFlux(nx,k,1)
    !set time-step size
    dtx=dx*cour/smax
 
    do j=1,nx   
      !update cell averages of conserved variables
      do l=1,5   
        u2(j,k,l)=u1(j,k,l)+(dtx/dx)*(flux(j-1,l)-flux(j,l)) 
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
    dty=dy*cour/smax
    
    do k=1,ny   
      !update cell averages of conserved variables  
      u2(j,k,1)=u1(j,k,1)+(dty/dy)*(flux(k-1,1)-flux(k,1)) 
      u2(j,k,2)=u1(j,k,2)+(dty/dy)*(flux(k-1,3)-flux(k,3)) 
      u2(j,k,3)=u1(j,k,3)+(dty/dy)*(flux(k-1,2)-flux(k,2)) 
      u2(j,k,4)=u1(j,k,4)+(dty/dy)*(flux(k-1,4)-flux(k,4))       
      u2(j,k,5)=u1(j,k,5)+(dty/dy)*(flux(k-1,5)-flux(k,5))       
    end do
  end do

  call bound()	
  u1=u2

  alt_flag=1

  else if(alt_flag==1)then

  !-----------------------------------------------
  !y-sweeps
  !-----------------------------------------------  
  do j=1,nx
    !compute numerical fluxes
    call computeFlux(ny,j,2)
    !set time-step size
    dty=dy*cour/smax
    
    do k=1,ny   
      !update cell averages of conserved variables  
      u2(j,k,1)=u1(j,k,1)+(dty/dy)*(flux(k-1,1)-flux(k,1)) 
      u2(j,k,2)=u1(j,k,2)+(dty/dy)*(flux(k-1,3)-flux(k,3)) 
      u2(j,k,3)=u1(j,k,3)+(dty/dy)*(flux(k-1,2)-flux(k,2)) 
      u2(j,k,4)=u1(j,k,4)+(dty/dy)*(flux(k-1,4)-flux(k,4))
      u2(j,k,5)=u1(j,k,5)+(dty/dy)*(flux(k-1,5)-flux(k,5))              
    end do
  end do

  call bound()	
  u1=u2

  !-----------------------------------------------
  !x-sweeps
  !-----------------------------------------------   
  do k=1,ny
    !compute numerical fluxes
    call computeFlux(nx,k,1)
    !set time-step size
    dtx=dx*cour/smax
 
    do j=1,nx   
      !update cell averages of conserved variables
      do l=1,5    
        u2(j,k,l)=u1(j,k,l)+(dtx/dx)*(flux(j-1,l)-flux(j,l)) 
      end do
    end do
  end do

  call bound()	
  u1=u2

  alt_flag=0

  end if  

  !print*,'dtx,dty=',dtx,dty
  
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

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start

end program EulerMH2D
