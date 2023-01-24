!Roe Upwind (1st Order) Method for 2d Ideal MHD Equations 


program MHD_Roe
use RoeSolver2D_mod
implicit none

integer::i,j,k,l,alt_flag
character(len=6)::uniti
real::start,finish

integer,parameter::tSkip=100

call cpu_time(start)
call init()


!Compute Solution 
do i=1,nt


  call timeStep() 
  print*,'timestep,dt',i,dt

  call x_sweep()
  call y_sweep()

  do j=1,nx
   do k=1,ny
    !update cell averages of conserved variables
    u2(j,k,1)=u1(j,k,1)+(dt/dx)*(flux_x(j-1,k,1)-flux_x(j,k,1)) &
              +(dt/dy)*(flux_y(j,k-1,1)-flux_y(j,k,1)) 
    u2(j,k,2)=u1(j,k,2)+(dt/dx)*(flux_x(j-1,k,2)-flux_x(j,k,2)) &
              +(dt/dy)*(flux_y(j,k-1,3)-flux_y(j,k,3)) 
    u2(j,k,3)=u1(j,k,3)+(dt/dx)*(flux_x(j-1,k,3)-flux_x(j,k,3)) &
              +(dt/dy)*(flux_y(j,k-1,2)-flux_y(j,k,2)) 
    u2(j,k,4)=u1(j,k,4)+(dt/dx)*(flux_x(j-1,k,4)-flux_x(j,k,4)) &
              +(dt/dy)*(flux_y(j,k-1,4)-flux_y(j,k,4)) 
    u2(j,k,5)=u1(j,k,5)+(dt/dy)*(flux_y(j,k-1,5)-flux_y(j,k,5))   
    u2(j,k,6)=u1(j,k,6)+(dt/dx)*(flux_x(j-1,k,5)-flux_x(j,k,5))
    u2(j,k,7)=u1(j,k,7)+(dt/dx)*(flux_x(j-1,k,6)-flux_x(j,k,6)) &
              +(dt/dy)*(flux_y(j,k-1,6)-flux_y(j,k,6)) 
    u2(j,k,8)=u1(j,k,8)+(dt/dx)*(flux_x(j-1,k,7)-flux_x(j,k,7)) &
              +(dt/dy)*(flux_y(j,k-1,7)-flux_y(j,k,7))       
   end do
  end do
  
  call bound()	
  u1=u2


  call horcut()
  call vercut()
  call diagcut()
  
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

  print*,'% Complete=',(real(i)/real(nt))*100.


end do 

close(unit=12)
close(unit=13)
close(unit=14)
close(unit=15)
close(unit=16)
close(unit=17)

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start


contains

!-----------------------------------------------
!x-sweeps
!----------------------------------------------- 
subroutine x_sweep() 
  
do k=1,ny
 !compute numerical fluxes
 call computeFlux(nx,k,1)
 
 !do j=1,nx   
   !update cell averages of conserved variables
  ! u2(j,k,1)=u1(j,k,1)+(dt/dx)*(flux(j-1,1)-flux(j,1)) 
  ! u2(j,k,2)=u1(j,k,2)+(dt/dx)*(flux(j-1,2)-flux(j,2))
  ! u2(j,k,3)=u1(j,k,3)+(dt/dx)*(flux(j-1,3)-flux(j,3))
  ! u2(j,k,4)=u1(j,k,4)+(dt/dx)*(flux(j-1,4)-flux(j,4))
  ! u2(j,k,5)=u1(j,k,5)
  ! u2(j,k,6)=u1(j,k,6)+(dt/dx)*(flux(j-1,5)-flux(j,5))
  ! u2(j,k,7)=u1(j,k,7)+(dt/dx)*(flux(j-1,6)-flux(j,6))
  ! u2(j,k,8)=u1(j,k,8)+(dt/dx)*(flux(j-1,7)-flux(j,7))   
 !end do
end do

!call protection()
!call bound()	
!u1=u2

end subroutine x_sweep 

!-----------------------------------------------
!y-sweeps
!-----------------------------------------------
subroutine y_sweep()

do j=1,nx
 !compute numerical fluxes

 call computeFlux(ny,j,2)
  
! do k=1,ny   
!   !update cell averages of conserved variables    
!   u2(j,k,1)=u1(j,k,1)+(dt/dy)*(flux(k-1,1)-flux(k,1)) 
!   u2(j,k,2)=u1(j,k,2)+(dt/dy)*(flux(k-1,3)-flux(k,3))
!   u2(j,k,3)=u1(j,k,3)+(dt/dy)*(flux(k-1,2)-flux(k,2))
!   u2(j,k,4)=u1(j,k,4)+(dt/dy)*(flux(k-1,4)-flux(k,4))
!   u2(j,k,5)=u1(j,k,5)+(dt/dy)*(flux(k-1,5)-flux(k,5))
!   u2(j,k,6)=u1(j,k,6)
!   u2(j,k,7)=u1(j,k,7)+(dt/dy)*(flux(k-1,6)-flux(k,6))
!   u2(j,k,8)=u1(j,k,8)+(dt/dy)*(flux(k-1,7)-flux(k,7))      
! end do
end do

!call protection()
!call bound()	
!u1=u2

end subroutine y_sweep

end program MHD_Roe
