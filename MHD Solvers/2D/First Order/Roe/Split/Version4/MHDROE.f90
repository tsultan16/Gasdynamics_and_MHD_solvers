!Roe Upwind (1st Order) Method for 2d Ideal MHD Equations 
program MHD_Roe

use vars_mod
use init_mod
use RoeSolverMHD_mod

implicit none

real*8::flux(-1:nx+2,7),u2(-1:nx+2,-1:ny+2,8),smax
integer::i,j,k,l,m,n,alt_flag
character(len=6)::uniti
real::start,finish

call cpu_time(start)
call init(u2)

!Reset smax
smax=0._8

alt_flag=1

!Compute Solution 
do i=1,nt

  call timeStep(u2,smax) 
  print*,'timestep,dt',i,dt

  if(alt_flag==0)then
    call x_sweep()
    call y_sweep()
    alt_flag=1
  else if(alt_flag==1)then
    call y_sweep()
    call x_sweep()
    alt_flag=0
  end if
   
  if(outOption==1)then
   call horcut(u2)
   call vercut(u2)
   call diagcut(u2)
  end if

  if(outOption==1)then
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
  call fileOutput(i+1000,u2)
  close(unit=i+1000)

  end if
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
  
do n=1,ny
 !compute numerical fluxes
 call computeFlux(u2,flux,nx,n,1)
 
 do m=1,nx   
   !update cell averages of conserved variables
   u2(m,n,1)=u2(m,n,1)+(dt/dx)*(flux(m-1,1)-flux(m,1)) 
   u2(m,n,2)=u2(m,n,2)+(dt/dx)*(flux(m-1,2)-flux(m,2))
   u2(m,n,3)=u2(m,n,3)+(dt/dx)*(flux(m-1,3)-flux(m,3))
   u2(m,n,4)=u2(m,n,4)+(dt/dx)*(flux(m-1,4)-flux(m,4))
   u2(m,n,5)=u2(m,n,5)
   u2(m,n,6)=u2(m,n,6)+(dt/dx)*(flux(m-1,5)-flux(m,5))
   u2(m,n,7)=u2(m,n,7)+(dt/dx)*(flux(m-1,6)-flux(m,6))
   u2(m,n,8)=u2(m,n,8)+(dt/dx)*(flux(m-1,7)-flux(m,7))   
 end do
end do

call bound()	

end subroutine x_sweep 

!-----------------------------------------------
!y-sweeps
!-----------------------------------------------
subroutine y_sweep()

do m=1,nx
 !compute numerical fluxes

 call computeFlux(u2,flux,ny,m,2)
  
 do n=1,ny   
   !update cell averages of conserved variables    
   u2(m,n,1)=u2(m,n,1)+(dt/dy)*(flux(n-1,1)-flux(n,1)) 
   u2(m,n,2)=u2(m,n,2)+(dt/dy)*(flux(n-1,3)-flux(n,3))
   u2(m,n,3)=u2(m,n,3)+(dt/dy)*(flux(n-1,2)-flux(n,2))
   u2(m,n,4)=u2(m,n,4)+(dt/dy)*(flux(n-1,4)-flux(n,4))
   u2(m,n,5)=u2(m,n,5)+(dt/dy)*(flux(n-1,5)-flux(n,5))
   u2(m,n,6)=u2(m,n,6)
   u2(m,n,7)=u2(m,n,7)+(dt/dy)*(flux(n-1,6)-flux(n,6))
   u2(m,n,8)=u2(m,n,8)+(dt/dy)*(flux(n-1,7)-flux(n,7))      
 end do
end do

call bound()	

end subroutine y_sweep


subroutine bound()
!Local Variables
integer::jj,kk,ll

if (boundaryType==1)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
end do

else if(boundaryType==2)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do

else if(boundaryType==3)then

do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
  u2(jj,0,3)=-u2(jj,1,3)	
  u2(jj,-1,3)=-u2(jj,1,3)
  u2(jj,ny+1,3)=-u2(jj,ny,3)
  u2(jj,ny+2,3)=-u2(jj,ny,3)
  u2(jj,0,6)=-u2(jj,1,6)	
  u2(jj,-1,6)=-u2(jj,1,6)
  u2(jj,ny+1,6)=-u2(jj,ny,6)
  u2(jj,ny+2,6)=-u2(jj,ny,6)

end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
  u2(0,kk,2)=-u2(1,kk,2)	
  u2(-1,kk,2)=-u2(1,kk,2)
  u2(nx+1,kk,2)=-u2(nx,kk,2)
  u2(nx+2,kk,2)=-u2(nx,kk,2)
  u2(0,kk,5)=-u2(1,kk,5)	
  u2(-1,kk,5)=-u2(1,kk,5)
  u2(nx+1,kk,5)=-u2(nx,kk,5)
  u2(nx+2,kk,5)=-u2(nx,kk,5)
end do

else if(boundaryType==4)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,ny,ll)	
    u2(jj,-1,ll)=u2(jj,ny-1,ll)
    u2(jj,ny+1,ll)=u2(jj,1,ll)
    u2(jj,ny+2,ll)=u2(jj,2,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do

end if

end subroutine bound

subroutine protection()



end subroutine protection

end program MHD_Roe
