!Van Leer MUSCL Scheme (2nd order) for Ieal MHD in 2D
!With CT for maintaining div(B)=0 on the grid

program MHD_MH
use RoeSolver2D_mod
implicit none

integer::i,j,k,l,N,N2,alt_flag
character(len=6)::uniti
real::start,finish
real*8::t,yx(nx),yx2(nx),xt0,yt0,yt2_0!,xt(0:nx),yt(0:ny),xt2(0:nx),yt2(0:ny)
t=0._8

open(unit=99,file='Output/energy.txt')


call cpu_time(start)
call init()


yx(1)=ymin+0.55*dy*(ny-0.5)
yx2(1)=ymin+0.4*dy*(ny-0.5)

!yt(1)=ymin+0.6*dy*(ny-0.5)
!xt(0)=xmin+0.5*dx
!xt2(0)=xt(0)
!yt2(1)=ymin+0.4*dy*(ny-0.5)
yt0=ymin+0.3*dy*(ny-0.5)
xt0=xmin
yt2_0=ymin+0.4*dy*(ny-0.5)

alt_flag=1

!Compute Solution 
do i=1,nt
  call timeStep() 
  print*,'timestep,dt',i,dt
  t=t+dt
  if(alt_flag==0)then

    call x_sweep()
    call y_sweep()
    alt_flag=1
  else if(alt_flag==1)then
    call y_sweep()
    call x_sweep()
    alt_flag=0
  end if
   
  !Update Cell interface Magnetic Field
  call magUpdate()

  !call horcut()
  !call vercut()
  !call diagcut()
  
  call computeEnergy()


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
  filename1=trim('Output/t=')//trim(uniti)//trim('.txt')
  filename2=trim('Output/b_t=')//trim(uniti)//trim('.txt')
  !filename3=trim('Output/fieldLine1_t=')//trim(uniti)//trim('.txt')

  open(unit=i+1000,file=filename1)
  open(unit=i+1000+7000,file=filename2)
  !open(unit=i+9000,file=filename3)

  call fileOutput(i+1000)
  !call computeFieldLinePar(i+9000)
  !  call computeFieldLinePar2(i+9000)
  close(unit=i+1000)
  close(unit=i+1000+7000)
  !close(unit=i+9000)
  
  else
   if(int(yt0/dy)<=ny .and. int(yt0/dy)>=0)then 
    yt0=yt0+(0.1*dx)*u2(int(xt0/dx),int(yt0/dy),6)
   end if
   !if(int(yt2_0/dy)<=ny .and. int(yt2_0/dy)>=0)then 
   ! yt2_0=yt2_0+(0.1*dx)*u2(int(xt0/dx),int(yt2_0/dy),6)
   !end if
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
  
do k=0,ny+1
 !compute numerical fluxes
 if(fluxType==1)then    
   call computeFluxRoe(nx,k,1)
 else if(fluxType==2)then
   call computeFluxHLL(nx,k,1)
 else if(fluxType==3)then
   call computeFluxHLLI(nx,k,1)
 end if
 do j=1,nx   
   !update cell averages of conserved variables
   u2(j,k,1)=u1(j,k,1)+(dt/dx)*(flux(j-1,1)-flux(j,1)) 
   u2(j,k,2)=u1(j,k,2)+(dt/dx)*(flux(j-1,2)-flux(j,2))
   u2(j,k,3)=u1(j,k,3)+(dt/dx)*(flux(j-1,3)-flux(j,3))
   u2(j,k,4)=u1(j,k,4)+(dt/dx)*(flux(j-1,4)-flux(j,4))
   !u2(j,k,5)=u1(j,k,5)
   !u2(j,k,6)=u1(j,k,6)+(dt/dx)*(flux(j-1,5)-flux(j,5))
   u2(j,k,7)=u1(j,k,7)+(dt/dx)*(flux(j-1,6)-flux(j,6))
   u2(j,k,8)=u1(j,k,8)+(dt/dx)*(flux(j-1,7)-flux(j,7))   
 end do
end do

call protection()
call bound()	
u1=u2

end subroutine x_sweep 

!-----------------------------------------------
!y-sweeps
!-----------------------------------------------
subroutine y_sweep()

do j=0,nx+1
 !compute numerical fluxes
 if(fluxType==1)then    
   call computeFluxRoe(ny,j,2)
 else if(fluxType==2)then    
   call computeFluxHLL(ny,j,2)
 else if(fluxType==3)then    
   call computeFluxHLLI(ny,j,2)
 end if 
 do k=1,ny   
   !update cell averages of conserved variables    
   u2(j,k,1)=u1(j,k,1)+(dt/dy)*(flux(k-1,1)-flux(k,1)) 
   u2(j,k,2)=u1(j,k,2)+(dt/dy)*(flux(k-1,3)-flux(k,3))
   u2(j,k,3)=u1(j,k,3)+(dt/dy)*(flux(k-1,2)-flux(k,2))
   u2(j,k,4)=u1(j,k,4)+(dt/dy)*(flux(k-1,4)-flux(k,4))
   !u2(j,k,5)=u1(j,k,5)+(dt/dy)*(flux(k-1,5)-flux(k,5))
   !u2(j,k,6)=u1(j,k,6)
   u2(j,k,7)=u1(j,k,7)+(dt/dy)*(flux(k-1,6)-flux(k,6))
   u2(j,k,8)=u1(j,k,8)+(dt/dy)*(flux(k-1,7)-flux(k,7))      
 end do
end do

call protection()
call bound()	
u1=u2

end subroutine y_sweep

subroutine magUpdate()

!------------------------------------------
!Cell Corner Electric Field: Ez[i,j]
!------------------------------------------
do j=0,nx
!$omp parallel do shared(Ez,fx,fy)
  do k=0,ny
   Ez(j,k)=0.5*(fy(j,k)+fy(j+1,k)-fx(j,k)-fx(j,k+1))
  end do
 !$omp end parallel do  
end do

!------------------------------------------------------------------------
!Update Cell interface magnetic field normal components: bx[i,j], by[i,j]
!-----------------------------------------------------------------------
do k=1,ny
  !$omp parallel do shared(Ez,bxInt)
  do j=0,nx
   bxInt(j,k)=bxInt(j,k)-(dt/dy)*(Ez(j,k)-Ez(j,k-1))
  end do
  !$omp end parallel do   
end do

do j=1,nx
  !$omp parallel do shared(Ez,byInt)
  do k=0,ny
   byInt(j,k)=byInt(j,k)+(dt/dx)*(Ez(j,k)-Ez(j-1,k))
  end do
   !$omp end parallel do    
end do


!-----------------------------------------------------------------
!Update Cell Center Magnetic Field Values (linear interpolation 
!of cell inteferface values)
!-----------------------------------------------------------------
do j=1,nx
  !$omp parallel do shared(u2,bxInt,byInt)
  do k=1,ny
    u2(j,k,5)=0.5*(bxInt(j-1,k)+bxInt(j,k))
    u2(j,k,6)=0.5*(byInt(j,k-1)+byInt(j,k))
  end do
  !$omp end parallel do
end do

call bound()
u1=u2


end subroutine magUpdate

subroutine computeEnergy()
real*8::eint,ekin,emag

eint=0.
ekin=0.
emag=0.

do j=1,nx 
 do k=1,ny
   ekin=ekin+dx*dy*0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1)
   emag=emag+dx*dy*0.5*(u2(j,k,5)**2.+u2(j,k,6)**2.+u2(j,k,7)**2.)
   eint=eint+dx*dy*( u2(j,k,8)&
   -0.5*(u2(j,k,2)**2.+u2(j,k,3)**2.+u2(j,k,4)**2.)/u2(j,k,1) &
   -0.5*(u2(j,k,5)**2.+u2(j,k,6)**2.+u2(j,k,7)**2.) )
 end do
end do

write(99,*) t,ekin,eint,emag 

end subroutine computeEnergy


subroutine computeFieldLine()

!Local Vars
real*8::x0,y0,y1
integer::k1,k2

y0=yx(1)
y1=yx2(1)

do j=1,nx
!Forward Euler for field line equation update
 x0=xmin+(j+0.5)*dx
 k1=(y0-ymin)/dy
 k2=(y1-ymin)/dy
 if(k1<=ny .and. k2>=1)then
 yx(j)=y0+dx*(u2(j,k1,6)/u2(j,k1,5)) !Field Line 1
 end if
 if(k2<=ny .and. k2>=1)then !Field Line 2
 yx2(j)=y1+dx*(u2(j,k2,6)/u2(j,k2,5))
 end if
 write(88,*) x0,yx(j)
 write(89,*) x0,yx2(j)
 y0=yx(j)
 y1=yx2(j)
end do


end subroutine computeFieldLine


subroutine computeFieldLinePar(fileUnit)

!Local Vars
real*8::x0,y0,x1,y1,xt,yt,dr
integer::N,flag,fileUnit

flag=0
N=0
dr=0.1*dx

y0=yt0
x0=xt0
!x0=xt(0)
!y0=yt(1)
!x1=xt2(0)
!y1=yt2(1)

do while(x0<=xmax .and. N<20000)
!Forward Euler for field line equation update

 if(int(y0/dy)<=ny .and. int(y0/dy)>=0 &
  .and. int(x0/dx)<=nx .and. int(x0/dx)>=0 )then
 xt=x0+dr*u2(int(x0/dx),int(y0/dy),5)
 yt=y0+dr*u2(int(x0/dx),int(y0/dy),6) !Field line 1
 end if


 write(fileUnit,*) xt,yt

 x0=xt
 y0=yt

 N=N+1

 if(flag==0)then
  !yt0=yt
  
  flag=1
 end if

end do
 print*,'N=',N

end subroutine computeFieldLinePar


subroutine computeFieldLinePar2(fileUnit)

!Local Vars
real*8::x0,y0,x1,y1,xt,yt,dr
integer::N,flag,fileUnit

flag=0
N=0
dr=0.1*dx

y0=yt2_0
x0=xt0

do while(x0<=xmax .and. N<20000)
!Forward Euler for field line equation update

 if(int(y0/dy)<=ny .and. int(y0/dy)>=0 &
  .and. int(x0/dx)<=nx .and. int(x0/dx)>=0 )then
 xt=x0+dr*u2(int(x0/dx),int(y0/dy),5)
 yt=y0+dr*u2(int(x0/dx),int(y0/dy),6) !Field line 1
 end if


 write(fileUnit,*) xt,yt

 x0=xt
 y0=yt

 N=N+1

 if(flag==0)then
  yt2_0=yt
  flag=1
 end if

end do
 print*,'N=',N

end subroutine computeFieldLinePar2





end program MHD_MH
