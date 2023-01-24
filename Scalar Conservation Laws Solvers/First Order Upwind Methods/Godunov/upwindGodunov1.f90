!Numerically solves linear advection equation using 
!Godunov's first order upwind scheme
!Flux function f(u)=u => f', no sonic points
!Wave Speed, a=1

program upwindGodunov1
implicit  none

integer,parameter::switch=1 !use this to check if godunov method gives equivalent result as FTBS
integer,parameter::option=1 ! 1:square wave 2:sinusoid 3:gaussian
integer,parameter::nt=100 !time steps
integer,parameter::nx=40 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2),uex(-1:nx+2)
real*8::xmin,xmax,dt,dx,f2,f1,lambda
integer::i,j

!initialization
open(unit=10,file='output_godunov1.txt')
open(unit=11,file='output_exact.txt')
xmin=-1.0
xmax=1.0
lambda=0.8
dx=(xmax-xmin)/(nx*1.)
dt=lambda*dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  uex(j)=u1(j)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+j*dx,u1(j)
	write(11,*) xmin+j*dx,u1(j)
  end if
end do  

!compute solution
do i=1,nt-1
  do j=1,nx 
    if(u1(j)<u1(j+1))then
      f2=min(u1(j),u1(j+1))
	else if(u1(j)>u1(j+1))then
      f2=max(u1(j),u1(j+1))
	else if(u1(j)==u1(j+1))then
      f2=u1(j)	  
    end if
    if(u1(j-1)<u1(j))then
      f1=min(u1(j-1),u1(j))
	else if(u1(j-1)>u1(j))then
      f1=max(u1(j-1),u1(j))
	else if(u1(j-1)==u1(j))then
      f1=u1(j)
    end if	
	if(switch==1)then
	  u2(j)=u1(j)-(dt/dx)*(f2-f1)
	else if(switch==0)then
	  u2(j)=u1(j)-(dt/dx)*(u1(j)-u1(j-1))
	end if
	!print*,'f2-f1, godunov=',f2-f1,', FTBS=',u1(j)-u1(j-1)

	
	uex(j)=u0(xmin+abs(mod(j*dx-dt*i,(xmax-xmin))))
	!uex(j)=u0(xmin+j*dx-dt*i)
	if(j>=1 .and. j<=nx)then
	  write(10,*) xmin+j*dx,u2(j)
	  write(11,*) xmin+j*dx,uex(j)
	end if
  end do
  !enforce peridic boundary conditions
  u2(-1)=u2(nx-1)
  u2(0)=u2(nx)
  u2(nx+1)=u2(1)
  u2(nx+2)=u2(2)
  uex(-1)=uex(nx-1)
  uex(0)=uex(nx)
  uex(nx+1)=uex(1)
  uex(nx+2)=uex(2)
  !!!!!!
  u1=u2       
end do  


print*,'Done...'  
close(unit=10)    
  
contains 

function u0(x) result(ux)
  real*8,intent(in)::x
  real*8::ux
  
  !square wave initial profile
  if(option==1)then
    if(x<=1./3. .and. x>=-1./3.) then
      ux=1.
    else
      ux=0.
    end if
  else if(option==2)then 
    !sinusoidal initial profile
    ux=-sin(PI*x)
  else if(option==3)then
    !Gaussian Initial Profile
    ux=(1./sqrt(0.1))*exp(-(x/0.1)**2)
  end if
end function u0



end program upwindGodunov1