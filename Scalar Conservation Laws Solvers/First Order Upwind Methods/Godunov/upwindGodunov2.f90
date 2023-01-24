!Numerically solves Burger's equation using 
!Godunov's first order upwind scheme
!Flux function f(u)=0.5*u^2 =>f'(u)=u =>sonic points: u*=0
!Wave Speed, a=u

program upwindGodunov2
implicit  none

integer,parameter::option=1 !1:square wave 2:sinusoide 3:gaussian
integer,parameter::switch=1 !set to 0 for FTBS
integer,parameter::nt=400!time steps
integer,parameter::nx=600 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2),uex(-1:nx+2)
real*8::xmin,xmax,dt,dx,f2,f1,us1,us2
integer::i,j


!initialization
open(unit=10,file='output_godunov2.txt')
!open(unit=11,file='output_exact.txt')
xmin=-1.0
xmax=1.0
dx=(xmax-xmin)/(nx*1.)
dt=0.8*dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  !uex(j)=u1(j)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+j*dx,u1(j)
	!write(11,*) xmin+j*dx,u1(j)
  end if
end do  

!compute solution
do i=1,nt-1
  do j=1,nx
    !check for sonic points
	!print*,'sgn u_',j,',sgn u_',j+1,'=',sign(1.0_8,u1(j)),sign(1.0_8,u1(j+1))
	if( sign(1.0_8,u1(j))/=sign(1.0_8,u1(j+1)) ) then 
	  us2=0.
	  !print*,'Sonic point detected between x=',xmin+dx*j,'..',xmin+dx*(j+1)
	else  
	  us2=u1(j)
	end if
	if( sign(1.0_8,u1(j-1))/=sign(1.0_8,u1(j)) ) then 
	  us1=0.
	  !print*,'Sonic point detected between x=',xmin+dx*(j-1),'..',xmin+dx*j
	else  
	  us1=u1(j)
	end if  
	!print*,'sgn(u_j-1),sgn(u_j),sgn(u_j+1)',sign(1,int(u1(j-1))),&
    !       sign(1,int(u1(j))),sign(1,int(u1(j+1)))	
    if(u1(j)<u1(j+1))then
      f2=min(f(u1(j)),f(u1(j+1)),f(us2))
	else if(u1(j)>u1(j+1))then
      f2=max(f(u1(j)),f(u1(j+1)),f(us2))
	else if(u1(j)==u1(j+1))then
	  f2=f(u1(j))
    end if
    if(u1(j-1)<u1(j))then
      f1=min(f(u1(j-1)),f(u1(j)),f(us1))
	else if(u1(j-1)>u1(j))then
      f1=max(f(u1(j-1)),f(u1(j)),f(us1))
	else if(u1(j-1)==u1(j))then
	  f1=f(u1(j))
    end if
	
    if(switch==1)then	
	  u2(j)=u1(j)-(dt/dx)*(f2-f1)
	else if(switch==0)then
	  u2(j)=u1(j)-(dt/dx)*(f(u1(j))-f(u1(j-1)))
	end if
	
	!print*,'u_j-1, u_j, u_j+1=',u1(j-1),u1(j),u1(j+1),', us1,us2=',us1,us2
	!print*,'f2-f1=',f2-f1
	if(j>=1 .and. j<=nx)then
	  write(10,*) xmin+j*dx,u2(j)
	  !write(11,*) xmin+j*dx,uex(j)
	end if
  end do
  !enforce peridic boundary conditions
  u2(-1)=u2(nx-1)
  u2(0)=u2(nx)
  u2(nx+1)=u2(1)
  u2(nx+2)=u2(2)
  !uex(-1)=uex(nx-1)
  !uex(0)=uex(nx)
  !uex(nx+1)=uex(1)
  !uex(nx+2)=uex(2)
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
   if(abs(x)<1./3.) then
    ux=1.
   else if(abs(x)>1./3. .and. abs(x)<=1.)then
    ux=-1.
   end if
  else if(option==2)then
  !sinusoidal initial profile
   ux=-sin(PI*x)
  else if(option==3)then
  !Gaussian Initial Profile
   ux=(1./sqrt(0.1))*exp(-(x/0.1)**2)
  end if
end function u0

function f(u) result(fx)
  real*8,intent(in)::u
  real*8::fx
  fx=0.5*u*u

end function f


end program upwindGodunov2