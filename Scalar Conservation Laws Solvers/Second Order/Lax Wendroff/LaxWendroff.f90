!Numerically solves Linear Advection equation/Burger's equation using 
!Lax-Wendroff scheme
!Flux function,f(u)=u (Liear advection eqn), 1/2 u^2(Burger's equation)
!Wave Speed, a=df/du

program LaxWendroff
implicit  none

integer,parameter::option=1 !1:square wave 2:sinusoid 3:square wave2
integer,parameter::option2=2 !1:f=u, 2:f=0.5*u^2
integer,parameter::nt=100!time steps
integer,parameter::nx=500 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2)
real*8::xmin,xmax,dt,dx
real*8::f2,f1,us1,us2,a2,a1,cour,lambda
integer::i,j


!initialization
open(unit=10,file='output_LW.txt')
xmin=-1.0
xmax=1.0
dx=(xmax-xmin)/(nx*1.)
cour=0.8 !CFL number
dt=cour*dx
lambda=dt/dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+(j-0.5)*dx,u1(j)
  end if
end do  

!compute solution
do i=1,nt-1
  do j=1,nx
    !a^n_i+1/2
	if(u1(j+1)/=u1(j))then
	  a2=(f(u1(j+1))-f(u1(j)))/(u1(j+1)-u1(j))
	else if(u1(j+1)==u1(j))then
      a2=a(u1(j))
	end if
	!a^n_i-1/2 
	if(u1(j)/=u1(j-1))then
	  a1=(f(u1(j))-f(u1(j-1)))/(u1(j)-u1(j-1))
	else if(u1(j)==u1(j-1))then
      a1=a(u1(j))	  
	end if
    !time-integral averaged fluxes 
	!f^n_1+1/2
	f2=0.5*( f(u1(j+1))+f(u1(j))-lambda*a2*a2*(u1(j+1)-u1(j)) )
    !f^n_1-1/2
	f1=0.5*( f(u1(j))+f(u1(j-1))-lambda*a1*a1*(u1(j)-u1(j-1)) )
	!update u in cell j
	u2(j)=u1(j)-lambda*(f2-f1)

	write(10,*) xmin+(j-0.5)*dx,u2(j)

  end do
  !enforce peridic boundary conditions
  u2(-1)=u2(nx-1)
  u2(0)=u2(nx)
  u2(nx+1)=u2(1)
  u2(nx+2)=u2(2)
  
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
    ux=0.
   end if
  else if(option==2)then
  !sinusoidal initial profile
   ux=-sin(PI*x)
  else if(option==3)then
  !Square Wave 2
   if(abs(x)<1./3.) then
    ux=1.
   else if(abs(x)>1./3. .and. abs(x)<=1.)then
    ux=-1.
   end if
  end if
end function u0

function f(u) result(fx)
  real*8,intent(in)::u
  real*8::fx
  if(option2==1)then
   fx=u
  else if(option2==2)then
   fx=0.5*u*u
  end if
end function f

function a(u) result(ax)
  real*8,intent(in)::u
  real*8::ax
  if(option2==1)then
   ax=1.
  else if(option2==2)then
   ax=u
  end if
end function a


end program LaxWendroff