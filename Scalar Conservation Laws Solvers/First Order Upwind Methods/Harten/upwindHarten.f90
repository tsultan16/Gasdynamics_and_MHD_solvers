!Numerically solves Linear Advection equation/Burger's equation using 
!Hartens's first order upwind scheme
!Flux function,f(u)=u (Liear advection eqn), 1/2 u^2(Burger's equation)
!Wave Speed, a=u

program upwindHarten
implicit  none

integer,parameter::option=3 !1:square wave 2:sinusoide 3:square wave2
integer,parameter::option2=2 !1:f=u, 2:f=0.5*u^2
integer,parameter::nt=400!time steps
integer,parameter::nx=600 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2)
real*8::xmin,xmax,dt,dx
real*8::f2,f1,us1,us2,del2,del1,del0,a2,a1,eps1,eps2
integer::i,j


!initialization
open(unit=10,file='output_harten.txt')
del0=1.0 !Note: If del0=0.5, then solution almost identical to Godunov upwind
xmin=-1.0
xmax=1.0
dx=(xmax-xmin)/(nx*1.)
dt=0.8*dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+j*dx,u1(j)
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
        !del^n_i+1/2
	del2=del0*(a(u1(j+1))-a(u1(j)))
	!del^n_i-1/2
	del1=del0*(a(u1(j))-a(u1(j-1)))
	!eps^n_i+1/2
	if(abs(a2)<del2)then
	  eps2=(a2*a2+del2*del2)/(2.*del2)
    else if(abs(a2)>del2)then
	  eps2=abs(a2)
	end if
	!eps^n_i-1/2
	if(abs(a1)<del1)then
	  eps1=(a1*a1+del1*del1)/(2.*del1)
    else if(abs(a1)>del1)then
	  eps1=abs(a1)
	end if
	!time-integral averaged fluxes 
	!f^n_1+1/2
	f2=0.5*(f(u1(j+1))+f(u1(j))-eps2*(u1(j+1)-u1(j)))
    !f^n_1-1/2
	f1=0.5*(f(u1(j))+f(u1(j-1))-eps1*(u1(j)-u1(j-1)))
	!update u in cell j
	u2(j)=u1(j)-(dt/dx)*(f2-f1)

	if(j>=1 .and. j<=nx)then
	  write(10,*) xmin+j*dx,u2(j)
	end if
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
   ax=1
  else if(option2==2)then
   ax=u
  end if
end function a


end program upwindHarten
