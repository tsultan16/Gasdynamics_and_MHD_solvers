!Numerically solves Linear Advection equation/Burger's equation using 
!Hartens-Osher-Engquist-Chakravarthy second order scheme (ENO)
!Flux function,f(u)=u (Linear advection eqn), 1/2 u^2(Burger's equation)
!Wave Speed, a=u

program ENO
implicit  none

integer,parameter::option=2 !1:square wave 2:sinusoide 3:square wave2
integer,parameter::option2=2 !1:f=u, 2:f=0.5*u^2
integer,parameter::option3=2 !1:godunov average, 2:harten average
integer,parameter::nt=300!time steps
integer,parameter::nx=100 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2)
real*8::a(-1:nx+2),a2(-1:nx+2),S(-1:nx+2),ut(-1:nx+2),f(-1:nx+2)
real*8::eps(-1:nx+2)
real*8::xmin,xmax,dt,dx,del
integer::i,j


!initialization
open(unit=10,file='output_scalar_ENO.txt')
del=1.
xmin=-1.0
xmax=1.0
dx=(xmax-xmin)/(nx*1.)
dt=0.5*dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+j*dx,u1(j)
  end if
end do  

!compute solution
do i=1,nt-1
 
  call computeFlux() 
  do j=1,nx
	!update u in cell j
	u2(j)=u1(j)-(dt/dx)*(f(j)-f(j-1))

	if(j>=1 .and. j<=nx)then
	  write(10,*) xmin+j*dx,u2(j)
	end if
  end do
  
  !enforce periodic boundary conditions
  u2(-1)=u2(nx-1)
  u2(0)=u2(nx)
  u2(nx+1)=u2(1)
  u2(nx+2)=u2(2)
  
  u1=u2       
end do  


print*,'Done...'  
close(unit=10)    
  
contains 

subroutine computeFlux()
  integer ::i
  real*8 ::us
  do i=0,nx+1
    !a_i
	a(i)=aspeed(u1(i))
	!a_i+1/2
	if(u1(i+1)/=u1(i))then
	  a2(i)=(flux(u1(i+1))-flux(u1(i)))/(u1(i+1)-u1(i))
	else
      a2(i)=a(i)
	end if
	!eps_i+1/2
	if(abs(a2(i))<del)then
	  eps(i)=(a2(i)*a2(i)+del*del)/(2.*del)
    else if(abs(a2(i))>del)then
	  eps(i)=abs(a2(i))
	end if
	!S_i
	if(abs(u1(i+1)-u1(i))<abs(u1(i)-u1(i-1)))then
	  S(i)=(u1(i+1)-u1(i))/dx
	else 
	  S(i)=(u1(i)-u1(i-1))/dx
	end if
	!~u_i
	ut(i)=u1(i)+0.5*(1.-(dt/dx)*a(i))*S(i)*dx
  end do

  do i=0,nx
    !^f_i+1/2
    if(option3==1)then
     !Godunov upwind scheme as the averaging operator
	 !check for sonic point
     if(sign(1._8,ut(i))/=sign(1._8,ut(i+1)))then
      us=0.
     else
      us=ut(i)
     end if
   	  
	 if(ut(i)<ut(i+1))then
       f(i)=min(flux(ut(i)),flux(ut(i+1)),flux(us))
	 else  
	   f(i)=max(flux(ut(i)),flux(ut(i+1)),flux(us))
	 end if
	else if(option3==2)then 
	  !Harten's upwind scheme as averaging operator
	  f(i)=0.5*(flux(ut(i))+flux(ut(i+1))-eps(i)*(u1(i+1)-u1(i)))
	end if
	 
  end do
end subroutine computeFlux

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

function flux(u) result(fx)
  real*8,intent(in)::u
  real*8::fx
  if(option2==1)then
   fx=u
  else if(option2==2)then
   fx=0.5*u*u
  end if
end function flux

function aspeed(u) result(ax)
  real*8,intent(in)::u
  real*8::ax
  if(option2==1)then
   ax=1
  else if(option2==2)then
   ax=u
  end if
end function aspeed


end program ENO
