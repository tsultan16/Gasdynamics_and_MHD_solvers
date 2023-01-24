!Numerically solves Linear Advection equation/Burger's equation using 
!Hartens's Second Order Flux Corrected TVD scheme
!Flux function,f(u)=u (Liear advection eqn), 1/2 u^2(Burger's equation)
!Wave Speed, a=u

program HartenFC
implicit  none

integer,parameter::option=3 !1:square wave 2:sinusoide 3:square wave2
integer,parameter::option2=2 !1:f=u, 2:f=0.5*u^2
integer,parameter::nt=100!time steps
integer,parameter::nx=500 !number of cells
real*8,parameter::PI=4*atan(1.0)
real*8::u1(-1:nx+2),u2(-1:nx+2)
real*8::xmin,xmax,dt,dx,cour,lambda
real*8::g2,g1,gt2,gt1,fc2,fc1,fc0,del0,a3,a2,a1,a0,gam2,gam1
integer::i,j


!initialization
open(unit=10,file='output_harten_FC.txt')
del0=1.0
cour=0.8 
xmin=-1.0
xmax=1.0
dx=(xmax-xmin)/(nx*1.)
dt=cour*dx
lambda=dt/dx
do j=-1,nx+2
  u1(j)=u0(xmin+j*dx)
  if(j>=1 .and. j<=nx) then
    write(10,*) xmin+j*dx,u1(j)
  end if
end do  

!compute solution
do i=1,nt-1
  do j=1,nx
    
	!a^n_i+3/2
	if(u1(j+2)/=u1(j+1))then
	  a3=(f(u1(j+2))-f(u1(j+1)))/(u1(j+2)-u1(j+1))
	else if(u1(j+2)==u1(j+1))then
      a3=a(u1(j+1))
	end if
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
	!a^n_i-3/2 
	if(u1(j-1)/=u1(j-2))then
	  a0=(f(u1(j-1))-f(u1(j-2)))/(u1(j-1)-u1(j-2))
	else if(u1(j-1)==u1(j-2))then
      a0=a(u1(j-1))	  
	end if
	
	!g_tilde^n_i+3/2 
	gt2=0.5*(u1(j+2)-u1(j+1))*(Q(a3)-lambda*a3*a3)
	!g_tilde^n_i+1/2
	gt1=0.5*(u1(j+1)-u1(j))*(Q(a2)-lambda*a2*a2)
	!f^(c)_i+1
	fc2=minmod(gt2,gt1)
	!g_tilde^n_i+1/2
	gt2=0.5*(u1(j+1)-u1(j))*(Q(a2)-lambda*a2*a2)
	!g_tilde^n_i-1/2
	gt1=0.5*(u1(j)-u1(j-1))*(Q(a1)-lambda*a1*a1)
	!f^(c)_i
	fc1=minmod(gt2,gt1)
	!g_tilde^n_i-1/2
	gt2=0.5*(u1(j)-u1(j-1))*(Q(a1)-lambda*a1*a1)
	!g_tilde^n_i-3/2
	gt1=0.5*(u1(j)-u1(j-1))*(Q(a0)-lambda*a0*a0)
	!f^(c)_i-1
	fc0=minmod(gt2,gt1)
    !gam^n_i+1/2
	gam2=(fc2-fc1)/(u1(j+1)-u1(j))
	!gam^n_i-1/2
	gam1=(fc1-fc0)/(u1(j)-u1(j-1))
	
	!compute numerical fluxes 
	!************************
	
	!g^n_1+1/2 	
	g2=0.5*( f(u1(j+1))+f(u1(j))+fc2+fc1-Q(a2+gam2)*(u1(j+1)-u1(j)) )
    !g^n_1-1/2
	g1=0.5*( f(u1(j))+f(u1(j-1))+fc1+fc0-Q(a1+gam1)*(u1(j)-u1(j-1)) )
	
	!update u in cell j
	u2(j)=u1(j)-lambda*(g2-g1)

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

function Q(x) result(qx)
  real*8,intent(in)::x
  real*8::qx
  
  if(abs(x)<del0)then
    qx=(x*x+del0*del0)/(2.*del0)
  else if(abs(x)>=del0)then
    qx=abs(x)
  end if

end function Q

function minmod(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z
  
  z=sign(1._8,x)*max(0.,min(abs(x),y*sign(1._8,x)))
  
end function minmod


end program HartenFC