!Lax-Wendroff Method for Euler's eqautions
!with Riemann Problem initial conditions

program EulerLW
implicit none

integer,parameter::option1=2  !1:Richtmyer 2:McCormack's
integer,parameter::option2=2  !1: w/o artificial viscosity 2: w/ artificial viscosity
integer,parameter::tSteps=100
integer,parameter::nx=400
real*8,parameter ::gam=5./3. !gam=cp/cv
real*8,parameter ::R=287. !gas constant in SI units

real*8::dt,dx
real*8::q1(-1:nx+2,3),q2(-1:nx+2,3) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::f(3),f1(3),f2(3),f3(3),Qh1(3),Qh2(3),temp(3),qbar1(3),qbar2(3)
real*8::xmin,xmax,cour,lambda,gamil,x,pressure,visc
real*8::eps!artificial viscosity parameter
integer::i,j,k
real*8::rhoL,uL,pL   !left state
real*8::rhoR,uR,pR   !right state

open(unit=10,file='input.txt')
open(unit=11,file='output_euler_LW.txt')

!initialization
eps=0.02  
xmin=-10.0
xmax=10.0
cour=0.6 !initial CFL number 
gamil=1./(gam-1.)
read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR
dx=(xmax-xmin)/(nx*1.)
dt=cour*dx/374.17
lambda=dt/dx

do i=-1,nx+2
  x=xmin+i*dx
  if(x<0) then
   q1(i,1)=rhoL
   q1(i,2)=rhoL*uL
   q1(i,3)=0.5*rhoL*uL*uL+gamil*pL
   pressure=pL
  else if(x>=0)then
   q1(i,1)=rhoR
   q1(i,2)=rhoR*uR
   q1(i,3)=0.5*rhoR*uR*uR+gamil*pR
   pressure=pR
  end if
  if(i>=1 .and. i<=nx) then  
    write(11,*) x,q1(i,1),q1(i,2)/q1(i,1),pressure
  end if    
end do



!compute solution
do i=1,tSteps-1
  do j=1,nx
    
    x=xmin+j*dx
    
	!***********Predictor Step*****************
	! Qh2=q^(n+1/2)_i+1/2 , Qh1=q^(n+1/2)_i-1/2 *
	!******************************************
	temp=(/q1(j+1,1),q1(j+1,2),q1(j+1,3)/)
	call flux(temp)
	f3=f
	temp=(/q1(j,1),q1(j,2),q1(j,3)/)
	call flux(temp)
	f2=f
	temp=(/q1(j-1,1),q1(j-1,2),q1(j-1,3)/)
	call flux(temp)
	f1=f
	if(option1==1)then
     do k=1,3
	  Qh2(k)=0.5*( q1(j+1,k)+q1(j,k)-lambda*(f3(k)-f2(k)) )
	  Qh1(k)=0.5*( q1(j,k)+q1(j-1,k)-lambda*(f2(k)-f1(k)) )
	 end do
	else if(option1==2)then	
	 do k=1,3
	  qbar2(k)=q1(j,k)-lambda*(f3(k)-f2(k))
	  qbar1(k)=q1(j-1,k)-lambda*(f2(k)-f1(k))
	 end do
    end if
	!******************Corrector Step******************************
	! q^(n+1)_i = q^n_i -lambda(f(Q2)-f(Q1))+artificial viscosity *
	!**************************************************************
    !do k=1,3	
	!  visc=eps*(q1(j+1,k)-2.*q1(j,k)+q1(j-1,k))
	!end do
	if(option1==1)then
	 call flux(Qh2)
	 f2=f
	 call flux(Qh1)
	 f1=f
	 if(option2==1)then 
	  do k=1,3
  	   q2(j,k)=q1(j,k)-lambda*(f2(k)-f1(k))
	  end do 
	 else if (option2==2)then
	  do k=1,3
	   q2(j,k)=q1(j,k)-lambda*(f2(k)-f1(k))&
	           +eps*(q1(j+1,k)-2.*q1(j,k)+q1(j-1,k))
	  end do
	 end if
	else if(option1==2)then
	 call flux(qbar2)
	 f2=f
	 call flux(qbar1)
	 f1=f
	 if(option2==1)then
	  do k=1,3
  	   q2(j,k)=0.5*(q1(j,k)+qbar2(k)-lambda*(f2(k)-f1(k)))
	  end do 
	 else if(option2==2)then
      do k=1,3
  	   q2(j,k)=0.5*(q1(j,k)+qbar2(k)-lambda*(f2(k)-f1(k)))&
	           +eps*(q1(j+1,k)-2.*q1(j,k)+q1(j-1,k))
	  end do
     end if	   
	end if
	
	
	pressure=(gam-1.)*(q2(j,3)-0.5*q2(j,2)*q2(j,2)/q2(j,1))
    if(j>=1 .and. j<=nx) then  
     write(11,*) x,q2(j,1),q2(j,2)/q2(j,1),pressure
    end if	


  end do
  !enforce open boundary conditions
  do k=1,3
    q2(0,k)=q2(1,k)	
    q2(-1,k)=q2(1,k)
    q2(nx+1,k)=q2(nx,k)
	q2(nx+2,k)=q2(nx,k)
  end do
  q1=q2
end do

print*,'Done...'

contains

subroutine flux(q)
  real*8::q(3)
  real*8::u(3)  !u=[rho,u,p]
  u(1)=q(1)
  u(2)=q(2)/q(1)
  u(3)=(gam-1.)*(q(3)-0.5*q(2)*q(2)/q(1))
  
  f(1)=u(1)*u(2)
  f(2)=u(1)*u(2)*u(2)+u(3)
  f(3)=u(1)*u(2)*(gam*gamil*u(3)/u(1)+0.5*u(2)*u(2))
 end 

end program EulerLW
