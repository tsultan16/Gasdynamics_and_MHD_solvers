!Roe's First Order Upwind Method for Euler's eqautions

program EulerRoe
use RiemannRoe_mod
implicit none

integer,parameter::tSteps=100
integer,parameter::nx=512

real*8::dt,dx
real*8::q1(-1:nx+2,3),q2(-1:nx+2,3) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::f1(3),f2(3)
real*8::xmin,xmax,cour,x,pressure
integer::i,j,k

open(unit=10,file='input.txt')
open(unit=11,file='output_euler_roe.txt')

!initialization
xmin=-10.0
xmax=10.0
cour=0.4 !initial CFL number 
gamil=1./(gam-1.)
read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR
dx=(xmax-xmin)/(nx*1.)
dt=cour*dx/374.17

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

   !solve Riemann Problem at both edge of jth cell
    rhoL=q1(j,1)
	uL=q1(j,2)/q1(j,1)
	pL=(gam-1.)*(q1(j,3)-0.5*q1(j,2)*q1(j,2)/q1(j,1))
	rhoR=q1(j+1,1)
	uR=q1(j+1,2)/q1(j+1,1)
	pR=(gam-1.)*(q1(j+1,3)-0.5*q1(j+1,2)*q1(j+1,2)/q1(j+1,1))
	call computeFlux()
    f2=f !flux at right edge

	!print*,'Cell#',j
	!print*,'left state_1=',rhoL,uL,pL
    !print*,'right state_1=',rhoR,uR,pR
	
    rhoL=q1(j-1,1)
	uL=q1(j-1,2)/q1(j-1,1)
	pL=(gam-1.)*(q1(j-1,3)-0.5*q1(j-1,2)*q1(j-1,2)/q1(j-1,1))
	rhoR=q1(j,1)
	uR=q1(j,2)/q1(j,1)
	pR=(gam-1.)*(q1(j,3)-0.5*q1(j,2)*q1(j,2)/q1(j,1))
	call computeFlux()
    f1=f!flux at left edge 
	
	!print*,'left state_2=',rhoL,uL,pL
    !print*,'right state_2=',rhoR,uR,pR
	!print*,'fL=',fL
	!print*,'fR=',fR
	
	!evolve fluid state
	do k=1,3
	 q2(j,k)=q1(j,k)-(dt/dx)*(f2(k)-f1(k))
	end do
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


end program EulerRoe