!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_entropyfixed_mod
implicit none

integer,parameter::debug=0
integer,parameter::shockSwitch=1 !0:off, 1:on
integer,parameter ::nt=100
integer,parameter ::nx=400 
real*8,parameter::del0=0.2
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL,HL !left state
real*8::rhoR,uR,pR,aR,HR !right state 
real*8::delrho,delp,delu,smax

real*8::gt(-1:nx+2,3),fc(-1:nx+2,3),gamma(-1:nx+2,3),sigma(-1:nx+2,3)
real*8::f(-1:nx+2,3),g(-1:nx+2,3),fl(3),u(3)

real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,lambda(-1:nx+2,3)
real*8::K1(-1:nx+2,3),K2(-1:nx+2,3),K3(-1:nx+2,3),alpha(-1:nx+2,3)	
real*8::u1(-1:nx+2,3),u2(-1:nx+2,3)

real*8 ::temp1,temp2,temp3,xc,xt,theta,w


contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_euler_harten_FC.txt')
open(unit=12,file='dt.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
cour=0.9
dx=(xmax-xmin)/nx
dt=0.

gam=1.4
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

smax=0.

read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR

do j=-1,nx+2
  if(j<nx/2) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=0.5*rhoL*uL*uL+pL*gamul
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=0.5*rhoR*uR*uR+pR*gamul
  end if
end do

print*,'Left State (rho,u,p) =',rhoL,uL,pL
print*,'Right State (rho,u,p) =',rhoR,uR,pR

end subroutine init

subroutine RoeSolver(jj)

!Input Variables
integer::jj

delrho=rhoR-rhoL
delu=uR-uL
delp=pR-pL

if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

!-------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^u
!-------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*ut*ut))

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(jj,1)=ut-at
lambda(jj,2)=ut
lambda(jj,3)=ut+at

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(jj,1)=1._8
K1(jj,2)=ut-at
K1(jj,3)=Ht-ut*at

K2(jj,1)=1
K2(jj,2)=ut
K2(jj,3)=0.5*ut*ut

K3(jj,1)=1._8
K3(jj,2)=ut+at
K3(jj,3)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(jj,1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(jj,2)=delrho-delp/(at*at)
alpha(jj,3)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
!smax=abs(ut)+abs(at)
smax=max(smax,abs(lambda(jj,1)),abs(lambda(jj,2)),abs(lambda(jj,3)))


end subroutine RoeSolver

subroutine computeFlux()

do j=-1,nx+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,u1(j,1))
uL=u1(j,2)/u1(j,1)
pL=max(premin,(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1)))
aL=sqrt(gam*pL/rhoL)
HL=0.5*uL*uL+gamul*aL*aL

rhoR=max(densmin,u1(j+1,1))
uR=u1(j+1,2)/u1(j+1,1)
pR=max(premin,(gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1)))
aR=sqrt(gam*pR/rhoR)
HR=0.5*uR*uR+gamul*aR*aR

call RoeSolver(j)

end do

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

!------------------------------------------------
!~g_i+1/2
!------------------------------------------------
do j=-1,nx+1
 do k=1,3
   gt(j,k)= ( 0.5*Q(lambda(j,k))-0.5*(dt/dx)*(lambda(j,k)*lambda(j,k)) )*alpha(j,k) 
   if(shockSwitch==1)then
     sigma(j,k)=0.5*(1.-Q(lambda(j,k)))
   end if	 
 end do  
end do

do j=0,nx+1
  !-----------------------------------------------
  !f(c)_j
  !-----------------------------------------------
   do k=1,3
     fc(j,k)=minmod(gt(j,k),gt(j-1,k))
   end do  
   
   if(shockSwitch==1)then
    do k=1,3
     if(abs(alpha(j,k))+abs(alpha(j-1,k)) .ne. 0._8)then
      theta=(abs(alpha(j,k)-alpha(j-1,k)))/(abs(alpha(j,k))+abs(alpha(j-1,k)))
     else
      theta=0.
     end if 
     w=minmod(alpha(j,k)*sigma(j,k),alpha(j-1,k)*sigma(j-1,k))
     
     !Modify the entropy wave only (will result in steeper contact discontinuity)
     if(k==2)then	 
       fc(j,k)=fc(j,k)+w*theta
     end if

    end do
   end if
   
 !-----------------------------------------------
 !f_j
 !-----------------------------------------------
 u(1)=u1(j,1)
 u(2)=u1(j,2)/u1(j,1)
 u(3)=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
 call flux(u)
   do k=1,3
     f(j,k)=fl(k)
   end do
end do

!------------------------------------------------
!gamma_i+1/2
!------------------------------------------------
do j=0,nx  
 do k=1,3
   if(alpha(j,k)>=1.d-30)then
     gamma(j,k)=(fc(j+1,k)-fc(j,k))/alpha(j,k)
   else
     gamma(j,k)=0.
   end if	  
 end do
  
!------------------------------------------------
!Modified Numerical Flux: ^g_i+1/2
!------------------------------------------------
 do k=1,3
!  g(j,k)=0.5*(f(j,k)+f(j+1,k)&
!        +K1(j,k)*(fc(j,1)+fc(j+1,1)-Q(lambda(j,1)+gamma(j,1))*alpha(j,1))&
!        +K2(j,k)*(fc(j,2)+fc(j+1,2)-Q(lambda(j,2)+gamma(j,2))*alpha(j,2))&
!	+K3(j,k)*(fc(j,3)+fc(j+1,3)-Q(lambda(j,3)+gamma(j,3))*alpha(j,3)) )

  g(j,k)=0.5*(f(j,k)+f(j+1,k))-0.5*&
        ( K1(j,k)*( Q(lambda(j,1)+gamma(j,1))*alpha(j,1)-fc(j,1)-fc(j+1,1) ) &
         +K2(j,k)*( Q(lambda(j,2)+gamma(j,2))*alpha(j,2)-fc(j,2)-fc(j+1,2) ) &
	 +K3(j,k)*( Q(lambda(j,3)+gamma(j,3))*alpha(j,3)-fc(j,3)-fc(j+1,3) ) ) 


 end do
end do
 
end subroutine computeFlux

subroutine flux(u)

!Input Variables    
real*8::u(3)

fl(1)=u(1)*u(2)
fl(2)=u(1)*u(2)*u(2)+u(3)
fl(3)=0.5*u(1)*u(2)*u(2)*u(2)+gam*gamul*u(2)*u(3)
  
end subroutine flux
 
function Q(x) result(qx)

!Input Variables  
real*8,intent(in)::x
!Output Variables 
real*8::qx
 
if(abs(x)<del0)then
  qx=(x*x+del0*del0)/(2.*del0)
else if(abs(x)>=del0)then
  qx=abs(x)
end if

end function Q
 
function minmod(x,y) result(z)

!Input Variables
real*8,intent(in)::x,y

!Output Variables
real*8::z
  
z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod

function limiter(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z,b

  b=1. !minmod limiter for b=1

  if(y>0._8)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter

subroutine bound()
  
do k=1,3
 u2(0,k)=u2(1,k)	
 u2(-1,k)=u2(1,k)
 u2(nx+1,k)=u2(nx,k)
 u2(nx+2,k)=u2(nx,k)
end do

end subroutine bound

end module RoeSolver_entropyfixed_mod
