!Roe's Approximate (Linear) Riemann Solver for the equations of Ideal MHD
!Primitive Variable Formulation
module RoeSolverMHD_mod
implicit none

integer,parameter::debug=0
integer,parameter ::nt=300
integer,parameter ::nx=200
real*8,parameter::pi=3.14159265359

real*8::tol=1.d-15
real*8::premin=1.d-15,densmin=1.d-15
integer ::i,j,k,l
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,byL,bzL,aL,FL(7) !left state
real*8::rhoR,uR,vR,wR,pR,byR,bzR,aR,FR(7) !right state 
real*8::smax,cs,cf,cA
real*8::dens,velx,vely,velz,pres,magy,magz
real*8::rhot,ut,vt,wt,pt,byt,bzt,at,lambda(7)
real*8::L1(7),L2(7),L3(7),L4(7),L5(7),L6(7),L7(7)	
real*8::R1(7),R2(7),R3(7),R4(7),R5(7),R6(7),R7(7)	
real*8::u1(-1:nx+2,7),u2(-1:nx+2,7),flux(-1:nx+2,7)
real*8::VLx(7),VRx(7),alpha(7)
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr,bxLR


contains 

subroutine init()
!------------------------------------------------------
!initialization of conserved variables
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_fluid.txt')
open(unit=12,file='output_mag.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
cour=0.8
dx=(xmax-xmin)/nx

gam=5./3.
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

smax=0.

!read(unit=10,fmt=*) rhoL,uL,vL,wL,byL,bzL,pL 
!read(unit=10,fmt=*) rhoR,uR,vR,wR,byR,bzR,pR
rhoL=1.
uL=0.
vL=0.
wL=0.
bxLR=3.
byL=5.
bzL=0.
pL=1.

rhoR=0.1
uR=0.
vR=0.
wR=0.
byR=2.
bzR=0.!/sqrt(4.*pi)
pR=10.


do j=-1,nx+2
  x=xmin+(j-0.5)*dx 
  if(j<nx/2) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=rhoL*vL
    u1(j,4)=rhoL*wL
    u1(j,5)=byL
    u1(j,6)=bzL  
    u1(j,7)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+(bxLR*bxLR+byL*byL+bzL*bzL)/(8.*pi)
    if(j>=0 .and. j<=nx)then
      write(11,*) x,rhoL,uL,vL,wL,pL
      write(12,*) x,byL,bzL
    end if
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=rhoR*vR
    u1(j,4)=rhoR*wR
    u1(j,5)=byR
    u1(j,6)=bzR  
    u1(j,7)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+(bxLR*bxLR+byR*byR+bzR*bzR)/(8.*pi)
    if(j>=0 .and. j<=nx)then
      write(11,*) x,rhoR,uR,vR,wR,pR
      write(12,*) x,byR,bzR
    end if
  end if

end do

print*,'Left State (rho,u,v,w,Bx,By,Bz,p) =',rhoL,uL,vL,wL,bxLR,byL,bzL,pL
print*,'Right State (rho,u,v,w,Bx,By,Bz,p) =',rhoR,uR,vR,wR,bxLR,byR,bzR,pR

end subroutine init

subroutine computeFlux()

real*8::BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
real*8::vxt(7),uxt(7),phi

do j=0,nx
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,u1(j,1))
uL=u1(j,2)/u1(j,1)
vL=u1(j,3)/u1(j,1)
wL=u1(j,4)/u1(j,1)
byL=u1(j,5)
bzL=u1(j,6)

BBL=bxLR*bxLR+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL

pL=max(premin,(gam-1.)*(u1(j,7)-0.5*rhoL*vvL-BBL/(8.*pi)))
aL=sqrt(gam*pL/rhoL)

!VLx(1)=u1(j,1)
!VLx(2)=u1(j,2)
!VLx(3)=u1(j,3)
!VLx(4)=u1(j,4)
!VLx(5)=u1(j,5)
!VLx(6)=u1(j,6)
!VLx(7)=u1(j,7)

VLx(1)=rhoL
VLx(2)=uL
VLx(3)=vL
VLx(4)=wL
VLx(5)=byL
VLx(6)=bzL
VLx(7)=pL

rhoR=max(densmin,u1(j+1,1))
uR=u1(j+1,2)/u1(j+1,1)
vR=u1(j+1,3)/u1(j+1,1)
wR=u1(j+1,4)/u1(j+1,1)
byR=u1(j+1,5)
bzR=u1(j+1,6)

BBR=bxLR*bxLR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR

pR=max(premin,(gam-1.)*(u1(j+1,7)-0.5*rhoR*vvR-BBR/(8.*pi)))
aR=sqrt(gam*pR/rhoR)

!VRx(1)=u1(j+1,1)
!VRx(2)=u1(j+1,2)
!VRx(3)=u1(j+1,3)
!VRx(4)=u1(j+1,4)
!VRx(5)=u1(j+1,5)
!VRx(6)=u1(j+1,6)
!VRx(7)=u1(j+1,7)

VRx(1)=rhoR
VRx(2)=uR
VRx(3)=vR
VRx(4)=wR
VRx(5)=byR
VRx(6)=bzR
VRx(7)=pR


FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL+(BBL-2.*bxLR*bxLR)/(8.*pi)
FL(3)=rhoL*uL*vL-bxLR*byL/(4.*pi)
FL(4)=rhoL*uL*wL-bxLR*bzL/(4.*pi)
FL(5)=uL*byL-vL*bxLR
FL(6)=uL*bzL-wL*bxLR
FL(7)=uL*(u1(j,7)+pL+BBL/(8.*pi)) &
      -bxLR*(bxLR*uL+byL*vL+bzL*wL)/(4.*pi)


FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR+(BBR-2.*bxLR*bxLR)/(8.*pi)
FR(3)=rhoR*uR*vR-bxLR*byR/(4.*pi)
FR(4)=rhoR*uR*wR-bxLR*bzR/(4.*pi)
FR(5)=uR*byR-vR*bxLR
FR(6)=uR*bzR-wR*bxLR
FR(7)=uR*(u1(j+1,7)+pR+BBR/(8.*pi)) &
      -bxLR*(bxLR*uR+byR*vR+bzR*wR)/(4.*pi)

if(debug==1)then
  print*,'Cell Inerface: ',j,'+1/2'
  print*,'Left State =',rhoL,uL,vL,wL,bxLR,byL,bzL,pL
  print*,'Right State =',rhoR,uR,vR,wR,bxLR,byR,bzR,pR
end if

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5*(rhoL+rhoR)
ut=0.5*(uL+uR)
vt=0.5*(vL+vR)
wt=0.5*(wL+wR)
byt=0.5*(byL+byR)
bzt=0.5*(bzL+bzR)
pt=0.5*(pL+pR)
at=sqrt(gam*pt/rhot)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bx=sqrt(bxLR**2./(4.*pi*rhot))
by=sqrt(byt**2./(4.*pi*rhot))
bz=sqrt(bzt**2./(4.*pi*rhot))
bsqr=bx*bx+by*by+bz*bz

cA=bx
cs=0.5*(at*at+bsqr-sqrt((at*at+bsqr)**2.-4.*at*at*bx*bx))
cf=0.5*(at*at+bsqr+sqrt((at*at+bsqr)**2.-4.*at*at*bx*bx))
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if

!---------------------------------------------
!Characteristic Eigenvalues: lambda_1..7
!---------------------------------------------
lambda(1)=ut-cf
lambda(2)=ut-cA
lambda(3)=ut-cs
lambda(4)=ut
lambda(5)=ut+cs
lambda(6)=ut+cA
lambda(7)=ut+cf

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(1),lambda(2),lambda(3) &
        ,lambda(4),lambda(5),lambda(6),lambda(7)
end if

!-------------------------------------------------
!Roe Jacobian Matrix Eigenvectors: R_1..7, L_1..7
!-------------------------------------------------
if(abs(by)<tol .and. abs(bz)<tol)then
 betay=1./sqrt(2.)
 betaz=1./sqrt(2.)
else
 betay=by/sqrt(by*by+bz*bz)
 betaz=bz/sqrt(by*by+bz*bz)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(by)<tol .and. abs(bz)<tol .and. abs(cA-at)<tol)then
 phi=atan(sqrt(by*by+bz*bz)/(abs(bx)-at))
 alphas=cos(0.5*phi)
 alphaf=sin(0.5*phi)
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt((cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt((at*at-cs*cs)/(cf*cf-cs*cs)) 
end if



R1(1)=alphaf*rhot
R1(2)=-alphaf*cf
R1(3)=alphas*cs*betay*sgn(bx)
R1(4)=alphas*cs*betaz*sgn(bx)
R1(5)=alphas*betay*at*sqrt(4.*pi*rhot)
R1(6)=alphas*betaz*at*sqrt(4.*pi*rhot)
R1(7)=alphaf*rhot*at*at

R2(1)=0.
R2(2)=0.
R2(3)=-betaz
R2(4)=betay
R2(5)=-betaz*sqrt(4.*pi*rhot)*sgn(bx)
R2(6)=betay*sqrt(4.*pi*rhot)*sgn(bx)
R2(7)=0.

R3(1)=alphas*rhot
R3(2)=-alphas*cs
R3(3)=-alphaf*betay*cf*sgn(bx)
R3(4)=-alphaf*betaz*cf*sgn(bx)
R3(5)=-alphaf*betay*at*sqrt(4.*pi*rhot)
R3(6)=-alphaf*betaz*at*sqrt(4.*pi*rhot)
R3(7)=alphas*rhot*at*at

R4(1)=1.
R4(2)=0.
R4(3)=0.
R4(4)=0.
R4(5)=0.
R4(6)=0.
R4(7)=0.

R5(1)=alphas*rhot
R5(2)=alphas*cs
R5(3)=alphaf*betay*cf*sgn(bx)
R5(4)=alphaf*betaz*cf*sgn(bx)
R5(5)=-alphaf*betay*at*sqrt(4.*pi*rhot)
R5(6)=-alphaf*betaz*at*sqrt(4.*pi*rhot)
R5(7)=alphas*rhot*at*at

R6(1)=0.
R6(2)=0.
R6(3)=betaz
R6(4)=-betay
R6(5)=-betaz*sqrt(4.*pi*rhot)*sgn(bx)
R6(6)=betay*sqrt(4.*pi*rhot)*sgn(bx)
R6(7)=0.

R7(1)=alphaf*rhot
R7(2)=alphaf*cf
R7(3)=-alphas*cs*betay*sgn(bx)
R7(4)=-alphas*cs*betaz*sgn(bx)
R7(5)=alphas*betay*at*sqrt(4.*pi*rhot)
R7(6)=alphas*betaz*at*sqrt(4.*pi*rhot)
R7(7)=alphaf*rhot*at*at


L1(1)=0.
L1(2)=-alphaf*cf/(2.*at*at)
L1(3)=alphas*cs*betay*sgn(bx)/(2.*at*at)
L1(4)=alphas*cs*betaz*sgn(bx)/(2.*at*at)
L1(5)=(alphas*at*betay/sqrt(4.*pi*rhot))/(2.*at*at)
L1(6)=(alphas*at*betaz/sqrt(4.*pi*rhot))/(2.*at*at)
L1(7)=(alphas/rhot)/(2.*at*at)

L2(1)=0.
L2(2)=0.
L2(3)=-0.5*betaz
L2(4)=0.5*betay
L2(5)=-0.5*betaz*sgn(bx)/sqrt(4.*pi*rhot)
L2(6)=0.5*betay*sgn(bx)/sqrt(4.*pi*rhot)
L2(7)=0.

L3(1)=0.
L3(2)=-alphas*cs/(2.*at*at)
L3(3)=-alphaf*cf*betay*sgn(bx)/(2.*at*at)
L3(4)=-alphaf*cf*betaz*sgn(bx)/(2.*at*at)
L3(5)=-(alphaf*at*betay/sqrt(4.*pi*rhot))/(2.*at*at)
L3(6)=-(alphaf*at*betaz/sqrt(4.*pi*rhot))/(2.*at*at)
L3(7)=(alphas/rhot)/(2.*at*at)


L4(1)=1.
L4(2)=0.
L4(3)=0.
L4(4)=0.
L4(5)=0.
L4(6)=0.
L4(7)=-1./(at*at)

L5(1)=0.
L5(2)=alphas*cs/(2.*at*at)
L5(3)=alphaf*cf*betay*sgn(bx)/(2.*at*at)
L5(4)=alphaf*cf*betaz*sgn(bx)/(2.*at*at)
L5(5)=-(alphaf*at*betay/sqrt(4.*pi*rhot))/(2.*at*at)
L5(6)=-(alphaf*at*betaz/sqrt(4.*pi*rhot))/(2.*at*at)
L5(7)=(alphas/rhot)/(2.*at*at)


L6(1)=0.
L6(2)=0.
L6(3)=0.5*betaz
L6(4)=-0.5*betay
L6(5)=-0.5*betaz*sgn(bx)/sqrt(4.*pi*rhot)
L6(6)=0.5*betay*sgn(bx)/sqrt(4.*pi*rhot)
L6(7)=0.


L7(1)=0.
L7(2)=alphaf*cf/(2.*at*at)
L7(3)=-alphas*cs*betay*sgn(bx)/(2.*at*at)
L7(4)=-alphas*cs*betaz*sgn(bx)/(2.*at*at)
L7(5)=(alphas*at*betay/sqrt(4.*pi*rhot))/(2.*at*at)
L7(6)=(alphas*at*betaz/sqrt(4.*pi*rhot))/(2.*at*at)
L7(7)=(alphaf/rhot)/(2.*at*at)

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------
alpha=0.

do k=1,7
  alpha(1)=alpha(1)+L1(k)*(VRx(k)-VLx(k))  
  alpha(2)=alpha(2)+L2(k)*(VRx(k)-VLx(k))
  alpha(3)=alpha(3)+L3(k)*(VRx(k)-VLx(k))
  alpha(4)=alpha(4)+L4(k)*(VRx(k)-VLx(k)) 
  alpha(5)=alpha(5)+L5(k)*(VRx(k)-VLx(k))
  alpha(6)=alpha(6)+L6(k)*(VRx(k)-VLx(k))
  alpha(7)=alpha(7)+L7(k)*(VRx(k)-VLx(k))
end do

!---------------------------------------
!Maximum Wave Speed
!---------------------------------------
do k=1,7
  smax=max(smax,abs(lambda(k)))
  !print*,'k,|lambda|,smax=',k,abs(lambda(k)),smax
end do


if(debug==1)then
  print*,'Max wave speed=',smax
end if

!-------------------------------------------------------
!Solution at Cell Interface: v((x/t)=0)
!-------------------------------------------------------
do k=1,7
  vxt(k)=0.5*(VLx(k)+VRx(k))-0.5*( &
         sign(1._8,lambda(1))*alpha(1)*R1(k) &
        +sign(1._8,lambda(2))*alpha(2)*R2(k) &
        +sign(1._8,lambda(3))*alpha(3)*R3(k) &  
        +sign(1._8,lambda(4))*alpha(4)*R4(k) &
        +sign(1._8,lambda(5))*alpha(5)*R5(k) &
        +sign(1._8,lambda(6))*alpha(6)*R6(k) &
        +sign(1._8,lambda(7))*alpha(7)*R7(k) )
end do

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!-------------------------------------------------------------
BBL=bxLR*bxLR+vxt(5)*vxt(5)+vxt(6)*vxt(6)
vvL=vxt(2)*vxt(2)+vxt(3)*vxt(3)+vxt(4)*vxt(4)

flux(j,1)=vxt(1)*vxt(2)
flux(j,2)=vxt(1)*vxt(2)*vxt(2)+vxt(7)+(BBL-2.*bxLR*bxLR)/(8.*pi)
flux(j,3)=vxt(1)*vxt(2)*vxt(3)-bxLR*vxt(5)/(4.*pi)
flux(j,4)=vxt(1)*vxt(2)*vxt(4)-bxLR*vxt(6)/(4.*pi)
flux(j,5)=vxt(2)*vxt(5)-vxt(3)*bxLR
flux(j,6)=vxt(2)*vxt(6)-vxt(4)*bxLR
flux(j,7)=vxt(2)*(0.5*vxt(1)*vvL+gam*gamul*vxt(7)+BBL/(4.*pi)) &
      -bxLR*(bxLR*vxt(2)+vxt(5)*vxt(3)+vxt(6)*vxt(4))/(4.*pi)
  
if(debug==1)then
  print*,'flux=',flux(j,1),flux(j,2),flux(j,3),flux(j,4),flux(j,5),flux(j,6),flux(j,7)
end if

end do

end subroutine computeFlux

subroutine bound()

!enforce outflow boundary conditions
do k=1,7
  u2(0,k)=u2(1,k)	
  u2(-1,k)=u2(1,k)
  u2(nx+1,k)=u2(nx,k)
  u2(nx+2,k)=u2(nx,k)
end do

end subroutine bound

function sgn(x) result(fx)
  real*8,intent(in)::x
  real*8::fx
  fx=sign(1._8,x)
end

function psi(x) result(fx)
  real*8,intent(in)::x
  real*8::fx,del

  del=0.
  if(abs(x)<del)then
    fx=(x*x+del*del)/(2.*del)
  else if(abs(x)>=del)then
    fx=abs(x)
  end if  
end

end module RoeSolverMHD_mod
