!Roe's Approximate (Linear) Riemann Solver for the equations of Ideal MHD
module RoeSolverMHD_mod
use vars_mod
implicit none


contains 


subroutine computeFlux(u1,flux)

!Input Variables
real*8::u1(-1:nx+2,8),flux(-1:nx+2,7)

!------------------------------------------------------------------------------
real*8::BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
integer::j,k
!------------------------------------------------------------------------------
real*8::rhoL,uL,vL,wL,pL,byL,bzL,ptotL,FL(7) !left state
real*8::rhoR,uR,vR,wR,pR,byR,bzR,ptotR,FR(7) !right state 
real*8::bxLR
real*8::smax,cs,cf,cA
real*8::dens,velx,vely,velz,pres,magy,magz
real*8::rhot,ut,vt,wt,pt,byt,bzt,at,ptot,lambda(7)
real*8::L1(7),L2(7),L3(7),L4(7),L5(7),L6(7),L7(7)	
real*8::R1(7),R2(7),R3(7),R4(7),R5(7),R6(7),R7(7)	
real*8::alpha(7)
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
!------------------------------------------------------------------------------

!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do j=0,nx
!set left and right states for local Riemann problem at j+1/2 cell interface

rhoL=u1(j,1)
uL=u1(j,2)/u1(j,1)
vL=u1(j,3)/u1(j,1)
wL=u1(j,4)/u1(j,1)
bxLR=u1(j,8)
byL=u1(j,5)
bzL=u1(j,6)


BBL=bxLR*bxLR+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(j,7)-0.5*rhoL*vvL-0.5*BBL)
ptotL=pL+0.5*BBL

rhoR=u1(j+1,1)
uR=u1(j+1,2)/u1(j+1,1)
vR=u1(j+1,3)/u1(j+1,1)
wR=u1(j+1,4)/u1(j+1,1)
byR=u1(j+1,5)
bzR=u1(j+1,6)

BBR=bxLR*bxLR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(j+1,7)-0.5*rhoR*vvR-0.5*BBR)
ptotR=pR+0.5*BBR


FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL+0.5*(BBL-2.*bxLR*bxLR)
FL(3)=rhoL*uL*vL-bxLR*byL
FL(4)=rhoL*uL*wL-bxLR*bzL
FL(5)=uL*byL-vL*bxLR
FL(6)=uL*bzL-wL*bxLR
FL(7)=uL*(u1(j,7)+pL+0.5*BBL) &
      -bxLR*(bxLR*uL+byL*vL+bzL*wL)


FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR+0.5*(BBR-2.*bxLR*bxLR)
FR(3)=rhoR*uR*vR-bxLR*byR
FR(4)=rhoR*uR*wR-bxLR*bzR
FR(5)=uR*byR-vR*bxLR
FR(6)=uR*bzR-wR*bxLR
FR(7)=uR*(u1(j+1,7)+pR+0.5*BBR) &
      -bxLR*(bxLR*uR+byR*vR+bzR*wR)



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
ptot=0.5*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5*(bxLR*bxLR+byt*byt+bzt*bzt)
at=sqrt(gam*pt/rhot)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bx=sqrt(bxLR**2./rhot)
by=sqrt(byt**2./rhot)
bz=sqrt(bzt**2./rhot)
bsqr=bx*bx+by*by+bz*bz
vsq=ut*ut+vt*vt+wt*wt

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
 betay=1._8/sqrt(2.)
 betaz=1._8/sqrt(2.)
else
 betay=byt/sqrt(byt*byt+bzt*bzt)
 betaz=bzt/sqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(by)<tol .and. abs(bz)<tol .and. abs(cA-at)<tol)then
 alphas=1._8
 alphaf=1._8
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
end if

gf7=gamul*alphaf*cf*cf+alphaf*cf*ut-alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gf1=gamul*alphaf*cf*cf-alphaf*cf*ut+alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphaf*(cf*cf-at*at)
gs5=gamul*alphas*cs*cs+alphas*cs*ut+alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)
gs3=gamul*alphas*cs*cs-alphas*cs*ut-alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2.)*alphas*(cs*cs-at*at)

theta1=alphaf*alphaf*at*at*(cf*cf+gamul*(2.-gam)*at*at) &
       +alphas*alphas*cf*cf*(cs*cs+gamul*(2.-gam)*at*at)

theta2=(alphaf*alphaf*cf*at+alphas*alphas*cs*cA)*sgn(bx)


R1(1)=alphaf
R1(2)=alphaf*(ut-cf)
R1(3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(5)=(alphas*betay*cf)/sqrt(rhot)
R1(6)=(alphas*betaz*cf)/sqrt(rhot)
R1(7)=alphaf*0.5*vsq+gf1

R2(1)=0.
R2(2)=0.
R2(3)=betaz*sgn(bxLR)
R2(4)=-betay*sgn(bxLR)
R2(5)=betaz/sqrt(rhot)
R2(6)=-betay/sqrt(rhot)
R2(7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(1)=alphas
R3(2)=alphas*(ut-cs)
R3(3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R3(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R3(7)=alphas*0.5*vsq+gs3

R4(1)=1.
R4(2)=ut
R4(3)=vt
R4(4)=wt
R4(5)=0.
R4(6)=0.
R4(7)=0.5*vsq

R5(1)=alphas
R5(2)=alphas*(ut+cs)
R5(3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R5(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R5(7)=alphas*0.5*vsq+gs5

R6(1)=0.
R6(2)=0.
R6(3)=-betaz*sgn(bxLR)
R6(4)=betay*sgn(bxLR)
R6(5)=betaz/sqrt(rhot)
R6(6)=-betay/sqrt(rhot)
R6(7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(1)=alphaf
R7(2)=alphaf*(ut+cf)
R7(3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(5)=(alphas*betay*cf)/sqrt(rhot)
R7(6)=(alphas*betaz*cf)/sqrt(rhot)
R7(7)=alphaf*0.5*vsq+gf7


L1(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L1(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L1(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

L1(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

L1(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(7)=(1./theta1)*0.5*alphaf*at*at

L2(1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L2(2)=0.
L2(3)=0.5*betaz*sgn(bxLR)
L2(4)=-0.5*betay*sgn(bxLR)
L2(5)=0.5*betaz*sqrt(rhot)
L2(6)=-0.5*betay*sqrt(rhot)
L2(7)=0.

L3(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L3(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L3(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

L3(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

L3(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(7)=(1./theta1)*0.5*alphas*cf*cf


L4(1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )


L4(5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L5(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L5(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

L5(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

L5(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(7)=(1./theta1)*0.5*alphas*cf*cf


L6(1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L6(2)=0.
L6(3)=-0.5*betaz*sgn(bxLR)
L6(4)=0.5*betay*sgn(bxLR)
L6(5)=0.5*betaz*sqrt(rhot)
L6(6)=-0.5*betay*sqrt(rhot)
L6(7)=0.


L7(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L7(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L7(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

L7(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

L7(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(7)=(1./theta1)*0.5*alphaf*at*at


if(at>cA)then
 if(byt<=0. .and. bzt<0.)then
  R3=-1.*R3
  R5=-1.*R5
  L3=-1.*L3
  L5=-1.*L5
 end if
else if(at<cA)then
 if(byt<=0. .and. bzt<0.)then
  R1=-1.*R1
  R7=-1.*R7
  L1=-1.*L1
  L7=-1.*L7
 end if
end if


!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------
alpha=0.

do k=1,7  
  alpha(1)=alpha(1)+L1(k)*(u1(j+1,k)-u1(j,k))  
  alpha(2)=alpha(2)+L2(k)*(u1(j+1,k)-u1(j,k))
  alpha(3)=alpha(3)+L3(k)*(u1(j+1,k)-u1(j,k))
  alpha(4)=alpha(4)+L4(k)*(u1(j+1,k)-u1(j,k)) 
  alpha(5)=alpha(5)+L5(k)*(u1(j+1,k)-u1(j,k))
  alpha(6)=alpha(6)+L6(k)*(u1(j+1,k)-u1(j,k))
  alpha(7)=alpha(7)+L7(k)*(u1(j+1,k)-u1(j,k))
end do

!---------------------------------------
!Maximum Wave Speed
!---------------------------------------
!smax=0.
!do k=1,7
!  smax=max(smax,abs(lambda(k)))
  !print*,'k,|lambda|,smax=',k,abs(lambda(k)),smax
!end do
smax=max(smax,abs(ut)+cf)

if(debug==1)then
  print*,'Max wave speed=',smax
end if

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do k=1,7
  flux(j,k)=0.5*(FL(k)+FR(k))-0.5*(     &
        abs(lambda(1))*alpha(1)*R1(k) &
       +abs(lambda(2))*alpha(2)*R2(k) &
       +abs(lambda(3))*alpha(3)*R3(k) &  
       +abs(lambda(4))*alpha(4)*R4(k) &
       +abs(lambda(5))*alpha(5)*R5(k) &
       +abs(lambda(6))*alpha(6)*R6(k) &
       +abs(lambda(7))*alpha(7)*R7(k) )

end do
 

if(debug==1)then
  print*,'flux=',flux(j,1),flux(j,2),flux(j,3),flux(j,4),flux(j,5),flux(j,6),flux(j,7)
end if

end do

end subroutine computeFlux

subroutine bound(u1)

!Input Variables
real*8::u1(-1:nx+2,7)

integer::k

!enforce outflow boundary conditions
do k=1,7
  u1(0,k)=u1(1,k)	
  u1(-1,k)=u1(1,k)
  u1(nx+1,k)=u1(nx,k)
  u1(nx+2,k)=u1(nx,k)
end do

end subroutine bound



function sgn(x) result(fx)
  real*8,intent(in)::x
  real*8::fx
  fx=sign(1._8,x)
end function sgn

function psi(x) result(fx)
  real*8,intent(in)::x
  real*8::fx,del

  del=0.2
  if(abs(x)<del)then
    fx=(x*x+del*del)/(2.*del)
  else if(abs(x)>=del)then
    fx=abs(x)
  end if  
end function psi

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


end module RoeSolverMHD_mod
