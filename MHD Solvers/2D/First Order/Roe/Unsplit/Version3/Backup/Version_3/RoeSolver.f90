!Roe's Approximate (Linear) Riemann Solver for the equations of Ideal MHD
module RoeSolverMHD_mod
use vars_mod
implicit none


contains 


subroutine computeFlux(u1,flux,g)

!Input Variables
real*8::u1(-1:nx+2,8),flux(-1:nx+2,7),g(-1:nx+2,7)

!------------------------------------------------------------------------------
real*8::BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
integer::j,k
!------------------------------------------------------------------------------
real*8::rhoL,uL,vL,wL,pL,byL,bzL,ptotL,FL(7) !left state
real*8::rhoR,uR,vR,wR,pR,byR,bzR,ptotR,FR(7) !right state 
real*8::bxLR
real*8::smax,cs,cf,cA
real*8::dens,velx,vely,velz,pres,magy,magz
real*8::rhot,ut,vt,wt,pt,byt,bzt,at,ptot
real*8::L1(-1:nx+2,7),L2(-1:nx+2,7),L3(-1:nx+2,7)&
,L4(-1:nx+2,7),L5(-1:nx+2,7),L6(-1:nx+2,7),L7(-1:nx+2,7)	
real*8::R1(-1:nx+2,7),R2(-1:nx+2,7),R3(-1:nx+2,7),R4(-1:nx+2,7)&
,R5(-1:nx+2,7),R6(-1:nx+2,7),R7(-1:nx+2,7)	
real*8::alpha(-1:nx+2,7),lambda(-1:nx+2,7),F(-1:nx+2,7)
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
!------------------------------------------------------------------------------


!Reset arrays
alpha=0.
F=0.
lambda=0.
R1=0.
R2=0.
R3=0.
R4=0.
R5=0.
R6=0.
R7=0.
L1=0.
L2=0.
L3=0.
L4=0.
L5=0.
L6=0.
L7=0.


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

F(j,1)=rhoL*uL
F(j,2)=rhoL*uL*uL+pL+0.5*(BBL-2.*bxLR*bxLR)
F(j,3)=rhoL*uL*vL-bxLR*byL
F(j,4)=rhoL*uL*wL-bxLR*bzL
F(j,5)=uL*byL-vL*bxLR
F(j,6)=uL*bzL-wL*bxLR
F(j,7)=uL*(u1(j,7)+pL+0.5*BBL) &
      -bxLR*(bxLR*uL+byL*vL+bzL*wL)

F(j+1,1)=rhoR*uR
F(j+1,2)=rhoR*uR*uR+pR+0.5*(BBR-2.*bxLR*bxLR)
F(j+1,3)=rhoR*uR*vR-bxLR*byR
F(j+1,4)=rhoR*uR*wR-bxLR*bzR
F(j+1,5)=uR*byR-vR*bxLR
F(j+1,6)=uR*bzR-wR*bxLR
F(j+1,7)=uR*(u1(j+1,7)+pR+0.5*BBR) &
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
lambda(j,1)=ut-cf
lambda(j,2)=ut-cA
lambda(j,3)=ut-cs
lambda(j,4)=ut
lambda(j,5)=ut+cs
lambda(j,6)=ut+cA
lambda(j,7)=ut+cf

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(j,1),lambda(j,2),lambda(j,3) &
        ,lambda(j,4),lambda(j,5),lambda(j,6),lambda(j,7)
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


R1(j,1)=alphaf
R1(j,2)=alphaf*(ut-cf)
R1(j,3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(j,4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(j,5)=(alphas*betay*cf)/sqrt(rhot)
R1(j,6)=(alphas*betaz*cf)/sqrt(rhot)
R1(j,7)=alphaf*0.5*vsq+gf1

R2(j,1)=0.
R2(j,2)=0.
R2(j,3)=betaz*sgn(bxLR)
R2(j,4)=-betay*sgn(bxLR)
R2(j,5)=betaz/sqrt(rhot)
R2(j,6)=-betay/sqrt(rhot)
R2(j,7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(j,1)=alphas
R3(j,2)=alphas*(ut-cs)
R3(j,3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(j,4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(j,5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R3(j,6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R3(j,7)=alphas*0.5*vsq+gs3

R4(j,1)=1.
R4(j,2)=ut
R4(j,3)=vt
R4(j,4)=wt
R4(j,5)=0.
R4(j,6)=0.
R4(j,7)=0.5*vsq

R5(j,1)=alphas
R5(j,2)=alphas*(ut+cs)
R5(j,3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(j,4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(j,5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R5(j,6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R5(j,7)=alphas*0.5*vsq+gs5

R6(j,1)=0.
R6(j,2)=0.
R6(j,3)=-betaz*sgn(bxLR)
R6(j,4)=betay*sgn(bxLR)
R6(j,5)=betaz/sqrt(rhot)
R6(j,6)=-betay/sqrt(rhot)
R6(j,7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(j,1)=alphaf
R7(j,2)=alphaf*(ut+cf)
R7(j,3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(j,4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(j,5)=(alphas*betay*cf)/sqrt(rhot)
R7(j,6)=(alphas*betaz*cf)/sqrt(rhot)
R7(j,7)=alphaf*0.5*vsq+gf7


L1(j,1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L1(j,2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L1(j,3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

L1(j,4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

L1(j,5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(j,6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(j,7)=(1./theta1)*0.5*alphaf*at*at

L2(j,1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L2(j,2)=0.
L2(j,3)=0.5*betaz*sgn(bxLR)
L2(j,4)=-0.5*betay*sgn(bxLR)
L2(j,5)=0.5*betaz*sqrt(rhot)
L2(j,6)=-0.5*betay*sqrt(rhot)
L2(j,7)=0.

L3(j,1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L3(j,2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L3(j,3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

L3(j,4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

L3(j,5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(j,6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(j,7)=(1./theta1)*0.5*alphas*cf*cf


L4(j,1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )


L4(j,5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(j,6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(j,7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(j,1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L5(j,2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L5(j,3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

L5(j,4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

L5(j,5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(j,6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(j,7)=(1./theta1)*0.5*alphas*cf*cf


L6(j,1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L6(j,2)=0.
L6(j,3)=-0.5*betaz*sgn(bxLR)
L6(j,4)=0.5*betay*sgn(bxLR)
L6(j,5)=0.5*betaz*sqrt(rhot)
L6(j,6)=-0.5*betay*sqrt(rhot)
L6(j,7)=0.


L7(j,1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L7(j,2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L7(j,3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

L7(j,4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

L7(j,5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(j,6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(j,7)=(1./theta1)*0.5*alphaf*at*at


if(at>cA)then
 if(byt<=0. .and. bzt<0.)then
  do k=1,7
  R3(j,k)=-1.*R3(j,k)
  R5(j,k)=-1.*R5(j,k)
  L3(j,k)=-1.*L3(j,k)
  L5(j,k)=-1.*L5(j,k)
  end do
 end if
else if(at<cA)then
 if(byt<=0. .and. bzt<0.)then
  do k=1,7
  R1(j,k)=-1.*R1(j,k)
  R7(j,k)=-1.*R7(j,k)
  L1(j,k)=-1.*L1(j,k)
  L7(j,k)=-1.*L7(j,k)
  end do
 end if
end if


!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------

do k=1,7  
  alpha(j,1)=alpha(j,1)+L1(j,k)*(u1(j+1,k)-u1(j,k))  
  alpha(j,2)=alpha(j,2)+L2(j,k)*(u1(j+1,k)-u1(j,k))
  alpha(j,3)=alpha(j,3)+L3(j,k)*(u1(j+1,k)-u1(j,k))
  alpha(j,4)=alpha(j,4)+L4(j,k)*(u1(j+1,k)-u1(j,k)) 
  alpha(j,5)=alpha(j,5)+L5(j,k)*(u1(j+1,k)-u1(j,k))
  alpha(j,6)=alpha(j,6)+L6(j,k)*(u1(j+1,k)-u1(j,k))
  alpha(j,7)=alpha(j,7)+L7(j,k)*(u1(j+1,k)-u1(j,k))
end do

!---------------------------------------
!Maximum Wave Speed
!---------------------------------------
smax=max(smax,abs(ut)+cf)

if(debug==1)then
  print*,'Max wave speed=',smax
end if

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

end do

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do j=0,nx
 do k=1,7
  flux(j,k)=0.5*(F(j,k)+F(j+1,k))-0.5*(     &
        abs(lambda(j,1))*alpha(j,1)*R1(j,k) &
       +abs(lambda(j,2))*alpha(j,2)*R2(j,k) &
       +abs(lambda(j,3))*alpha(j,3)*R3(j,k) &  
       +abs(lambda(j,4))*alpha(j,4)*R4(j,k) &
       +abs(lambda(j,5))*alpha(j,5)*R5(j,k) &
       +abs(lambda(j,6))*alpha(j,6)*R6(j,k) &
       +abs(lambda(j,7))*alpha(j,7)*R7(j,k) )

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
