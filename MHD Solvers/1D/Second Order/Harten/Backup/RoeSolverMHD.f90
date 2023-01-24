!Roe's Approximate (Linear) Riemann Solver for the equations of Ideal MHD
module RoeSolverMHD_mod
use vars_mod
implicit none


integer ::i,j,k
integer,parameter::flux_correction=0 !0:off 1:on

contains 


subroutine computeFlux(u1,flux)

!Input Variables
real*8::u1(-1:nx+2,7),flux(-1:nx+2,7)

!Local Variables
!------------------------------------------------------------------------------
real*8 ::temp1,temp2,temp3,xc,xt,theta,w
!------------------------------------------------------------------------------
real*8::gt(-1:nx+2,7),fc(-1:nx+2,7),gamma(-1:nx+2,7),sigma(-1:nx+2,7)
real*8::f(-1:nx+2,7)
!------------------------------------------------------------------------------
real*8::R1(-1:nx+2,7),R2(-1:nx+2,7),R3(-1:nx+2,7),R4(-1:nx+2,7)&
,R5(-1:nx+2,7),R6(-1:nx+2,7),R7(-1:nx+2,7)
real*8::L1(-1:nx+2,7),L2(-1:nx+2,7),L3(-1:nx+2,7),L4(-1:nx+2,7)&
,L5(-1:nx+2,7),L6(-1:nx+2,7),L7(-1:nx+2,7)	
real*8::alpha(-1:nx+2,7),lambda(-1:nx+2,7),smax	
real*8::uleft(-1:nx+2,7),uright(-1:nx+2,7)

real*8::gf1,gf7,gs3,gs5,theta1,theta2,vsq
real*8::rhoL,uL,vL,wL,pL,byL,bzL,ptotL !left state
real*8::rhoR,uR,vR,wR,pR,byR,bzR,ptotR !right state 
real*8::BBR,BBL,vvR,vvL
real*8::rhot,ut,vt,wt,bxt,byt,bzt,ptot,pt,at
real*8::cs,cf,cA
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr

!------------------------------------------------------------------------------
!Solve Local Riemann Problems: w_j+1/2(0)

!--------------------------------------------
!Vector of Left and Right States: uL_j, uR_j
!--------------------------------------------
do j=-1,nx+1
 do k=1,7
  uleft(j,k)=u1(j,k)
  uright(j,k)=u1(j+1,k)
 end do
end do

!-----------------------------------------------------------------------------------------
!Roe Matrix Eigenvalues, Right and Left Eigenvectors: lambda^k_j+1/2, R^k_j+1/2, L^k_j+1/2
!-----------------------------------------------------------------------------------------
do j=-1,nx+1

 rhoL=uleft(j,1)
 uL=uleft(j,2)/uleft(j,1)
 vL=uleft(j,3)/uleft(j,1)
 wL=uleft(j,4)/uleft(j,1)
 byL=uleft(j,5)
 bzL=uleft(j,6)
 BBL=bxLR*bxLR+byL*byL+bzL*bzL
 vvL=uL*uL+vL*vL+wL*wL
 pL=(gam-1._8)*(uleft(j,7)-0.5_8*rhoL*vvL-0.5_8*BBL)
 ptotL=pL+0.5_8*BBL

 rhoR=uright(j,1)
 uR=uright(j,2)/uright(j,1)
 vR=uright(j,3)/uright(j,1)
 wR=uright(j,4)/uright(j,1)
 byR=uright(j,5)
 bzR=uright(j,6)
 BBR=bxLR*bxLR+byR*byR+bzR*bzR
 vvR=uR*uR+vR*vR+wR*wR
 pR=(gam-1._8)*(uright(j,7)-0.5_8*rhoR*vvR-0.5_8*BBR)
 ptotR=pR+0.5_8*BBR

 if(debug==1)then
  print*,'Cell Inerface: ',j,'+1/2'
  print*,'Left State =',rhoL,uL,vL,wL,bxLR,byL,bzL,pL
  print*,'Right State =',rhoR,uR,vR,wR,bxLR,byR,bzR,pR
 end if

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
bxt=bxLR
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
at=dsqrt(gam*pt/rhot)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bx=bxt*bxt/rhot
by=byt*byt/rhot
bz=bzt*bzt/rhot
bsqr=bx+by+bz
vsq=ut*ut+vt*vt+wt*wt

cA=dsqrt(bx)
cs=0.5_8*(at*at+bsqr-dsqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx))
cf=0.5_8*(at*at+bsqr+dsqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx))
if(debug==1)then
  print*,'cs^2,cf^2=',cs,cf
end if
cs=dsqrt(cs)
cf=dsqrt(cf)

if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if

!---------------------------------------
!Maximum Wave Speed
!---------------------------------------
smax=max(smax,abs(ut)+cf)
if(debug==1)then
  print*,'Max wave speed=',smax
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
if(abs(byt)<tol .and. abs(bzt)<tol)then
 betay=1._8/dsqrt(2._8)
 betaz=1._8/dsqrt(2._8)
else
 betay=byt/dsqrt(byt*byt+bzt*bzt)
 betaz=bzt/dsqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
 alphas=1._8
 alphaf=1._8
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=dsqrt((cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=dsqrt((cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
end if

gf7=gamul*alphaf*cf*cf+alphaf*cf*ut-alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2._8)*alphaf*(cf*cf-at*at)
gf1=gamul*alphaf*cf*cf-alphaf*cf*ut+alphas*cA*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2._8)*alphaf*(cf*cf-at*at)
gs5=gamul*alphas*cs*cs+alphas*cs*ut+alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2._8)*alphas*(cs*cs-at*at)
gs3=gamul*alphas*cs*cs-alphas*cs*ut-alphaf*at*sgn(bxLR)*(betay*vt+betaz*wt) &
    +gamul*(gam-2._8)*alphas*(cs*cs-at*at)

theta1=alphaf*alphaf*at*at*(cf*cf+gamul*(2._8-gam)*at*at) &
       +alphas*alphas*cf*cf*(cs*cs+gamul*(2._8-gam)*at*at)

theta2=(alphaf*alphaf*cf*at+alphas*alphas*cs*cA)*sgn(bxLR)


R1(j,1)=alphaf
R1(j,2)=alphaf*(ut-cf)
R1(j,3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(j,4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(j,5)=(alphas*betay*cf)/dsqrt(rhot)
R1(j,6)=(alphas*betaz*cf)/dsqrt(rhot)
R1(j,7)=alphaf*0.5_8*vsq+gf1

R2(j,1)=0.
R2(j,2)=0.
R2(j,3)=betaz*sgn(bxLR)
R2(j,4)=-betay*sgn(bxLR)
R2(j,5)=betaz/dsqrt(rhot)
R2(j,6)=-betay/dsqrt(rhot)
R2(j,7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(j,1)=alphas
R3(j,2)=alphas*(ut-cs)
R3(j,3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(j,4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(j,5)=-(alphaf*betay*at*at)/(cf*dsqrt(rhot))
R3(j,6)=-(alphaf*betaz*at*at)/(cf*dsqrt(rhot))
R3(j,7)=alphas*0.5_8*vsq+gs3

R4(j,1)=1._8
R4(j,2)=ut
R4(j,3)=vt
R4(j,4)=wt
R4(j,5)=0._8
R4(j,6)=0._8
R4(j,7)=0.5_8*vsq

R5(j,1)=alphas
R5(j,2)=alphas*(ut+cs)
R5(j,3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(j,4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(j,5)=-(alphaf*betay*at*at)/(cf*dsqrt(rhot))
R5(j,6)=-(alphaf*betaz*at*at)/(cf*dsqrt(rhot))
R5(j,7)=alphas*0.5_8*vsq+gs5

R6(j,1)=0._8
R6(j,2)=0._8
R6(j,3)=-betaz*sgn(bxLR)
R6(j,4)=betay*sgn(bxLR)
R6(j,5)=betaz/dsqrt(rhot)
R6(j,6)=-betay/dsqrt(rhot)
R6(j,7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(j,1)=alphaf
R7(j,2)=alphaf*(ut+cf)
R7(j,3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(j,4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(j,5)=(alphas*betay*cf)/dsqrt(rhot)
R7(j,6)=(alphas*betaz*cf)/dsqrt(rhot)
R7(j,7)=alphaf*0.5_8*vsq+gf7


L1(j,1)=(1._8/theta1)*0.25_8*alphaf*at*at*vsq & 
     +(1._8/theta2)*( 0.5_8*alphaf*at*ut*sgn(bxLR) &
     -0.5_8*alphas*cs*(betay*vt+betaz*wt) )

L1(j,2)=-(1._8/theta1)*0.5_8*alphaf*at*at*ut & 
      -(1._8/theta2)*0.5_8*alphaf*at*sgn(bxLR)

L1(j,3)=-(1._8/theta1)*0.5_8*alphaf*at*at*vt & 
      +(1._8/theta2)*0.5_8*alphas*betay*cs

L1(j,4)=-(1._8/theta1)*0.5_8*alphaf*at*at*wt & 
      +(1._8/theta2)*0.5_8*alphas*betaz*cs

L1(j,5)=(1._8/theta1)*0.5_8*alphas*betay*cf* &
      (cs*cs+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L1(j,6)=(1._8/theta1)*0.5_8*alphas*betaz*cf* &
      (cs*cs+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L1(j,7)=(1._8/theta1)*0.5_8*alphaf*at*at

L2(j,1)=-0.5_8*(betaz*vt-betay*wt)*sgn(bxLR)
L2(j,2)=0._8
L2(j,3)=0.5_8*betaz*sgn(bxLR)
L2(j,4)=-0.5_8*betay*sgn(bxLR)
L2(j,5)=0.5_8*betaz*dsqrt(rhot)
L2(j,6)=-0.5_8*betay*dsqrt(rhot)
L2(j,7)=0._8

L3(j,1)=(1._8/theta1)*0.25_8*alphas*cf*cf*vsq & 
     +(1._8/theta2)*( 0.5_8*alphas*cA*ut*sgn(bxLR) &
     +0.5_8*alphaf*cf*(betay*vt+betaz*wt) )

L3(j,2)=-(1._8/theta1)*0.5_8*alphas*cf*cf*ut & 
      -(1._8/theta2)*0.5_8*alphas*cA*sgn(bxLR)

L3(j,3)=-(1._8/theta1)*0.5_8*alphas*cf*cf*vt & 
      -(1._8/theta2)*0.5_8*alphaf*betay*cf

L3(j,4)=-(1._8/theta1)*0.5_8*alphas*cf*cf*wt & 
      -(1._8/theta2)*0.5_8*alphaf*betaz*cf

L3(j,5)=-(1._8/theta1)*0.5_8*alphaf*betay*cf* &
      (cf*cf+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L3(j,6)=-(1._8/theta1)*0.5_8*alphaf*betaz*cf* &
      (cf*cf+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L3(j,7)=(1._8/theta1)*0.5_8*alphas*cf*cf


L4(j,1)=1._8-(1._8/theta1)*0.5_8*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,2)=(1._8/theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,3)=(1._8/theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,4)=(1._8/theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(j,5)=(1._8/theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*dsqrt(rhot)

L4(j,6)=(1._8/theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*dsqrt(rhot)

L4(j,7)=-(1._8/theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(j,1)=(1._8/theta1)*0.25_8*alphas*cf*cf*vsq & 
     -(1._8/theta2)*( 0.5_8*alphas*cA*ut*sgn(bxLR) &
     +0.5_8*alphaf*cf*(betay*vt+betaz*wt) )

L5(j,2)=-(1._8/theta1)*0.5_8*alphas*cf*cf*ut & 
      +(1._8/theta2)*0.5_8*alphas*cA*sgn(bxLR)

L5(j,3)=-(1._8/theta1)*0.5_8*alphas*cf*cf*vt & 
      +(1._8/theta2)*0.5_8*alphaf*betay*cf

L5(j,4)=-(1._8/theta1)*0.5_8*alphas*cf*cf*wt & 
      +(1._8/theta2)*0.5_8*alphaf*betaz*cf

L5(j,5)=-(1._8/theta1)*0.5_8*alphaf*betay*cf* &
      (cf*cf+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L5(j,6)=-(1._8/theta1)*0.5_8*alphaf*betaz*cf* &
      (cf*cf+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L5(j,7)=(1._8/theta1)*0.5_8*alphas*cf*cf


L6(j,1)=0.5_8*(betaz*vt-betay*wt)*sgn(bxLR)
L6(j,2)=0._8
L6(j,3)=-0.5_8*betaz*sgn(bxLR)
L6(j,4)=0.5_8*betay*sgn(bxLR)
L6(j,5)=0.5_8*betaz*dsqrt(rhot)
L6(j,6)=-0.5_8*betay*dsqrt(rhot)
L6(j,7)=0._8


L7(j,1)=(1._8/theta1)*0.25_8*alphaf*at*at*vsq & 
     -(1._8/theta2)*( 0.5_8*alphaf*at*ut*sgn(bxLR) &
     -0.5_8*alphas*cs*(betay*vt+betaz*wt) )

L7(j,2)=-(1._8/theta1)*0.5_8*alphaf*at*at*ut & 
      +(1._8/theta2)*0.5_8*alphaf*at*sgn(bxLR)

L7(j,3)=-(1._8/theta1)*0.5_8*alphaf*at*at*vt & 
      -(1._8/theta2)*0.5_8*alphas*betay*cs

L7(j,4)=-(1._8/theta1)*0.5_8*alphaf*at*at*wt & 
      -(1._8/theta2)*0.5_8*alphas*betaz*cs

L7(j,5)=(1._8/theta1)*0.5_8*alphas*betay*cf* &
      (cs*cs+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L7(j,6)=(1._8/theta1)*0.5_8*alphas*betaz*cf* &
      (cs*cs+gamul*(2._8-gam)*at*at)*dsqrt(rhot)

L7(j,7)=(1._8/theta1)*0.5_8*alphaf*at*at


if(at>cA)then
 if(byt<=tol .and. bzt<tol)then 
  do k=1,7
    R3(j,k)=-1._8*R3(j,k)
    R5(j,k)=-1._8*R5(j,k)
    L3(j,k)=-1._8*L3(j,k)
    L5(j,k)=-1._8*L5(j,k)
  end do
 end if
else if(at<cA)then
 if(byt<=tol .and. bzt<tol)then
  do k=1,7
   R1(j,k)=-1._8*R1(j,k)
   R7(j,k)=-1._8*R7(j,k)
   L1(j,k)=-1._8*R1(j,k)
   L7(j,k)=-1._8*L7(j,k)
  end do
 end if
end if


end do

!---------------------------------------------
!Roe Averaged Wave Strengths:alpha^k_j
!---------------------------------------------
alpha=0.
do j=-1,nx+1
 do k=1,7  
  alpha(j,1)=alpha(j,1)+L1(j,k)*(uright(j,k)-uleft(j,k))  
  alpha(j,2)=alpha(j,2)+L2(j,k)*(uright(j,k)-uleft(j,k))
  alpha(j,3)=alpha(j,3)+L3(j,k)*(uright(j,k)-uleft(j,k))
  alpha(j,4)=alpha(j,4)+L4(j,k)*(uright(j,k)-uleft(j,k)) 
  alpha(j,5)=alpha(j,5)+L5(j,k)*(uright(j,k)-uleft(j,k))
  alpha(j,6)=alpha(j,6)+L6(j,k)*(uright(j,k)-uleft(j,k))
  alpha(j,7)=alpha(j,7)+L7(j,k)*(uright(j,k)-uleft(j,k))
 end do


end do 

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

!------------------------------------------------
!~g_i+1/2
!------------------------------------------------
do j=-1,nx+1
 do k=1,7
   !gt(j,k)=(0.5_8*Q((dt/dx)*lambda(j,k),del(k))-0.5_8*(dt/dx)*(dt/dx) &
   !          *(lambda(j,k)*lambda(j,k)))*alpha(j,k)
   gt(j,k)=0.5_8*(Q(lambda(j,k),del(k))-(dt/dx) &
             *lambda(j,k)*lambda(j,k))*alpha(j,k) 
    
   if(shockSwitch==1)then
     !sigma(j,k)=0.5_8*(1._8-abs(lambda(j,k)) )
     sigma(j,k)=0.5_8*(1._8-abs(lambda(j,k)))
   end if	 
 end do  
end do

!-----------------------------------------------
!f(c)_j
!-----------------------------------------------
do j=0,nx+1
   do k=1,7
     fc(j,k)=minmod(gt(j-1,k),gt(j,k))
   end do  
   
   if(shockSwitch==1)then
    do k=1,7
     if(abs(alpha(j,k))+abs(alpha(j-1,k))>tol)then
      theta=(abs(alpha(j,k)-alpha(j-1,k)))/(abs(alpha(j,k))+abs(alpha(j-1,k)))
     else
      theta=0._8
     end if 
     w=minmod(alpha(j-1,k)*sigma(j-1,k),alpha(j,k)*sigma(j,k))
     
     !Modify the Alfven and entropy modes only (will result in steeper contact and rotational discontinuities)
     if(k==2 .or. k==4 .or.k==6)then 	 
       fc(j,k)=fc(j,k)+w*theta
     end if

    end do
   end if

end do

!----------------------------------------
!Cell Center Fluxes: F_j
!---------------------------------------
do j=0,nx+1 
 rhoL=uleft(j,1)
 uL=uleft(j,2)/uleft(j,1)
 vL=uleft(j,3)/uleft(j,1)
 wL=uleft(j,4)/uleft(j,1)
 byL=uleft(j,5)
 bzL=uleft(j,6)
 BBL=bxLR*bxLR+byL*byL+bzL*bzL
 vvL=uL*uL+vL*vL+wL*wL
 pL=(gam-1._8)*(uleft(j,7)-0.5_8*rhoL*vvL-0.5_8*BBL)


 f(j,1)=rhoL*uL
 f(j,2)=rhoL*uL*uL+pL+0.5_8*(BBL-2._8*bxLR*bxLR)
 f(j,3)=rhoL*uL*vL-bxLR*byL
 f(j,4)=rhoL*uL*wL-bxLR*bzL
 f(j,5)=uL*byL-vL*bxLR
 f(j,6)=uL*bzL-wL*bxLR
 f(j,7)=uL*(uleft(j,7)+pL+0.5_8*BBL) &
      -bxLR*(bxLR*uL+byL*vL+bzL*wL)

end do

!------------------------------------------------
!gamma_i+1/2
!------------------------------------------------
do j=0,nx  
 do k=1,7
   !if(abs(alpha(j,k))>tol)then
   if(alpha(j,k)>tol)then
     gamma(j,k)=(fc(j+1,k)-fc(j,k))/alpha(j,k)
   else
     gamma(j,k)=0._8
   end if	  
 end do

!------------------------------------------------
!Modified Numerical Flux: ^g_i+1/2
!------------------------------------------------
 do k=1,7
  if(flux_correction==1)then
   flux(j,k)=0.5_8*(f(j,k)+f(j+1,k))-0.5_8* &
        ( R1(j,k)*( Q(lambda(j,1)+gamma(j,1),del(k))*alpha(j,1)-fc(j,1)-fc(j+1,1) ) &
         +R2(j,k)*( Q(lambda(j,2)+gamma(j,2),del(k))*alpha(j,2)-fc(j,2)-fc(j+1,2) ) &
	 +R3(j,k)*( Q(lambda(j,3)+gamma(j,3),del(k))*alpha(j,3)-fc(j,3)-fc(j+1,3) ) &
         +R4(j,k)*( Q(lambda(j,4)+gamma(j,4),del(k))*alpha(j,4)-fc(j,4)-fc(j+1,4) ) &  
         +R5(j,k)*( Q(lambda(j,5)+gamma(j,5),del(k))*alpha(j,5)-fc(j,5)-fc(j+1,5) ) &
         +R6(j,k)*( Q(lambda(j,6)+gamma(j,6),del(k))*alpha(j,6)-fc(j,6)-fc(j+1,6) ) &  
         +R7(j,k)*( Q(lambda(j,7)+gamma(j,7),del(k))*alpha(j,7)-fc(j,7)-fc(j+1,7) ) )
  else if(flux_correction==0)then
   flux(j,k)=0.5_8*(f(j,k)+f(j+1,k))-0.5_8*&
        ( R1(j,k)*( abs(lambda(j,1))*alpha(j,1) ) &
         +R2(j,k)*( abs(lambda(j,2))*alpha(j,2) ) &
	 +R3(j,k)*( abs(lambda(j,3))*alpha(j,3) ) &
         +R4(j,k)*( abs(lambda(j,4))*alpha(j,4) ) &  
         +R5(j,k)*( abs(lambda(j,5))*alpha(j,5) ) &
         +R6(j,k)*( abs(lambda(j,6))*alpha(j,6) ) &  
         +R7(j,k)*( abs(lambda(j,7))*alpha(j,7)) ) 
  end if
 end do

end do

end subroutine computeFlux


subroutine bound(u2)

!Input Variables
real*8::u2(-1:nx+2,7)

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

function Q(x,del0) result(qx)

!Input Variables  
real*8,intent(in)::x,del0
!Output Variables 
real*8::qx
 
if(abs(x)<del0)then
  qx=(x*x+del0*del0)/(2.*del0)
else if(abs(x)>=del0)then
  qx=abs(x)
end if

end function Q

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

function minmod(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z
  
  z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod


end module RoeSolverMHD_mod
