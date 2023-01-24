module WENOflux_mod
use globalVariables_mod
implicit none

contains

subroutine compute_flux(q,flux)
!Subroutine I/O variables
real::q(-3:nx+4,8),flux(0:nx,8)

!Subroutine local variables
real::lambda(-3:nx+4,7)
real::R(-3:nx+4,7,7),L(-3:nx+4,7,7)

real::smax,cs,cf,cA
real::rhot,ut,vt,wt,pt,byt,bzt,at,ptotL,ptotR,ptot,bxx,byy,bzz,att
real::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
real::bxLR,BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
real::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state

real::F(-3:nx+4,7)
real::Fs_plus(-3:nx+4,7,-2:3),Fs_minus(-3:nx+4,7,-2:3)!,Fs(-3:nx+4,7,-2:3)

real::a1,b1,c1,d1
real::a2,b2,c2,d2

real::temp1,temp2
integer::i,j,k,s


do i=-2,nx+3
  rhoL=q(i,1)
  uL=q(i,2)/q(i,1)
  vL=q(i,3)/q(i,1)
  wL=q(i,4)/q(i,1) 
  bxL=q(i,5)
  byL=q(i,6)
  bzL=q(i,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(q(i,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  rhoR=q(i+1,1)
  uR=q(i+1,2)/q(i+1,1)
  vR=q(i+1,3)/q(i+1,1)
  wR=q(i+1,4)/q(i+1,1)
  bxR=q(i+1,5)
  byR=q(i+1,6)
  bzR=q(i+1,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(q(i+1,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR

  !-----------------------------
  bxLR=bxL
  !-----------------------------

  !----------------------------------------------
  !Averaged State Variables (arithmetic mean):
  !----------------------------------------------
  rhot=0.5*(rhoL+rhoR)
  ut=0.5*(uL+uR)
  vt=0.5*(vL+vR)
  wt=0.5*(wL+wR)
  byt=0.5*(byL+byR)
  bzt=0.5*(bzL+bzR)
  ptot=0.5*(ptotL+ptotR) !total pressure: fluid+magnetic field
  pt=ptot-0.5*(bxLR*bxLR+byt*byt+bzt*bzt)
  att=gam*pt/rhot
  at=sqrt(att)

  !-------------------------------------------------
  !Alfven and magnetosonic wave speeds: c_s,c_A,c_f
  !-------------------------------------------------
  bxx=(bxLR*bxLR)/rhot
  byy=(byt*byt)/rhot
  bzz=(bzt*bzt)/rhot
  bsqr=bxx+byy+bzz
  vsq=ut*ut+vt*vt+wt*wt

  cA=bxx
  cA=sqrt(bxx)
  cs=0.5*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4.*att*bxx))
  !cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
  cf=att+bsqr-cs
  cs=sqrt(cs)
  cf=sqrt(cf)

  !---------------------------------------------
  !Characteristic Eigenvalues: lambda_1..7
  !---------------------------------------------
  lambda(i,1)=ut-cf
  lambda(i,2)=ut-cA
  lambda(i,3)=ut-cs
  lambda(i,4)=ut
  lambda(i,5)=ut+cs
  lambda(i,6)=ut+cA
  lambda(i,7)=ut+cf

  !-------------------------------------------------
  !Average Jacobian Matrix Eigenvectors: R_1..7, L_1..7
  !-------------------------------------------------
  !if(abs(byt)<tol .and. abs(bzt)<tol)then
  if(byt**2+bzt**2<tol)then
    betay=1./sqrt(2.)
    betaz=1./sqrt(2.)
  else
    betay=byt/sqrt(byt*byt+bzt*bzt)
    betaz=bzt/sqrt(byt*byt+bzt*bzt)
  end if

  !if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
  if(byt**2+bzt**2<tol .or. abs(att-cA**2)<tol)then
    alphas=1.!1./sqrt(2.)!1.
    alphaf=1.!1./sqrt(2.)!1.
  else
  !!***NOTE**** Inside the square roots in renormalization constants, 
  !! can take abs of numerator to avoid negative values caused by numerical precision errors
    alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
    alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs))
    !alphaf=sqrt(abs(at*at-cs*cs)/(cf*cf-cs*cs)) 
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
  theta2=(alphaf*alphaf*cf*at+alphas*alphas*cs*cA)*sgn(bxLR)

  R(i,1,1)=alphaf
  R(i,1,2)=alphaf*(ut-cf)
  R(i,1,3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
  R(i,1,4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
  R(i,1,5)=(alphas*betay*cf)/sqrt(rhot)
  R(i,1,6)=(alphas*betaz*cf)/sqrt(rhot)
  R(i,1,7)=alphaf*0.5*vsq+gf1

  R(i,2,1)=0.
  R(i,2,2)=0.
  R(i,2,3)=betaz*sgn(bxLR)
  R(i,2,4)=-betay*sgn(bxLR)
  R(i,2,5)=betaz/sqrt(rhot)
  R(i,2,6)=-betay/sqrt(rhot)
  R(i,2,7)=(betaz*vt-betay*wt)*sgn(bxLR)

  R(i,3,1)=alphas
  R(i,3,2)=alphas*(ut-cs)
  R(i,3,3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
  R(i,3,4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
  R(i,3,5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
  R(i,3,6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
  R(i,3,7)=alphas*0.5*vsq+gs3

  R(i,4,1)=1.
  R(i,4,2)=ut
  R(i,4,3)=vt
  R(i,4,4)=wt
  R(i,4,5)=0.
  R(i,4,6)=0.
  R(i,4,7)=0.5*vsq

  R(i,5,1)=alphas
  R(i,5,2)=alphas*(ut+cs)
  R(i,5,3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
  R(i,5,4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
  R(i,5,5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
  R(i,5,6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
  R(i,5,7)=alphas*0.5*vsq+gs5

  R(i,6,1)=0.
  R(i,6,2)=0.
  R(i,6,3)=-betaz*sgn(bxLR)
  R(i,6,4)=betay*sgn(bxLR)
  R(i,6,5)=betaz/sqrt(rhot)
  R(i,6,6)=-betay/sqrt(rhot)
  R(i,6,7)=-(betaz*vt-betay*wt)*sgn(bxLR)

  R(i,7,1)=alphaf
  R(i,7,2)=alphaf*(ut+cf)
  R(i,7,3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
  R(i,7,4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
  R(i,7,5)=(alphas*betay*cf)/sqrt(rhot)
  R(i,7,6)=(alphas*betaz*cf)/sqrt(rhot)
  R(i,7,7)=alphaf*0.5*vsq+gf7


  L(i,1,1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

  L(i,1,2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

  L(i,1,3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

  L(i,1,4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

  L(i,1,5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,1,6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,1,7)=(1./theta1)*0.5*alphaf*at*at

  L(i,2,1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
  L(i,2,2)=0.
  L(i,2,3)=0.5*betaz*sgn(bxLR)
  L(i,2,4)=-0.5*betay*sgn(bxLR)
  L(i,2,5)=0.5*betaz*sqrt(rhot)
  L(i,2,6)=-0.5*betay*sqrt(rhot)
  L(i,2,7)=0.

  L(i,3,1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

  L(i,3,2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

  L(i,3,3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

  L(i,3,4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

  L(i,3,5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,3,6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,3,7)=(1./theta1)*0.5*alphas*cf*cf


  L(i,4,1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

  L(i,4,2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

  L(i,4,3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

  L(i,4,4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

  L(i,4,5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

  L(i,4,6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

  L(i,4,7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


  L(i,5,1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

  L(i,5,2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

  L(i,5,3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

  L(i,5,4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

  L(i,5,5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,5,6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,5,7)=(1./theta1)*0.5*alphas*cf*cf


  L(i,6,1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
  L(i,6,2)=0.
  L(i,6,3)=-0.5*betaz*sgn(bxLR)
  L(i,6,4)=0.5*betay*sgn(bxLR)
  L(i,6,5)=0.5*betaz*sqrt(rhot)
  L(i,6,6)=-0.5*betay*sqrt(rhot)
  L(i,6,7)=0.


  L(i,7,1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

  L(i,7,2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

  L(i,7,3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

  L(i,7,4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

  L(i,7,5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,7,6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

  L(i,7,7)=(1./theta1)*0.5*alphaf*at*at

  if((at*at-cA*cA)>tol)then
   if(byt**2+bzt**2<tol)then
    do j=1,7
      R(i,3,j)=-1.*R(i,3,j)
      R(i,5,j)=-1.*R(i,5,j)
      L(i,3,j)=-1.*L(i,3,j)
      L(i,5,j)=-1.*L(i,5,j)
    end do
   end if
  else 
   if(byt**2+bzt**2<tol)then
    do j=1,7
      R(i,1,j)=-1.*R(i,1,j)
      R(i,7,j)=-1.*R(i,7,j)
      L(i,1,j)=-1.*L(i,1,j)
      L(i,7,j)=-1.*L(i,7,j)
    end do
   end if
  end if

end do

!----------------------------------------------
!Physical Flux: F_i
!----------------------------------------------
do i=-3,nx+4
  rhoL=q(i,1)
  uL=q(i,2)/q(i,1)
  vL=q(i,3)/q(i,1)
  wL=q(i,4)/q(i,1) 
  bxL=q(i,5)
  byL=q(i,6)
  bzL=q(i,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(q(i,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  F(i,1)=rhoL*uL
  F(i,2)=rhoL*uL*uL+ptotL-bxL*bxL
  F(i,3)=rhoL*uL*vL-bxL*byL
  F(i,4)=rhoL*uL*wL-bxL*bzL
  F(i,5)=uL*byL-vL*bxL
  F(i,6)=uL*bzL-wL*bxL
  F(i,7)=uL*(q(i,8)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)

end do

!-----------------------------------------------------------------
!Physical Flux Projection onto Right eigenvectors: F^s_k
!-----------------------------------------------------------------
do i=0,nx
  do s=1,7
    !Max wave speed for local lax-friedrichs split flux
    alphas=0.
    do k=-2,3
      alphas=max(alphas,abs(lambda(i+k,s)))
    end do

    do k=-2,3
      temp1=0.
      temp2=0. 
      do j=1,7
        temp1=temp1+L(i,s,j)*F(i+k,j)
        temp2=temp2+L(i,s,j)*q(i+k,j)
      end do
      !Fs(i,s,k)=temp1
      Fs_plus(i,s,k)=0.5*(temp1+alphas*temp2)
      Fs_minus(i,s,k)=0.5*(temp1-alphas*temp2)
    end do
  end do
end do

!-----------------------------------------------------------------
!WENO Numerical Flux: F^_i+1/2
!-----------------------------------------------------------------
do i=0,nx
  do j=1,7
    flux(i,j)=(1./12.)*(-F(i-1,j)+7.*F(i,j)+7.*F(i+1,j)-F(i+2,j))
    do s=1,7
      a1=Fs_plus(i,s,-1)-Fs_plus(i,s,-2)
      b1=Fs_plus(i,s,0)-Fs_plus(i,s,-1)
      c1=Fs_plus(i,s,1)-Fs_plus(i,s,0)
      d1=Fs_plus(i,s,2)-Fs_plus(i,s,1)
    
      a2=Fs_minus(i,s,3)-Fs_minus(i,s,2)
      b2=Fs_minus(i,s,2)-Fs_minus(i,s,1)
      c2=Fs_minus(i,s,1)-Fs_minus(i,s,0)
      d2=Fs_minus(i,s,0)-Fs_minus(i,s,-1)
    
      flux(i,j)=flux(i,j)+(-phi(a1,b1,c1,d1)+phi(a2,b2,c2,d2))*R(i,s,j)
    end do
  end do
  flux(i,8)=flux(i,7)
  flux(i,7)=flux(i,6)
  flux(i,6)=flux(i,5)
  flux(i,5)=0.
end do

end subroutine compute_flux

subroutine timeStep(q,dt)
!Subroutine I/O variables
real::q(-3:nx+4,8),dt
!Subroutine local variables
integer::i,j

!Subroutine local variables
real::smax,cs,cf,cA
real::rhot,ut,vt,wt,pt,byt,bzt,at,ptotL,ptotR,ptot,bxx,byy,bzz,att

real::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
real::bxLR,BBR,BBL,vvR,vvL,vsq
real::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state

smax=0.

do i=0,nx+1
  rhoL=q(i,1)
  uL=q(i,2)/q(i,1)
  vL=q(i,3)/q(i,1)
  wL=q(i,4)/q(i,1) 
  bxL=q(i,5)
  byL=q(i,6)
  bzL=q(i,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=(gam-1.)*(q(i,8)-0.5*rhoL*vvL-0.5*BBL)
  aL=sqrt(gam*pL/rhoL)
  ptotL=pL+0.5*BBL

  rhoR=q(i+1,1)
  uR=q(i+1,2)/q(i+1,1)
  vR=q(i+1,3)/q(i+1,1)
  wR=q(i+1,4)/q(i+1,1)
  bxR=q(i+1,5)
  byR=q(i+1,6)
  bzR=q(i+1,7)
  BBR=bxR*bxR+byR*byR+bzR*bzR
  vvR=uR*uR+vR*vR+wR*wR
  pR=(gam-1.)*(q(i+1,8)-0.5*rhoR*vvR-0.5*BBR)
  aR=sqrt(gam*pR/rhoR)
  ptotR=pR+0.5*BBR

  !-----------------------------
  bxLR=bxL
  !-----------------------------

  !----------------------------------------------
  !Averaged State Variables (arithmetic mean):
  !----------------------------------------------
  rhot=0.5*(rhoL+rhoR)
  ut=0.5*(uL+uR)
  vt=0.5*(vL+vR)
  wt=0.5*(wL+wR)
  byt=0.5*(byL+byR)
  bzt=0.5*(bzL+bzR)
  ptot=0.5*(ptotL+ptotR) !total pressure: fluid+magnetic field
  pt=ptot-0.5*(bxLR*bxLR+byt*byt+bzt*bzt)
  att=gam*pt/rhot
  at=sqrt(att)

  !-------------------------------------------------
  !Alfven and magnetosonic wave speeds: c_s,c_A,c_f
  !-------------------------------------------------
  bxx=(bxLR*bxLR)/rhot
  byy=(byt*byt)/rhot
  bzz=(bzt*bzt)/rhot
  bsqr=bxx+byy+bzz
  vsq=ut*ut+vt*vt+wt*wt

  cA=bxx
  cA=sqrt(bxx)
  cs=0.5*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4.*att*bxx))
  !cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
  cf=att+bsqr-cs
  cs=sqrt(cs)
  cf=sqrt(cf)

  !---------------------------------------------
  !Characteristic Eigenvalues: lambda_1..7
  !---------------------------------------------
  !lambda(i,1)=ut-cf
  !lambda(i,2)=ut-cA
  !lambda(i,3)=ut-cs
  !lambda(i,4)=ut
  !lambda(i,5)=ut+cs
  !lambda(i,6)=ut+cA
  !lambda(i,7)=ut+cf

  smax=max(smax,abs(ut)+abs(cf))
  
end do


dt=cour*dx/smax


end subroutine timeStep


function phi(a,b,c,d) result(phix)
!Inputs
real,intent(in)::a,b,c,d
!Outputs
real::phix
!Local Variables
real::w0,w2,alpha0,alpha1,alpha2,alphasum,IS0,IS1,IS2

IS0=13.*(a-b)**2+3.*(a-3.*b)**2
IS1=13.*(b-c)**2+3.*(b+c)**2
IS2=13.*(c-d)**2+3.*(3.*c-d)**2.

alpha0=1./((eps+IS0)**2)
alpha1=6./((eps+IS1)**2)
alpha2=3./((eps+IS2)**2)

alphasum=alpha0+alpha1+alpha2

w0=alpha0/alphasum
w2=alpha2/alphasum

phix=(1./3.)*w0*(a-2.*b+c)+(1./6.)*(w2-0.5)*(b-2.*c+d)

end function phi

function sgn(x) result(fx)
  real,intent(in)::x
  real::fx
  if(x>tol)then
    fx=1.0!sign(1.0,x)
  else
    fx=-1.0
  end if
end

end module WENOflux_mod
