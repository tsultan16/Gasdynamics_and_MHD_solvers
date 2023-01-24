module WENOflux_mod
use globalVariables_mod
implicit none

contains

subroutine compute_flux(q,flux)
!Subroutine I/O variables
real::q(-3:nx+4,8),flux(0:nx,8)

!Subroutine local variables
real,allocatable::lambda(:,:)!lambda(-3:nx+4,7)
real,allocatable::R(:,:,:),L(:,:,:)!R(-3:nx+4,7,7),L(-3:nx+4,7,7)

real::smax,cs,cf,cA
real::rhot,ut,vt,wt,pt,byt,bzt,at,ptotL,ptotR,ptot,bxx,byy,bzz,att
real::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
real::bxLR,BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq,gama,gamf,gams,tau
real::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state

real,allocatable::F(:,:)!F(-3:nx+4,7)
real,allocatable::Fs_plus(:,:,:),Fs_minus(:,:,:),Fs(:,:,:)!Fs_plus(-3:nx+4,7,-2:3),Fs_minus(-3:nx+4,7,-2:3),Fs(-3:nx+4,7,-2:3)

real::a1,b1,c1,d1
real::a2,b2,c2,d2

real::temp1,temp2
integer::i,j,k,s


allocate(lambda(-3:nx+4,7),R(-3:nx+4,7,7),L(-3:nx+4,7,7),F(-3:nx+4,7))
allocate(Fs_plus(-3:nx+4,7,-2:3),Fs_minus(-3:nx+4,7,-2:3),Fs(-3:nx+4,7,-2:3))

do i=-3,nx+3

  if(debug==1)then
    print*,'Cell =',i
  end if
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
  
  !pt=ptot-0.5*(bxLR*bxLR+byt*byt+bzt*bzt) !first option for computing p_avg
  pt=0.5*(pL+pR) !second option for computing p_avg
  att=max(tol,gam*pt/rhot)
  at=sqrt(att)


  if(debug==1)then
    print*,'Roe Average(rho,vx,vy,vz,bx,by,bz,p) =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
  end if

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
  temp1=max(tol,(att+bsqr)**2-4.*att*bxx)
  cs=0.5*(att+bsqr-sqrt(temp1))
  cf=0.5*(att+bsqr+sqrt(temp1))
  !cf=att+bsqr-cs
  cs=sqrt(max(tol,cs))
  cf=sqrt(max(tol,cf))

  if(debug==1)then
    print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
  end if

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

  if(debug==1)then
    print*,'betay,betaz=',betay,betaz
  end if

  !if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
  if(sqrt(cf**2-cs**2)<tol)then
  !if(byt**2+bzt**2<tol .or. abs(att-cA**2)<tol)then
    alphas=1.!/sqrt(2.)!1.
    alphaf=1.!/sqrt(2.)!1.
  else
  !!***NOTE**** 
  !!For zero magnetic field, set alphas/f to zero
    if(abs(cf-at)<tol)then
      alphas=0.0
    else 
      alphas=sqrt(max(0.0,cf**2-att)/(cf**2-cs**2))
    end if
   if(abs(cs-at)<tol)then
      alphaf=0.0  
    else 
      alphaf=sqrt(max(0.0,att-cs**2)/(cf**2-cs**2))
    end if 
  end if

  if(debug==1)then
    print*,'alphas,alphaf=',alphas,alphaf
  end if

  gama=sgn(bxLR)*(betaz*vt-betay*wt)
  gamf=alphaf*cf*ut-alphas*cs*sgn(bxLR)*(betay*vt+betaz*wt)
  gams=alphas*cs*ut+alphaf*cf*sgn(bxLR)*(betay*vt+betaz*wt)
  tau=(gam-1.)/att

  R(i,1,1)=alphaf
  R(i,1,2)=alphaf*(ut-cf)
  R(i,1,3)=alphaf*vt+cs*alphas*betay*sgn(bxLR)
  R(i,1,4)=alphaf*wt+cs*alphas*betaz*sgn(bxLR)
  R(i,1,5)=(alphas*betay*at)/sqrt(rhot)
  R(i,1,6)=(alphas*betaz*at)/sqrt(rhot)
  R(i,1,7)=alphaf*(0.5*vsq+cf**2-gam2*att)-gamf

  R(i,2,1)=0.0
  R(i,2,2)=0.0
  R(i,2,3)=-betaz*sgn(bxLR)
  R(i,2,4)=betay*sgn(bxLR)
  R(i,2,5)=-betaz/sqrt(rhot)
  R(i,2,6)=betay/sqrt(rhot)
  R(i,2,7)=-gama

  R(i,3,1)=alphas
  R(i,3,2)=alphas*(ut-cs)
  R(i,3,3)=alphas*vt-cf*alphaf*betay*sgn(bxLR)
  R(i,3,4)=alphas*wt-cf*alphaf*betaz*sgn(bxLR)
  R(i,3,5)=-(at*alphaf*betay)/sqrt(rhot)
  R(i,3,6)=-(at*alphaf*betaz)/sqrt(rhot)
  R(i,3,7)=alphas*(0.5*vsq+cs**2-gam2*att)-gams

  R(i,4,1)=1.
  R(i,4,2)=ut
  R(i,4,3)=vt
  R(i,4,4)=wt
  R(i,4,5)=0.0
  R(i,4,6)=0.0
  R(i,4,7)=0.5*vsq

  R(i,5,1)=alphas
  R(i,5,2)=alphas*(ut+cs)
  R(i,5,3)=alphas*vt+cf*alphaf*betay*sgn(bxLR)
  R(i,5,4)=alphas*wt+cf*alphaf*betaz*sgn(bxLR)
  R(i,5,5)=-(at*alphaf*betay)/sqrt(rhot)
  R(i,5,6)=-(at*alphaf*betaz)/sqrt(rhot)
  R(i,5,7)=alphas*(0.5*vsq+cs**2-gam2*att)+gams

  R(i,6,1)=0.0
  R(i,6,2)=0.0
  R(i,6,3)=-betaz*sgn(bxLR)
  R(i,6,4)=betay*sgn(bxLR)
  R(i,6,5)=betaz/sqrt(rhot)
  R(i,6,6)=-betay/sqrt(rhot)
  R(i,6,7)=-gama

  R(i,7,1)=alphaf
  R(i,7,2)=alphaf*(ut+cf)
  R(i,7,3)=alphaf*vt-cs*alphas*betay*sgn(bxLR)
  R(i,7,4)=alphaf*wt-cs*alphas*betaz*sgn(bxLR)
  R(i,7,5)=(alphas*betay*at)/sqrt(rhot)
  R(i,7,6)=(alphas*betaz*at)/sqrt(rhot)
  R(i,7,7)=alphaf*(0.5*vsq+cf**2-gam2*att)+gamf


  L(i,1,1)=(1./(2.*att))*(gam1*alphaf*vsq+gamf)

  L(i,1,2)=(1./(2.*att))*((1.-gam)*alphaf*ut-alphaf*cf)

  L(i,1,3)=(1./(2.*att))*((1.-gam)*alphaf*vt+cs*alphas*betay*sgn(bxLR))

  L(i,1,4)=(1./(2.*att))*((1.-gam)*alphaf*wt+cs*alphas*betaz*sgn(bxLR))

  !L(i,1,5)=(1./(2.*att))*((1.-gam)*alphaf*byt-sqrt(rhot)*at*alphas*betay)
  L(i,1,5)=(1./(2.*att))*((1.-gam)*alphaf*byt+sqrt(rhot)*at*alphas*betay)

  !L(i,1,6)=(1./(2.*att))*((1.-gam)*alphaf*bzt-sqrt(rhot)*at*alphas*betaz)
  L(i,1,6)=(1./(2.*att))*((1.-gam)*alphaf*bzt+sqrt(rhot)*at*alphas*betaz)


  L(i,1,7)=(1./(2.*att))*(gam-1.)*alphaf

  L(i,2,1)=0.5*gama
  L(i,2,2)=0.0
  L(i,2,3)=-0.5*betaz*sgn(bxLR)
  L(i,2,4)=0.5*betay*sgn(bxLR)
  L(i,2,5)=-0.5*betaz*sqrt(rhot)
  L(i,2,6)=0.5*betay*sqrt(rhot)
  L(i,2,7)=0.0


  L(i,3,1)=(1./(2.*att))*(gam1*alphas*vsq+gams)

  L(i,3,2)=(1./(2.*att))*((1.-gam)*alphas*ut-alphas*cs)

  L(i,3,3)=(1./(2.*att))*((1.-gam)*alphas*vt-cf*alphaf*betay*sgn(bxLR))

  L(i,3,4)=(1./(2.*att))*((1.-gam)*alphas*wt-cf*alphaf*betaz*sgn(bxLR))

  L(i,3,5)=(1./(2.*att))*((1.-gam)*alphas*byt-sqrt(rhot)*at*alphaf*betay)

  L(i,3,6)=(1./(2.*att))*((1.-gam)*alphas*bzt-sqrt(rhot)*at*alphaf*betaz)

  L(i,3,7)=(1./(2.*att))*(gam-1.)*alphas


  L(i,4,1)=1.-0.5*tau*vsq
  L(i,4,2)=tau*ut
  L(i,4,3)=tau*vt
  L(i,4,4)=tau*wt
  L(i,4,5)=tau*byt
  L(i,4,6)=tau*bzt
  L(i,4,7)=-tau


  L(i,5,1)=(1./(2.*att))*(gam1*alphas*vsq-gams)

  L(i,5,2)=(1./(2.*att))*((1.-gam)*alphas*ut+alphas*cs)

  L(i,5,3)=(1./(2.*att))*((1.-gam)*alphas*vt+cf*alphaf*betay*sgn(bxLR))

  L(i,5,4)=(1./(2.*att))*((1.-gam)*alphas*wt+cf*alphaf*betaz*sgn(bxLR))

  L(i,5,5)=(1./(2.*att))*((1.-gam)*alphas*byt-sqrt(rhot)*at*alphaf*betay)

  L(i,5,6)=(1./(2.*att))*((1.-gam)*alphas*bzt-sqrt(rhot)*at*alphaf*betaz)

  L(i,5,7)=(1./(2.*att))*(gam-1.)*alphas


  L(i,6,1)=0.5*gama
  L(i,6,2)=0.0
  L(i,6,3)=-0.5*betaz*sgn(bxLR)
  L(i,6,4)=0.5*betay*sgn(bxLR)
  L(i,6,5)=0.5*betaz*sqrt(rhot)
  L(i,6,6)=-0.5*betay*sqrt(rhot)
  L(i,6,7)=0.0


  L(i,7,1)=(1./(2.*att))*(gam1*alphaf*vsq-gamf)

  L(i,7,2)=(1./(2.*att))*((1.-gam)*alphaf*ut+alphaf*cf)

  L(i,7,3)=(1./(2.*att))*((1.-gam)*alphaf*vt-cs*alphas*betay*sgn(bxLR))

  L(i,7,4)=(1./(2.*att))*((1.-gam)*alphaf*wt-cs*alphas*betaz*sgn(bxLR))

  !L(i,7,5)=(1./(2.*att))*((1.-gam)*alphaf*byt-sqrt(rhot)*at*alphas*betay) !typo in chiang-wu paper
  L(i,7,5)=(1./(2.*att))*((1.-gam)*alphaf*byt+sqrt(rhot)*at*alphas*betay)

  !L(i,7,6)=(1./(2.*att))*((1.-gam)*alphaf*bzt-sqrt(rhot)*at*alphas*betaz) !typo in chiang-wu paper
  L(i,7,6)=(1./(2.*att))*((1.-gam)*alphaf*bzt+sqrt(rhot)*at*alphas*betaz)

  L(i,7,7)=(1./(2.*att))*(gam-1.)*alphaf

  !go to 111
  if((att-cA*cA)>tol)then
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
  !111 continue


end do

!----------------------------------------------
!Physical Flux: F_i
!----------------------------------------------
do i=-3,nx+4
  rhoL=max(tol,q(i,1))
  uL=q(i,2)/q(i,1)
  vL=q(i,3)/q(i,1)
  wL=q(i,4)/q(i,1) 
  bxL=q(i,5)
  byL=q(i,6)
  bzL=q(i,7)
  BBL=bxL*bxL+byL*byL+bzL*bzL
  vvL=uL*uL+vL*vL+wL*wL
  pL=max(tol,(gam-1.)*(q(i,8)-0.5*rhoL*vvL-0.5*BBL))
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
  do k=-2,3 
    do s=1,7
      !Max wave speed for local lax-friedrichs split flux
      !alphas=0.
      !do k=-2,3
      !  alphas=max(alphas,abs(lambda(i+k,s)))
      !end do
      alphas=max(abs(lambda(i,s)),abs(lambda(i+1,s)))

      temp1=L(i,s,1)*F(i+k,1)+L(i,s,2)*F(i+k,2)+L(i,s,3)*F(i+k,3)&
           +L(i,s,4)*F(i+k,4)+L(i,s,5)*F(i+k,5)+L(i,s,6)*F(i+k,6)+L(i,s,7)*F(i+k,7)
      temp2=L(i,s,1)*q(i+k,1)+L(i,s,2)*q(i+k,2)+L(i,s,3)*q(i+k,3)+L(i,s,4)*q(i+k,4)&
           +L(i,s,5)*q(i+k,6)+L(i,s,6)*q(i+k,7)+L(i,s,7)*q(i+k,8)
      
      Fs(i,s,k)=temp1
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
    !flux(i,j)=(1./12.)*(-F(i-1,j)+7.*F(i,j)+7.*F(i+1,j)-F(i+2,j))
    flux(i,j)=0.0
    do s=1,7
      a1=Fs_plus(i,s,-1)-Fs_plus(i,s,-2)
      b1=Fs_plus(i,s,0)-Fs_plus(i,s,-1)
      c1=Fs_plus(i,s,1)-Fs_plus(i,s,0)
      d1=Fs_plus(i,s,2)-Fs_plus(i,s,1)
    
      a2=Fs_minus(i,s,3)-Fs_minus(i,s,2)
      b2=Fs_minus(i,s,2)-Fs_minus(i,s,1)
      c2=Fs_minus(i,s,1)-Fs_minus(i,s,0)
      d2=Fs_minus(i,s,0)-Fs_minus(i,s,-1)
    
      !flux(i,j)=flux(i,j)+(-phi(a1,b1,c1,d1)+phi(a2,b2,c2,d2))*R(i,s,j)
      flux(i,j)=flux(i,j)+((1./12.)*(-Fs(i,s,-1)+7.*Fs(i,s,0)+7.*Fs(i,s,1)-Fs(i,s,2))&
               -phi(a1,b1,c1,d1)+phi(a2,b2,c2,d2))*R(i,s,j)
    end do
  end do
  flux(i,8)=flux(i,7)
  flux(i,7)=flux(i,6)
  flux(i,6)=flux(i,5)
  flux(i,5)=0.
end do

deallocate(lambda,R,L,F)
deallocate(Fs_plus,Fs_minus,Fs)

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
real::temp1

smax=0.

do i=0,nx
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
  
  !pt=ptot-0.5*(bxLR*bxLR+byt*byt+bzt*bzt) !first option for computing p_avg
  pt=0.5*(pL+pR) !second option for computing p_avg
  att=max(tol,gam*pt/rhot)
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
  temp1=max(0.0,(att+bsqr)**2-4.*att*bxx)
  cs=0.5*(att+bsqr-sqrt(temp1))
  cf=0.5*(att+bsqr+sqrt(temp1))
  !cf=att+bsqr-cs
  cs=sqrt(max(tol,cs))
  cf=sqrt(max(tol,cf))

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
IS2=13.*(c-d)**2+3.*(3.*c-d)**2

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
