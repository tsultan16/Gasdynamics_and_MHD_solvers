!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_mod
implicit none

integer,parameter::flux_option=3 !1:mh 2:harten 3:ENO

integer,parameter::debug=0
integer,parameter::shockSwitch=1 !0:off, 1:on
integer,parameter ::nt=100
integer,parameter ::nx=400 
real*8,parameter::del0=0.5
real*8,parameter::Q_user=2.0 
integer,parameter::characteristicLimiting=0 !0:off, 1:on 
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,aL,HL,FLx(5) !left state
real*8::rhoR,uR,vR,wR,pR,aR,HR,FRx(5) !right state 
real*8::delrho,delp,delu,delv,delw,smax

!----------------------------------------------------------------------
real*8::gt(-1:nx+2,5),fc(-1:nx+2,5),gamma(-1:nx+2,5),sigma(-1:nx+2,5)
real*8::f(-1:nx+2,5),g(-1:nx+2,5),fl(5),u(5)

real*8::dens,velx,vely,velz,pres,a,H
real*8::rhot,ut,vt,wt,Ht,at,lambda(-1:nx+2,5),alpha(-1:nx+2,5)
real*8::K1(-1:nx+2,5),K2(-1:nx+2,5),K3(-1:nx+2,5)
real*8::K4(-1:nx+2,5),K5(-1:nx+2,5)	
real*8::u1(-1:nx+2,5),u2(-1:nx+2,5)

real*8 ::temp1,temp2,temp3,xc,xt,theta,w
!------------------------------------------------------------------------------
real*8::lambda1L(-1:nx+2),lambda1R(-1:nx+2),lambda5L(-1:nx+2),lambda5R(-1:nx+2)
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::fright(-1:nx+2,5),fleft(-1:nx+2,5)
real*8::qL(-1:nx+2,5),qR(-1:nx+2,5),slope(-1:nx+2,5)
real*8::delta1(-1:nx+2,5),delta2(-1:nx+2,5),delta3(-1:nx+2,5)
real*8::delta4(-1:nx+2,5),delta5(-1:nx+2,5)
real*8::ubL(-1:nx+2,5),ubR(-1:nx+2,5)

!------------------------------------------------------------------------------
real*8::qh(-1:nx+2,5) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::v(-1:nx+2,5) !characteristic variables
real*8::Qx(-1:nx+2,5,5),Qinv(-1:nx+2,5,5) !characeteristic matrix
real*8::S(-1:nx+2,5)
real*8::dqdx(-1:nx+2,5),dqdt(-1:nx+2,5)
real*8::dpdx,dudx,dpdt,dudt

!------------------------------------------------------------------------------


contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
if(flux_option==1)then
  open(unit=11,file='output_EulerMH.txt')
else if(flux_option==2)then
  open(unit=11,file='output_EulerHarten.txt')
else if(flux_option==3)then
  open(unit=11,file='output_EulerENO.txt')
end if

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

read(unit=10,fmt=*) rhoL,uL,vL,wL,pL 
read(unit=10,fmt=*) rhoR,uR,vR,wR,pR

do j=-1,nx+2
  if(j<nx/2) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=rhoL*vL
    u1(j,4)=rhoL*wL
    u1(j,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=rhoR*vR
    u1(j,4)=rhoR*wR
    u1(j,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul
  end if
end do

print*,'Left State (rho,u,v,w,p) =',rhoL,uL,vL,wL,pL
print*,'Right State (rho,u,v,w,p) =',rhoR,uR,vR,wR,pR

end subroutine init

subroutine RoeSolver(jj)

!Input Variables
integer::jj

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delw=wR-wL
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
vt=(sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
wt=(sqrt(rhoL)*wL+sqrt(rhoR)*wR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*(ut**2.+vt**2.+wt**2.)) )

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(jj,1)=ut-at
lambda(jj,2)=ut
lambda(jj,3)=ut
lambda(jj,4)=ut
lambda(jj,5)=ut+at

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(jj,1)=1._8
K1(jj,2)=ut-at
K1(jj,3)=vt
K1(jj,4)=wt
K1(jj,3)=Ht-ut*at

K2(jj,1)=1
K2(jj,2)=ut
K2(jj,3)=vt
K2(jj,4)=wt
K2(jj,3)=0.5*(ut**2.+vt**2.+wt**2.)

K3(jj,1)=0.
K3(jj,2)=0.
K3(jj,3)=1.
K3(jj,4)=0.
K3(jj,5)=vt

K4(jj,1)=0.
K4(jj,2)=0.
K4(jj,3)=0.
K4(jj,4)=1.
K4(jj,5)=wt

K5(jj,1)=1._8
K5(jj,2)=ut+at
K5(jj,3)=vt
K5(jj,4)=wt
K5(jj,5)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(jj,1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(jj,2)=delrho-delp/(at*at)
alpha(jj,3)=rhot*delv
alpha(jj,4)=rhot*delw
alpha(jj,5)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
!smax=abs(ut)+abs(at)
smax=max(smax,abs(lambda(jj,1)),abs(lambda(jj,2)),abs(lambda(jj,5)))


end subroutine RoeSolver

subroutine computeFluxHarten()

do j=-1,nx+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,u1(j,1))
uL=u1(j,2)/u1(j,1)
vL=u1(j,3)/u1(j,1)
wL=u1(j,4)/u1(j,1)
pL=max(premin,(gam-1.)*( u1(j,5)-0.5*rhoL &
  *(uL**2.+vL**2.+wL**2.) ) )
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

rhoR=max(densmin,u1(j+1,1))
uR=u1(j+1,2)/u1(j+1,1)
vR=u1(j+1,3)/u1(j+1,1)
wR=u1(j+1,4)/u1(j+1,1)
pR=max(premin,(gam-1.)*( u1(j+1,5)-0.5*rhoR &
  *(uR**2.+vR**2.+wR**2.) ) )
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR

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
 do k=1,5
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
   do k=1,5
     fc(j,k)=minmod(gt(j,k),gt(j-1,k))
   end do  
   
   if(shockSwitch==1)then
    do k=1,5
     if(abs(alpha(j,k))+abs(alpha(j-1,k)) .ne. 0._8)then
      theta=(abs(alpha(j,k)-alpha(j-1,k)))/(abs(alpha(j,k))+abs(alpha(j-1,k)))
     else
      theta=0.
     end if 
     w=minmod(alpha(j,k)*sigma(j,k),alpha(j-1,k)*sigma(j-1,k))
     
     !Modify the entropy wave only (will result in steeper contact discontinuity)
     if(k==2)then ! .or. k==3 .or. k==4)then	 
       fc(j,k)=fc(j,k)+w*theta
     end if

    end do
   end if
   
 !-----------------------------------------------
 !f_j
 !-----------------------------------------------
 u(1)=u1(j,1)
 u(2)=u1(j,2)/u1(j,1)
 u(3)=u1(j,3)/u1(j,1)
 u(4)=u1(j,4)/u1(j,1)
 u(5)=(gam-1.)*(u1(j,5)-0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.))
 call flux(u)
   do k=1,5
     f(j,k)=fl(k)
   end do
end do

!------------------------------------------------
!gamma_i+1/2
!------------------------------------------------
do j=0,nx  
 do k=1,5
   if(alpha(j,k)>=1.d-30)then
     gamma(j,k)=(fc(j+1,k)-fc(j,k))/alpha(j,k)
   else
     gamma(j,k)=0.
   end if	  
 end do
  
!------------------------------------------------
!Modified Numerical Flux: ^g_i+1/2
!------------------------------------------------
 do k=1,5
!  g(j,k)=0.5*(f(j,k)+f(j+1,k)&
!        +K1(j,k)*(fc(j,1)+fc(j+1,1)-Q(lambda(j,1)+gamma(j,1))*alpha(j,1))&
!        +K2(j,k)*(fc(j,2)+fc(j+1,2)-Q(lambda(j,2)+gamma(j,2))*alpha(j,2))&
!	+K3(j,k)*(fc(j,3)+fc(j+1,3)-Q(lambda(j,3)+gamma(j,3))*alpha(j,3)) )

  g(j,k)=0.5*(f(j,k)+f(j+1,k))-0.5*&
        ( K1(j,k)*( Q(lambda(j,1)+gamma(j,1))*alpha(j,1)-fc(j,1)-fc(j+1,1) ) &
         +K2(j,k)*( Q(lambda(j,2)+gamma(j,2))*alpha(j,2)-fc(j,2)-fc(j+1,2) ) &
	 +K3(j,k)*( Q(lambda(j,3)+gamma(j,3))*alpha(j,3)-fc(j,3)-fc(j+1,3) ) &
         +K4(j,k)*( Q(lambda(j,4)+gamma(j,4))*alpha(j,4)-fc(j,4)-fc(j+1,4) ) &  
         +K5(j,k)*( Q(lambda(j,5)+gamma(j,5))*alpha(j,5)-fc(j,5)-fc(j+1,5) ) ) 


 end do
end do
 
end subroutine computeFluxHarten


subroutine computeFluxMH()

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(characteristicLimiting==1)then

do j=-1,nx+2 
  !set left and right states for local Riemann problem at j+1/2 cell interface 
  rhoL=max(densmin,u1(j,1))
  uL=u1(j,2)/u1(j,1)
  vL=u1(j,3)/u1(j,1)
  wL=u1(j,4)/u1(j,1)
  pL=max(premin,(gam-1.)*( u1(j,5)-0.5*rhoL &
    *(uL**2.+vL**2.+wL**2.) ) )
  aL=sqrt(gam*pL/rhoL)
  HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

  rhoR=max(densmin,u1(j+1,1))
  uR=u1(j+1,2)/u1(j+1,1)
  vR=u1(j+1,3)/u1(j+1,1)
  wR=u1(j+1,4)/u1(j+1,1)
  pR=max(premin,(gam-1.)*( u1(j+1,5)-0.5*rhoR &
    *(uR**2.+vR**2.+wR**2.) ) )
  aR=sqrt(gam*pR/rhoR)
  HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR

  call RoeSolver(j)
  
  do k=1,5
    delta1(j,k)=alpha(j,1)*K1(j,k)
    delta2(j,k)=alpha(j,2)*K2(j,k)
    delta3(j,k)=alpha(j,3)*K3(j,k)
    delta4(j,k)=alpha(j,4)*K4(j,k)
    delta5(j,k)=alpha(j,5)*K5(j,k)
  end do
end do

do j=0,nx+1
  do k=1,5
   slope(j,k)=limiter(delta1(j,k)-delta1(j-1,k) &
               ,delta1(j+1,k)-delta1(j,k)) &
             +limiter(delta2(j,k)-delta2(j-1,k) &
               ,delta2(j+1,k)-delta2(j,k)) &
             +limiter(delta3(j,k)-delta3(j-1,k) &
               ,delta3(j+1,k)-delta3(j,k)) &
             +limiter(delta4(j,k)-delta4(j-1,k) &
               ,delta4(j+1,k)-delta4(j,k)) &
             +limiter(delta5(j,k)-delta5(j-1,k) &
               ,delta5(j+1,k)-delta5(j,k)) 
  end do
end do

else if(characteristicLimiting==0)then

do j=0,nx+1
  do k=1,5
   slope(j,k)=limiter(u1(j,k)-u1(j-1,k),u1(j+1,k)-u1(j,k))
  end do
end do

end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do j=0,nx+1
  do k=1,5
   qL(j,k)=u1(j,k)-0.5*slope(j,k)
   qR(j,k)=u1(j,k)+0.5*slope(j,k)
  end do
end do

!-------------------------------------------------------------------
!Evolve Boudary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------
do j=0,nx+1
  FLx(1)=qL(j,2)
  FLx(2)=(gam-1.)*qL(j,5)+0.5*(3.-gam)*((qL(j,2)**2. &
         +qL(j,3)**2.+qL(j,4)**2.))/qL(j,1)
  FLx(3)=qL(j,2)*qL(j,3)/qL(j,1)
  FLx(4)=qL(j,2)*qL(j,4)/qL(j,1)
  FLx(5)=gam*(qL(j,2)*qL(j,5)/qL(j,1))-0.5*(gam-1.)*qL(j,2)* &
        (qL(j,2)**2.+qL(j,3)**2.+qL(j,4)**2.)/(qL(j,1)**2.)



  FRx(1)=qR(j,2)
  FRx(2)=(gam-1.)*qR(j,5)+0.5*(3.-gam)*((qR(j,2)**2. &
         +qR(j,3)**2.+qR(j,4)**2.))/qR(j,1)
  FRx(3)=qR(j,2)*qR(j,3)/qR(j,1)
  FRx(4)=qR(j,2)*qR(j,4)/qR(j,1)
  FRx(5)=gam*(qR(j,2)*qR(j,5)/qR(j,1))-0.5*(gam-1.)*qR(j,2)* &
        (qR(j,2)**2.+qR(j,3)**2.+qR(j,4)**2.)/(qR(j,1)**2.)

  do k=1,5
    ubL(j,k)=qL(j,k)+0.5*(dt/dx)*(FLx(k)-FRx(k))
    ubR(j,k)=qR(j,k)+0.5*(dt/dx)*(FLx(k)-FRx(k))
  end do
end do

!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do j=0,nx
!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,ubR(j,1))
uL=ubR(j,2)/ubR(j,1)
vL=ubR(j,3)/ubR(j,1)
wL=ubR(j,4)/ubR(j,1)
pL=max(premin,(gam-1.)*( ubR(j,5)-0.5*rhoL &
   *(uL**2.+vL**2.+wL**2.) ) )
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

rhoR=max(densmin,ubL(j+1,1))
uR=ubL(j+1,2)/ubL(j+1,1)
vR=ubL(j+1,3)/ubL(j+1,1)
wR=ubL(j+1,4)/ubL(j+1,1)
pR=max(premin,(gam-1.)*( ubL(j+1,5)-0.5*rhoR &
   *(uR**2.+vR**2.+wR**2.) ) )
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR


fleft(j,1)=rhoL*uL
fleft(j,2)=rhoL*uL*uL+pL
fleft(j,3)=rhoL*uL*vL
fleft(j,4)=rhoL*uL*wL
fleft(j,5)=0.5*rhoL*(uL**2.+vL**2.+wL**2.)*uL+gam*gamul*pL*uL

fright(j,1)=rhoR*uR
fright(j,2)=rhoR*uR*uR+pR
fright(j,3)=rhoR*uR*vR
fright(j,4)=rhoR*uR*wR
fright(j,5)=0.5*rhoR*(uR**2.+vR**2.+wR**2.)*uR+gam*gamul*pR*uR

if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

call RoeSolver(j)


if(debug==1)then
  print*,'Max wave speed=',smax
end if

end do

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------
do j=0,nx
!Check for transonic rarefaction waves
if(lambda1L(j)<0. .and. lambda1R(j)>0.)then
  do k=1,5
!    flux(j,k)=FL(k)+(lambda1L(j)*(lambda1R(j)-lambda(j,1))/(lambda1R(j)-lambda1L(j)))&
!            *alpha(j,1)*K1(j,k)  
     g(j,k)= fleft(j,k)+(lambda1L(j)*(lambda1R(j)-lambda(j,1))/(lambda1R(j)-lambda1L(j)))&
            *alpha(j,1)*K1(j,k)  
  end do
else if(lambda5L(j)<0. .and. lambda5R(j)>0.)then
  do k=1,5
!    flux(j,k)=FR(k)-(lambda3R(j)*(lambda(j,3)-lambda3L(j))/(lambda3R(j)-lambda3L(j)))&
!            *alpha(j,3)*K3(j,k)
     g(j,k)=fright(j,k)-(lambda5R(j)*(lambda(j,5)-lambda5L(j))/(lambda5R(j)-lambda5L(j)))&
            *alpha(j,5)*K5(j,k)
  end do
else
  do k=1,5
!    flux(j,k)=0.5*(FL(k)+FR(k))-0.5*(alpha(j,1)*abs(lambda(j,1))*K1(j,k)&
!          +alpha(j,2)*abs(lambda(j,2))*K2(j,k)+alpha(j,3)*abs(lambda(j,3))*K3(j,k)) 
    g(j,k)=0.5*(fleft(j,k)+fright(j,k))-0.5*(alpha(j,1)*abs(lambda(j,1))*K1(j,k) &
          +alpha(j,2)*abs(lambda(j,2))*K2(j,k)+alpha(j,3)*abs(lambda(j,3))*K3(j,k) &
          +alpha(j,4)*abs(lambda(j,4))*K4(j,k)+alpha(j,5)*abs(lambda(j,5))*K5(j,k) ) 
  end do
end if

end do

end subroutine computeFluxMH

subroutine computeFluxENO()

!-----------------------------------------------------------------
!Characteristic Variables: Q_i, Qinv_i v_i
!----------------------------------------------------------------- 
do j=-1,nx+2
  dens=max(u1(j,1),densmin)   
  velx=u1(j,2)/u1(j,1)
  vely=u1(j,3)/u1(j,1)
  velz=u1(j,4)/u1(j,1)
  pres=max(premin,(gam-1.)*( u1(j,5)-0.5*dens &
     *(velx**2.+vely**2.+velz**2.) ) )
  a=sqrt(gam*pres/dens)
  H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a   
 
  Qx(j,1,1)=1.
  Qx(j,2,1)=velx-a
  Qx(j,3,1)=vely
  Qx(j,4,1)=velz
  Qx(i,5,1)=H-velx*a

  Qx(j,1,2)=1.
  Qx(j,2,2)=velx
  Qx(j,3,2)=vely
  Qx(j,4,2)=velz
  Qx(j,5,2)=0.5*(velx**2.+vely**2.+velz**2.)

  Qx(j,1,3)=0.
  Qx(j,2,3)=0.
  Qx(j,3,3)=1.
  Qx(j,4,3)=0.
  Qx(j,5,3)=vely

  Qx(j,1,4)=0.
  Qx(j,2,4)=0.
  Qx(j,3,4)=0.
  Qx(j,4,4)=1.
  Qx(j,5,4)=velz

  Qx(j,1,5)=1.
  Qx(j,2,5)=velx+a
  Qx(j,3,5)=vely
  Qx(j,4,5)=velz
  Qx(i,5,5)=H+velx*a

  Qinv(j,1,1)=((gam-1.)/(2.*a*a))*(H &
              +a*(velx-a)/(gam-1.))
  Qinv(j,1,2)=-((gam-1.)/(2.*a*a))*(velx+a/(gam-1.))
  Qinv(j,1,3)=-((gam-1.)/(2.*a*a))*vely
  Qinv(j,1,4)=-((gam-1.)/(2.*a*a))*velz
  Qinv(j,1,5)=((gam-1.)/(2.*a*a))
  
  Qinv(j,2,1)=((gam-1.)/(2*a*a))*(-2.*H+4.*a*a/(gam-1.))
  Qinv(j,2,2)=((gam-1.)/(a*a))*velx
  Qinv(j,2,3)=((gam-1.)/(a*a))*vely
  Qinv(j,2,4)=((gam-1.)/(a*a))*velz
  Qinv(j,2,5)=-((gam-1.)/(a*a))

  Qinv(j,3,1)=-vely
  Qinv(j,3,2)=0.
  Qinv(j,3,3)=1.
  Qinv(j,3,4)=0.
  Qinv(j,3,5)=0.

  Qinv(j,4,1)=-velz
  Qinv(j,4,2)=0.
  Qinv(j,4,3)=0.
  Qinv(j,4,4)=1.
  Qinv(j,4,5)=0.

  Qinv(j,5,1)=((gam-1.)/(2.*a*a))*(H &
              -a*(velx+a)/(gam-1.))
  Qinv(j,5,2)=((gam-1.)/(2.*a*a))*(-velx+a/(gam-1.))
  Qinv(j,5,3)=-((gam-1.)/(2.*a*a))*vely
  Qinv(j,5,4)=-((gam-1.)/(2.*a*a))*velz
  Qinv(j,5,5)=((gam-1.)/(2.*a*a))
     
  do k=1,5
    v(j,k)=Qinv(j,k,1)*u1(j,1)+Qinv(j,k,2)*u1(j,2)+Qinv(j,k,3)*u1(i,3) &
          +Qinv(j,k,4)*u1(j,4)+Qinv(j,k,5)*u1(j,5)
  end do   
end do

!------------------------------------------------------------------
!Slope Limiter: S_j
!------------------------------------------------------------------  
do j=0,nx+1  
 do k=1,5 
   S(j,k)=minmod(v(j+1,k)-v(j,k),v(j,k)-v(j-1,k))/dx	
 end do
end do

!------------------------------------------------------------------
!(dq/dx)_j
!------------------------------------------------------------------ 
do j=0,nx+1  
 do k=1,5
   dqdx(j,k)=Qx(j,k,1)*S(j,1)+Qx(j,k,2)*S(j,2)+Qx(j,k,3)*S(j,3) &
            +Qx(j,k,4)*S(j,4)+Qx(j,k,5)*S(j,5)
 end do
end do

!------------------------------------------------------------------
!(dq/dt)_j
!------------------------------------------------------------------ 
do j=0,nx+1  
  dens=max(u1(j,1),densmin)   
  velx=u1(j,2)/u1(j,1)
  vely=u1(j,3)/u1(j,1)
  velz=u1(j,4)/u1(j,1)
  pres=max(premin,(gam-1.)*( u1(j,5)-0.5*dens &
     *(velx**2.+vely**2.+velz**2.) ) )
  a=sqrt(gam*pres/dens)
  H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a 
   
  dqdt(j,1)=-dqdx(j,2)
  dqdt(j,2)=-((gam-1.)*H-velx**2.-a**2.)*dqdx(j,1) &
            -(3.-gam)*velx*dqdx(j,2)+(gam-1.)*vely*dqdx(j,3) &
            +(gam-1.)*velz*dqdx(j,4)-(gam-1.)*dqdx(j,5)
  dqdt(j,3)=velx*vely*dqdx(j,1)-vely*dqdx(j,2)-velx*dqdx(j,3)
  dqdt(j,4)=velx*velz*dqdx(j,1)-velz*dqdx(j,2)-velx*dqdx(j,4)
  dqdt(j,5)=-0.5*velx*((gam-3.)*H-a**2.)*dqdx(j,1) &
            -(H-(gam-1.)*velx**2.)*dqdx(j,2) &
            +(gam-1.)*velx*vely*dqdx(j,3) &
            +(gam-1.)*velx*velz*dqdx(j,4)-gam*velx*dqdx(j,5)

end do
 
!------------------------------------------------------------------
!^q_j
!------------------------------------------------------------------
do j=0,nx+1  
 do k=1,5
   qh(j,k)=u1(j,k)+0.5*(dt*dqdt(j,k)+dx*dqdx(j,k))
 end do 
end do
 
!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do j=0,nx
 !set left and right states for local Riemann problem at j+1/2 cell interface
 rhoL=max(densmin,qh(j,1))
 uL=qh(j,2)/qh(j,1)
 vL=qh(j,3)/qh(j,1)
 wL=qh(j,4)/qh(j,1)
 pL=max(premin,(gam-1.)*( qh(j,5)-0.5*rhoL &
    *(uL**2.+vL**2.+wL**2.) ) )
 aL=sqrt(gam*pL/rhoL)
 HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

 rhoR=max(densmin,qh(j+1,1))
 uR=qh(j+1,2)/qh(j+1,1)
 vR=qh(j+1,3)/qh(j+1,1)
 wR=qh(j+1,4)/qh(j+1,1)
 pR=max(premin,(gam-1.)*( qh(j+1,5)-0.5*rhoR &
    *(uR**2.+vR**2.+wR**2.) ) )
 aR=sqrt(gam*pR/rhoR)
 HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR
   
 if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
 end if

 call RoeSolver(j)

 if(debug==1)then
  print*,'Max wave speed=',smax
 end if

end do
   
!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=dx*cour/smax 

!-----------------------------------------------
!f_j
!-----------------------------------------------
do j=0,nx+1
 u(1)=qh(j,1)
 u(2)=qh(j,2)/qh(j,1)
 u(3)=qh(j,3)/qh(j,1)
 u(4)=qh(j,4)/qh(j,1)
 u(5)=(gam-1.)*(qh(j,5)-0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.))
 call flux(u)
   do k=1,5
     f(j,k)=fl(k)
   end do    
end do

!------------------------------------------------
!^g_i+1/2
!------------------------------------------------
do j=0,nx 
 do k=1,5
   g(j,k)=0.5*(f(j,k)+f(j+1,k))-0.5*&
        ( K1(j,k)*( Q(lambda(j,1))*alpha(j,1) ) &
         +K2(j,k)*( Q(lambda(j,2))*alpha(j,2) ) &
	 +K3(j,k)*( Q(lambda(j,3))*alpha(j,3) ) &
         +K4(j,k)*( Q(lambda(j,4))*alpha(j,4) ) &  
         +K5(j,k)*( Q(lambda(j,5))*alpha(j,5) ) ) 

 end do
end do
 
end subroutine computeFluxENO


subroutine starRegion()

!Local Variables
real*8::rhoB,aB,pMin,pMax,p0

rhoB=0.5*(rhoL+rhoR)
aB=0.5*(aL+aR)

pMin=min(pL,pR)
pMax=max(pL,pR)

pS=0.5*(pL+pR)+0.5*(uL-uR)*rhoB*aB

if(pS<pMax .and. pS>pMin .and. (pMax/pMin)<Q_user)then

!PV solver(solution of linearized primitive variable Euler eqns)
!u_star
uS=0.5*(uL+uR)+0.5*(pL-pR)/(rhoB*aB)
!rho_star
rhoSL=rhoL+(uL-uS)*(rhoB/aB)
rhoSR=rhoR+(uS-uR)*(rhoB/aB)

else if(pS<=pMin) then

!Two-Rarefaction Riemann Solver
!p_star
pS=aL+aR-0.5*(gam-1.)*(uR-uL)
pS=pS/(aL/(pL**gamee)+aR/(pR**gamee))
pS=pS**(1./gamee)
!u_star
uS=uL-2.*aL*gamul*(-1.+(pS/pL)**gamee)
!rho_star
rhoSL=rhoL*(pS/pL)**(1./gam)
rhoSR=rhoR*(pS/pR)**(1./gam)

else if(pS>pMin) then

!Two-Shock Riemann Solver
!p_star
p0=max(0._8,pS)
pS=gL(p0)*pL+gR(p0)*pR-(uR-uL)
pS=pS/(gL(p0)+gR(p0))
!u_star
uS=0.5*(uL+uR)+0.5*((pS-pR)*gR(p0)-(pS-pL)*gL(p0))
!rho_star
rhoSL=rhoL*((pS/pL)+gamel)/((pS/pL)*gamel+1.)
rhoSR=rhoR*((pS/pR)+gamel)/((pS/pR)*gamel+1.)
end if

aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)

end subroutine starRegion

function gL(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(x+gamel*pL))
end function gL

function gR(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(x+gamel*pR))
end function gR


subroutine flux(u)

!Input Variables    
real*8::u(5)

fl(1)=u(1)*u(2)
fl(2)=u(1)*u(2)*u(2)+u(5)
fl(3)=u(1)*u(2)*u(3)
fl(4)=u(1)*u(2)*u(4)
fl(5)=0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.)*u(2)+gam*gamul*u(2)*u(5)
  
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
  
do k=1,5
 u2(0,k)=u2(1,k)	
 u2(-1,k)=u2(1,k)
 u2(nx+1,k)=u2(nx,k)
 u2(nx+2,k)=u2(nx,k)
end do

end subroutine bound

end module RoeSolver_mod
