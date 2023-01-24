!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_entropyfixed_mod
implicit none

integer,parameter::debug=0
integer,parameter ::nt=200
integer,parameter ::nx=400
real*8,parameter::Q_user=2.0 
integer,parameter::characteristicLimiting=0 !0:off, 1:on 
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL,HL,FL(3) !left state
real*8::rhoR,uR,pR,aR,HR,FR(3) !right state 
real*8::delrho,delp,delu,smax
real*8::lambda1L(-1:nx+2),lambda1R(-1:nx+2),lambda3L(-1:nx+2),lambda3R(-1:nx+2)
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,lambda(-1:nx+2,3)
real*8::K1(-1:nx+2,3),K2(-1:nx+2,3),K3(-1:nx+2,3),alpha(-1:nx+2,3)	
real*8::u1(-1:nx+2,3),u2(-1:nx+2,3),flux(-1:nx+2,3)
real*8::fright(-1:nx+2,3),fleft(-1:nx+2,3)
real*8::qL(-1:nx+2,3),qR(-1:nx+2,3),slope(-1:nx+2,3)
real*8::delta1(-1:nx+2,3),delta2(-1:nx+2,3),delta3(-1:nx+2,3)
real*8::ubL(-1:nx+2,3),ubR(-1:nx+2,3)

contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
if(characteristicLimiting==0)then
 open(unit=11,file='output_EulerMH.txt')
else if(characteristicLimiting==1)then
 open(unit=11,file='output_EulerMH_lim.txt')
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

read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR

rhoL=1.
rhoR=rhoL
pL=gam*(1/30.)*(1/30.)
pR=pL
uL=-1.
uR=uL

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

aL=sqrt(gam*pL/rhoL)
HL=0.5*uL*uL+gamul*aL*aL

aR=sqrt(gam*pR/rhoR)
HR=0.5*uR*uR+gamul*aR*aR

delrho=rhoR-rhoL
delu=uR-uL
delp=pR-pL

!Compute star region state
call starRegion()

!-------------------------------------------
!Wave Speeds for Harten-Hyman Entropy Fix
!-------------------------------------------
lambda1L(jj)=uL-aL
lambda1R(jj)=uS-aSL

lambda3L(jj)=uS+aSR
lambda3R(jj)=uR+aR

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

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(jj,1),lambda(jj,2),lambda(jj,3)
end if
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


end subroutine RoeSolver

subroutine computeFlux()

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(characteristicLimiting==1)then

do j=-1,nx+2 

  !set left and right states for local Riemann problem at j+1/2 cell interface 
  rhoL=max(densmin,u1(j,1))
  uL=u1(j,2)/u1(j,1)
  pL=max(premin,(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1)))

  rhoR=max(densmin,u1(j+1,1))
  uR=u1(j+1,2)/u1(j+1,1)
  pR=max(premin,(gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1)))

  call RoeSolver(j)
  
  do k=1,3
    delta1(j,k)=alpha(j,1)*K1(j,k)
    delta2(j,k)=alpha(j,2)*K2(j,k)
    delta3(j,k)=alpha(j,3)*K3(j,k)
  end do
end do

do j=0,nx+1
  do k=1,3
   slope(j,k)=limiter(delta1(j,k)-delta1(j-1,k) &
               ,delta1(j+1,k)-delta1(j,k)) &
             +limiter(delta2(j,k)-delta2(j-1,k) &
               ,delta2(j+1,k)-delta2(j,k)) &
             +limiter(delta3(j,k)-delta3(j-1,k) &
               ,delta3(j+1,k)-delta3(j,k)) 
  end do
end do

else if(characteristicLimiting==0)then

do j=0,nx+1
  do k=1,3
   slope(j,k)=limiter(u1(j,k)-u1(j-1,k),u1(j+1,k)-u1(j,k))
  end do
end do

end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do j=0,nx+1
  do k=1,3
   qL(j,k)=u1(j,k)-0.5*slope(j,k)
   qR(j,k)=u1(j,k)+0.5*slope(j,k)
  end do
end do

!-------------------------------------------------------------------
!Evolve Boudary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------
do j=0,nx+1
  FL(1)=qL(j,2)
  FL(2)=(gam-1.)*qL(j,3)+0.5*(3.-gam)*(qL(j,2)**2.)/qL(j,1)
  FL(3)=gam*(qL(j,2)*qL(j,3)/qL(j,1))-0.5*(gam-1.)*(qL(j,2)**3.)/(qL(j,1)**2.)

  FR(1)=qR(j,2)
  FR(2)=(gam-1.)*qR(j,3)+0.5*(3.-gam)*(qR(j,2)**2.)/qR(j,1)
  FR(3)=gam*(qR(j,2)*qR(j,3)/qR(j,1))-0.5*(gam-1.)*(qR(j,2)**3.)/(qR(j,1)**2.)

  do k=1,3
    ubL(j,k)=qL(j,k)+0.5*(dt/dx)*(FL(k)-FR(k))
    ubR(j,k)=qR(j,k)+0.5*(dt/dx)*(FL(k)-FR(k))
  end do
end do

!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do j=0,nx
!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,ubR(j,1))
uL=ubR(j,2)/ubR(j,1)
pL=max(premin,(gam-1.)*(ubR(j,3)-0.5*ubR(j,2)*ubR(j,2)/ubR(j,1)))

rhoR=max(densmin,ubL(j+1,1))
uR=ubL(j+1,2)/ubL(j+1,1)
pR=max(premin,(gam-1.)*(ubL(j+1,3)-0.5*ubL(j+1,2)*ubL(j+1,2)/ubL(j+1,1)))


FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL
FL(3)=0.5*rhoL*uL*uL*uL+gam*gamul*pL*uL

FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR
FR(3)=0.5*rhoR*uR*uR*uR+gam*gamul*pR*uR

fleft(j,1)=rhoL*uL
fleft(j,2)=rhoL*uL*uL+pL
fleft(j,3)=0.5*rhoL*uL*uL*uL+gam*gamul*pL*uL
fright(j,1)=rhoR*uR
fright(j,2)=rhoR*uR*uR+pR
fright(j,3)=0.5*rhoR*uR*uR*uR+gam*gamul*pR*uR



if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

call RoeSolver(j)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(j,1)),abs(lambda(j,2)),abs(lambda(j,3)))

if(debug==1)then
  print*,'Max wave speed=',smax
end if

end do

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------
do j=0,nx
!Check for transonic rarefaction waves
if(lambda1L(j)<0. .and. lambda1R(j)>0.)then
  do k=1,3
!    flux(j,k)=FL(k)+(lambda1L(j)*(lambda1R(j)-lambda(j,1))/(lambda1R(j)-lambda1L(j)))&
!            *alpha(j,1)*K1(j,k)  
     flux(j,k)= fleft(j,k)+(lambda1L(j)*(lambda1R(j)-lambda(j,1))/(lambda1R(j)-lambda1L(j)))&
            *alpha(j,1)*K1(j,k)  
  end do
else if(lambda3L(j)<0. .and. lambda3R(j)>0.)then
  do k=1,3
!    flux(j,k)=FR(k)-(lambda3R(j)*(lambda(j,3)-lambda3L(j))/(lambda3R(j)-lambda3L(j)))&
!            *alpha(j,3)*K3(j,k)
     flux(j,k)=fright(j,k)-(lambda3R(j)*(lambda(j,3)-lambda3L(j))/(lambda3R(j)-lambda3L(j)))&
            *alpha(j,3)*K3(j,k)
  end do
else
  do k=1,3
!    flux(j,k)=0.5*(FL(k)+FR(k))-0.5*(alpha(j,1)*abs(lambda(j,1))*K1(j,k)&
!          +alpha(j,2)*abs(lambda(j,2))*K2(j,k)+alpha(j,3)*abs(lambda(j,3))*K3(j,k)) 
    flux(j,k)=0.5*(fleft(j,k)+fright(j,k))-0.5*(alpha(j,1)*abs(lambda(j,1))*K1(j,k)&
          +alpha(j,2)*abs(lambda(j,2))*K2(j,k)+alpha(j,3)*abs(lambda(j,3))*K3(j,k)) 
  end do
end if

end do

end subroutine computeFlux


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

function limiter(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z,b

  b=1.

  if(y>0._8)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter

end module RoeSolver_entropyfixed_mod
