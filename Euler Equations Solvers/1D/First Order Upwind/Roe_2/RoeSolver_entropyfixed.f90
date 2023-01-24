!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_entropyfixed_mod
implicit none

integer,parameter::debug=0
integer,parameter ::nt=100
integer,parameter ::nx=200
real*8,parameter::Q_user=2.0 
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL,HL,FL(3) !left state
real*8::rhoR,uR,pR,aR,HR,FR(3) !right state 
real*8::delrho,delp,delu,smax
real*8::lambda1L,lambda1R,lambda3L,lambda3R
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,lambda(3),K1(3),K2(3),K3(3),alpha(3)	
real*8::u1(-1:nx+2,3),u2(-1:nx+2,3),flux(-1:nx+2,3)

contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_Roe_fix.txt')
open(unit=12,file='dt.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
cour=0.9
dx=(xmax-xmin)/nx

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

subroutine computeFlux()

do j=0,nx
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

delrho=rhoR-rhoL
delu=uR-uL
delp=pR-pL

FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL
FL(3)=0.5*rhoL*uL*uL*uL+gam*gamul*pL*uL

FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR
FR(3)=0.5*rhoR*uR*uR*uR+gam*gamul*pR*uR

if(debug==1)then
  print*,'Cell: ',j
  print*,'Left State=',rhoL,uL,pL
  print*,'Right State=',rhoR,uR,pR
end if

!Compute star region state
call starRegion()

!-------------------------------------------
!Wave Speeds for Harten-Hyman Entropy Fix
!-------------------------------------------
lambda1L=uL-aL
lambda1R=uS-aSL

lambda3L=uS+aSR
lambda3R=uR+aR

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
lambda(1)=ut-at
lambda(2)=ut
lambda(3)=ut+at

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(1),lambda(2),lambda(3)
end if
!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(1)=1._8
K1(2)=ut-at
K1(3)=Ht-ut*at

K2(1)=1
K2(2)=ut
K2(3)=0.5*ut*ut

K3(1)=1._8
K3(2)=ut+at
K3(3)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(2)=delrho-delp/(at*at)
alpha(3)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(1)),abs(lambda(2)),abs(lambda(3)))

if(debug==1)then
  print*,'Max wave speed=',smax
end if
!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------

!Check for transonic rarefaction waves
if(lambda1L<0. .and. lambda1R>0.)then
  do k=1,3
    flux(j,k)=FL(k)+(lambda1L*(lambda1R-lambda(1))/(lambda1R-lambda1L))&
            *alpha(1)*K1(k)  
  end do
else if(lambda3L<0. .and. lambda3R>0.)then
  do k=1,3
    flux(j,k)=FR(k)-(lambda3R*(lambda(3)-lambda3L)/(lambda3R-lambda3L))&
            *alpha(3)*K3(k)
  end do
else
  do k=1,3
    flux(j,k)=0.5*(FL(k)+FR(k))-0.5*(alpha(1)*abs(lambda(1))*K1(k)&
          +alpha(2)*abs(lambda(2))*K2(k)+alpha(3)*abs(lambda(3))*K3(k)) 
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


end module RoeSolver_entropyfixed_mod
