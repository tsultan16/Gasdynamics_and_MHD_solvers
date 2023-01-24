!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_mod
implicit none

integer,parameter::debug=0
integer,parameter ::nt=100
integer,parameter ::nx=200
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL,HL,FL(3) !left state
real*8::rhoR,uR,pR,aR,HR,FR(3) !right state 
real*8::delrho,delp,delu,smax
real*8::dens,vel,pres
real*8::rhot,ut,Ht,at,lambda(3),K1(3),K2(3),K3(3),alpha(3)	
real*8::u1(-1:nx+2,3),u2(-1:nx+2,3),flux(-1:nx+2,3)

contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_Roe.txt')
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
flux(j,1)=0.5*(FL(1)+FR(1))-0.5*(alpha(1)*abs(lambda(1))*K1(1)&
          +alpha(2)*abs(lambda(2))*K2(1)+alpha(3)*abs(lambda(3))*K3(1))
flux(j,2)=0.5*(FL(2)+FR(2))-0.5*(alpha(1)*abs(lambda(1))*K1(2)&
          +alpha(2)*abs(lambda(2))*K2(2)+alpha(3)*abs(lambda(3))*K3(2))
flux(j,3)=0.5*(FL(3)+FR(3))-0.5*(alpha(1)*abs(lambda(1))*K1(3)&
          +alpha(2)*abs(lambda(2))*K2(3)+alpha(3)*abs(lambda(3))*K3(3))

end do

end subroutine computeFlux

end module RoeSolver_mod
