!Adaptive Non-iterative Approximate State Riemann Solver
module HLLRiemannSolver_mod
implicit none

integer,parameter::rootSolverOption=2 !1:Secant 2:Newton
integer,parameter ::nt=100
integer,parameter ::nx=500
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL,FL(3) !left state
real*8::rhoR,uR,pR,aR,FR(3) !right state
real*8::rhoSL,rhoSR,uS,pS,aSL,aSR !star region state
real*8::rhoB,aB,p0,pMax,pMin
real*8,parameter::Q_user=2.0 
real*8::sL,sR,smax,smax1
real*8::w1,w2,w3,dens,vel,pres

real*8::u1(-1:nx+2,3),u2(-1:nx+2,3),flux(-1:nx+2,3)

contains 

subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_Godunov_appr.txt')
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

smax1=0.
do j=0,nx
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=u1(j,1)
uL=u1(j,2)/u1(j,1)
pL=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))

rhoR=u1(j+1,1)
uR=u1(j+1,2)/u1(j+1,1)
pR=(gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1))

FL(1)=rhoL*uL
FL(2)=rhoL*uL*uL+pL
FL(3)=0.5*rhoL*uL*uL*uL+gamul*pL*uL

FR(1)=rhoR*uR
FR(2)=rhoR*uR*uR+pR
FR(3)=0.5*rhoR*uR*uR*uR+gamul*pR*uR


call computeWaveSpeeds()

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax1=max(smax1,smax)

!---------------------------------------
!F_j+1/2
!---------------------------------------
 if (sL>=0._8)then
   do k=1,3
     flux(j,k)=FL(k)
   end do
 else if(sL<0._8 .and. sR>0._8) then
   do k=1,3
     flux(j,k)=sR*FL(k)-sL*FR(k)+sL*sR*(u1(j+1,k)-u1(j,k))
     flux(j,k)=flux(j,k)/(sR-sL)
   end do
 else if(sR<=0._8)then
   do k=1,3
     flux(j,k)=FR(k)
   end do 
 end  if

end do

end subroutine computeFlux

subroutine computeWaveSpeeds()

!sound speeds for left and right states
aL=sqrt(gam*pL/rhoL)
aR=sqrt(gam*pR/rhoR)

!------------------------------------------------------
!compute pressure and velocity in star region
!------------------------------------------------------
rhoB=0.5*(rhoL+rhoR)
aB=0.5*(aL+aR)

pMin=min(pL,pR)
pMax=max(pL,pR)

pS=0.5*(pL+pR)+0.5*(uL-uR)*rhoB*aB

if(pS<pMax .and. pS>pMin .and. (pMax/pMin)<Q_user)then

!PV solver(solution of linearized primitive variable Euler eqns)

else if(pS<=pMin) then

!Two-Rarefaction Riemann Solver
!p_star
pS=aL+aR-0.5*(gam-1.)*(uR-uL)
pS=pS/(aL/(pL**gamee)+aR/(pR**gamee))
pS=pS**(1./gamee)
else if(pS>pMin) then

!Two-Shock Riemann Solver
!p_star
p0=max(0._8,pS)
pS=gL(p0)*pL+gR(p0)*pR-(uR-uL)
pS=pS/(gL(p0)+gR(p0))
end if

!Compute Wave Speeds

!------------------------------------------------------
!Two-Shocks
!------------------------------------------------------
if(pS>pR .and. pS>pL)then

!HLL Wave speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)
!------------------------------------------------------
!Two-Rarefactions
!------------------------------------------------------
else if(pS<=pL .and. pS<=pR)then

!HLL Wave speeds
sL=uL-aL
sR=uR+aR
!------------------------------------------------------
!Left Shock, Right Rarefaction
!------------------------------------------------------
else if(pS>pL .and. pS<=pR)then

!HLL Wave speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR
!------------------------------------------------------
!Left Rarefaction, Right Shock
!------------------------------------------------------
else if(pS<=pL .and. pS>pR)then

!HLL Wave speeds
sL=uL-aL
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)
end if

smax=abs(sR)

end subroutine computeWaveSpeeds


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


end module HLLRiemannSolver_mod
