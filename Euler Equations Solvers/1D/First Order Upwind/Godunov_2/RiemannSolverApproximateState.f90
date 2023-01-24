!Adaptive Non-iterative Approximate State Riemann Solver
module ApproximateStateRiemannSolver_mod
implicit none

integer,parameter::rootSolverOption=2 !1:Secant 2:Newton
integer,parameter ::nt=100
integer,parameter ::nx=500
integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL !left state
real*8::rhoR,uR,pR,aR !right state
real*8::rhoSL,rhoSR,uS,pS,aSL,aSR,fR !star region state
real*8::rhoB,aB,p0,pMax,pMin
real*8,parameter::Q_user=2.0 
real*8::sL,sR,sHL,sTL,sHR,sTR,smax,smax1
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

call computeSolution(0._8)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax1=max(smax1,smax)

!---------------------------------------
!F_j+1/2
!---------------------------------------
flux(j,1)=w1*w2
flux(j,2)=w1*w2*w2+w3
flux(j,3)=0.5*w1*w2*w2*w2+gam*gamul*w3*w2

end do

end subroutine computeFlux

subroutine computeSolution(xt)
real*8::xt
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

!Construct solution profile

!------------------------------------------------------
!Two-Shocks
!------------------------------------------------------
if(pS>pR .and. pS>pL)then
!Compute shock speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)

smax=max(abs(sL),abs(sR))

call solutionSample1(xt)

!------------------------------------------------------
!Two-Rarefactions
!------------------------------------------------------
else if(pS<=pL .and. pS<=pR)then
!Rarefaction head,tail
aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)
sHL=uL-aL
sTL=uS-aSL
sHR=uR+aR
sTR=uS+aSR

smax=max(abs(sHR),abs(sHL))

call solutionSample2(xt)

!------------------------------------------------------
!Left Shock, Right Rarefaction
!------------------------------------------------------
else if(pS>pL .and. pS<=pR)then
!left shock speed
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
!Right rarefaction head,tail
aSR=sqrt(gam*pS/rhoSR)
sHR=uR+aR
sTR=uS+aSR

smax=max(abs(sL),abs(sHR))
 
call solutionSample3(xt)

!------------------------------------------------------
!Left Rarefaction, Right Shock
!------------------------------------------------------
else if(pS<=pL .and. pS>pR)then
!Left rarefaction head,tail
asL=sqrt(gam*pS/rhoSL)
sHL=uL-aL
sTL=uS-aSL
!Compute right shock speed
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)

smax=max(abs(sHL),abs(sR))

call solutionSample4(xt)

end if

end subroutine computeSolution

!-----------------------------------------------------------------

!pressure function
function f(p) result(fx)
  real*8,intent(in)::p
  real*8 ::fx,fL,fR

  !fL
  if(p>pL)then
    fL=(p-pL)*sqrt((2.*gamil/rhoL)/(p+gamel*pL))
  else
    fL=2*aL*gamul*(-1.+(p/pL)**(gamee))
  end if
  !fR
  if(p>pR)then
    fR=(p-pR)*sqrt((2.*gamil/rhoR)/(p+gamel*pR))
  else
    fR=2*aR*gamul*(-1.+(p/pR)**(gamee))
  end if
  !f
  fx=fL+fR+uR-uL

end function f

function fprime(p) result(fxp)
  real*8,intent(in)::p
  real*8 ::fxp,fLp,fRp

  !fL
  if(p>pL)then
    fLp=(1-0.5*(p-pL)/(p+gamel*pL))*sqrt((2.*gamil/rhoL)/(p+gamel*pL))
  else
    fLp=((p/pL)**(-gamuu))/(rhoL*aL)
  end if
  !fR
  if(p>pR)then
    fRp=(1-0.5*(p-pR)/(p+gamel*pR))*sqrt((2.*gamil/rhoR)/(p+gamel*pR))
  else
    fRp=((p/pR)**(-gamuu))/(rhoR*aR)
  end if
  !f
  fxP=fLP+fRP

end function fprime

subroutine findroot() !secant method root solver
   !use pL and pR as initial guesses for pS
   integer ::Niter
   integer ::i, iterCount

   real*8 ::p0,p1 !initial guesses
   real*8 ::tol

   Niter = 30
   iterCount=0
   !set tolerance
   tol=0.00000001
   p0=0.75*min(pL,pR)
   p1=1.2*max(pL,pR)

   do i=1,Niter
     !print*,'Iteration# ',i,', p0,p1,ps=',p0,p1,pS

     if(abs((p1-p0)/p0)<tol) then  !stop interations after convergernce criterion satisfied
       !print*,'# of iterations:',i
       exit
     else if(abs((p1-p0)/p0)>tol .and. i<Niter-1) then
       Niter=Niter+10 !add more iterations if solution does not converge
       iterCount=iterCount+1
     else if(abs((p1-p0)/p0)>tol .and. i<Niter-1 .and. iterCount >20) then
       print*,'findRoot() unable to converge to the root...'
       stop
     end if
     pS=max(p1-(f(p1)*(p1-p0)/(f(p1)-f(p0))),0.00001) !avoid negative pressure
     p0=max(p1,0.00002)
     p1=pS

   end do

end subroutine findRoot

subroutine findRoot2() !Newton-Raphson root solver
   integer ::Niter
   integer ::i, iterCount

   real*8 ::p0 !initial guess
   real*8 ::tol,p1

   Niter = 30

   iterCount=0
   !set tolerance
   tol=0.00000001
   p0=0.5*(pL+pR)
   p1=2.*p0

   do i=1,Niter
     !print*,'Iteration# ',i,', p0,ps=',p0,pS

     if(abs((p1-p0)/p0)<tol) then  !stop interations after convergernce criterion satisfied
       !print*,'# of iterations:',i
       exit
     else if(abs((p1-p0)/p0)>tol .and. i<Niter-1) then
       Niter=Niter+10 !add more iterations if solution does not converge
       iterCount=iterCount+1
     else if(abs((p1-p0)/p0)>tol .and. i<Niter-1 .and. iterCount >20) then
       print*,'findRoot() unable to converge to the root...'
       stop
     end if
     pS=max(p0-(f(p0)/fprime(p0)),0.0001) !avoid negative pressure
     p1=p0
     p0=pS
     
   end do

end subroutine findRoot2

!Two shock
subroutine solutionSample1(xt)
real*8::xt

if(xt<=sL)then
  w1=rhoL
  w2=uL
  w3=pL
else if(xt>sL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=pS
else if(xt>=uS .and. xt<sR)then
  w1=rhoSR
  w2=uS
  w3=pS
else if(xt>=sR)then
  w1=rhoR
  w2=uR
  w3=pR
end if

end subroutine solutionSample1

!Two Rarefaction
subroutine solutionSample2(xt)
real*8::xt

if(xt<=sHL)then
  w1=rhoL
  w2=uL
  w3=pL
else if(xt>sHL .and. xt<sTL)then
  w1=rhoL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gamul)
  w2=2.*gamil*(aL+0.5*(gam-1.)*uL+xt)
  w3=pL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gam*gamul)
else if(xt>=sTL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=pS
else if(xt>=uS .and. xt<sTR)then
  w1=rhoSR
  w2=uS
  w3=pS
else if(xt>=sTR .and. xt<sHR)then
  w1=rhoR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gamul)
  w2=2.*gamil*(-aR+0.5*(gam-1.)*uR+xt)
  w3=pR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gam*gamul)
else if(xt>=sHR)then
  w1=rhoR
  w2=uR
  w3=pR
end if

end subroutine solutionSample2

!left shock, right rarefaction
subroutine solutionSample3(xt)
real*8::xt

if(xt<=sL)then
  w1=rhoL
  w2=uL
  w3=pL
else if(xt>sL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=pS
else if(xt>=uS .and. xt<sTR)then
  w1=rhoSR
  w2=uS
  w3=pS
else if(xt>=sTR .and. xt<sHR)then
  w1=rhoR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gamul)
  w2=2.*gamil*(-aR+0.5*(gam-1.)*uR+xt)
  w3=pR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gam*gamul)
else if(xt>=sHR)then
  w1=rhoR
  w2=uR
  w3=pR
end if

end subroutine solutionSample3

!left rarefaction, right shock
subroutine solutionSample4(xt)
real*8::xt

if(xt<=sHL)then
  w1=rhoL
  w2=uL
  w3=pL
else if(xt>sHL .and. xt<sTL)then
  w1=rhoL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gamul)
  w2=2.*gamil*(aL+0.5*(gam-1.)*uL+xt)
  w3=pL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gam*gamul)
else if(xt>=sTL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=pS
else if(xt>=uS .and. xt<sR)then
  w1=rhoSR
  w2=uS
  w3=pS
else if(xt>=sR)then
  w1=rhoR
  w2=uR
  w3=pR
end if

end subroutine solutionSample4

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


end module ApproximateStateRiemannSolver_mod
