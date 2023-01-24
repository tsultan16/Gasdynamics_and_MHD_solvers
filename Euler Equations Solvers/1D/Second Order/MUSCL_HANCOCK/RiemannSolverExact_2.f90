!Exact Riemann Solver of 1D Eulers Equiations (Improved version, nearly perfect)
!Two root solvers included: Secant Method and Newton
!Newton solver performs better

program RiemannSolver
implicit none

integer,parameter::rootSolverOption=2 !1:Secant 2:Newton
integer,parameter ::nt=100
integer,parameter ::nx=400
integer ::i,j
real*8 ::dt(nt),dx !time step and spatial resolution
real*8::xmin,xmax,xt,x,t
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,pL,aL !left state
real*8::rhoR,uR,pR,aR !right state
real*8::rhoSL,rhoSR,uS,pS,aSL,aSR,fR; !star region
real*8::sL,sR,sHL,sTL,sHR,sTR
real*8::w1,w2,w3

!------------------------------------------------------
!initialization
!------------------------------------------------------
open(unit=10,file='input.txt')
open(unit=11,file='output.txt')
open(unit=12,file='dt.txt')

read(unit=10,fmt=*) rhoL,uL,pL
read(unit=10,fmt=*) rhoR,uR,pR
print*,'Left State=',rhoL,uL,pL
print*,'Right State=',rhoR,uR,pR

do i=1,nt
 read(12,*) dt(i)
end do

!Discountinuity in initial condition at the middle of [xmin,xmax]
xmin=0.0
xmax=1.0

print*,'xmin,xmax=',xmin,xmax

t=0._8
dx=(xmax-xmin)/nx

gam=1.4
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

!sound speeds for left and right states
aL=sqrt(gam*pL/rhoL)
aR=sqrt(gam*pR/rhoR)

!------------------------------------------------------
!compute pressure and velocity in star region
!------------------------------------------------------
if (rootSolverOption==1)then
  call findroot()
else if(rootSolverOption==2)then
  call findRoot2()
end if

if(pS>pR)then
    fR=(pS-pR)*sqrt((2.*gamil/rhoR)/(pS+gamel*pR))
  else
    fR=2*aR*gamul*(-1.+(pS/pR)**(gamee))
  end if
uS=uR+fR


print*,'p_star=',pS
print*,'u_star=',uS

!go to 1000

!Construct solution profile

!------------------------------------------------------
!Two-Shocks
!------------------------------------------------------
if(pS>pR .and. pS>pL)then

!rho_star
rhoSL=rhoL*(((pS/pL)+gamel)/(gamel*(pS/pL)+1))
rhoSR=rhoR*(((pS/pR)+gamel)/(gamel*(pS/pR)+1))

!Compute shock speeds
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)

 do i=1,nt
  t=t+dt(i)
  do j=1,nx
    x=(j-0.5)*dx
    xt=(-0.5*(xmax-xmin)+(j-0.5)*dx)/t
    call solutionSample1()
    write(11,*) x,w1,w2,w3
  end do
 end do

!------------------------------------------------------
!Two-Rarefactions
!------------------------------------------------------
else if(pS<=pL .and. pS<=pR)then

!rho_star
rhoSL=rhoL*(pS/pL)**(1./gam)
rhoSR=rhoR*(pS/pR)**(1./gam)
!Rarefaction head,tail
aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)
sHL=uL-aL
sTL=uS-aSL
sHR=uR+aR
sTR=uS+aSR

 do i=1,nt
  t=t+dt(i)
  do j=1,nx
    x=(j-0.5)*dx
    xt=(-0.5*(xmax-xmin)+(j-0.5)*dx)/t
    call solutionSample2()
    write(11,*) x,w1,w2,w3
  end do
 end do

!------------------------------------------------------
!Left Shock, Right Rarefaction
!------------------------------------------------------
else if(pS>pL .and. pS<=pR)then

!rho_star
rhoSL=rhoL*(((pS/pL)+gamel)/(gamel*(pS/pL)+1))
rhoSR=rhoR*(pS/pR)**(1./gam)
!left shock speed
sL=uL-aL*sqrt(gamuu*(pS/pL)+gamee)
!Right rarefaction head,tail
aSR=sqrt(gam*pS/rhoSR)
sHR=uR+aR
sTR=uS+aSR
 do i=1,nt
  t=t+dt(i)
  do j=1,nx
    x=(j-0.5)*dx
    xt=(-0.5*(xmax-xmin)+(j-0.5)*dx)/t
    call solutionSample3()
    write(11,*) x,w1,w2,w3
  end do
 end do

!------------------------------------------------------
!Left Rarefaction, Right Shock
!------------------------------------------------------
else if(pS<=pL .and. pS>pR)then

!rho_star
rhoSL=rhoL*(pS/pL)**(1./gam)
rhoSR=rhoR*(((pS/pR)+gamel)/(gamel*(pS/pR)+1))
!Left rarefaction head,tail
asL=sqrt(gam*pS/rhoSL)
sHL=uL-aL
sTL=uS-aSL
!Compute right shock speed
sR=uR+aR*sqrt(gamuu*(pS/pR)+gamee)

 do i=1,nt
  t=t+dt(i)
  do j=1,nx
    x=(j-0.5)*dx
    xt=(-0.5*(xmax-xmin)+(j-0.5)*dx)/t
    call solutionSample4()
    write(11,*) x,w1,w2,w3
  end do
 end do
end if

!1000 continue

print*,'Done...'

close(unit=10)
close(unit=11)
close(unit=12)

!-----------------------------------------------------------------s
contains

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
   integer ::Niter = 30
   integer ::i, iterCount

   real*8 ::p0,p1 !initial guesses
   real*8 ::tol

   iterCount=0
   !set tolerance
   tol=0.00000001
   p0=0.75*min(pL,pR)
   p1=1.2*max(pL,pR)

   do i=1,Niter
     !print*,'Iteration# ',i,', p0,p1,ps=',p0,p1,pS

     if(abs((p1-p0)/p0)<tol) then  !stop interations after convergernce criterion satisfied
       print*,'# of iterations:',i
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
   integer ::Niter = 30
   integer ::i, iterCount

   real*8 ::p0 !initial guess
   real*8 ::tol,p1

   iterCount=0
   !set tolerance
   tol=0.00000001
   p0=0.5*(pL+pR)
   p1=2.*p0

   do i=1,Niter
     !print*,'Iteration# ',i,', p0,ps=',p0,pS

     if(abs((p1-p0)/p0)<tol) then  !stop interations after convergernce criterion satisfied
       print*,'# of iterations:',i
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
subroutine solutionSample1()

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
subroutine solutionSample2()

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
subroutine solutionSample3()

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
subroutine solutionSample4()

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

!-----------------------------------------------------------
end program RiemannSolver
