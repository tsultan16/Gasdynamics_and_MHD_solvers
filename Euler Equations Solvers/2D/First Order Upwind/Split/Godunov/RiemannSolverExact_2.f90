!Exact Riemann Solver of 1D Eulers Equiations (Improved version, nearly perfect)
!Two root solvers included: Secant Method and Newton
!Newton solver performs better

module RiemannSolver_mod
implicit none

integer,parameter::debug=0 !0:off 1:on
integer,parameter::rootSolverOption=2 !1:Secant 2:Newton
integer,parameter::RStype=2 !1:Exact Riemann Solver, 2: Approximate State(adaptive)
integer,parameter::nt=100
integer,parameter::nx=200
integer,parameter::ny=200
real*8,parameter::s_tol=1.d-25
real,parameter::min_pres=1.d-30
real,parameter::min_dens=1.d-30
integer,parameter::vshockInit=0 !0:off :1:on 

real*8 ::dt,dx,dy !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,pL,aL !left state
real*8::rhoR,uR,vR,pR,aR !right state
real*8::rhoSL,rhoSR,uS,vS,pS,aSL,aSR,fR !star region
real*8::sL,sR,sHL,sTL,sHR,sTR,smax,smax1
real*8::w1,w2,w3,w4,dens,velx,vely,pres
real*8::u1(-1:nx+2,-1:ny+2,4),u2(-1:nx+2,-1:ny+2,4),xflux(-1:nx+2,4),yflux(-1:ny+2,4)

character(len=20) :: filename

contains 

!------------------------------------------------------
!initialization
!------------------------------------------------------
subroutine init()

!Local Variables
integer::jj,kk
!-----------------------------------------------------
real*8,parameter::pi=3.14159265359
real*8::x1,y1
real*8 velx0,velx1,vely0,vely1,velz0,velz1,rho0,rho1
real*8 velx2,vely2,velx3,vely3
real*8 pre0,pre1,M0,temp1,temp2,temp3,temp4
real*8 theta !shock angle w.r.t. x axis
real*8 alpha
!-----------------------------------------------------
open(unit=10,file='input.txt')
open(unit=0+1000,file='Output/t=0.txt')
open(unit=12,file='Output/horcut.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
cour=0.5

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny

gam=1.4
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

read(unit=10,fmt=*) rhoL,uL,vL,pL 
read(unit=10,fmt=*) rhoR,uR,vR,pR


if(vshockInit==1)then

!Intersecting shocks 	 
         
!pre-shock conditions
theta=pi/4.   !shock angle set to 45 degrees		
rho0=1.
pre0=1.	
velx0=10.!1.
vely0=0.
M0=sin(theta)*velx0/sqrt(gam*pre0/rho0) !20.	  
print*,'M0=',M0       
!using Rankine Hugoniot relations to determine post shock fluid state	  
rho1=rho0*(((gam+1.)*M0*M0)/(2.+(gam-1.)*M0*M0))
pre1=pre0*(2.*gam*M0*M0-gam+1.)/(gam+1.)
velx1=velx0*(cos(theta)**2+(rho0/rho1)*sin(theta)**2)
vely1=-velx0*cos(theta)*sin(theta)*(1-(rho0/rho1))
velx2=velx1
vely2=-vely1
		
alpha=atan(abs(vely2)/abs(velx2))
print*,'alpha=',alpha*180./pi,'degrees'			

do kk=-1,ny+2
  y1 = (kk - 0.5)*dy
do jj=-1,nx+2
  x1 = (jj - 0.5)*dx
				
  temp1=0.5-tan(theta)*(x1-0.5)
  temp2=0.5+tan(theta)*(x1-0.5)
  temp3=0.5+tan(alpha)*(x1-0.5)
  temp4=0.5-tan(alpha)*(x1-0.5)

  !pre-shock	
  if (y1 .le. temp1 .or. y1 .ge. temp2) then 
    u1(jj,kk,1) = rho0
    u1(jj,kk,2) = velx0*u1(jj,kk,1)
    u1(jj,kk,3) = vely0*u1(jj,kk,1) 
    u1(jj,kk,4) = (pre0/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3))/u1(jj,kk,1)
   !post-shock
		
   !region 2
  elseif (y1 .lt. temp2 .and. y1 .gt. temp3) then
   !print*,'Region 2 checkpoint...'
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx2*u1(jj,kk,1)
    u1(jj,kk,3) = vely2*u1(jj,kk,1) 
    u1(jj,kk,4) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3))/u1(jj,kk,1)
    !region 1
  elseif (y1 .lt. temp4 .and. y1 .gt. temp1) then
   !print*,'Region 1 checkpoint...'
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx1*u1(jj,kk,1)
    u1(jj,kk,3) = vely1*u1(jj,kk,1) 
    u1(jj,kk,4) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3))/u1(jj,kk,1)
   !region 3(deadzone)
  elseif (y1 .gt. temp4 .and. y1 .lt. temp3) then
    !print*,'Deadzone checkpoint...'
    u1(jj,kk,1) = rho1
    velx3=velx1
    !linear interpolation to get vely3, continuous across dead zone boundary  
    vely3=vely1+((vely1-vely2)/(temp4-temp3))*(y1-temp4)
    u1(jj,kk,2) = velx3*u1(jj,kk,1)
    u1(jj,kk,3) = vely3*u1(jj,kk,1) 
    u1(jj,kk,4) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3))/u1(jj,kk,1)
   end if			
end do
end do

else

do jj=-1,nx+2
 do kk=-1,ny+2
  if(jj<nx/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=0.5*rhoL*(uL*uL+vL*vL)+pL*gamul
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=0.5*rhoR*(uR*uR+vR*vR)+pR*gamul 
  end if
 end do
end do

end if
u2=u1

print*,'Left State (rho,u,v,p) =',rhoL,uL,vL,pL
print*,'Right State (rho,u,v,p) =',rhoR,uR,vR,pR

call fileOutput(0+1000)
call horcut()
close(unit=0+1000)

end subroutine init


subroutine fileOutput(iunit)
!Input variables
integer::iunit
!local variables
integer::jj,kk

do kk=1,ny
  y=ymin+(kk-0.5)*dy
  do jj=1,nx
   x=xmin+(jj-0.5)*dx
   dens=max(u2(jj,kk,1),min_dens)
   velx=u2(jj,kk,2)/u2(jj,kk,1)
   vely=u2(jj,kk,3)/u2(jj,kk,1)
   pres=max((gam-1.)*(u2(jj,kk,4)-0.5*(u2(jj,kk,2)*u2(jj,kk,2)/u2(jj,kk,1)&
        +u2(jj,kk,3)*u2(jj,kk,3)/u2(jj,kk,1))),min_pres)
   write(iunit,*) x,y,dens,velx,vely,pres
  end do
end do

end subroutine fileOutput

subroutine horcut()
!Local Variables
integer,parameter::ycell=ny*0.65
integer::jj

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=max(u2(jj,ycell,1),min_dens)
   velx=u2(jj,ycell,2)/u2(jj,ycell,1)
   vely=u2(jj,ycell,3)/u2(jj,ycell,1)
   pres=max((gam-1.)*(u2(jj,ycell,4)-0.5*(u2(jj,ycell,2)*u2(jj,ycell,2)/u2(jj,ycell,1)&
        +u2(jj,ycell,3)*u2(jj,ycell,3)/u2(jj,ycell,1))),min_pres)
   write(12,*) x,dens,velx,vely,pres
end do

end subroutine horcut

subroutine computeFlux(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj

smax1=0.
do jj=0,n
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
if(dir==1)then

rhoL=max(u1(jj,ind,1),min_dens)
uL=u1(jj,ind,2)/u1(jj,ind,1)
vL=u1(jj,ind,3)/u1(jj,ind,1)
pL=max((gam-1.)*(u1(jj,ind,4)-0.5*(u1(jj,ind,2)*u1(jj,ind,2)/u1(jj,ind,1)&
        +u1(jj,ind,3)*u1(jj,ind,3)/u1(jj,ind,1))),min_pres)

rhoR=max(u1(jj+1,ind,1),min_dens)
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
pR=max((gam-1.)*(u1(jj+1,ind,4)-0.5*(u1(jj+1,ind,2)*u1(jj+1,ind,2)/u1(jj+1,ind,1)&
        +u1(jj+1,ind,3)*u1(jj+1,ind,3)/u1(jj+1,ind,1))),min_pres)

!y-sweeps
else if(dir==2)then

rhoL=max(u1(ind,jj,1),min_dens)
uL=u1(ind,jj,2)/u1(ind,jj,1)
vL=u1(ind,jj,3)/u1(ind,jj,1)
pL=max((gam-1.)*(u1(ind,jj,4)-0.5*(u1(ind,jj,2)*u1(ind,jj,2)/u1(ind,jj,1)&
        +u1(ind,jj,3)*u1(ind,jj,3)/u1(ind,jj,1))),min_pres)

rhoR=max(u1(ind,jj+1,1),min_dens)
uR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
pR=max((gam-1.)*(u1(ind,jj+1,4)-0.5*(u1(ind,jj+1,2)*u1(ind,jj+1,2)/u1(ind,jj+1,1)&
        +u1(ind,jj+1,3)*u1(ind,jj+1,3)/u1(ind,jj+1,1))),min_pres)
end if

if(debug==1)then
  if(dir==1)then
    print*,'Cell: ',jj,ind
  end if
  if(dir==2)then
    print*,'Cell: ',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,p) =',rhoL,uL,vL,pL
  print*,'Local RP Right State (rho,u,v,p) =',rhoR,uR,vR,pR
end if

if(RStype==1)then
  call computeSolutionExact(0._8)
else if(RStype==2)then
  call computeSolutionApprox(0._8)
end if
!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax1=max(smax1,smax)


if(dir==1)then
!---------------------------------------
!F_j+1/2
!---------------------------------------
xflux(jj,1)=w1*w2
xflux(jj,2)=w1*w2*w2+w4
xflux(jj,3)=w1*w2*w3
xflux(jj,4)=0.5*w1*(w2*w2+w3*w3)*w2+gam*gamul*w2*w4
else if(dir==2)then
!---------------------------------------
!G_j+1/2
!---------------------------------------
yflux(jj,1)=w1*w3
yflux(jj,2)=w1*w2*w3
yflux(jj,3)=w1*w3*w3+w4
yflux(jj,4)=0.5*w1*(w2*w2+w3*w3)*w3+gam*gamul*w3*w4
end if

end do

end subroutine computeFlux

subroutine computeSolutionExact(xt)
!Input Variabeles
real*8::xt

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
    fR=2.*aR*gamul*(-1.+(pS/pR)**(gamee))
  end if
uS=uR+fR


!print*,'p_star=',pS
!print*,'u_star=',uS

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

smax=max(abs(sL),abs(sR))

call solutionSample1(xt)

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

smax=max(abs(sHR),abs(sHL))

call solutionSample2(xt)

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

smax=max(abs(sL),abs(sHR))
 
call solutionSample3(xt)

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

smax=max(abs(sHL),abs(sR))

call solutionSample4(xt)

end if

!check for zero wave speed
if(smax<s_tol)then
  smax=abs(w2)+sqrt(gam*w4/w1)   
end if

!check for negative pressure
if(w4<0.)then
  print*,'Negative pressure, p=',w4
end if
w4=max(w4,min_pres)

!print*,'smax=',smax

end subroutine computeSolutionExact

subroutine computeSolutionApprox(xt)
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

end subroutine computeSolutionApprox
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
  w3=vL
  w4=pL
else if(xt>sL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=vL
  w4=pS
else if(xt>=uS .and. xt<sR)then
  w1=rhoSR
  w2=uS
  w3=vR
  w4=pS
else if(xt>=sR)then
  w1=rhoR
  w2=uR
  w3=vR
  w4=pR
end if

end subroutine solutionSample1

!Two Rarefaction
subroutine solutionSample2(xt)
real*8::xt

if(xt<=sHL)then
  w1=rhoL
  w2=uL
  w3=vL
  w4=pL
else if(xt>sHL .and. xt<sTL)then
  w1=rhoL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gamul)
  w2=2.*gamil*(aL+0.5*(gam-1.)*uL+xt)
  w3=vL
  w4=pL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gam*gamul)
else if(xt>=sTL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=vL
  w4=pS
else if(xt>=uS .and. xt<sTR)then
  w1=rhoSR
  w2=uS
  w3=vR
  w4=pS
else if(xt>=sTR .and. xt<sHR)then
  w1=rhoR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gamul)
  w2=2.*gamil*(-aR+0.5*(gam-1.)*uR+xt)
  w3=vR
  w4=pR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gam*gamul)
else if(xt>=sHR)then
  w1=rhoR
  w2=uR
  w3=vR
  w4=pR
end if

end subroutine solutionSample2

!left shock, right rarefaction
subroutine solutionSample3(xt)
real*8::xt

if(xt<=sL)then
  w1=rhoL
  w2=uL
  w3=vL
  w4=pL
else if(xt>sL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=vL
  w4=pS
else if(xt>=uS .and. xt<sTR)then
  w1=rhoSR
  w2=uS
  w3=vR
  w4=pS
else if(xt>=sTR .and. xt<sHR)then
  w1=rhoR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gamul)
  w2=2.*gamil*(-aR+0.5*(gam-1.)*uR+xt)
  w3=vR
  w4=pR*(2*gamil-(gamel/aR)*(uR-xt))**(2.*gam*gamul)
else if(xt>=sHR)then
  w1=rhoR
  w2=uR
  w3=vR
  w4=pR
end if

end subroutine solutionSample3

!left rarefaction, right shock
subroutine solutionSample4(xt)
real*8::xt

if(xt<=sHL)then
  w1=rhoL
  w2=uL
  w3=vL
  w4=pL
else if(xt>sHL .and. xt<sTL)then
  w1=rhoL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gamul)
  w2=2.*gamil*(aL+0.5*(gam-1.)*uL+xt)
  w3=vL
  w4=pL*(2*gamil+(gamel/aL)*(uL-xt))**(2.*gam*gamul)
else if(xt>=sTL .and. xt<uS)then
  w1=rhoSL
  w2=uS
  w3=vL
  w4=pS
else if(xt>=uS .and. xt<sR)then
  w1=rhoSR
  w2=uS
  w3=vR
  w4=pS
else if(xt>=sR)then
  w1=rhoR
  w2=uR
  w3=vR
  w4=pR
end if

end subroutine solutionSample4

subroutine bound()
!Local Variables
integer::jj,kk,ll
do jj=1,nx
  do ll=1,4
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=1,ny
  do ll=1,4
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
end do
end subroutine bound

function gL(px) result(pxx)
!Input Variables
real*8,intent(in)::px
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(px+gamel*pL))
end function gL

function gR(px) result(pxx)
!Input Variables
real*8,intent(in)::px
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(px+gamel*pR))
end function gR

!-----------------------------------------------------------
end module ExactRiemannSolver_mod
