!Roe first order Upwind Scheme for Ieal MHD in 2D
!(Dimensionally split)  
!No Div(B) constraint

module RoeSolver2D_mod
implicit none


integer,parameter::debug=0 !0:off 1:on
integer,parameter::discontinuityDir=1!1:x, 2:y, 3:oblique, 4:partial grid circle 
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic, 5:peridoc in x reflecting in y
integer,parameter::initOption=0
integer,parameter::fluxType=3 !1:Roe 2:HLL 3:HLLI
integer,parameter::HLLOption=0

integer,parameter::KH_test=0 !0:off 1:on
integer,parameter::perturbationType=1!1:Sinusoidal 2:Random


real,parameter::eps=0.00001

integer,parameter::tSkip=50

integer,parameter::nt=200
integer,parameter::nx=150
integer,parameter::ny=150

real*8,parameter::tol=1.d-20
real*8,parameter::min_pres=1.d-5
real*8,parameter::min_dens=1.d-5
real*8,parameter::pi=3.14159265359


real*8::dt,dx,dy !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real*8::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state

!-------------------------------------------------------------
real*8::dens,velx,vely,velz,pres,magx,magy,magz
!-------------------------------------------------------------
real*8::smax,cs,cf,cA
real*8::rhot,ut,vt,wt,pt,byt,bzt,at,ptotL,ptotR,ptot
real*8::lambda(7)
real*8::L1(7),L2(7),L3(7),L4(7),L5(7),L6(7),L7(7)	
real*8::R1(7),R2(7),R3(7),R4(7),R5(7),R6(7),R7(7)
real*8::alpha(7),uright(7),uleft(7)
real*8::alphas,alphaf,betay,betaz,bx,by,bz,bsqr
real*8::BBR,BBL,vvR,vvL,gf1,gf7,gs3,gs5,theta1,theta2,vsq
!-------------------------------------------------------------
real*8::u1(-1:nx+2,-1:ny+2,8),u2(-1:nx+2,-1:ny+2,8),flux(-1:nx+2,8)
real*8::fleft(7),fright(7)
!-------------------------------------------------------------


character(len=20) :: filename

contains 

!------------------------------------------------------
!initialization
!------------------------------------------------------
subroutine init()

!Local Variables
integer::jj,kk,mode
real*8::vparL,vparR,vperpL,vperpR,bparL,bparR,bperpL,bperpR,kz

open(unit=12,file='Output/horcut_fluid.txt')
open(unit=13,file='Output/horcut_B.txt')
open(unit=14,file='Output/vercut_fluid.txt')
open(unit=15,file='Output/vercut_B.txt')
open(unit=16,file='Output/diagcut_fluid.txt')
open(unit=17,file='Output/diagcut_B.txt')
open(unit=0+1000,file='Output/t=0.txt')

print*,'First Order Roe Upwind Method...'

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
cour=0.7

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dt=0.

gam=5._8/3._8
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)


!------------------------------------------------------
!Set left and right states 
!------------------------------------------------------
if(initOption==0)then

bxL=0.0!0.1
bxR=bxL

rhoL=1.0
uL=0.
vL=0.
wL=0.
byL=0.
bzL=0.
pL=20.

rhoR=0.125
uR=0.
vR=0.
wR=0.
byR=0.
bzR=0.
pR=5.
end if

if(initOption==1)then

bxL=5./sqrt(4.*pi)
bxR=bxL

rhoL=1.
uL=10.
vL=0.
wL=0.
byL=5./sqrt(4.*pi)
bzL=0.
pL=20.

rhoR=1.
uR=-10.
vR=0.
wR=0.
byR=5./sqrt(4.*pi)
bzR=0.
pR=1.
end if

if(initOption==2)then
bxL=3./sqrt(4.*pi)
bxR=bxL

rhoL=1.0
uL=0.
vL=0.
wL=0.
byL=5./sqrt(4.*pi)
bzL=0.
pL=1.

rhoR=0.1
uR=0.
vR=0.
wR=0.
byR=2./sqrt(4.*pi)
bzR=0.
pR=10.
end if

if(initOption==3)then
bxL=2./sqrt(4.*pi)
bxR=bxL

rhoL=1.08
uL=1.2
vL=0.01
wL=0.5
byL=3.6/sqrt(4.*pi)
bzL=2./sqrt(4.*pi)
pL=0.95

rhoR=1.
uR=0.
vR=0.
wR=0.
byR=4./sqrt(4.*pi)
bzR=2./sqrt(4.*pi)
pR=1.
end if

if(initOption==4)then
bxL=3./sqrt(4.*pi)
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=6./sqrt(4.*pi)
bzL=0.
pL=1.

rhoR=0.1
uR=0.
vR=2.
wR=1.
byR=1./sqrt(4.*pi)
bzR=0.
pR=10.
end if


if(initOption==5)then
bxL=0.
bxR=bxL

rhoL=0.1
uL=50.
vL=0.
wL=0.
byL=-1./sqrt(4.*pi)
bzL=-2./sqrt(4.*pi)
pL=0.4

rhoR=0.1
uR=0.
vR=0.
wR=0.
byR=1./sqrt(4.*pi)
bzR=2./sqrt(4.*pi)
pR=0.2
end if

if(initOption==6)then
bxL=0.
bxR=bxL

rhoL=1.
uL=-1.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=1.
uR=1.
vR=0.
wR=0.
byR=1.
bzR=0.
pR=1.
end if

if(initOption==7)then
bxL=1.
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=0.2
uR=0.
vR=0.
wR=0.
byR=0.
bzR=0.
pR=0.1
end if


if(initOption==8)then
bxL=1.3
bxR=bxL

rhoL=0.4
uL=-0.66991
vL=0.98263
wL=0.
byL=0.0025293
bzL=0.
pL=0.52467

rhoR=1.
uR=0.
vR=0.
wR=0.
byR=1.
bzR=0.
pR=1.
end if

if(initOption==9)then
bxL=0.75
bxR=bxL

rhoL=0.65
uL=0.667
vL=-0.257
wL=0.
byL=0.55
bzL=0.
pL=0.5

rhoR=1.
uR=0.4
vR=-0.94
wR=0.
byR=0.
bzR=0.
pR=0.75
end if

if(initOption==10)then
bxL=0.7
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=0.
bzL=0.
pL=1.

rhoR=0.3
uR=0.
vR=0.
wR=1.
byR=1.
bzR=0.
pR=0.2
end if


if(initOption==11)then
bxL=0.75
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=0.125
uR=0.
vR=0.
wR=0.
byR=-1.
bzR=0.
pR=0.1
end if


if(initOption==12)then
bxL=1.3
bxR=bxL

rhoL=1.
uL=0.
vL=0.
wL=0.
byL=1.
bzL=0.
pL=1.

rhoR=0.4
uR=0.
vR=0.
wR=0.
byR=-1.
bzR=0.
pR=0.4

end if



if(discontinuityDir==1)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if(jj<nx/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
  end if
 end do
 end do

else if(discontinuityDir==2)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if(kk<ny/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*vL
    u1(jj,kk,3)=rhoL*uL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=byL
    u1(jj,kk,6)=bxL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*vR
    u1(jj,kk,3)=rhoR*uR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=byR
    u1(jj,kk,6)=bxR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
  end if
 end do
 end do
else if(discontinuityDir==3)then

 vparL=uL
 vperpL=vL
 bparL=bxL
 bperpL=byL
 vparR=uR
 vperpR=vR
 bparR=bxR
 bperpR=byR

 uL=(vparL-vperpL)/sqrt(2.)
 vL=(vparL+vperpL)/sqrt(2.)
 bxL=(bparL-bperpL)/sqrt(2.)
 byL=(bparL+bperpL)/sqrt(2.)
 uR=(vparR-vperpR)/sqrt(2.)
 vR=(vparR+vperpR)/sqrt(2.)
 bxR=(bparR-bperpR)/sqrt(2.)
 byR=(bparR+bperpR)/sqrt(2.)


 do jj=-1,nx+2
 do kk=-1,ny+2
  if(kk<=nx/2-jj) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
  end if
 end do
 end do

else if(discontinuityDir==4)then
 do jj=-1,nx+2
 do kk=-1,ny+2
  if((jj-nx/2)*(jj-nx/2)+(kk-ny/2)*(kk-ny/2)<=(nx/4)*(nx/4)) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
  end if
 end do
 end do

end if

if(KH_test==1)then
!Kelvin-Helmholtz Instability 
rhoL=1.0
uL=-1.0
vL=0.0
bxL=0.0!0.001!0./sqrt(4.*pi)
byL=0.!/sqrt(4.*pi)
bzL=0.
wL=0.
pL=10.0

rhoR=2.0
uR=1.0
vR=0.0
wR=0.
bxR=0.0!0.001!0./sqrt(4.*pi)
byR=0.!/sqrt(4.*pi)
bzR=0.
pR=10.0

mode=10

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2
 x=xmin+(jj-0.5)*dx
 do kk=-1,ny+2
  y=ymin+(kk-0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=0.05*(ymax-ymin)*sin(kz*x)*exp(-300.*(kk-ny/2)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.05*(ymax-ymin)*rand(0)*exp(-300.*(kk-ny/2)**2.)
  end if
  vR=vL
  if(kk<ny/2) then !Below Interface
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=bxL
    u1(jj,kk,6)=byL
    u1(jj,kk,7)=bzL  
    u1(jj,kk,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=bxR
    u1(jj,kk,6)=byR
    u1(jj,kk,7)=bzR  
    u1(jj,kk,8)=0.5*rhoR*(uR*uR+vL*vL+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
  end if
 end do   
end do
end if

u2=u1

!-------------------------------
call bound()
u1=u2
!-------------------------------


print*,'Riemann Problem'
print*,'Left State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoL,uL,vL,wL,bxL,byL,bzL,pL
print*,'Right State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoR,uR,vR,wR,bxR,byR,bzR,pR

call vercut()
call horcut()
call diagcut()
call fileOutput(0+1000)

print*,'Initialization done...'


end subroutine init


subroutine computeFluxRoe(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk

!Local Variables
real*8::bxLR,LdotR(7),bxx,byy,bzz,att


do jj=-1,n+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=u1(jj,ind,1)
uL=u1(jj,ind,2)/u1(jj,ind,1)
vL=u1(jj,ind,3)/u1(jj,ind,1)
wL=u1(jj,ind,4)/u1(jj,ind,1)
bxL=u1(jj,ind,5)
byL=u1(jj,ind,6)
bzL=u1(jj,ind,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(jj,ind,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=u1(jj,ind,1)
uleft(2)=u1(jj,ind,2)
uleft(3)=u1(jj,ind,3)
uleft(4)=u1(jj,ind,4)
uleft(5)=u1(jj,ind,6)
uleft(6)=u1(jj,ind,7)
uleft(7)=u1(jj,ind,8)

rhoR=u1(jj+1,ind,1)
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
wR=u1(jj+1,ind,4)/u1(jj+1,ind,1)
bxR=u1(jj+1,ind,5)
byR=u1(jj+1,ind,6)
bzR=u1(jj+1,ind,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(jj+1,ind,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(jj+1,ind,1)
uright(2)=u1(jj+1,ind,2)
uright(3)=u1(jj+1,ind,3)
uright(4)=u1(jj+1,ind,4)
uright(5)=u1(jj+1,ind,6)
uright(6)=u1(jj+1,ind,7)
uright(7)=u1(jj+1,ind,8)

!--------------------------------
else if(dir==2)then
!--------------------------------


rhoL=u1(ind,jj,1)
uL=u1(ind,jj,3)/u1(ind,jj,1)
vL=u1(ind,jj,2)/u1(ind,jj,1)
wL=u1(ind,jj,4)/u1(ind,jj,1)
bxL=u1(ind,jj,6)
byL=u1(ind,jj,5)
bzL=u1(ind,jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(ind,jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL


uleft(1)=u1(ind,jj,1)
uleft(2)=u1(ind,jj,3)
uleft(3)=u1(ind,jj,2)
uleft(4)=u1(ind,jj,4)
uleft(5)=u1(ind,jj,5)
uleft(6)=u1(ind,jj,7)
uleft(7)=u1(ind,jj,8)

rhoR=u1(ind,jj+1,1)
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
wR=u1(ind,jj+1,4)/u1(ind,jj+1,1)
bxR=u1(ind,jj+1,6)
byR=u1(ind,jj+1,5)
bzR=u1(ind,jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(ind,jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(ind,jj+1,1)
uright(2)=u1(ind,jj+1,3)
uright(3)=u1(ind,jj+1,2)
uright(4)=u1(ind,jj+1,4)
uright(5)=u1(ind,jj+1,5)
uright(6)=u1(ind,jj+1,7)
uright(7)=u1(ind,jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL-vL*bxL
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR-vR*bxR
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if

!---------------------------------------------
!Characteristic Eigenvalues: lambda_1..7
!---------------------------------------------
lambda(1)=ut-cf
lambda(2)=ut-cA
lambda(3)=ut-cs
lambda(4)=ut
lambda(5)=ut+cs
lambda(6)=ut+cA
lambda(7)=ut+cf

if(debug==1)then
  print*,'Roe average Eigenvalues=',lambda(1),lambda(2),lambda(3) &
        ,lambda(4),lambda(5),lambda(6),lambda(7)
end if

!-------------------------------------------------
!Roe Jacobian Matrix Eigenvectors: R_1..7, L_1..7
!-------------------------------------------------
if(abs(byt)<tol .and. abs(bzt)<tol)then
 betay=1./sqrt(2.)
 betaz=1./sqrt(2.)
else
 betay=byt/sqrt(byt*byt+bzt*bzt)
 betaz=bzt/sqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
 alphas=1.
 alphaf=1.
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
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


R1(1)=alphaf
R1(2)=alphaf*(ut-cf)
R1(3)=alphaf*vt+alphas*betay*cA*sgn(bxLR)
R1(4)=alphaf*wt+alphas*betaz*cA*sgn(bxLR)
R1(5)=(alphas*betay*cf)/sqrt(rhot)
R1(6)=(alphas*betaz*cf)/sqrt(rhot)
R1(7)=alphaf*0.5*vsq+gf1

R2(1)=0.
R2(2)=0.
R2(3)=betaz*sgn(bxLR)
R2(4)=-betay*sgn(bxLR)
R2(5)=betaz/sqrt(rhot)
R2(6)=-betay/sqrt(rhot)
R2(7)=(betaz*vt-betay*wt)*sgn(bxLR)

R3(1)=alphas
R3(2)=alphas*(ut-cs)
R3(3)=alphas*vt-alphaf*betay*at*sgn(bxLR)
R3(4)=alphas*wt-alphaf*betaz*at*sgn(bxLR)
R3(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R3(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R3(7)=alphas*0.5*vsq+gs3

R4(1)=1.
R4(2)=ut
R4(3)=vt
R4(4)=wt
R4(5)=0.
R4(6)=0.
R4(7)=0.5*vsq

R5(1)=alphas
R5(2)=alphas*(ut+cs)
R5(3)=alphas*vt+alphaf*betay*at*sgn(bxLR)
R5(4)=alphas*wt+alphaf*betaz*at*sgn(bxLR)
R5(5)=-(alphaf*betay*at*at)/(cf*sqrt(rhot))
R5(6)=-(alphaf*betaz*at*at)/(cf*sqrt(rhot))
R5(7)=alphas*0.5*vsq+gs5

R6(1)=0.
R6(2)=0.
R6(3)=-betaz*sgn(bxLR)
R6(4)=betay*sgn(bxLR)
R6(5)=betaz/sqrt(rhot)
R6(6)=-betay/sqrt(rhot)
R6(7)=-(betaz*vt-betay*wt)*sgn(bxLR)

R7(1)=alphaf
R7(2)=alphaf*(ut+cf)
R7(3)=alphaf*vt-alphas*betay*cA*sgn(bxLR)
R7(4)=alphaf*wt-alphas*betaz*cA*sgn(bxLR)
R7(5)=(alphas*betay*cf)/sqrt(rhot)
R7(6)=(alphas*betaz*cf)/sqrt(rhot)
R7(7)=alphaf*0.5*vsq+gf7


L1(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     +(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L1(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      -(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L1(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      +(1./theta2)*0.5*alphas*betay*cs

L1(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      +(1./theta2)*0.5*alphas*betaz*cs

L1(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L1(7)=(1./theta1)*0.5*alphaf*at*at

L2(1)=-0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L2(2)=0.
L2(3)=0.5*betaz*sgn(bxLR)
L2(4)=-0.5*betay*sgn(bxLR)
L2(5)=0.5*betaz*sqrt(rhot)
L2(6)=-0.5*betay*sqrt(rhot)
L2(7)=0.

L3(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     +(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L3(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      -(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L3(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      -(1./theta2)*0.5*alphaf*betay*cf

L3(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      -(1./theta2)*0.5*alphaf*betaz*cf

L3(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L3(7)=(1./theta1)*0.5*alphas*cf*cf


L4(1)=1.-(1./theta1)*0.5*vsq*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(2)=(1./theta1)*ut*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(3)=(1./theta1)*vt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )

L4(4)=(1./theta1)*wt*(alphaf*alphaf*at*at &
      +alphas*alphas*cf*cf )


L4(5)=(1./theta1)*alphaf*alphas*betay*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(6)=(1./theta1)*alphaf*alphas*betaz*cf*(cf*cf-cs*cs)*sqrt(rhot)

L4(7)=-(1./theta1)*(alphaf*alphaf*at*at+alphas*alphas*cf*cf)


L5(1)=(1./theta1)*0.25*alphas*cf*cf*vsq & 
     -(1./theta2)*( 0.5*alphas*cA*ut*sgn(bxLR) &
     +0.5*alphaf*cf*(betay*vt+betaz*wt) )

L5(2)=-(1./theta1)*0.5*alphas*cf*cf*ut & 
      +(1./theta2)*0.5*alphas*cA*sgn(bxLR)

L5(3)=-(1./theta1)*0.5*alphas*cf*cf*vt & 
      +(1./theta2)*0.5*alphaf*betay*cf

L5(4)=-(1./theta1)*0.5*alphas*cf*cf*wt & 
      +(1./theta2)*0.5*alphaf*betaz*cf

L5(5)=-(1./theta1)*0.5*alphaf*betay*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(6)=-(1./theta1)*0.5*alphaf*betaz*cf* &
      (cf*cf+gamul*(2.-gam)*at*at)*sqrt(rhot)

L5(7)=(1./theta1)*0.5*alphas*cf*cf


L6(1)=0.5*(betaz*vt-betay*wt)*sgn(bxLR)
L6(2)=0.
L6(3)=-0.5*betaz*sgn(bxLR)
L6(4)=0.5*betay*sgn(bxLR)
L6(5)=0.5*betaz*sqrt(rhot)
L6(6)=-0.5*betay*sqrt(rhot)
L6(7)=0.


L7(1)=(1./theta1)*0.25*alphaf*at*at*vsq & 
     -(1./theta2)*( 0.5*alphaf*at*ut*sgn(bxLR) &
     -0.5*alphas*cs*(betay*vt+betaz*wt) )

L7(2)=-(1./theta1)*0.5*alphaf*at*at*ut & 
      +(1./theta2)*0.5*alphaf*at*sgn(bxLR)

L7(3)=-(1./theta1)*0.5*alphaf*at*at*vt & 
      -(1./theta2)*0.5*alphas*betay*cs

L7(4)=-(1./theta1)*0.5*alphaf*at*at*wt & 
      -(1./theta2)*0.5*alphas*betaz*cs

L7(5)=(1./theta1)*0.5*alphas*betay*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(6)=(1./theta1)*0.5*alphas*betaz*cf* &
      (cs*cs+gamul*(2.-gam)*at*at)*sqrt(rhot)

L7(7)=(1./theta1)*0.5*alphaf*at*at

if((at*at-cA*cA)>tol)then
 if(byt<=tol .and. bzt<tol)then
  R3=-1.*R3
  R5=-1.*R5
  L3=-1.*L3
  L5=-1.*L5
 end if
else 
 if(byt<=tol .and. bzt<tol)then
  R1=-1.*R1
  R7=-1.*R7
  L1=-1.*L1
  L7=-1.*L7
 end if
end if

LdotR=0.
do kk=1,7
LdotR(1)=LdotR(1)+L1(kk)*R1(kk)
LdotR(2)=LdotR(2)+L2(kk)*R2(kk)
LdotR(3)=LdotR(3)+L3(kk)*R3(kk)
LdotR(4)=LdotR(4)+L4(kk)*R4(kk)
LdotR(5)=LdotR(5)+L5(kk)*R5(kk)
LdotR(6)=LdotR(6)+L6(kk)*R6(kk)
LdotR(7)=LdotR(7)+L7(kk)*R7(kk)
end do

if(dir==1)then
 !print*,'Cell:',jj,ind
else
!print*,'Cell:',ind,jj
end if
!print*,'||L,R||=',LdotR

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1..7
!---------------------------------------------
alpha=0.

do kk=1,7  
   alpha(1)=alpha(1)+L1(kk)*(uright(kk)-uleft(kk))  
   alpha(2)=alpha(2)+L2(kk)*(uright(kk)-uleft(kk))
   alpha(3)=alpha(3)+L3(kk)*(uright(kk)-uleft(kk))
   alpha(4)=alpha(4)+L4(kk)*(uright(kk)-uleft(kk))
   alpha(5)=alpha(5)+L5(kk)*(uright(kk)-uleft(kk))
   alpha(6)=alpha(6)+L6(kk)*(uright(kk)-uleft(kk))
   alpha(7)=alpha(7)+L7(kk)*(uright(kk)-uleft(kk))
end do


!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7
  flux(jj,kk)=0.5*(fleft(kk)+fright(kk))-0.5*(     &
         abs(lambda(1))*alpha(1)*R1(kk) &
        +abs(lambda(2))*alpha(2)*R2(kk) &
        +abs(lambda(3))*alpha(3)*R3(kk) &  
        +abs(lambda(4))*alpha(4)*R4(kk) &
        +abs(lambda(5))*alpha(5)*R5(kk) &
        +abs(lambda(6))*alpha(6)*R6(kk) &
        +abs(lambda(7))*alpha(7)*R7(kk) )

end do


if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 


end do

end subroutine computeFluxRoe

subroutine computeFluxHLL(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk

!Local Variables
real*8::bxLR,bxx,byy,bzz,att,SL,SR,lambdaL,lambdaR


do jj=-1,n+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=u1(jj,ind,1)
uL=u1(jj,ind,2)/u1(jj,ind,1)
vL=u1(jj,ind,3)/u1(jj,ind,1)
wL=u1(jj,ind,4)/u1(jj,ind,1)
bxL=u1(jj,ind,5)
byL=u1(jj,ind,6)
bzL=u1(jj,ind,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(jj,ind,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=u1(jj,ind,1)
uleft(2)=u1(jj,ind,2)
uleft(3)=u1(jj,ind,3)
uleft(4)=u1(jj,ind,4)
uleft(5)=u1(jj,ind,6)
uleft(6)=u1(jj,ind,7)
uleft(7)=u1(jj,ind,8)

rhoR=u1(jj+1,ind,1)
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
wR=u1(jj+1,ind,4)/u1(jj+1,ind,1)
bxR=u1(jj+1,ind,5)
byR=u1(jj+1,ind,6)
bzR=u1(jj+1,ind,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(jj+1,ind,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(jj+1,ind,1)
uright(2)=u1(jj+1,ind,2)
uright(3)=u1(jj+1,ind,3)
uright(4)=u1(jj+1,ind,4)
uright(5)=u1(jj+1,ind,6)
uright(6)=u1(jj+1,ind,7)
uright(7)=u1(jj+1,ind,8)

!--------------------------------
else if(dir==2)then
!--------------------------------

rhoL=u1(ind,jj,1)
uL=u1(ind,jj,3)/u1(ind,jj,1)
vL=u1(ind,jj,2)/u1(ind,jj,1)
wL=u1(ind,jj,4)/u1(ind,jj,1)
bxL=u1(ind,jj,6)
byL=u1(ind,jj,5)
bzL=u1(ind,jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(ind,jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL


uleft(1)=u1(ind,jj,1)
uleft(2)=u1(ind,jj,3)
uleft(3)=u1(ind,jj,2)
uleft(4)=u1(ind,jj,4)
uleft(5)=u1(ind,jj,5)
uleft(6)=u1(ind,jj,7)
uleft(7)=u1(ind,jj,8)

rhoR=u1(ind,jj+1,1)
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
wR=u1(ind,jj+1,4)/u1(ind,jj+1,1)
bxR=u1(ind,jj+1,6)
byR=u1(ind,jj+1,5)
bzR=u1(ind,jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(ind,jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(ind,jj+1,1)
uright(2)=u1(ind,jj+1,3)
uright(3)=u1(ind,jj+1,2)
uright(4)=u1(ind,jj+1,4)
uright(5)=u1(ind,jj+1,5)
uright(6)=u1(ind,jj+1,7)
uright(7)=u1(ind,jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL-vL*bxL
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR-vR*bxR
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if


!----------------------------------------------
!Extremal Wave Speeds
!----------------------------------------------

lambdaL=0.5*(aL*aL+(BBL/rhoL)) &
+0.5*sqrt( (aL*aL+(BBL/rhoL))**2.-4._8*aL*aL*bxL*bxL/rhoL)
lambdaL=uL-sqrt(lambdaL)

lambdaR=0.5*(aR*aR+(BBR/rhoR)) &
+0.5*sqrt( (aR*aR+(BBR/rhoR))**2.-4._8*aR*aR*bxR*bxR/rhoR)
lambdaR=uR+sqrt(lambdaR)

if(HLLOption==0)then
 SL=min(lambdaL,ut-cf)
 SR=max(lambdaR,ut+cf)
else if(HLLOption==1)then
 SL=min(lambdaL,ut-cf,-eps)
 SR=max(lambdaR,ut+cf,eps)
end if

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7 
  if(SL>=0.)then
    flux(jj,kk)=fleft(kk)
  else if(SR<=0.)then
    flux(jj,kk)=fright(kk)
  else if(SL<0. .and. SR>0.)then
    flux(jj,kk)=SR*fleft(kk)-SL*fright(kk)+SL*SR*(uright(kk)-uleft(kk))
    flux(jj,kk)=flux(jj,kk)/(SR-SL)
  end if
end do


if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 

end do

end subroutine computeFluxHLL

subroutine computeFluxHLLI(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk,ll,mm

!Local Variables
real*8::bxLR,bxx,byy,bzz,att,SL,SR,lambdaL,lambdaR
real*8::VLx(7),VRx(7),DUDV(7,7),Rx(7,7),Lx(7,7),Rxu(7,7),Lxu(7,7)
real*8::delta(7)

do jj=-1,n+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

!set left and right states for local Riemann problem at cell interface
!x-sweeps
!--------------------------------
if(dir==1)then
!--------------------------------
rhoL=u1(jj,ind,1)
uL=u1(jj,ind,2)/u1(jj,ind,1)
vL=u1(jj,ind,3)/u1(jj,ind,1)
wL=u1(jj,ind,4)/u1(jj,ind,1)
bxL=u1(jj,ind,5)
byL=u1(jj,ind,6)
bzL=u1(jj,ind,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(jj,ind,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL

uleft(1)=u1(jj,ind,1)
uleft(2)=u1(jj,ind,2)
uleft(3)=u1(jj,ind,3)
uleft(4)=u1(jj,ind,4)
uleft(5)=u1(jj,ind,6)
uleft(6)=u1(jj,ind,7)
uleft(7)=u1(jj,ind,8)

rhoR=u1(jj+1,ind,1)
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
wR=u1(jj+1,ind,4)/u1(jj+1,ind,1)
bxR=u1(jj+1,ind,5)
byR=u1(jj+1,ind,6)
bzR=u1(jj+1,ind,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(jj+1,ind,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(jj+1,ind,1)
uright(2)=u1(jj+1,ind,2)
uright(3)=u1(jj+1,ind,3)
uright(4)=u1(jj+1,ind,4)
uright(5)=u1(jj+1,ind,6)
uright(6)=u1(jj+1,ind,7)
uright(7)=u1(jj+1,ind,8)

!--------------------------------
else if(dir==2)then
!--------------------------------

rhoL=u1(ind,jj,1)
uL=u1(ind,jj,3)/u1(ind,jj,1)
vL=u1(ind,jj,2)/u1(ind,jj,1)
wL=u1(ind,jj,4)/u1(ind,jj,1)
bxL=u1(ind,jj,6)
byL=u1(ind,jj,5)
bzL=u1(ind,jj,7)

BBL=bxL*bxL+byL*byL+bzL*bzL
vvL=uL*uL+vL*vL+wL*wL
pL=(gam-1.)*(u1(ind,jj,8)-0.5*rhoL*vvL-0.5*BBL)
aL=sqrt(gam*pL/rhoL)
ptotL=pL+0.5*BBL


uleft(1)=u1(ind,jj,1)
uleft(2)=u1(ind,jj,3)
uleft(3)=u1(ind,jj,2)
uleft(4)=u1(ind,jj,4)
uleft(5)=u1(ind,jj,5)
uleft(6)=u1(ind,jj,7)
uleft(7)=u1(ind,jj,8)

rhoR=u1(ind,jj+1,1)
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
wR=u1(ind,jj+1,4)/u1(ind,jj+1,1)
bxR=u1(ind,jj+1,6)
byR=u1(ind,jj+1,5)
bzR=u1(ind,jj+1,7)

BBR=bxR*bxR+byR*byR+bzR*bzR
vvR=uR*uR+vR*vR+wR*wR
pR=(gam-1.)*(u1(ind,jj+1,8)-0.5*rhoR*vvR-0.5*BBR)
aR=sqrt(gam*pR/rhoR)
ptotR=pR+0.5*BBR

uright(1)=u1(ind,jj+1,1)
uright(2)=u1(ind,jj+1,3)
uright(3)=u1(ind,jj+1,2)
uright(4)=u1(ind,jj+1,4)
uright(5)=u1(ind,jj+1,5)
uright(6)=u1(ind,jj+1,7)
uright(7)=u1(ind,jj+1,8)

end if

fleft(1)=rhoL*uL
fleft(2)=rhoL*uL*uL+ptotL-bxL*bxL
fleft(3)=rhoL*uL*vL-bxL*byL
fleft(4)=rhoL*uL*wL-bxL*bzL
fleft(5)=uL*byL-vL*bxL
fleft(6)=uL*bzL-wL*bxL
fleft(7)=uL*(uleft(7)+ptotL) &
      -bxL*(bxL*uL+byL*vL+bzL*wL)


fright(1)=rhoR*uR
fright(2)=rhoR*uR*uR+ptotR-bxR*bxR
fright(3)=rhoR*uR*vR-bxR*byR
fright(4)=rhoR*uR*wR-bxR*bzR
fright(5)=uR*byR-vR*bxR
fright(6)=uR*bzR-wR*bxR
fright(7)=uR*(uright(7)+ptotR) &
      -bxR*(bxR*uR+byR*vR+bzR*wR)

VLx(1)=rhoL
VLx(2)=uL
VLx(3)=vL
VLx(4)=wL
VLx(5)=byL
VLx(6)=bzL
VLx(7)=pL

VRx(1)=rhoR
VRx(2)=uR
VRx(3)=vR
VRx(4)=wR
VRx(5)=byR
VRx(6)=bzR
VRx(7)=pR

if(debug==1)then
  if(dir==1)then
   print*,'Cell:',jj,ind
  else if(dir==2)then
   print*,'Cell:',ind,jj
  end if
  print*,'Local RP Left State (rho,u,v,w,bx,by,bz,p) =',rhoL,uL,vL,wL,bxL,byL,bzL,pL
  print*,'Local RP Right State (rho,u,v,w,bx,by,bz,p) =',rhoR,uR,vR,wR,bxR,byR,bzR,pR
end if

!-----------------------------
bxLR=bxL
!-----------------------------

!----------------------------------------------
!Roe-Averaged State Variables (arithmetic mean):
!----------------------------------------------
rhot=0.5_8*(rhoL+rhoR)
ut=0.5_8*(uL+uR)
vt=0.5_8*(vL+vR)
wt=0.5_8*(wL+wR)
byt=0.5_8*(byL+byR)
bzt=0.5_8*(bzL+bzR)
ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
pt=ptot-0.5_8*(bxLR*bxLR+byt*byt+bzt*bzt)
att=gam*pt/rhot
at=sqrt(att)

if(debug==1)then
  print*,'Roe Average State =',rhot,ut,vt,wt,bxLR,byt,bzt,pt
end if

!-------------------------------------------------
!Alfven and magnetosonic wave speeds: c_s,c_A,c_f
!-------------------------------------------------
bxx=(bxLR*bxLR)/rhot
byy=(byt*byt)/rhot
bzz=(bzt*bzt)/rhot
bsqr=bxx+byy+bzz
vsq=ut*ut+vt*vt+wt*wt

if(debug==1)then
 print*,'bx^2,by^2,bz^2=',bxx,byy,bzz
end if

cA=bxx
cA=sqrt(bxx)
cs=0.5_8*(att+bsqr-sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
cf=0.5_8*(att+bsqr+sqrt((att+bsqr)*(att+bsqr)-4._8*att*bxx))
!cf=att+bsqr-cs
cs=sqrt(cs)
cf=sqrt(cf)


if(debug==1)then
  print*,'Slow,Alfven,Fast,Sound wave speeds=',cs,cA,cf,at
end if


!---------------------------------------------
!Roe Averaged Characteristic Eigenvalues: lambda_1..7
!---------------------------------------------
lambda(1)=ut-cf
lambda(2)=ut-cA
lambda(3)=ut-cs
lambda(4)=ut
lambda(5)=ut+cs
lambda(6)=ut+cA
lambda(7)=ut+cf


!-------------------------------------------------
!Transformation Matrix for V->U
!-------------------------------------------------
DUDV(1,1)=1.
DUDV(1,2)=0.
DUDV(1,3)=0.
DUDV(1,4)=0.
DUDV(1,5)=0.
DUDV(1,6)=0.
DUDV(1,7)=0.

DUDV(2,1)=ut
DUDV(2,2)=rhot
DUDV(2,3)=0.
DUDV(2,4)=0.
DUDV(2,5)=0.
DUDV(2,6)=0.
DUDV(2,7)=0.

DUDV(3,1)=vt
DUDV(3,2)=0.
DUDV(3,3)=rhot
DUDV(3,4)=0.
DUDV(3,5)=0.
DUDV(3,6)=0.
DUDV(3,7)=0.

DUDV(4,1)=wt
DUDV(4,2)=0.
DUDV(4,3)=0.
DUDV(4,4)=rhot
DUDV(4,5)=0.
DUDV(4,6)=0.
DUDV(4,7)=0.

DUDV(5,1)=0.
DUDV(5,2)=0.
DUDV(5,3)=0.
DUDV(5,4)=0.
DUDV(5,5)=1.
DUDV(5,6)=0.
DUDV(5,7)=0.

DUDV(6,1)=0.
DUDV(6,2)=0.
DUDV(6,3)=0.
DUDV(6,4)=0.
DUDV(6,5)=0.
DUDV(6,6)=1.
DUDV(6,7)=0.

DUDV(7,1)=0.5*vsq
DUDV(7,2)=rhot*ut
DUDV(7,3)=rhot*vt
DUDV(7,4)=rhot*wt
DUDV(7,5)=byt
DUDV(7,6)=bzt
DUDV(7,7)=gamul


!----------------------------------------------------
!Right and Left Eigenvectors for primitive variables:
!----------------------------------------------------
if(abs(byt)<tol .and. abs(bzt)<tol)then
 betay=1./sqrt(2.)
 betaz=1./sqrt(2.)
else
 betay=byt/sqrt(byt*byt+bzt*bzt)
 betaz=bzt/sqrt(byt*byt+bzt*bzt)
end if

if(debug==1)then
  print*,'betay,betaz=',betay,betaz
end if

if(abs(byt)<tol .and. abs(bzt)<tol .and. abs(cA-at)<tol)then
 alphas=1.
 alphaf=1.
else
 !!***NOTE**** Inside the square roots in renormalization constants, 
 !! can take abs of numerator to avoid negative values caused by numerical precision errors
 alphas=sqrt(abs(cf*cf-at*at)/(cf*cf-cs*cs))
 alphaf=sqrt(abs(cf*cf-cA*cA)/(cf*cf-cs*cs)) 
end if

if(debug==1)then
  print*,'alphas,alphaf=',alphas,alphaf
end if

Rx(1,1)=alphaf*rhot
Rx(2,1)=-alphaf*cf
Rx(3,1)=alphas*cs*betay*sgn(bx)
Rx(4,1)=alphas*cs*betaz*sgn(bx)
Rx(5,1)=alphas*betay*at*sqrt(rhot)
Rx(6,1)=alphas*betaz*at*sqrt(rhot)
Rx(7,1)=alphaf*rhot*at*at

Rx(1,2)=0.
Rx(2,2)=0.
Rx(3,2)=-betaz
Rx(4,2)=betay
Rx(5,2)=-betaz*sqrt(rhot)*sgn(bxLR)
Rx(6,2)=betay*sqrt(rhot)*sgn(bxLR)
Rx(7,2)=0.

Rx(1,3)=alphas*rhot
Rx(2,3)=-alphas*cs
Rx(3,3)=-alphaf*betay*cf*sgn(bxLR)
Rx(4,3)=-alphaf*betaz*cf*sgn(bxLR)
Rx(5,3)=-alphaf*betay*at*sqrt(rhot)
Rx(6,3)=-alphaf*betaz*at*sqrt(rhot)
Rx(7,3)=alphas*rhot*at*at

Rx(1,4)=1.
Rx(2,4)=0.
Rx(3,4)=0.
Rx(4,4)=0.
Rx(5,4)=0.
Rx(6,4)=0.
Rx(7,4)=0.

Rx(1,5)=alphas*rhot
Rx(2,5)=alphas*cs
Rx(3,5)=alphaf*betay*cf*sgn(bxLR)
Rx(4,5)=alphaf*betaz*cf*sgn(bxLR)
Rx(5,5)=-alphaf*betay*at*sqrt(rhot)
Rx(6,5)=-alphaf*betaz*at*sqrt(rhot)
Rx(7,5)=alphas*rhot*at*at

Rx(1,6)=0.
Rx(2,6)=0.
Rx(3,6)=betaz
Rx(4,6)=-betay
Rx(5,6)=-betaz*sqrt(rhot)*sgn(bxLR)
Rx(6,6)=betay*sqrt(rhot)*sgn(bxLR)
Rx(7,6)=0.

Rx(1,7)=alphaf*rhot
Rx(2,7)=alphaf*cf
Rx(3,7)=-alphas*cs*betay*sgn(bxLR)
Rx(4,7)=-alphas*cs*betaz*sgn(bxLR)
Rx(5,7)=alphas*betay*at*sqrt(rhot)
Rx(6,7)=alphas*betaz*at*sqrt(rhot)
Rx(7,7)=alphaf*rhot*at*at


Lx(1,1)=0.
Lx(1,2)=-alphaf*cf/(2.*at*at)
Lx(1,3)=alphas*cs*betay*sgn(bxLR)/(2.*at*at)
Lx(1,4)=alphas*cs*betaz*sgn(bxLR)/(2.*at*at)
Lx(1,5)=(alphas*at*betay/sqrt(rhot))/(2.*at*at)
Lx(1,6)=(alphas*at*betaz/sqrt(rhot))/(2.*at*at)
Lx(1,7)=(alphas/rhot)/(2.*at*at)

Lx(2,1)=0.
Lx(2,2)=0.
Lx(2,3)=-0.5*betaz
Lx(2,4)=0.5*betay
Lx(2,5)=-0.5*betaz*sgn(bxLR)/sqrt(rhot)
Lx(2,6)=0.5*betay*sgn(bxLR)/sqrt(rhot)
Lx(2,7)=0.

Lx(3,1)=0.
Lx(3,2)=-alphas*cs/(2.*at*at)
Lx(3,3)=-alphaf*cf*betay*sgn(bxLR)/(2.*at*at)
Lx(3,4)=-alphaf*cf*betaz*sgn(bxLR)/(2.*at*at)
Lx(3,5)=-(alphaf*at*betay/sqrt(rhot))/(2.*at*at)
Lx(3,6)=-(alphaf*at*betaz/sqrt(rhot))/(2.*at*at)
Lx(3,7)=(alphas/rhot)/(2.*at*at)

Lx(4,1)=1.
Lx(4,2)=0.
Lx(4,3)=0.
Lx(4,4)=0.
Lx(4,5)=0.
Lx(4,6)=0.
Lx(4,7)=-1./(at*at)

Lx(5,1)=0.
Lx(5,2)=alphas*cs/(2.*at*at)
Lx(5,3)=alphaf*cf*betay*sgn(bxLR)/(2.*at*at)
Lx(5,4)=alphaf*cf*betaz*sgn(bxLR)/(2.*at*at)
Lx(5,5)=-(alphaf*at*betay/sqrt(rhot))/(2.*at*at)
Lx(5,6)=-(alphaf*at*betaz/sqrt(rhot))/(2.*at*at)
Lx(5,7)=(alphas/rhot)/(2.*at*at)


Lx(6,1)=0.
Lx(6,2)=0.
Lx(6,3)=0.5*betaz
Lx(6,4)=-0.5*betay
Lx(6,5)=-0.5*betaz*sgn(bxLR)/sqrt(rhot)
Lx(6,6)=0.5*betay*sgn(bxLR)/sqrt(rhot)
Lx(6,7)=0.


Lx(7,1)=0.
Lx(7,2)=alphaf*cf/(2.*at*at)
Lx(7,3)=-alphas*cs*betay*sgn(bxLR)/(2.*at*at)
Lx(7,4)=-alphas*cs*betaz*sgn(bxLR)/(2.*at*at)
Lx(7,5)=(alphas*at*betay/sqrt(rhot))/(2.*at*at)
Lx(7,6)=(alphas*at*betaz/sqrt(rhot))/(2.*at*at)
Lx(7,7)=(alphaf/rhot)/(2.*at*at)

Rxu=0.
Lxu=0.
!----------------------------------------------
!Eigenvectors for Conserved Variables:
!----------------------------------------------
do kk=1,7
  do ll=1,7
    do mm=1,7
      Rxu(kk,ll)=Rxu(kk,ll)+DUDV(kk,mm)*Rx(mm,ll)
      Lxu(kk,ll)=Lxu(kk,ll)+Lx(kk,mm)*DUDV(mm,ll)
    end do
  end do
end do


!----------------------------------------------
!Jump in U accross kth wave: alpha_k
!----------------------------------------------
do kk=1,7  
   alpha(1)=alpha(1)+Lxu(1,kk)*(uright(kk)-uleft(kk))  
   alpha(2)=alpha(2)+Lxu(2,kk)*(uright(kk)-uleft(kk))
   alpha(3)=alpha(3)+Lxu(3,kk)*(uright(kk)-uleft(kk))
   alpha(4)=alpha(4)+Lxu(4,kk)*(uright(kk)-uleft(kk))
   alpha(5)=alpha(5)+Lxu(5,kk)*(uright(kk)-uleft(kk))
   alpha(6)=alpha(6)+Lxu(6,kk)*(uright(kk)-uleft(kk))
   alpha(7)=alpha(7)+Lxu(7,kk)*(uright(kk)-uleft(kk))
end do


!----------------------------------------------
!Extremal Wave Speeds
!----------------------------------------------

lambdaL=0.5*(aL*aL+(BBL/rhoL)) &
+0.5*sqrt( (aL*aL+(BBL/rhoL))**2.-4._8*aL*aL*bxL*bxL/rhoL)
lambdaL=uL-sqrt(lambdaL)

lambdaR=0.5*(aR*aR+(BBR/rhoR)) &
+0.5*sqrt( (aR*aR+(BBR/rhoR))**2.-4._8*aR*aR*bxR*bxR/rhoR)
lambdaR=uR+sqrt(lambdaR)

if(HLLOption==0)then
 SL=min(lambdaL,ut-cf)
 SR=max(lambdaR,ut+cf)
else if(HLLOption==1)then
 SL=min(lambdaL,ut-cf,-eps)
 SR=max(lambdaR,ut+cf,eps)
end if

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2 = F(v((x/t)=0))
!------------------------------------------------------------- 
do kk=1,7
  delta(kk)=1.-(min(lambda(kk),0.)/SL)-(max(lambda(kk),0.)/SR)
end do

do kk=1,7 
  if(SL>=0.)then
    flux(jj,kk)=fleft(kk)
  else if(SR<=0.)then
    flux(jj,kk)=fright(kk)
  else if(SL<0. .and. SR>0.)then
    flux(jj,kk)=SR*fleft(kk)-SL*fright(kk)+SL*SR*(uright(kk)-uleft(kk))
    !HLLI correction to entropy wave only
    flux(jj,kk)=flux(jj,kk)-(SR*SL)*delta(4)*Rxu(kk,4)*alpha(4) 
    flux(jj,kk)=flux(jj,kk)/(SR-SL)
  end if
end do


if(debug==1)then
  print*,'flux=',flux(jj,1),flux(jj,2),flux(jj,3), &
                 flux(jj,4),flux(jj,5),flux(jj,6),flux(jj,7)
end if 

end do

end subroutine computeFluxHLLI



subroutine timeStep()
integer::ii,jj
real*8::att

!Compute max wave speed. Determine Roe averaged states 
!at each cell interface in both x and y-directions.
!Compute maximum wave propagation speed from among these.

!reset smax
!smax=0._8
print*,'Start computing time-step.'

!x-sweeps
do jj=1,ny
  do ii=1,nx

   rhoL=u1(ii,jj,1)
   uL=u1(ii,jj,2)/u1(ii,jj,1)
   vL=u1(ii,jj,3)/u1(ii,jj,1)
   wL=u1(ii,jj,4)/u1(ii,jj,1)
   bxL=u1(ii,jj,5)
   byL=u1(ii,jj,6)
   bzL=u1(ii,jj,7)
   BBL=bxL*bxL+byL*byL+bzL*bzL
   vvL=uL*uL+vL*vL+wL*wL
   pL=(gam-1.)*(u1(ii,jj,8)-0.5_8*rhoL*vvL-0.5_8*BBL)
   ptotL=pL+0.5_8*BBL


   rhoR=u1(ii+1,jj,1)
   uR=u1(ii+1,jj,2)/u1(ii+1,jj,1)
   vR=u1(ii+1,jj,3)/u1(ii+1,jj,1)
   wR=u1(ii+1,jj,4)/u1(ii+1,jj,1)
   bxR=u1(ii+1,jj,5)
   byR=u1(ii+1,jj,6)
   bzR=u1(ii+1,jj,7)
   BBR=bxR*bxR+byR*byR+bzR*bzR
   vvR=uR*uR+vR*vR+wR*wR
   pR=(gam-1.)*(u1(ii+1,jj,8)-0.5_8*rhoR*vvR-0.5_8*BBR)
   ptotR=pR+0.5_8*BBR
 
   !print*,'i,j=',ii,jj
   !print*,'rhot,pt',rhot,pt
 

   !----------------------------------------------
   !Roe-Averaged State Variables (arithmetic mean):
   !----------------------------------------------
   rhot=0.5_8*(rhoL+rhoR)
   ut=0.5_8*(uL+uR)
   byt=0.5_8*(byL+byR)
   bzt=0.5_8*(bzL+bzR)
   ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
   pt=ptot-0.5_8*(bxL*bxL+byt*byt+bzt*bzt)
   
   if(pt>0._8 .and. rhot>0._8)then
   at=sqrt(gam*pt/rhot)


   !-------------------------------------------------
   !Fast Mode wave speedc_f
   !-------------------------------------------------
   bx=sqrt(bxL*bxL/rhot)
   by=sqrt(byt*byt/rhot)
   bz=sqrt(bzt*bzt/rhot)
   bsqr=(bxL*bxL+byt*byt+bzt*bzt)/rhot

   cf=0.5_8*(at*at+bsqr+sqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx*bx)) 
   cf=sqrt(cf)

   !find maximum speed
   smax=max(smax,max(abs(ut),abs(vt),abs(wt))+cf)

   end if

  
  end do 
end do


!y-sweeps
do ii=1,nx
  do jj=1,ny
  
   rhoL=u1(ii,jj,1)
   uL=u1(ii,jj,3)/u1(ii,jj,1)
   vL=u1(ii,jj,2)/u1(ii,jj,1)
   wL=u1(ii,jj,4)/u1(ii,jj,1)
   bxL=u1(ii,jj,6)
   byL=u1(ii,jj,5)
   bzL=u1(ii,jj,7)
   BBL=bxL*bxL+byL*byL+bzL*bzL
   vvL=uL*uL+vL*vL+wL*wL
   pL=(gam-1.)*(u1(ii,jj,8)-0.5_8*rhoL*vvL-0.5*BBL)
   ptotL=pL+0.5_8*BBL

   rhoR=u1(ii,jj+1,1)
   uR=u1(ii,jj+1,3)/u1(ii,jj+1,1)
   vR=u1(ii,jj+1,2)/u1(ii,jj+1,1)
   wR=u1(ii,jj+1,4)/u1(ii,jj+1,1)
   bxR=u1(ii,jj+1,6)
   byR=u1(ii,jj+1,5)
   bzR=u1(ii,jj+1,7)
   BBR=bxR*bxR+byR*byR+bzR*bzR
   vvR=uR*uR+vR*vR+wR*wR
   pR=(gam-1.)*(u1(ii,jj+1,8)-0.5_8*rhoR*vvR-0.5*BBR)
   ptotR=pR+0.5_8*BBR
 
   !----------------------------------------------
   !Roe-Averaged State Variables (arithmetic mean):
   !----------------------------------------------
   rhot=0.5_8*(rhoL+rhoR)
   ut=0.5_8*(uL+uR)
   byt=0.5_8*(byL+byR)
   bzt=0.5_8*(bzL+bzR)
   ptot=0.5_8*(ptotL+ptotR) !total pressure: fluid+magnetic field
   pt=ptot-0.5_8*(bxL*bxL+byt*byt+bzt*bzt)
   
   if(pt>0._8 .and. rhot>0._8)then

   at=sqrt(gam*pt/rhot)

   
   !-------------------------------------------------
   !Fast Mode wave speedc_f
   !-------------------------------------------------
   bx=sqrt(bxL*bxL/rhot)
   by=sqrt(byt*byt/rhot)
   bz=sqrt(bzt*bzt/rhot)
   bsqr=(bxL*bxL+byt*byt+bzt*bzt)/rhot

   cf=0.5_8*(at*at+bsqr+sqrt((at*at+bsqr)*(at*at+bsqr)-4._8*at*at*bx*bx))
   cf=sqrt(cf)

   !find maximum speed
   smax=max(smax,max(abs(ut),abs(vt),abs(wt))+cf)
   end if

  end do
end do 
 
!--------------------------------------------
!Set time-step size
!--------------------------------------------
dt=dx*cour/smax 

print*,'Done computing time-step.'

end subroutine timeStep


subroutine horcut()
!Local Variables
integer,parameter::ycell=ny/2
integer::jj

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=u2(jj,ycell,1)
  velx=u1(jj,ycell,2)/u1(jj,ycell,1)
  vely=u1(jj,ycell,3)/u1(jj,ycell,1)
  velz=u1(jj,ycell,4)/u1(jj,ycell,1)
  magx=u1(jj,ycell,5)
  magy=u1(jj,ycell,6)
  magz=u1(jj,ycell,7)
  pres=(gam-1.)*(u1(jj,ycell,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(12,*) x,dens,velx,vely,velz,pres
  write(13,*) x,magx,magy,magz
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx/2
integer::jj

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(xcell,jj,1)
  velx=u1(xcell,jj,2)/u1(xcell,jj,1)
  vely=u1(xcell,jj,3)/u1(xcell,jj,1)
  velz=u1(xcell,jj,4)/u1(xcell,jj,1)
  magx=u1(xcell,jj,5)
  magy=u1(xcell,jj,6)
  magz=u1(xcell,jj,7)
  pres=(gam-1.)*(u1(xcell,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(14,*) y,dens,velx,vely,velz,pres
  write(15,*) y,magx,magy,magz
end do

end subroutine vercut

subroutine diagcut()
integer::jj
real*8::vperp,vpar,bpar,bperp

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(jj,jj,1)
  velx=u1(jj,jj,2)/u1(jj,jj,1)
  vely=u1(jj,jj,3)/u1(jj,jj,1)
  velz=u1(jj,jj,4)/u1(jj,jj,1)
  vpar=(vely+velx)/sqrt(2.)
  vperp=(vely-velx)/sqrt(2.)
  magx=u1(jj,jj,5)
  magy=u1(jj,jj,6)
  magz=u1(jj,jj,7)
  bpar=(magy+magx)/sqrt(2.)
  bperp=(magy-magx)/sqrt(2.)
  pres=(gam-1.)*(u1(jj,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(16,*) y,dens,vpar,vperp,velz,pres
  write(17,*) y,bpar,bperp,magz
end do

end subroutine diagcut

subroutine fileOutput(iunit)
!Input variables
integer::iunit
!local variables
integer::jj,kk
real*8::magpres

do kk=1,ny
  y=ymin+(kk-0.5)*dy
  do jj=1,nx
   x=xmin+(jj-0.5)*dx
   dens=u1(jj,kk,1)
   velx=u1(jj,kk,2)/u1(jj,kk,1)
   vely=u1(jj,kk,3)/u1(jj,kk,1)
   velz=u1(jj,kk,4)/u1(jj,kk,1)
   magx=u1(jj,kk,5)
   magy=u1(jj,kk,6)
   magz=u1(jj,kk,7)
   magpres=0.5*(magx**2.+magy**2.+magz**2.)
   pres=(gam-1.)*(u1(jj,kk,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))    

   write(iunit,*) x,y,dens,velx,vely,magx,magy,pres,magpres
  end do
end do

end subroutine fileOutput


subroutine bound()
!Local Variables
integer::jj,kk,ll

if (boundaryType==1)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
end do

else if(boundaryType==2)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do

else if(boundaryType==3)then

do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
  u2(jj,0,3)=-u2(jj,1,3)	
  u2(jj,-1,3)=-u2(jj,1,3)
  u2(jj,ny+1,3)=-u2(jj,ny,3)
  u2(jj,ny+2,3)=-u2(jj,ny,3)
  u2(jj,0,6)=-u2(jj,1,6)	
  u2(jj,-1,6)=-u2(jj,1,6)
  u2(jj,ny+1,6)=-u2(jj,ny,6)
  u2(jj,ny+2,6)=-u2(jj,ny,6)

end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
  u2(0,kk,2)=-u2(1,kk,2)	
  u2(-1,kk,2)=-u2(1,kk,2)
  u2(nx+1,kk,2)=-u2(nx,kk,2)
  u2(nx+2,kk,2)=-u2(nx,kk,2)
  u2(0,kk,5)=-u2(1,kk,5)	
  u2(-1,kk,5)=-u2(1,kk,5)
  u2(nx+1,kk,5)=-u2(nx,kk,5)
  u2(nx+2,kk,5)=-u2(nx,kk,5)
end do

else if(boundaryType==4)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,ny,ll)	
    u2(jj,-1,ll)=u2(jj,ny-1,ll)
    u2(jj,ny+1,ll)=u2(jj,1,ll)
    u2(jj,ny+2,ll)=u2(jj,2,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do

else if(boundaryType==5)then
do jj=-1,nx+2
  do ll=1,8
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
  u2(jj,0,3)=-u2(jj,1,3)	
  u2(jj,-1,3)=-u2(jj,1,3)
  u2(jj,ny+1,3)=-u2(jj,ny,3)
  u2(jj,ny+2,3)=-u2(jj,ny,3)
  u2(jj,0,6)=-u2(jj,1,6)	
  u2(jj,-1,6)=-u2(jj,1,6)
  u2(jj,ny+1,6)=-u2(jj,ny,6)
  u2(jj,ny+2,6)=-u2(jj,ny,6)

end do

do kk=-1,ny+2
  do ll=1,8
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do


end if

end subroutine bound


subroutine protection()
!Local Variables
integer::jj,kk,ll
real*8::pres0

do jj=1,nx
  do kk=1,ny
   
    dens=u2(jj,kk,1)
    if(dens<0._8)then
      !if(debug==1)then
        print*,'Negative density detected in cell:',jj,kk,'. Using protection.'
      !end if
      dens=max(u1(jj,kk,1),min_dens)
      u2(jj,kk,1)=dens
    end if

    pres=(gam-1.)*( u2(jj,kk,8)- &
         0.5*(u2(jj,kk,2)**2.+u2(jj,kk,3)**2.+u2(jj,kk,4)**2.)/u2(jj,kk,1) &
        -0.5*(u2(jj,kk,5)**2.+u2(jj,kk,6)**2.+u2(jj,kk,7)**2.) )    
    pres0=(gam-1.)*( u1(jj,kk,8)- &
         0.5*(u1(jj,kk,2)**2.+u1(jj,kk,3)**2.+u1(jj,kk,4)**2.)/u1(jj,kk,1) &
        -0.5*(u1(jj,kk,5)**2.+u1(jj,kk,6)**2.+u1(jj,kk,7)**2.) ) 

    if(pres<0._8)then
    !if(debug==1)then
      print*,'Negative pressure detected in cell:',jj,kk,'. Using protection.'
      print*,'pres=',pres
     ! end if
      pres=max(min_pres,pres0)
      print*,'new value=',pres
      u2(jj,kk,8)=gamul*pres &
        +0.5*(u1(jj,kk,2)**2.+u1(jj,kk,3)**2.+u1(jj,kk,4)**2.)/u1(jj,kk,1) &
        +0.5*(u1(jj,kk,5)**2.+u1(jj,kk,6)**2.+u1(jj,kk,7)**2.) 
    end if
   end do
end do

end subroutine protection

function sgn(x) result(fx)
  real*8,intent(in)::x
  real*8::fx
  fx=sign(1._8,x)
end


end module RoeSolver2D_mod
