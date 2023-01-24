module init_mod
use vars_mod

implicit none


contains

subroutine init(u1)
!------------------------------------------------------
!initialization of conserved variables
!------------------------------------------------------

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
real*8::rhoL,uL,vL,wL,pL,bxL,byL,bzL !left state
real*8::rhoR,uR,vR,wR,pR,bxR,byR,bzR !right state 
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

gam=5./3.

gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

if(inputOption==1)then

if(initOption==0)then

bxL=0.
bxR=bxL

rhoL=1.0
uL=0.
vL=0.
wL=0.
byL=0.
bzL=0.
pL=1.

rhoR=0.125
uR=0.
vR=0.
wR=0.
byR=0.
bzR=0.
pR=0.1
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



if(initOption==1)then

bxL=5./sqrt(4.*pi)
bxR=bxL

rhoL=1.0
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


else if(inputOption==2)then
  read(unit=10,fmt=*) rhoL,uL,vL,wL,byL,bzL,pL 
  read(unit=10,fmt=*) rhoR,uR,vR,wR,byR,bzR,pR
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
uL=-0.5
vL=0.0
bxL=0./sqrt(4.*pi)
byL=0./sqrt(4.*pi)
bzL=0.
wL=0.
pL=1.0

rhoR=1.5
uR=0.5
vR=0.0
wR=0.
bxR=0./sqrt(4.*pi)
byR=0./sqrt(4.*pi)
bzR=0.
pR=1.0

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
    vL=1.0*(ymax-ymin)*sin(kz*x)*exp(10*(kk-ny/2)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.5*(ymax-ymin)*rand(0)*exp(-300.*(kk-ny/2)**2.)
  end if

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
    u1(jj,kk,3)=rhoR*vL
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


print*,'Riemann Problem'
print*,'Left State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoL,uL,vL,wL,bxL,byL,bzL,pL
print*,'Right State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoR,uR,vR,wR,bxR,byR,bzR,pR

print*,'nx,ny=',nx,ny

if (outOption==1)then
call vercut(u1)
call horcut(u1)
call diagcut(u1)
call fileOutput(0+1000,u1)
end if

print*,'Initialization done...'


end subroutine init

subroutine horcut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
integer,parameter::ycell=38
integer::jj
real*8::dens,velx,vely,velz,pres,magx,magy,magz

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=u1(jj,ycell,1)
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

subroutine vercut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
integer,parameter::xcell=nx/2
integer::jj
real*8::dens,velx,vely,velz,pres,magx,magy,magz

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

subroutine diagcut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

integer::jj
real*8::vperp,vpar,bpar,bperp
real*8::dens,velx,vely,velz,pres,magx,magy,magz

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

subroutine fileOutput(iunit,u1)

!Input variables
integer::iunit
real*8::u1(-1:nx+2,-1:ny+2,8)
!local variables
integer::jj,kk
real*8::dens,velx,vely,velz,pres,magx,magy,magz

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
   pres=(gam-1.)*(u1(jj,kk,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))    

   write(iunit,*) x,y,dens,velx,vely,magx,magy,pres
  end do
end do

end subroutine fileOutput

end module init_mod
