module init_mod
use vars_mod
implicit none


contains

subroutine init(u1)
!------------------------------------------------------
!initialization of conserved variables
!------------------------------------------------------

!Local Variables
integer::j
real*8::rhoL,uL,vL,wL,pL,bxLR,byL,bzL,aL !left state
real*8::rhoR,uR,vR,wR,pR,byR,bzR,aR !right state 

!Input Variables
real*8::u1(-1:nx+2,8)

open(unit=10,file='input.txt')
open(unit=11,file='output_fluid.txt')
open(unit=12,file='output_mag.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
cour=0.8
dx=(xmax-xmin)/nx

gam=5./3.

gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)


if(inputOption==1)then

if(initOption==1)then
bxLR=5./sqrt(4.*pi)

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
bxLR=3./sqrt(4.*pi)

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
bxLR=2./sqrt(4.*pi)

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
bxLR=3./sqrt(4.*pi)

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
bxLR=0.

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
bxLR=0.

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
bxLR=1.

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
bxLR=1.3

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
bxLR=0.75

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
bxLR=0.7

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
bxLR=0.75

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
bxLR=1.3

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

do j=-1,nx+2
  x=xmin+(j-0.5)*dx 
  if(j<nx/2) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=rhoL*vL
    u1(j,4)=rhoL*wL
    u1(j,5)=byL
    u1(j,6)=bzL  
    u1(j,7)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxLR*bxLR+byL*byL+bzL*bzL)/(4.*pi)
    u1(j,8)=bxLR    

    if(j>=0 .and. j<=nx)then
      write(11,*) x,rhoL,uL,vL,wL,pL
      write(12,*) x,byL,bzL
    end if
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=rhoR*vR
    u1(j,4)=rhoR*wR
    u1(j,5)=byR
    u1(j,6)=bzR  
    u1(j,7)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxLR*bxLR+byR*byR+bzR*bzR)/(4.*pi)
    u1(j,8)=bxLR    
    
    if(j>=0 .and. j<=nx)then
      write(11,*) x,rhoR,uR,vR,wR,pR
      write(12,*) x,byR,bzR
    end if
  end if

end do

print*,'Left State (rho,u,v,w,Bx,By,Bz,p) =',rhoL,uL,vL,wL,bxLR,byL,bzL,pL
print*,'Right State (rho,u,v,w,Bx,By,Bz,p) =',rhoR,uR,vR,wR,bxLR,byR,bzR,pR

end subroutine init


end module init_mod
