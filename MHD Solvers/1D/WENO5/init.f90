module init_mod
use globalVariables_mod
implicit none


contains

subroutine initialize(q)

!Subroutine I/O variables
real::q(-3:nx+4,8)

!Subroutine local variables
real::rhoL,uL,vL,wL,bxL,byL,bzL,pL,aL !left state
real::rhoR,uR,vR,wR,bxR,byR,bzR,pR,aR !right state
integer::i,j,k


xmin=0.0
xmax=1.0
!ymin=0.0
!ymax=1.0

dx=(xmax-xmin)/nx
!dy=(ymax-ymin)/ny
dt=0.

gam=5./3.
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)
gam1=(gam-1.)/2.
gam2=(gam-2.)/(gam-1.)


if(initOption==0)then

bxL=tol
bxR=bxL
byL=tol
byR=byL
bzL=tol
bzR=bzL

rhoL=1.0
uL=0.
vL=0.
wL=0.

pL=20.

rhoR=0.125
uR=0.
vR=0.
wR=0.
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

rhoL=1.
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


do i=-3,nx+4
  if(i<nx/2) then
    q(i,1)=rhoL
    q(i,2)=rhoL*uL
    q(i,3)=rhoL*vL
    q(i,4)=rhoL*wL
    q(i,5)=bxL
    q(i,6)=byL
    q(i,7)=bzL  
    q(i,8)=0.5*rhoL*(uL*uL+vL*vL+wL*wL) &
           +gamul*pL+0.5*(bxL*bxL+byL*byL+bzL*bzL)

    !bxInt(jj,kk)=bxL
    !byInt(jj,kk)=byL
  else
    q(i,1)=rhoR
    q(i,2)=rhoR*uR
    q(i,3)=rhoR*vR
    q(i,4)=rhoR*wR
    q(i,5)=bxR
    q(i,6)=byR
    q(i,7)=bzR  
    q(i,8)=0.5*rhoR*(uR*uR+vR*vR+wR*wR) &
           +gamul*pR+0.5*(bxR*bxR+byR*byR+bzR*bzR) 
    !bxInt(jj,kk)=bxR 
    !byInt(jj,kk)=byR 
  end if
end do


print*,'Shock-Tube Test'
print*,'Left State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoL,uL,vL,wL,bxL,byL,bzL,pL
print*,'Right State(rho,vx,vy,vz,Bx,By,Bz,p)=',rhoR,uR,vR,wR,bxR,byR,bzR,pR


end subroutine initialize 



end module init_mod
