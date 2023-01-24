!Exact Riemann Solver of 1D Eulers Equiations (Improved version, nearly perfect)
!Two root solvers included: Secant Method and Newton
!Newton solver performs better

module HLLRiemannSolver_mod
implicit none

integer,parameter::discontinuityDir=2 !1:x, 2:y
integer,parameter::debug=0 !0:off 1:on
integer,parameter::rootSolverOption=2 !1:Secant 2:Newton
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

real*8::rhoL,uL,vL,pL,aL,FL(4),GL(4) !left state
real*8::rhoR,uR,vR,pR,aR,FR(4),GR(4) !right state
real*8::pS !star region state
real*8::rhoB,aB,p0,pMax,pMin
real*8,parameter::Q_user=2.0 
real*8::sL,sR,smax,smax1
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
open(unit=13,file='Output/vercut.txt')


!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
cour=0.9

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

if(discontinuityDir==1)then
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
else if(discontinuityDir==2)then
do jj=-1,nx+2
 do kk=-1,ny+2
  if(kk<ny/2) then
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

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx*0.55
integer::jj

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=max(u2(xcell,jj,1),min_dens)
   velx=u2(xcell,jj,2)/u2(xcell,jj,1)
   vely=u2(xcell,jj,3)/u2(xcell,jj,1)
   pres=max((gam-1.)*(u2(xcell,jj,4)-0.5*(u2(xcell,jj,2)*u2(xcell,jj,2)/u2(xcell,jj,1)&
        +u2(xcell,jj,3)*u2(xcell,jj,3)/u2(xcell,jj,1))),min_pres)
   write(13,*) y,dens,velx,vely,pres
end do

end subroutine vercut


subroutine computeFlux(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk

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
 
FL(1)=rhoL*uL 
FL(2)=rhoL*uL*uL+pL
FL(3)=rhoL*uL*vL
FL(4)=0.5*rhoL*(uL*uL+vL*vL)*uL+gam*gamul*pL*uL

FR(1)=rhoR*uR  
FR(2)=rhoR*uR*uR+pR
FR(3)=rhoR*uR*vR
FR(4)=0.5*rhoR*(uR*uR+vR*vR)*uR+gam*gamul*pR*uR
!y-sweeps
else if(dir==2)then

rhoL=max(u1(ind,jj,1),min_dens)
uL=u1(ind,jj,3)/u1(ind,jj,1)
vL=u1(ind,jj,2)/u1(ind,jj,1)
pL=max((gam-1.)*(u1(ind,jj,4)-0.5*(u1(ind,jj,2)*u1(ind,jj,2)/u1(ind,jj,1)&
        +u1(ind,jj,3)*u1(ind,jj,3)/u1(ind,jj,1))),min_pres)

rhoR=max(u1(ind,jj+1,1),min_dens)
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
pR=max((gam-1.)*(u1(ind,jj+1,4)-0.5*(u1(ind,jj+1,2)*u1(ind,jj+1,2)/u1(ind,jj+1,1)&
        +u1(ind,jj+1,3)*u1(ind,jj+1,3)/u1(ind,jj+1,1))),min_pres)

GL(1)=rhoL*uL
GL(2)=rhoL*uL*vL
GL(3)=rhoL*uL*uL+pL
GL(4)=0.5*rhoL*(uL*uL+vL*vL)*uL+gam*gamul*pL*uL

GR(1)=rhoR*uR
GR(2)=rhoR*uR*vR   
GR(3)=rhoR*uR*uR+pR
GR(4)=0.5*rhoR*(uR*uR+vR*vR)*uR+gam*gamul*pR*uR
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

call computeWaveSpeeds()

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax1=max(smax1,smax)

if(dir==1)then
!---------------------------------------
!F_j+1/2
!---------------------------------------
 if (sL>=0._8)then
   do kk=1,4
     xflux(jj,kk)=FL(kk)
   end do
 else if(sL<0._8 .and. sR>0._8) then
   do kk=1,4
     xflux(jj,kk)=sR*FL(kk)-sL*FR(kk)+sL*sR*(u1(jj+1,ind,kk)-u1(jj,ind,kk))
     xflux(jj,kk)=xflux(jj,kk)/(sR-sL)
   end do
 else if(sR<=0._8)then
   do kk=1,4
     xflux(jj,kk)=FR(kk)
   end do 
 end  if

else if(dir==2)then
!---------------------------------------
!G_j+1/2
!---------------------------------------
 if (sL>=0._8)then
   do kk=1,4
     yflux(jj,kk)=GL(kk)
   end do
 else if(sL<0._8 .and. sR>0._8) then
   do kk=1,4
     yflux(jj,kk)=sR*GL(kk)-sL*GR(kk)+sL*sR*(u1(ind,jj+1,kk)-u1(ind,jj,kk))
     yflux(jj,kk)=yflux(jj,kk)/(sR-sL)
   end do
 else if(sR<=0._8)then
   do kk=1,4
     yflux(jj,kk)=GR(kk)
   end do 
 end  if

end if

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
pS=gLx(p0)*pL+gRx(p0)*pR-(uR-uL)
pS=pS/(gLx(p0)+gRx(p0))
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

function gLx(px) result(pxx)
!Input Variables
real*8,intent(in)::px
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(px+gamel*pL))
end function gLx

function gRx(px) result(pxx)
!Input Variables
real*8,intent(in)::px
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(px+gamel*pR))
end function gRx

!-----------------------------------------------------------
end module HLLRiemannSolver_mod
