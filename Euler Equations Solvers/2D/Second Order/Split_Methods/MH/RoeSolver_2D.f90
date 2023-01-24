!MUSCL-Hancock Method for Euler's Equation
!Riemann Solver: Roe

module RoeSolver2D_mod
implicit none

integer,parameter::discontinuityDir=2 !1:x, 2:y
integer,parameter::debug=0 !0:off 1:on
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y
real*8,parameter::Q_user=3.0 
integer,parameter::nt=200
integer,parameter::nx=40
integer,parameter::ny=40
real*8,parameter::s_tol=1.d-25
real*8,parameter::min_pres=1.d-30
real*8,parameter::min_dens=1.d-30
integer,parameter::fluidInit=6!1:V Shock, 2:Oblique Riemann Problem, 3:Riemann Problem File input 4:Kelvin-Helmholtz Instability, 5: 3 shear layers,  6:Square Sod 
integer,parameter::perturbationType=2 !1:Sinusoidal 2:Random

real*8 ::dt,dtx,dty,dx,dy !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,aL,HL,FL(5) !left state
real*8::rhoR,uR,vR,wR,pR,aR,HR,FR(5) !right state
real*8::rhot,ut,vt,wt,Ht,at,lambda(-1:nx+2,5),alpha(-1:nx+2,5)
real*8::K1(-1:nx+2,5),K2(-1:nx+2,5),K3(-1:nx+2,5),K4(-1:nx+2,5),K5(-1:nx+2,5)
	
real*8::delrho,delp,delu,delv,delw,smax
real*8::lambda1L(-1:nx+2),lambda1R(-1:nx+2),lambda5L(-1:nx+2),lambda5R(-1:nx+2)
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::dens,velx,vely,velz,pres
real*8::u1(-1:nx+2,-1:ny+2,5),u2(-1:nx+2,-1:ny+2,5),flux(-1:nx+2,5)
real*8::fleft(-1:nx+2,5),fright(-1:nx+2,5)

character(len=20) :: filename

contains 

!------------------------------------------------------
!initialization
!------------------------------------------------------
subroutine init()

!Local Variables
integer::jj,kk,mode
!-----------------------------------------------------
real*8,parameter::pi=3.14159265359
real*8::x1,y1
real*8::velx0,velx1,vely0,vely1,velz0,velz1,rho0,rho1
real*8::velx2,vely2,velx3,vely3
real*8::pre0,pre1,M0,temp1,temp2,temp3,temp4
real*8::theta !shock angle w.r.t. x axis
real*8::alpha
real*8::Sh !shock speed 
real*8::zeta,kz
real*8::vL1,vL2,vL3
!-----------------------------------------------------
open(unit=10,file='input.txt')
open(unit=0+1000,file='Output/t=0.txt')
open(unit=12,file='Output/horcut.txt')
open(unit=13,file='Output/vercut.txt')

print*,'Second Order MUSCL-Hancock Method, using Roe Solver...'

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
ymin=0.0
ymax=1.0
cour=0.9

dx=(xmax-xmin)/nx
dy=(ymax-ymin)/ny
dt=0.
dtx=0.
dty=0.

gam=1.4
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

smax=0.

read(unit=10,fmt=*) rhoL,uL,vL,wL,pL 
read(unit=10,fmt=*) rhoR,uR,vR,wR,pR

print*,'Riemann Problem'
print*,'Left State=',rhoL,uL,vL,wL,pL
print*,'Right State=',rhoR,uR,vR,wR,pR

if(fluidInit==1)then

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
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx2*u1(jj,kk,1)
    u1(jj,kk,3) = vely2*u1(jj,kk,1) 
    u1(jj,kk,4) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3))/u1(jj,kk,1)
    !region 1
  elseif (y1 .lt. temp4 .and. y1 .gt. temp1) then
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

else if(fluidInit==2)then
!Oblique Discontinuity 	 

do kk=-1,ny+2
 y1=(kk-0.5)*dy
 do jj=-1,nx+2
  x1=(jj-0.5)*dx
				
  if (y1 .le. x1) then 
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR**2.)+pR*gamul 
  else
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL  
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL**2.)+pL*gamul
  end if			
 end do
end do


else if(fluidInit==3)then

if(discontinuityDir==1)then
do jj=-1,nx+2
 do kk=-1,ny+2
  if(jj<nx/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL  
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL**2.)+pL*gamul
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR**2.)+pR*gamul 
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
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL**2.)+pL*gamul
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR**2.)+pR*gamul 
  end if
 end do
end do

end if

else if(fluidInit==4)then
!Kelvin-Helmholtz Instability 
rhoL=1.0
uL=-0.6
vL=0.0
wL=0.
pL=1.0

rhoR=3.0
uR=0.6
vR=0.0
wR=0.
pR=1.0

mode=4

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2

 x1 =xmin+(jj - 0.5)*dx
 
 do kk=-1,ny+2
  y1 =ymin+(kk - 0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=1.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-ny/2)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.75*(ymax-ymin)*rand(0)*exp(-500.*(kk-ny/2)**2.)
  end if
  


  if(kk<ny/2) then !Below Interface
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  else !Above Interface
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vL
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul 
  end if
 end do   
end do

else if(fluidInit==5)then
!3 Shear Layers 
rhoL=1.0
uL=-0.6
vL=0.0
wL=0.
pL=1.0

rhoR=3.0
uR=0.6
vR=0.0
wR=0.
pR=1.0

mode=4

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2

 x1 =xmin+(jj - 0.5)*dx
 
 do kk=-1,ny+2
  y1 =ymin+(kk - 0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
  !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=1.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-ny/3.)**2.)
    vL1=1.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-2.*ny/3.)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.75*(ymax-ymin)*rand(0)*exp(-500.*(kk-ny/3.)**2.)
    vL1=0.75*(ymax-ymin)*rand(0)*exp(-500.*(kk-2.*ny/3.)**2.)
  end if
  

  if(kk<ny/3.) then 
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  else if (kk<2.*ny/3. .and. kk>ny/3. )then 
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vL1
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul 
  else
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  end if
 end do   
end do

else if(fluidInit==6)then

rhoL=1.0
uL=0.
vL=0.
wL=0.
pL=1.0

rhoR=0.125
uR=0.
vR=0.
wR=0.
pR=0.1

do kk=-1,ny+2
 do jj=-1,nx+2				
  if(kk==ny/2 .and. jj==nx/2) then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL  
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL**2.)+pL*gamul 
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR**2.)+pR*gamul
  end if			
 end do
end do


end if

u2=u1

!print*,'Left State (rho,u,v,p) =',rhoL,uL,vL,pL
!print*,'Right State (rho,u,v,p) =',rhoR,uR,vR,pR

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
   velz=u2(jj,kk,4)/u2(jj,kk,1)
   pres=(gam-1.)*(u2(jj,kk,5)-0.5*dens*(velx**2.+vely**2.+velz**2.))
   write(iunit,*) x,y,dens,velx,vely,pres
  end do
end do

end subroutine fileOutput

subroutine horcut()
!Local Variables
integer,parameter::ycell=ny*0.4
integer::jj

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=max(u2(jj,ycell,1),min_dens)
   velx=u2(jj,ycell,2)/u2(jj,ycell,1)
   vely=u2(jj,ycell,3)/u2(jj,ycell,1)
   velz=u2(jj,ycell,4)/u2(jj,ycell,1)
   pres=(gam-1.)*(u2(jj,ycell,5)-0.5*dens*(velx**2.+vely**2.+velz**2.))   
   write(12,*) x,dens,velx,vely,pres,velz
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx*0.4
integer::jj

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=max(u2(xcell,jj,1),min_dens)
   velx=u2(xcell,jj,2)/u2(xcell,jj,1)
   vely=u2(xcell,jj,3)/u2(xcell,jj,1)
   velz=u2(xcell,jj,4)/u2(xcell,jj,1)
   pres=(gam-1.)*(u2(xcell,jj,5)-0.5*dens*(velx**2.+vely**2.+velz**2.))   
   write(13,*) y,dens,velx,vely,pres,velz
end do

end subroutine vercut


subroutine computeFlux(n,ind,dir)

!Input Variables
integer::n,ind,dir,jj,kk
!Local Variables
real*8::qL(-1:nx+2,5),qR(-1:nx+2,5),slope(-1:nx+2,5)
real*8::ubL(-1:nx+2,5),ubR(-1:nx+2,5)

!--------------------------------------------------
!Limited Slopes: delta_j
!--------------------------------------------------
if(dir==1)then

do jj=0,n+1
  do kk=1,5
    slope(jj,kk)=limiter(u1(jj,ind,kk)-u1(jj-1,ind,kk),u1(jj+1,ind,kk)-u1(jj,ind,kk))
  end do
end do

else if(dir==2)then
do jj=0,n+1
  do kk=1,5
    slope(jj,kk)=limiter(u1(ind,jj,kk)-u1(ind,jj-1,kk),u1(ind,jj+1,kk)-u1(ind,jj,kk))
  end do
end do

end if

!----------------------------------------------
!Boundary Extrapolated Values: uL_i, uR_i
!----------------------------------------------
do jj=0,n+1
  do kk=1,5
   if(dir==1)then
     qL(jj,kk)=u1(jj,ind,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(jj,ind,kk)+0.5*slope(jj,kk)
   else if(dir==2)then
     qL(jj,kk)=u1(ind,jj,kk)-0.5*slope(jj,kk)
     qR(jj,kk)=u1(ind,jj,kk)+0.5*slope(jj,kk)
   end if
  end do
end do

!-------------------------------------------------------------------
!Evolve Boundary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------
if(dir==1)then
do jj=0,n+1 
  dens=max(qL(jj,1),min_dens)
  velx=qL(jj,2)/qL(jj,1)
  vely=qL(jj,3)/qL(jj,1)
  velz=qL(jj,4)/qL(jj,1)
  pres=max((gam-1.)*(qL(jj,5)-0.5*dens*(velx**2.+vely**2.+velz**2.)),min_pres)


  FL(1)=dens*velx 
  FL(2)=dens*velx*velx+pres
  FL(3)=dens*velx*vely
  FL(4)=dens*velx*velz
  FL(5)=0.5*dens*(velx*velx+vely*vely+velz*velz)*velx+gam*gamul*pres*velx

  dens=max(qR(jj,1),min_dens)
  velx=qR(jj,2)/qR(jj,1)
  vely=qR(jj,3)/qR(jj,1)
  velz=qR(jj,4)/qR(jj,1)
  pres=max((gam-1.)*(qR(jj,5)-0.5*dens*(velx**2.+vely**2.+velz**2.)),min_pres)

  FR(1)=dens*velx 
  FR(2)=dens*velx*velx+pres
  FR(3)=dens*velx*vely
  FR(4)=dens*velx*velz
  FR(5)=0.5*dens*(velx*velx+vely*vely+velz*velz)*velx+gam*gamul*pres*velx

  do kk=1,5
    ubL(jj,kk)=qL(jj,kk)+0.5*(dtx/dx)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dtx/dx)*(FL(kk)-FR(kk))
  end do

end do

else if(dir==2)then
do jj=0,n+1 
  dens=max(qL(jj,1),min_dens)
  velx=qL(jj,2)/qL(jj,1)
  vely=qL(jj,3)/qL(jj,1)
  velz=qL(jj,4)/qL(jj,1)
  pres=max((gam-1.)*(qL(jj,5)-0.5*dens*(velx**2.+vely**2.+velz**2.)),min_pres)


  FL(1)=dens*vely 
  FL(2)=dens*velx*vely
  FL(3)=dens*vely*vely+pres
  FL(4)=dens*velz*vely
  FL(5)=0.5*dens*(velx*velx+vely*vely+velz*velz)*vely+gam*gamul*pres*vely


  dens=max(qR(jj,1),min_dens)
  velx=qR(jj,2)/qR(jj,1)
  vely=qR(jj,3)/qR(jj,1)
  velz=qR(jj,4)/qR(jj,1)
  pres=max((gam-1.)*(qR(jj,5)-0.5*dens*(velx**2.+vely**2.+velz**2.)),min_pres)


  FR(1)=dens*vely 
  FR(2)=dens*velx*vely
  FR(3)=dens*vely*vely+pres
  FR(4)=dens*velz*vely
  FR(5)=0.5*dens*(velx*velx+vely*vely+velz*velz)*vely+gam*gamul*pres*vely


  do kk=1,5
    ubL(jj,kk)=qL(jj,kk)+0.5*(dtx/dx)*(FL(kk)-FR(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dtx/dx)*(FL(kk)-FR(kk))
  end do 

  
end do

end if

!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do jj=0,n
!set left and right states for local Riemann problem at cell interface

if(dir==1)then
rhoL=max(min_dens,ubR(jj,1))
uL=ubR(jj,2)/ubR(jj,1)
vL=ubR(jj,3)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
pL=max((gam-1.)*(ubR(jj,5)-0.5*rhoL*(uL*uL+vL*vL+wL*wL)),min_pres)

rhoR=max(min_dens,ubL(jj+1,1))
uR=ubL(jj+1,2)/ubL(jj+1,1)
vR=ubL(jj+1,3)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
pR=max((gam-1.)*(ubL(jj+1,5)-0.5*rhoR*(uR*uR+vR*vR+wR*wR)),min_pres)

else if(dir==2)then
rhoL=max(min_dens,ubR(jj,1))
uL=ubR(jj,3)/ubR(jj,1)
vL=ubR(jj,2)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
pL=max((gam-1.)*(ubR(jj,5)-0.5*rhoL*(uL*uL+vL*vL+wL*wL)),min_pres)

rhoR=max(min_dens,ubL(jj+1,1))
uR=ubL(jj+1,3)/ubL(jj+1,1)
vR=ubL(jj+1,2)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
pR=max((gam-1.)*(ubL(jj+1,5)-0.5*rhoR*(uR*uR+vR*vR+wR*wR)),min_pres)

end if

aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL*uL+vL*vL+wL*wL)+gamul*aL*aL

aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR*uR+vR*vR+wR*wR)+gamul*aR*aR


fleft(jj,1)=rhoL*uL 
fleft(jj,2)=rhoL*uL*uL+pL
fleft(jj,3)=rhoL*uL*vL
fleft(jj,4)=rhoL*uL*wL
fleft(jj,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)*uL+gam*gamul*pL*uL

fright(jj,1)=rhoR*uR  
fright(jj,2)=rhoR*uR*uR+pR
fright(jj,3)=rhoR*uR*vR
fright(jj,4)=rhoR*uR*wR
fright(jj,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)*uR+gam*gamul*pR*uR

call RoeSolver(jj)

end do


if(debug==1)then
  print*,'Max wave speed=',smax
end if


!---------------------------------------
!F_j+1/2
!---------------------------------------
do jj=0,n
!check for transonic rarefaction waves
if(lambda1L(jj)<0. .and. lambda1R(jj)>0.)then
  do kk=1,5
    flux(jj,kk)=fleft(jj,kk)+(lambda1L(jj)*(lambda1R(jj)&
       -lambda(jj,1))/(lambda1R(jj)-lambda1L(jj)))*alpha(jj,1)*K1(jj,kk) 
  end do
else if(lambda5L(jj)<0. .and. lambda5R(jj)>0.)then
  do kk=1,5
    flux(jj,kk)=fright(jj,kk)-(lambda5R(jj)*(lambda(jj,5)&
       -lambda5L(jj))/(lambda5R(jj)-lambda5L(jj)))*alpha(jj,5)*K5(jj,kk)
  end do
else
  do kk=1,5
    flux(jj,kk)=0.5*(fleft(jj,kk)+fright(jj,kk))&
               -0.5*( alpha(jj,1)*abs(lambda(jj,1))*K1(jj,kk)&
               +alpha(jj,2)*abs(lambda(jj,2))*K2(jj,kk)&
     	       +alpha(jj,3)*abs(lambda(jj,3))*K3(jj,kk)&
     	       +alpha(jj,4)*abs(lambda(jj,4))*K4(jj,kk)&
     	       +alpha(jj,5)*abs(lambda(jj,5))*K5(jj,kk) )   
  end do
end if

end do

end subroutine computeFlux

subroutine RoeSolver(jj)

!Input Variables
integer::jj

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delw=wR-wL
delp=pR-pL
 


if(debug==1)then
  print*,'Cell:',jj
  print*,'Local RP Left State (rho,u,v,w,p) =',rhoL,uL,vL,wL,pL
  print*,'Local RP Right State (rho,u,v,w,p) =',rhoR,uR,vR,wR,pR
end if


!-------------------------------------------
!Wave Speeds for Harten-Hyman Entropy Fix
!-------------------------------------------
!Compute star region state
call starRegion()

lambda1L(jj)=uL-aL
lambda1R(jj)=uS-aSL

lambda5L(jj)=uS+aSR
lambda5R(jj)=uR+aR

!----------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^v,^w
!----------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
vt=(sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
wt=(sqrt(rhoL)*wL+sqrt(rhoR)*wR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*(ut*ut+vt*vt+wt*wt)))

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(jj,1)=ut-at
lambda(jj,2)=ut
lambda(jj,3)=ut
lambda(jj,4)=ut
lambda(jj,5)=ut+at

if(debug==1)then
  print*,'Characetristic Speeds(1,3,5)=',lambda(jj,1),lambda(jj,3),lambda(jj,5)
end if

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(jj,1)=1._8
K1(jj,2)=ut-at
K1(jj,3)=vt
K1(jj,4)=wt
K1(jj,5)=Ht-ut*at

K2(jj,1)=1._8
K2(jj,2)=ut
K2(jj,3)=vt
K2(jj,4)=wt
K2(jj,5)=0.5*(ut**2.+vt**2.+wt**2.)

K3(jj,1)=0._8
K3(jj,2)=0._8
K3(jj,3)=1._8
K3(jj,4)=0._8
K3(jj,5)=vt

K4(jj,1)=0._8
K4(jj,2)=0._8
K4(jj,3)=0._8
K4(jj,4)=1._8
K4(jj,5)=wt

K5(jj,1)=1._8
K5(jj,2)=ut+at
K5(jj,3)=vt
K5(jj,4)=wt
K5(jj,5)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(jj,1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(jj,2)=delrho-delp/(at*at)
alpha(jj,3)=rhot*delv
alpha(jj,4)=rhot*delw
alpha(jj,5)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(jj,1)),abs(lambda(jj,3)),abs(lambda(jj,5)))


end subroutine RoeSolver


subroutine starRegion()

!Input Variables

!Local Variables
real*8::rhoB,aB,pMin,pMax,p0


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
pS=gLx(p0)*pL+gRx(p0)*pR-(uR-uL)
pS=pS/(gLx(p0)+gRx(p0))

!u_star
uS=0.5*(uL+uR)+0.5*((pS-pR)*gRx(p0)-(pS-pL)*gLx(p0))
!rho_star
rhoSL=rhoL*((pS/pL)+gamel)/((pS/pL)*gamel+1.)
rhoSR=rhoR*((pS/pR)+gamel)/((pS/pR)*gamel+1.)
end if

aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)

if(debug==1)then
  print*,'Star Region, pS,uS=',pS,uS
end if

end subroutine starRegion

subroutine bound()
!Local Variables
integer::jj,kk,ll

if (boundaryType==1)then
do jj=-1,nx+2
  do ll=1,5
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,5
    u2(0,kk,ll)=u2(1,kk,ll)	
    u2(-1,kk,ll)=u2(1,kk,ll)
    u2(nx+1,kk,ll)=u2(nx,kk,ll)
    u2(nx+2,kk,ll)=u2(nx,kk,ll)
  end do
end do
else if(boundaryType==2)then
do jj=-1,nx+2
  do ll=1,5
    u2(jj,0,ll)=u2(jj,1,ll)	
    u2(jj,-1,ll)=u2(jj,1,ll)
    u2(jj,ny+1,ll)=u2(jj,ny,ll)
    u2(jj,ny+2,ll)=u2(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,5
    u2(0,kk,ll)=u2(nx,kk,ll)	
    u2(-1,kk,ll)=u2(nx-1,kk,ll)
    u2(nx+1,kk,ll)=u2(1,kk,ll)
    u2(nx+2,kk,ll)=u2(2,kk,ll)
  end do
end do

end if

end subroutine bound

function gLx(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(x+gamel*pL))
end function gLx

function gRx(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(x+gamel*pR))
end function gRx

function limiter(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z,b

  b=1.

  if(y>0._8)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter


end module RoeSolver2D_mod
