!MUSCL-Hancock Method for Euler's Equation
!2D-Unsplit
!Riemann Solver: Roe

module RoeSolver2D_mod
implicit none

integer,parameter::discontinuityDir=2 !1:x, 2:y
integer,parameter::debug=1 !0:off 1:on
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y
real*8,parameter::Q_user=3.0 
integer,parameter::nt=200
integer,parameter::nx=100
integer,parameter::ny=100
integer,parameter::tSkip=5
real*8::xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,cour=0.9

real*8,parameter::s_tol=1.d-25
real*8,parameter::min_pres=1.d-30
real*8,parameter::min_dens=1.d-30
integer,parameter::fluidInit=2 !1:V Shock, 2:Oblique Riemann Problem, 3:Riemann Problem File input 4:KH Instability 5: 3-Shear-Layer KH 

integer,parameter::perturbationType=1 !1:Sinusoidal 2:Random

real*8::dt,dx,dy !time step and spatial resolution
real*8::x,y
real*8::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,aL,HL,FL(5),GL(5) !left state
real*8::rhoR,uR,vR,wR,pR,aR,HR,FR(5),GR(5) !right state
real*8::rhot,ut,vt,wt,Ht,at,lambda(5),alpha(5)
real*8::K1(5),K2(5),K3(5),K4(5),K5(5)
	
real*8::delrho,delp,delu,delv,delw,smax
real*8::lambda1L,lambda1R,lambda5L,lambda5R
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::dens,velx,vely,velz,pres
real*8::u1(-1:nx+2,-1:ny+2,5),u2(-1:nx+2,-1:ny+2,5)
real*8::f(-1:nx+2,-1:ny+2,5),g(-1:nx+2,-1:ny+2,5)



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
real*8::vL1,vL2
!-----------------------------------------------------
open(unit=10,file='input.txt')
open(unit=0+1000,file='Output/t=0.txt')
open(unit=12,file='Output/horcut.txt')
open(unit=13,file='Output/vercut.txt')

print*,'Second Order MUSCL-Hancock Method, using Roe Solver...'

!Discountinuity in initial condition at the middle of [xmin,xmax]

dx=(xmax-xmin)/real(nx)
dy=(ymax-ymin)/real(ny)

print*,'dx,dy=',dx,dy

dt=0.

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
				
  if (y1 .le. ymax-x1) then 
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

rhoR=1.5
uR=0.6
vR=0.0
wR=0.
pR=1.0

mode=6

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
    vL=1.25*(ymax-ymin)*rand(0)*exp(-500.*(kk-ny/2)**2.)
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
!3-Shear Layer Kelvin-Helmholtz Instability 
rhoL=1.0
uL=-0.6
vL=0.0
wL=0.
pL=1.0

rhoR=1.5
uR=0.6
vR=0.0
wR=0.
pR=1.0

mode=3

kz=2.*pi*mode/(xmax-xmin)

do jj=-1,nx+2

 x1 =xmin+(jj - 0.5)*dx
 
 do kk=-1,ny+2
  y1 =ymin+(kk - 0.5)*dy

  !Perturbation in velocity y-component only
  !Exponential term localizes the perturbation near the interface
  
   !Sinusoidal Perturbation
  if(perturbationType==1)then
    vL=1.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-ny/3)**2.)
    vL1=1.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-2*ny/3)**2.)
  !Random Perturbation
  else if(perturbationType==2)then
    vL=0.75*(ymax-ymin)*rand(0)*exp(-500.*(kk-ny/3)**2.)
    vL1=0.75*(ymax-ymin)*rand(0)*exp(-500.*(kk-2*ny/3)**2.)
  end if
  

  if(kk<ny/3) then 
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  else if(kk>=ny/3 .and. kk<2*ny/3)then 
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vL1
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul 
  else if(kk>=2*ny/3)then
    u1(jj,kk,1)=rhoL
    u1(jj,kk,2)=rhoL*uL
    u1(jj,kk,3)=rhoL*vL
    u1(jj,kk,4)=rhoL*wL
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  end if
  end do   
end do

end if

u2=u1

!print*,'Left State (rho,u,v,p) =',rhoL,uL,vL,pL
!print*,'Right State (rho,u,v,p) =',rhoR,uR,vR,pR

call fileOutput(0+1000)
call horcut()
call vercut()

print*,'Completed fluid initialization...'

end subroutine init


subroutine computeFlux()

!Input Variables
integer::ii,jj,kk
!Local Variables
!real*8::qLx(-1:nx+2,-1:ny+2,5),qRx(-1:nx+2,-1:ny+2,5)
!real*8::qLy(-1:nx+2,-1:ny+2,5),qRy(-1:nx+2,-1:ny+2,5)
!real*8::slopex(-1:nx+2,-1:ny+2,5),slopey(-1:nx+2,-1:ny+2,5)
!real*8::ubLx(-1:nx+2,-1:ny+2,5),ubRx(-1:nx+2,-1:ny+2,5)
!real*8::ubLy(-1:nx+2,-1:ny+2,5),ubRy(-1:nx+2,-1:ny+2,5)


!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do ii=0,nx
 do jj=0,ny
 !set left and right states for local Riemann problem at x,y cell interface
 !x_interface
 rhoL=max(min_dens,u1(ii,jj,1))
 uL=u1(ii,jj,2)/u1(ii,jj,1)
 vL=u1(ii,jj,3)/u1(ii,jj,1)
 wL=u1(ii,jj,4)/u1(ii,jj,1)
 pL=max((gam-1.)*(u1(ii,jj,5)-0.5*rhoL*(uL*uL+vL*vL+wL*wL)),min_pres)
 aL=sqrt(gam*pL/rhoL)
 HL=0.5*(uL*uL+vL*vL+wL*wL)+gamul*aL*aL

 rhoR=max(min_dens,u1(ii+1,jj,1))
 uR=u1(ii+1,jj,2)/u1(ii+1,jj,1)
 vR=u1(ii+1,jj,3)/u1(ii+1,jj,1)
 wR=u1(ii+1,jj,4)/u1(ii+1,jj,1)
 pR=max((gam-1.)*(u1(ii+1,jj,5)-0.5*rhoR*(uR*uR+vR*vR+wR*wR)),min_pres)
 aR=sqrt(gam*pR/rhoR)
 HR=0.5*(uR*uR+vR*vR+wR*wR)+gamul*aR*aR

 FL(1)=rhoL*uL 
 FL(2)=rhoL*uL*uL+pL
 FL(3)=rhoL*uL*vL
 FL(4)=rhoL*uL*wL
 FL(5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)*uL+gam*gamul*pL*uL

 FR(1)=rhoR*uR  
 FR(2)=rhoR*uR*uR+pR
 FR(3)=rhoR*uR*vR
 FR(4)=rhoR*uR*wR
 FR(5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)*uR+gam*gamul*pR*uR


 if(debug==1)then
  print*,'Cell Iterface:',ii,'+1/2',jj
  print*,'Local RP Left State (rho,u,v,w,p) =',rhoL,uL,vL,wL,pL
  print*,'Local RP Right State (rho,u,v,w,p) =',rhoR,uR,vR,wR,pR
 end if


 call RoeSolver()

 !---------------------------------------
 !F_i+1/2,j
 !---------------------------------------
 !check for transonic rarefaction waves
 if(lambda1L<0. .and. lambda1R>0.)then
  do kk=1,5
    f(ii,jj,kk)=FL(kk)+(lambda1L*(lambda1R-lambda(1)) &
               /(lambda1R-lambda1L))*alpha(1)*K1(kk) 
  end do
 else if(lambda5L<0. .and. lambda5R>0.)then
  do kk=1,5
    f(ii,jj,kk)=FR(kk)-(lambda5R*(lambda(5)&
       -lambda5L)/(lambda5R-lambda5L))*alpha(5)*K5(kk)
  end do
 else
  do kk=1,5
    f(ii,jj,kk)=0.5*(FL(kk)+FR(kk))&
               -0.5*( alpha(1)*abs(lambda(1))*K1(kk)&
               +alpha(2)*abs(lambda(2))*K2(kk)&
     	       +alpha(3)*abs(lambda(3))*K3(kk)&
     	       +alpha(4)*abs(lambda(4))*K4(kk)&
     	       +alpha(5)*abs(lambda(5))*K5(kk) )   
  end do
 end if

 !y-interface
 rhoL=max(min_dens,u1(ii,jj,1))
 uL=u1(ii,jj,3)/u1(ii,jj,1)
 vL=u1(ii,jj,2)/u1(ii,jj,1)
 wL=u1(ii,jj,4)/u1(ii,jj,1)
 pL=max((gam-1.)*(u1(ii,jj,5)-0.5*rhoL*(uL*uL+vL*vL+wL*wL)),min_pres)
 aL=sqrt(gam*pL/rhoL)
 HL=0.5*(uL*uL+vL*vL+wL*wL)+gamul*aL*aL

 rhoR=max(min_dens,u1(ii,jj+1,1))
 uR=u1(ii,jj+1,3)/u1(ii,jj+1,1)
 vR=u1(ii,jj+1,2)/u1(ii,jj+1,1)
 wR=u1(ii,jj+1,4)/u1(ii,jj+1,1)
 pR=max((gam-1.)*(u1(ii,jj+1,5)-0.5*rhoR*(uR*uR+vR*vR+wR*wR)),min_pres)
 aR=sqrt(gam*pR/rhoR)
 HR=0.5*(uR*uR+vR*vR+wR*wR)+gamul*aR*aR


 GL(1)=rhoL*uL 
 GL(2)=rhoL*uL*uL+pL
 GL(3)=rhoL*uL*vL
 GL(4)=rhoL*uL*wL
 GL(5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)*uL+gam*gamul*pL*uL

 GR(1)=rhoR*uR  
 GR(2)=rhoR*uR*uR+pR
 GR(3)=rhoR*uR*vR
 GR(4)=rhoR*uR*wR
 GR(5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)*uR+gam*gamul*pR*uR

 if(debug==1)then
  print*,'Cell Iterface:',ii,jj,'+1/2'
  print*,'Local RP Left State (rho,u,v,w,p) =',rhoL,uL,vL,wL,pL
  print*,'Local RP Right State (rho,u,v,w,p) =',rhoR,uR,vR,wR,pR
 end if

 call RoeSolver()

 !---------------------------------------
 !G_i,j+1/2
 !---------------------------------------
 !check for transonic rarefaction waves
 if(lambda1L<0. .and. lambda1R>0.)then
  do kk=1,5
    g(ii,jj,kk)=GL(kk)+(lambda1L*(lambda1R-lambda(1)) &
               /(lambda1R-lambda1L))*alpha(1)*K1(kk) 
  end do
 else if(lambda5L<0. .and. lambda5R>0.)then
  do kk=1,5
    g(ii,jj,kk)=GR(kk)-(lambda5R*(lambda(5)&
       -lambda5L)/(lambda5R-lambda5L))*alpha(5)*K5(kk)
  end do
 else
  do kk=1,5
    g(ii,jj,kk)=0.5*(GL(kk)+GR(kk))&
               -0.5*( alpha(1)*abs(lambda(1))*K1(kk)&
               +alpha(2)*abs(lambda(2))*K2(kk)&
     	       +alpha(3)*abs(lambda(3))*K3(kk)&
     	       +alpha(4)*abs(lambda(4))*K4(kk)&
     	       +alpha(5)*abs(lambda(5))*K5(kk) )   
  end do
 end if


 end do
end do

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
dt=min(dx*cour/smax,dy*cour/smax) 



if(debug==1)then
  print*,'Max wave speed=',smax
end if


end subroutine computeFlux

subroutine RoeSolver()

!Input Variables

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delw=wR-wL
delp=pR-pL
 

!-------------------------------------------
!Wave Speeds for Harten-Hyman Entropy Fix
!-------------------------------------------
!Compute star region state
call starRegion()

lambda1L=uL-aL
lambda1R=uS-aSL

lambda5L=uS+aSR
lambda5R=uR+aR

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
lambda(1)=ut-at
lambda(2)=ut
lambda(3)=ut
lambda(4)=ut
lambda(5)=ut+at

if(debug==1)then
 print*,'Characetristic Speeds(1,3,5)=',lambda(1),lambda(3),lambda(5)
end if

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(1)=1._8
K1(2)=ut-at
K1(3)=vt
K1(4)=wt
K1(5)=Ht-ut*at

K2(1)=1._8
K2(2)=ut
K2(3)=vt
K2(4)=wt
K2(5)=0.5*(ut**2.+vt**2.+wt**2.)

K3(1)=0._8
K3(2)=0._8
K3(3)=1._8
K3(4)=0._8
K3(5)=vt

K4(1)=0._8
K4(2)=0._8
K4(3)=0._8
K4(4)=1._8
K4(5)=wt

K5(1)=1._8
K5(2)=ut+at
K5(3)=vt
K5(4)=wt
K5(5)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(2)=delrho-delp/(at*at)
alpha(3)=rhot*delv
alpha(4)=rhot*delw
alpha(5)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,abs(lambda(1)),abs(lambda(3)),abs(lambda(5)))


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

function minmod(x,y) result(z)

!Input Variables
real*8,intent(in)::x,y

!Output Variables
real*8::z
  
z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod


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
integer,parameter::ycell=ny*0.5
integer::jj

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=max(u2(jj,ycell,1),min_dens)
   velx=u2(jj,ycell,2)/u2(jj,ycell,1)
   vely=u2(jj,ycell,3)/u2(jj,ycell,1)
   velz=u2(jj,ycell,4)/u2(jj,ycell,1)
   pres=(gam-1.)*(u2(jj,ycell,5)-0.5*dens*(velx**2.+vely**2.+velz**2.))   
   write(12,*) x,dens,velx,vely,pres,velz
   write(22,*) x,f(jj,ycell,1),f(jj,ycell,2),f(jj,ycell,3),f(jj,ycell,4),f(jj,ycell,5)
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx*1.0
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


end module RoeSolver2D_mod
