!Exact Riemann Solver of 1D Eulers Equiations (Improved version, nearly perfect)
!Two root solvers included: Secant Method and Newton
!Newton solver performs better

module RoeSolver2D_mod
implicit none

integer,parameter::discontinuityDir=1 !1:x, 2:y
integer,parameter::debug=0 !0:off 1:on
integer,parameter::nt=100
integer,parameter::nx=200
integer,parameter::ny=200
real*8,parameter::s_tol=1.d-25
real,parameter::min_pres=1.d-30
real,parameter::min_dens=1.d-30
integer,parameter::vshockInit=1 !0:off :1:on 

real*8 ::dt,dx,dy !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,pL,aL,HL,FL(4),GL(4) !left state
real*8::rhoR,uR,vR,pR,aR,HR,FR(4),GR(4) !right state
real*8::rhot,ut,vt,Ht,at,lambda(4),K1(4),K2(4),K3(4),K4(4),alpha(4)	
real*8::delrho,delp,delu,delv,smax
real*8::dens,velx,vely,pres
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

print*,'First Order Godunov Method, using Roe Solver...'

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

smax=0.

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
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL*uL+vL*vL)+gamul*aL*aL


rhoR=max(u1(jj+1,ind,1),min_dens)
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
pR=max((gam-1.)*(u1(jj+1,ind,4)-0.5*(u1(jj+1,ind,2)*u1(jj+1,ind,2)/u1(jj+1,ind,1)&
        +u1(jj+1,ind,3)*u1(jj+1,ind,3)/u1(jj+1,ind,1))),min_pres)
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR*uR+vR*vR)+gamul*aR*aR

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delp=pR-pL
 
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
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL*uL+vL*vL)+gamul*aL*aL

rhoR=max(u1(ind,jj+1,1),min_dens)
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
pR=max((gam-1.)*(u1(ind,jj+1,4)-0.5*(u1(ind,jj+1,2)*u1(ind,jj+1,2)/u1(ind,jj+1,1)&
        +u1(ind,jj+1,3)*u1(ind,jj+1,3)/u1(ind,jj+1,1))),min_pres)
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR*uR+vR*vR)+gamul*aR*aR

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delp=pR-pL

 
GL(1)=rhoL*uL 
GL(2)=rhoL*uL*uL+pL
GL(3)=rhoL*uL*vL
GL(4)=0.5*rhoL*(uL*uL+vL*vL)*uL+gam*gamul*pL*uL

GR(1)=rhoR*uR  
GR(2)=rhoR*uR*uR+pR
GR(3)=rhoR*uR*vR
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

!-------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^u
!-------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
vt=(sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*(ut*ut+vt*vt)))

!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(1)=ut-at
lambda(2)=ut
lambda(3)=ut
lambda(4)=ut+at

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(1)=1._8
K1(2)=ut-at
K1(3)=vt
K1(4)=Ht-ut*at

K2(1)=1
K2(2)=ut
K2(3)=vt
K2(4)=0.5*ut*ut

K3(1)=0.
K3(2)=0.
K3(3)=1.
K3(4)=vt


K4(1)=1._8
K4(2)=ut+at
K4(3)=vt
K4(4)=Ht+ut*at

!---------------------------------------------
!Roe Averaged Wave Strengths:^alpha_1,2,3
!---------------------------------------------
alpha(1)=(0.5/(at*at))*(delp-rhot*at*delu)
alpha(2)=delrho-delp/(at*at)
alpha(3)=rhot*delv
alpha(4)=(0.5/(at*at))*(delp+rhot*at*delu)

!---------------------------------------
!Max Wave Speed
!---------------------------------------
smax=max(smax,lambda(1),lambda(3),lambda(4))

if(dir==1)then
!---------------------------------------
!F_j+1/2
!---------------------------------------
 do kk=1,4
   xflux(jj,kk)=0.5*(FL(kk)+FR(kk))-0.5*( alpha(1)*abs(lambda(1))*K1(kk)&
          +alpha(2)*abs(lambda(2))*K2(kk)+alpha(3)*abs(lambda(3))*K3(kk)&
          +alpha(4)*abs(lambda(4))*K4(kk) ) 
 end do 
else if(dir==2)then
!---------------------------------------
!G_j+1/2
!---------------------------------------
 do kk=1,4
   yflux(jj,kk)=0.5*(GL(kk)+GR(kk))-0.5*( alpha(1)*abs(lambda(1))*K1(kk)&
          +alpha(2)*abs(lambda(2))*K2(kk)+alpha(3)*abs(lambda(3))*K3(kk)&
          +alpha(4)*abs(lambda(4))*K4(kk) ) 
 end do 
end if

end do

end subroutine computeFlux

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

end module RoeSolver2D_mod
