!Roe's Approximate (Linear) Riemann Solver
module RoeSolver_mod
implicit none

!-------------------------------------------------------------------------
integer,parameter::flux_option=1 !1:mh 2:harten 3:ENO
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y
integer,parameter::discontinuityDir=2 !1:x, 2:y
integer,parameter::debug=1 !0:off 1:on
integer,parameter::shockSwitch=1 !0:off, 1:on
integer,parameter::fluidInit=3 !1:V Shock, 2:Oblique Shock, 3:Riemann Problem File input
!-------------------------------------------------------------------------
integer,parameter::nt=100
integer,parameter::nx=50 
integer,parameter::ny=50
real*8,parameter::del0=0.5
real*8,parameter::Q_user=2.0 
integer,parameter::characteristicLimiting=0 !0:off, 1:on 
real*8::premin=1.d-25,densmin=1.d-25
integer ::i,j,k,l


real*8 ::dt,dtx,dty,dx,dy,dr !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,aL,HL,FLx(5) !left state
real*8::rhoR,uR,vR,wR,pR,aR,HR,FRx(5) !right state 
real*8::delrho,delp,delu,delv,delw,smax

real*8::u1(-1:nx+2,-1:ny+2,5),u2(-1:nx+2,-1:ny+2,5)

!----------------------------------------------------------------------
real*8::gt(-1:nx+2,5),fc(-1:nx+2,5),gamma(-1:nx+2,5),sigma(-1:nx+2,5)
real*8::f(-1:nx+2,5),g(-1:nx+2,5),fl(5),u(5)

real*8::dens,velx,vely,velz,pres,a,H
real*8::rhot,ut,vt,wt,Ht,at,lambda(-1:nx+2,5),alpha(-1:nx+2,5)
real*8::K1(-1:nx+2,5),K2(-1:nx+2,5),K3(-1:nx+2,5)
real*8::K4(-1:nx+2,5),K5(-1:nx+2,5)	

real*8 ::temp1,temp2,temp3,xc,xt,theta,w
!------------------------------------------------------------------------------
real*8::lambda1L(-1:nx+2),lambda1R(-1:nx+2),lambda5L(-1:nx+2),lambda5R(-1:nx+2)
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR
real*8::fright(-1:nx+2,5),fleft(-1:nx+2,5)
real*8::qL(-1:nx+2,5),qR(-1:nx+2,5),slope(-1:nx+2,5)
real*8::delta1(-1:nx+2,5),delta2(-1:nx+2,5),delta3(-1:nx+2,5)
real*8::delta4(-1:nx+2,5),delta5(-1:nx+2,5)
real*8::ubL(-1:nx+2,5),ubR(-1:nx+2,5)

!------------------------------------------------------------------------------
real*8::qh(-1:nx+2,5) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::v(-1:nx+2,5) !characteristic variables
real*8::Qx(-1:nx+2,5,5),Qinv(-1:nx+2,5,5) !characeteristic matrix
real*8::S(-1:nx+2,5)
real*8::dqdx(-1:nx+2,5),dqdt(-1:nx+2,5)
real*8::dpdx,dudx,dpdt,dudt

!------------------------------------------------------------------------------

character(len=20) :: filename

!------------------------------------------------------------------------------

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
real*8::velx0,velx1,velx2,velx3,vely0,vely1,vely2,vely3
real*8::velz0,velz1,velz2,velz3,rho0,rho1
real*8::pre0,pre1,M0,temp1,temp2,temp3,temp4
real*8::theta !shock angle w.r.t. x axis
real*8::alpha
real*8::Sh !shock speed 
real*8::zeta,kz
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
read(unit=10,fmt=*) rhoR,uR,vR,wL,pR


if(fluidInit==1)then

!Intersecting shocks 	 
         
!pre-shock conditions
theta=pi/4.   !shock angle set to 45 degrees		
rho0=1.
pre0=1.	
velx0=10.!1.
vely0=0.
velz0=0.
M0=sin(theta)*velx0/sqrt(gam*pre0/rho0) !20.	  
print*,'M0=',M0       
!using Rankine Hugoniot relations to determine post shock fluid state	  
rho1=rho0*(((gam+1.)*M0*M0)/(2.+(gam-1.)*M0*M0))
pre1=pre0*(2.*gam*M0*M0-gam+1.)/(gam+1.)
velx1=velx0*(cos(theta)**2+(rho0/rho1)*sin(theta)**2)
vely1=-velx0*cos(theta)*sin(theta)*(1-(rho0/rho1))
velz1=0.
velx2=velx1
vely2=-vely1
velz2=0.
velz3=0.
		
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
    u1(jj,kk,4) = velz0*u1(jj,kk,1)
    u1(jj,kk,5) = (pre0/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
   !post-shock
		
   !region 2
  elseif (y1 .lt. temp2 .and. y1 .gt. temp3) then
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx2*u1(jj,kk,1)
    u1(jj,kk,3) = vely2*u1(jj,kk,1)
    u1(jj,kk,4) = velz2*u1(jj,kk,1)
    u1(jj,kk,5) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
    !region 1
  elseif (y1 .lt. temp4 .and. y1 .gt. temp1) then
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx1*u1(jj,kk,1)
    u1(jj,kk,3) = vely1*u1(jj,kk,1)
    u1(jj,kk,4) = velz1*u1(jj,kk,1)
    u1(jj,kk,5) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
   !region 3(deadzone)
  elseif (y1 .gt. temp4 .and. y1 .lt. temp3) then
    !print*,'Deadzone checkpoint...'
    u1(jj,kk,1) = rho1
    velx3=velx1
    !linear interpolation to get vely3, continuous across dead zone boundary  
    vely3=vely1+((vely1-vely2)/(temp4-temp3))*(y1-temp4)
    u1(jj,kk,2) = velx3*u1(jj,kk,1)
    u1(jj,kk,3) = vely3*u1(jj,kk,1)
    u1(jj,kk,4) = velz3*u1(jj,kk,1)
    u1(jj,kk,5) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
   end if			
end do
end do

else if(fluidInit==2)then
!Oblique shocks 	 
         
!pre-shock conditions
theta=pi/4.   !shock angle set to 45 degrees		
rho0=1.
pre0=1.	
velx0=10.!1.
vely0=0.
velz0=0.

!Define Shock Speed
Sh=5.*velx0

M0=sin(theta)*(velx0-Sh)/sqrt(gam*pre0/rho0)	  
print*,'M0=',M0       

!using Rankine Hugoniot relations to determine post shock fluid state	  
rho1=rho0*(((gam+1.)*M0*M0)/(2.+(gam-1.)*M0*M0))
pre1=pre0*(2.*gam*M0*M0-gam+1.)/(gam+1.)
velx1=velx0*(cos(theta)**2+(rho0/rho1)*sin(theta)**2)
vely1=velx0*cos(theta)*sin(theta)*(1-(rho0/rho1))
velz1=0.		
		

do kk=-1,ny+2
 y1 = (kk - 0.5)*dy
 do jj=-1,nx+2
  x1 = (jj - 0.5)*dx
				
  !pre-shock	
  if (y1 .ge. x1) then 
    u1(jj,kk,1) = rho0
    u1(jj,kk,2) = velx0*u1(jj,kk,1)
    u1(jj,kk,3) = vely0*u1(jj,kk,1) 
    u1(jj,kk,4) = velz0*u1(jj,kk,1)
    u1(jj,kk,5) = (pre0/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
  !post-shock  
  else
    u1(jj,kk,1) = rho1
    u1(jj,kk,2) = velx1*u1(jj,kk,1)
    u1(jj,kk,3) = vely1*u1(jj,kk,1) 
    u1(jj,kk,4) = velz1*u1(jj,kk,1)
    u1(jj,kk,5) = (pre1/(gam - 1.0d0)) + 5.0d-1&
                 *(u1(jj,kk,2)*u1(jj,kk,2)+u1(jj,kk,3)*u1(jj,kk,3) &
                  +u1(jj,kk,4)*u1(jj,kk,4))/u1(jj,kk,1)
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
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*wL*wL)+pL*gamul
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul 
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
    u1(jj,kk,5)=0.5*rhoL*(uL*uL+vL*vL+wL*wL)+pL*gamul
  else
    u1(jj,kk,1)=rhoR
    u1(jj,kk,2)=rhoR*uR
    u1(jj,kk,3)=rhoR*vR
    u1(jj,kk,4)=rhoR*wR
    u1(jj,kk,5)=0.5*rhoR*(uR*uR+vR*vR+wR*wR)+pR*gamul 
  end if
 end do
end do

end if

else if(fluidInit==4)then
!Kelvin-Helmholtz Instability 
rhoL=1.0
uL=-1.
vL=0.0
wL=0.
pL=1.0

rhoR=1.0
uR=1.
vR=0.0
wR=0.
pR=1.0

kz=2.*pi*3./(xmax-xmin)

do jj=-1,nx+2

 x1 =xmin+(jj - 0.5)*dx
 
 do kk=-1,ny+2
  y1 =ymin+(kk - 0.5)*dy
  
  vL=1.*(ymax-ymin)*sin(kz*x1)*exp(-abs(y1-(ymax-ymin)))
 
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


end if

u2=u1



!call fileOutput(0+1000)
call horcut()
call vercut()
close(unit=0+1000)

end subroutine init

subroutine RoeSolver(jj)

!Input Variables
integer::jj

delrho=rhoR-rhoL
delu=uR-uL
delv=vR-vL
delw=wR-wL
delp=pR-pL

if(debug==1)then
  print*,'Cell: ',jj
  print*,'Left State(rho,u,v,w,p)=',rhoL,uL,vL,wL,pL
  print*,'Right State(rho,u,v,w,p)=',rhoR,uR,vR,wR,pR
end if



!-------------------------------------------
!Roe-Averaged State Variables: ^rho,^u,^H,^u
!-------------------------------------------
rhot=sqrt(rhoL*rhoR)
ut=(sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR))
vt=(sqrt(rhoL)*vL+sqrt(rhoR)*vR)/(sqrt(rhoL)+sqrt(rhoR))
wt=(sqrt(rhoL)*wL+sqrt(rhoR)*wR)/(sqrt(rhoL)+sqrt(rhoR))
Ht=(sqrt(rhoL)*HL+sqrt(rhoR)*HR)/(sqrt(rhoL)+sqrt(rhoR))
at=sqrt((gam-1.)*(Ht-0.5*(ut**2.+vt**2.+wt**2.)) )



!---------------------------------------------
!Roe Jacobian Matrix Eigenvalues:^lambda_1,2,3
!---------------------------------------------
lambda(jj,1)=ut-at
lambda(jj,2)=ut
lambda(jj,3)=ut
lambda(jj,4)=ut
lambda(jj,5)=ut+at

!---------------------------------------------
!Roe Jacobian Matrix Eigenvectors:^K_1,2,3
!---------------------------------------------
K1(jj,1)=1._8
K1(jj,2)=ut-at
K1(jj,3)=vt
K1(jj,4)=wt
K1(jj,3)=Ht-ut*at

K2(jj,1)=1
K2(jj,2)=ut
K2(jj,3)=vt
K2(jj,4)=wt
K2(jj,3)=0.5*(ut**2.+vt**2.+wt**2.)

K3(jj,1)=0.
K3(jj,2)=0.
K3(jj,3)=1.
K3(jj,4)=0.
K3(jj,5)=vt

K4(jj,1)=0.
K4(jj,2)=0.
K4(jj,3)=0.
K4(jj,4)=1.
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

!print*,'Checkpoint 5.1'

!---------------------------------------
!Max Wave Speed
!---------------------------------------
!smax=abs(ut)+abs(at)
smax=max(smax,abs(lambda(jj,1)),abs(lambda(jj,2)),abs(lambda(jj,5)))

end subroutine RoeSolver

subroutine computeFluxHarten(n,ind,dir)

!Input Variables
integer::n,ind,dir
!Local Variables
integer::jj,kk


do jj=-1,n+1
!---------------------------------
!w_j+1/2(0)
!---------------------------------

if(dir==1)then

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,u1(jj,ind,1))
uL=u1(jj,ind,2)/u1(jj,ind,1)
vL=u1(jj,ind,3)/u1(jj,ind,1)
wL=u1(jj,ind,4)/u1(jj,ind,1)
pL=max(premin,(gam-1.)*( u1(jj,ind,5)-0.5*rhoL &
  *(uL**2.+vL**2.+wL**2.) ) )
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

rhoR=max(densmin,u1(jj+1,ind,1))
uR=u1(jj+1,ind,2)/u1(jj+1,ind,1)
vR=u1(jj+1,ind,3)/u1(jj+1,ind,1)
wR=u1(jj+1,ind,4)/u1(jj+1,ind,1)
pR=max(premin,(gam-1.)*( u1(jj+1,ind,5)-0.5*rhoR &
  *(uR**2.+vR**2.+wR**2.) ) )
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR

else if(dir==2)then

!set left and right states for local Riemann problem at j+1/2 cell interface
rhoL=max(densmin,u1(ind,jj,1))
uL=u1(ind,jj,3)/u1(ind,jj,1)
vL=u1(ind,jj,2)/u1(ind,jj,1)
wL=u1(ind,jj,4)/u1(ind,jj,1)
pL=max(premin,(gam-1.)*( u1(ind,jj,5)-0.5*rhoL &
  *(uL**2.+vL**2.+wL**2.) ) )
aL=sqrt(gam*pL/rhoL)
HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

rhoR=max(densmin,u1(ind,jj+1,1))
uR=u1(ind,jj+1,3)/u1(ind,jj+1,1)
vR=u1(ind,jj+1,2)/u1(ind,jj+1,1)
wR=u1(ind,jj+1,4)/u1(ind,jj+1,1)
pR=max(premin,(gam-1.)*( u1(ind,jj+1,5)-0.5*rhoR &
  *(uR**2.+vR**2.+wR**2.) ) )
aR=sqrt(gam*pR/rhoR)
HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR


end if

!-----------------------------------------------
!f_j
!-----------------------------------------------
u(1)=rhoL
u(2)=uL
u(3)=vL
u(4)=wL
u(5)=pL
call flux(u)
do kk=1,5
  f(jj,kk)=fl(kk)
end do

call RoeSolver(jj)

end do

!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
if(dir==1)then
 dtx=dx*cour/smax 
 dt=dtx
 dr=dx
else if(dir==2)then
 dty=dy*cour/smax
 dt=dty
 dr=dx 
end if

!------------------------------------------------
!~g_i+1/2
!------------------------------------------------
do jj=-1,nx+1
 do kk=1,5
   gt(jj,kk)=(0.5*Q(lambda(jj,kk))-0.5*(dt/dr) &
             *(lambda(jj,kk)*lambda(jj,kk)))*alpha(jj,kk) 
   if(shockSwitch==1)then
     sigma(jj,kk)=0.5*(1.-Q(lambda(jj,kk)))
   end if	 
 end do  
end do

do jj=0,nx+1
  !-----------------------------------------------
  !f(c)_j
  !-----------------------------------------------
   do kk=1,5
     fc(jj,kk)=minmod(gt(jj,kk),gt(jj-1,kk))
   end do  
   
   if(shockSwitch==1)then
    do kk=1,5
     if(abs(alpha(jj,kk))+abs(alpha(jj-1,kk)) .ne. 0._8)then
      theta=(abs(alpha(jj,kk)-alpha(jj-1,kk)))/(abs(alpha(jj,kk))+abs(alpha(jj-1,kk)))
     else
      theta=0.
     end if 
     w=minmod(alpha(jj,kk)*sigma(jj,kk),alpha(jj-1,kk)*sigma(jj-1,kk))
     
     !Modify the entropy wave only (will result in steeper contact discontinuity)
     if(kk==2)then ! .or. k==3 .or. k==4)then	 
       fc(jj,kk)=fc(jj,kk)+w*theta
     end if

    end do
   end if
   
 !-----------------------------------------------
 !f_j
 !-----------------------------------------------
! u(1)=u1(j,1)
! u(2)=u1(j,2)/u1(j,1)
! u(3)=u1(j,3)/u1(j,1)
! u(4)=u1(j,4)/u1(j,1)
! u(5)=(gam-1.)*(u1(j,5)-0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.))
! call flux(u)
!   do k=1,5
!     f(j,k)=fl(k)
!   end do

end do

!------------------------------------------------
!gamma_i+1/2
!------------------------------------------------
do jj=0,nx  
 do kk=1,5
   if(alpha(jj,kk)>=1.d-30)then
     gamma(jj,kk)=(fc(jj+1,kk)-fc(jj,kk))/alpha(jj,kk)
   else
     gamma(jj,kk)=0.
   end if	  
 end do
  
!------------------------------------------------
!Modified Numerical Flux: ^g_i+1/2
!------------------------------------------------
 do kk=1,5
  g(jj,kk)=0.5*(f(jj,kk)+f(jj+1,kk))-0.5*&
        ( K1(jj,kk)*( Q(lambda(jj,1)+gamma(jj,1))*alpha(jj,1)-fc(jj,1)-fc(jj+1,1) ) &
         +K2(jj,kk)*( Q(lambda(jj,2)+gamma(jj,2))*alpha(jj,2)-fc(jj,2)-fc(jj+1,2) ) &
	 +K3(jj,kk)*( Q(lambda(jj,3)+gamma(jj,3))*alpha(jj,3)-fc(jj,3)-fc(jj+1,3) ) &
         +K4(jj,kk)*( Q(lambda(jj,4)+gamma(jj,4))*alpha(jj,4)-fc(jj,4)-fc(jj+1,4) ) &  
         +K5(jj,kk)*( Q(lambda(jj,5)+gamma(jj,5))*alpha(jj,5)-fc(jj,5)-fc(jj+1,5) ) ) 
 end do
end do
 
end subroutine computeFluxHarten

subroutine computeFluxMH(n,ind,dir)

!Input Variables
integer::n,ind,dir
!Local Variables
integer::jj,kk
real*8::vsq

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
if(dir==1)then
 do jj=0,n+1
  do kk=1,5  
    qL(jj,kk)=u1(jj,ind,kk)-0.5*slope(jj,kk)
    qR(jj,kk)=u1(jj,ind,kk)+0.5*slope(jj,kk)
  end do
 end do
else if(dir==2)then
 do jj=0,nx+1
  do kk=1,5    
    qL(jj,kk)=u1(ind,jj,kk)-0.5*slope(jj,kk)
    qR(jj,kk)=u1(ind,jj,kk)+0.5*slope(jj,kk)
  end do
 end do
end if

!-------------------------------------------------------------------
!Evolve Boudary Extrapolated values by a half-step: ubarL_i, ubarR_i
!-------------------------------------------------------------------
if(dir==1)then
 !dt=dtx
 !dr=dx

 do jj=0,n+1
  dens=max(qL(jj,1),densmin)
  velx=qL(jj,2)/qL(jj,1)
  vely=qL(jj,3)/qL(jj,1)
  velz=qL(jj,4)/qL(jj,1)
  vsq=velx**2.+vely**2.+velz**2.
  pres=max(premin,(gam-1.)*(qL(jj,5)-0.5*dens*vsq))
 !pres=max( (gam-1.)*( qL(jj,5)-0.5*dens*(velx**2. & 
 !      +vely**2.+velz**2.) ),premin )

  FLx(1)=dens*velx
  FLx(2)=dens*velx*velx+pres
  FLx(3)=dens*velx*vely
  FLx(4)=dens*velx*velz
  FLx(5)=0.5*dens*vsq*velx+gam*gamul*pres*velx

  dens=max(qR(jj,1),densmin)
  velx=qR(jj,2)/qR(jj,1)
  vely=qR(jj,3)/qR(jj,1)
  velz=qR(jj,4)/qR(jj,1)
  vsq=velx**2.+vely**2.+velz**2.
  pres=max(premin,(gam-1.)*(qR(jj,5)-0.5*dens*vsq))
  !pres=max( (gam-1.)*( qR(jj,5)-0.5*dens*(velx**2. & 
  !     +vely**2.+velz**2.) ),premin )

  FRx(1)=dens*velx
  FRx(2)=dens*velx*velx+pres
  FRx(3)=dens*velx*vely
  FRx(4)=dens*velx*velz
  FRx(5)=0.5*dens*vsq*velx+gam*gamul*pres*velx

  do kk=1,5
   ubL(jj,kk)=qL(jj,kk)+0.5*(dtx/dx)*(FLx(kk)-FRx(kk))
   ubR(jj,kk)=qR(jj,kk)+0.5*(dtx/dx)*(FLx(kk)-FRx(kk))
  end do
end do

else if(dir==2)then
 !dt=dty
 !dr=dy

do jj=0,n+1
  dens=max(qL(jj,1),densmin)
  velx=qL(jj,2)/qL(jj,1)
  vely=qL(jj,3)/qL(jj,1)
  velz=qL(jj,4)/qL(jj,1)
  vsq=velx**2.+vely**2.+velz**2.
  pres=max(premin,(gam-1.)*(qL(jj,5)-0.5*dens*vsq))

  FLx(1)=dens*vely
  FLx(2)=dens*velx*vely
  FLx(3)=dens*vely*vely+pres
  FLx(4)=dens*velz*vely
  FLx(5)=0.5*dens*vsq*vely+gam*gamul*pres*vely

  dens=max(qR(jj,1),densmin)
  velx=qR(jj,2)/qR(jj,1)
  vely=qR(jj,3)/qR(jj,1)
  vely=qR(jj,4)/qR(jj,1)
  vsq=velx**2.+vely**2.+velz**2.
  pres=max(premin,(gam-1.)*(qR(jj,5)-0.5*dens*vsq))

  FRx(1)=dens*vely
  FRx(2)=dens*velx*vely
  FRx(3)=dens*vely*vely+pres
  FRx(4)=dens*velz*vely
  FRx(5)=0.5*dens*vsq*vely+gam*gamul*pres*vely

   do kk=1,5
    ubL(jj,kk)=qL(jj,kk)+0.5*(dty/dy)*(FLx(kk)-FRx(kk))
    ubR(jj,kk)=qR(jj,kk)+0.5*(dty/dy)*(FLx(kk)-FRx(kk))
   end do
end do

end if


!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do jj=0,n
!set left and right states for local Riemann problem at j+1/2 cell interface

if(dir==1)then
rhoL=max(densmin,ubR(jj,1))
uL=ubR(jj,2)/ubR(jj,1)
vL=ubR(jj,3)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
vsq=uL**2.+vL**2.+wL**2.
pL=max( premin,(gam-1.)*(ubR(jj,5)-0.5*rhoL*vsq) )
  
aL=sqrt(gam*pL/rhoL)
HL=0.5*vsq+gamul*aL*aL

rhoR=max(densmin,ubL(jj+1,1))
uR=ubL(jj+1,2)/ubL(jj+1,1)
vR=ubL(jj+1,3)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
vsq=uR**2.+vR**2.+wR**2.
pR=max(premin,(gam-1.)*(ubR(jj+1,5)-0.5*rhoR*vsq) )

aR=sqrt(gam*pR/rhoR)
HR=0.5*vsq+gamul*aR*aR

else if(dir==2)then
rhoL=max(densmin,ubR(jj,1))
uL=ubR(jj,3)/ubR(jj,1)
vL=ubR(jj,2)/ubR(jj,1)
wL=ubR(jj,4)/ubR(jj,1)
vsq=uL**2.+vL**2.+wL**2.
pL=max(premin,(gam-1.)*(ubR(jj,5)-0.5*rhoL*vsq) )

aL=sqrt(gam*pL/rhoL)
HL=0.5*vsq+gamul*aL*aL

rhoR=max(densmin,ubL(jj+1,1))
uR=ubL(jj+1,3)/ubL(jj+1,1)
vR=ubL(jj+1,2)/ubL(jj+1,1)
wR=ubL(jj+1,4)/ubL(jj+1,1)
vsq=uR**2.+vR**2.+wR**2.
pR=max(premin,(gam-1.)*(ubL(jj+1,5)-0.5*rhoR*vsq) )

aR=sqrt(gam*pR/rhoR)
HR=0.5*vsq+gamul*aR*aR

end if

fleft(jj,1)=rhoL*uL
fleft(jj,2)=rhoL*uL*uL+pL
fleft(jj,3)=rhoL*uL*vL
fleft(jj,4)=rhoL*uL*wL
fleft(jj,5)=0.5*rhoL*(uL**2.+vL**2.+wL**2.)*uL+gam*gamul*pL*uL

fright(jj,1)=rhoR*uR
fright(jj,2)=rhoR*uR*uR+pR
fright(jj,3)=rhoR*uR*vR
fright(jj,4)=rhoR*uR*wR
fright(jj,5)=0.5*rhoR*(uR**2.+vR**2.+wR**2.)*uR+gam*gamul*pR*uR

!if(dir==1)then
!  print*,'Cell=',jj,ind
!else
!  print*,'Cell=',ind,jj
!end if
!print*,'pL,pR=',pL,pR

call RoeSolver(jj) 

!print*,'Checkpoint 5.5...jj=',jj

if(debug==1)then
  print*,'Max wave speed=',smax
end if

end do


!-----------------------------------------------
!Time-Step Size: dt
!-----------------------------------------------
if(dir==1)then
 dtx=dx*cour/smax 
 dt=dtx
else if(dir==2)then
 dty=dy*cour/smax
 dt=dty
end if

!-------------------------------------------------------------
!Numerical Flux: F_j+1/2
!-------------------------------------------------------------
do jj=0,n
!Check for transonic rarefaction waves
if(lambda1L(jj)<0. .and. lambda1R(jj)>0.)then
  do kk=1,5
     g(jj,kk)= fleft(jj,kk)+(lambda1L(jj)*(lambda1R(jj)-lambda(jj,1))/(lambda1R(jj)-lambda1L(jj)))&
            *alpha(jj,1)*K1(jj,kk)  
  end do
else if(lambda5L(jj)<0. .and. lambda5R(jj)>0.)then
  do kk=1,5
     g(jj,kk)=fright(jj,kk)-(lambda5R(jj)*(lambda(jj,5)-lambda5L(jj))/(lambda5R(jj)-lambda5L(jj)))&
            *alpha(jj,5)*K5(jj,kk)
  end do
else
  do kk=1,5
     g(jj,kk)=0.5*(fleft(jj,kk)+fright(jj,kk))-0.5*(alpha(jj,1)*abs(lambda(jj,1))*K1(jj,kk) &
          +alpha(jj,2)*abs(lambda(jj,2))*K2(jj,kk)+alpha(jj,3)*abs(lambda(jj,3))*K3(jj,kk) &
          +alpha(jj,4)*abs(lambda(jj,4))*K4(jj,kk)+alpha(jj,5)*abs(lambda(jj,5))*K5(jj,kk) ) 
  end do
end if

end do

end subroutine computeFluxMH


subroutine starRegion()

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
pS=gL(p0)*pL+gR(p0)*pR-(uR-uL)
pS=pS/(gL(p0)+gR(p0))
!u_star
uS=0.5*(uL+uR)+0.5*((pS-pR)*gR(p0)-(pS-pL)*gL(p0))
!rho_star
rhoSL=rhoL*((pS/pL)+gamel)/((pS/pL)*gamel+1.)
rhoSR=rhoR*((pS/pR)+gamel)/((pS/pR)*gamel+1.)
end if

aSL=sqrt(gam*pS/rhoSL)
aSR=sqrt(gam*pS/rhoSR)

end subroutine starRegion

function gL(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoL)/(x+gamel*pL))
end function gL

function gR(x) result(pxx)
!Input Variables
real*8,intent(in)::x
!Output Variables
real*8::pxx

  pxx=sqrt((2.*gamil/rhoR)/(x+gamel*pR))
end function gR


subroutine flux(u)

!Input Variables    
real*8::u(5)

fl(1)=u(1)*u(2)
fl(2)=u(1)*u(2)*u(2)+u(5)
fl(3)=u(1)*u(2)*u(3)
fl(4)=u(1)*u(2)*u(4)
fl(5)=0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.)*u(2)+gam*gamul*u(2)*u(5)
  
end subroutine flux
 
function Q(x) result(qx)

!Input Variables  
real*8,intent(in)::x
!Output Variables 
real*8::qx
 
if(abs(x)<del0)then
  qx=(x*x+del0*del0)/(2.*del0)
else if(abs(x)>=del0)then
  qx=abs(x)
end if

end function Q
 
function minmod(x,y) result(z)

!Input Variables
real*8,intent(in)::x,y

!Output Variables
real*8::z
  
z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod

function limiter(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z,b

  b=1. !minmod limiter for b=1

  if(y>0._8)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter

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

subroutine fileOutput(iunit)
!Input variables
integer::iunit
!local variables
integer::jj,kk

do kk=1,ny
  y=ymin+(kk-0.5)*dy
  do jj=1,nx
   x=xmin+(jj-0.5)*dx
   dens=u2(jj,kk,1)
   velx=u2(jj,kk,2)/u2(jj,kk,1)
   vely=u2(jj,kk,3)/u2(jj,kk,1)
   pres=(gam-1.)*(u2(jj,kk,4)-0.5*(u2(jj,kk,2)*u2(jj,kk,2)/u2(jj,kk,1)&
        +u2(jj,kk,3)*u2(jj,kk,3)/u2(jj,kk,1)))
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
  dens=u2(jj,ycell,1)
  velx=u2(jj,ycell,2)/u2(jj,ycell,1)
  vely=u2(jj,ycell,3)/u2(jj,ycell,1)
  velz=u2(jj,ycell,4)/u2(jj,ycell,1)
  pres=(gam-1.)*( u2(jj,ycell,5)-0.5*dens &
  *(velx**2.+vely**2.+velz**2.) )  
   write(12,*) x,dens,velx,vely,pres
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx*0.55
integer::jj

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u2(xcell,jj,1)
   velx=u2(xcell,jj,2)/u2(xcell,jj,1)
   vely=u2(xcell,jj,3)/u2(xcell,jj,1)
   velz=u2(xcell,jj,4)/u2(xcell,jj,1)
   pres=(gam-1.)*( u2(xcell,jj,5)-0.5*dens &
  *(velx**2.+vely**2.+velz**2.) ) 
   write(13,*) y,dens,velx,vely,pres
end do

end subroutine vercut


end module RoeSolver_mod
