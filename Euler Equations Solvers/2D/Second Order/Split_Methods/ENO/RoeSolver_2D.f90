!ENO Scheme Harten-Osher-Engquist-Chakravathy for Euler's Equation
!Riemann Solver: Roe w/ Harten-Hyman Entropy Fix


module RoeSolver2D_mod
implicit none

integer,parameter::discontinuityDir=1 !1:x, 2:y
integer,parameter::debug=0 !0:off 1:on
integer,parameter::boundaryType=2 !1:outflow 2:periodic in x, outflow in y
real*8,parameter::Q_user=2.0 
integer,parameter::perturbationType=2 !1:Sinusoidal 2:Random

integer,parameter::nt=2000
integer,parameter::nx=250
integer,parameter::ny=250
real*8,parameter::del0=0.15
real*8,parameter::s_tol=1.d-25
real*8,parameter::min_pres=1.d-30
real*8,parameter::min_dens=1.d-30
integer,parameter::fluidInit=4 !1:V Shock, 2:Oblique Riemann Problem, 3:Riemann Problem File input 4:Kelvin-Helmholtz Instability 

real*8 ::dt,dtx,dty,dx,dy,dr !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

real*8::rhoL,uL,vL,wL,pL,aL,HL !left state
real*8::rhoR,uR,vR,wR,pR,aR,HR !right state
real*8::rhot,ut,vt,wt,Ht,at,lambda(-1:nx+2,5),alpha(-1:nx+2,5)
real*8::K1(-1:nx+2,5),K2(-1:nx+2,5),K3(-1:nx+2,5),K4(-1:nx+2,5),K5(-1:nx+2,5)
	
real*8::delrho,delp,delu,delv,delw,smax
real*8::dens,velx,vely,velz,pres,a,H,u(5),fl(5)
real*8::u1(-1:nx+2,-1:ny+2,5),u2(-1:nx+2,-1:ny+2,5)
real*8::f(-1:nx+2,5),g(-1:nx+2,5)

real*8::lambda1L(-1:nx+2),lambda1R(-1:nx+2),lambda5L(-1:nx+2),lambda5R(-1:nx+2)
real*8::pS,rhoSL,rhoSR,uS,aSL,aSR

!------------------------------------------------------------------------------


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
real*8::velx0,velx1,vely0,vely1,velz0,velz1,rho0,rho1
real*8::velx2,vely2,velx3,vely3
real*8::pre0,pre1,M0,temp1,temp2,temp3,temp4
real*8::theta !shock angle w.r.t. x axis
real*8::alpha
real*8::Sh !shock speed 
real*8::zeta,kz
integer::mode
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
    vL=0.5*(ymax-ymin)*sin(kz*x1)*exp(-100.*(kk-ny/2)**2.)
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
end do

end subroutine horcut

subroutine vercut()
!Local Variables
integer,parameter::xcell=nx*0.25
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
real*8::qh(-1:nx+2,5) 
real*8::v(-1:nx+2,5) 
real*8::Qx(-1:nx+2,5,5),Qinv(-1:nx+2,5,5) 
real*8::S(-1:nx+2,5)
real*8::dqdx(-1:nx+2,5),dqdt(-1:nx+2,5)
real*8::dpdx,dudx,dpdt,dudt

 

!-----------------------------------------------------------------
!Characteristic Variables: Q_i, Qinv_i v_i
!----------------------------------------------------------------- 
do jj=-1,n+2
  if(dir==1)then
    dens=max(min_dens,u1(jj,ind,1))
    velx=u1(jj,ind,2)/u1(jj,ind,1)
    vely=u1(jj,ind,3)/u1(jj,ind,1)
    velz=u1(jj,ind,4)/u1(jj,ind,1)
    pres=max(min_pres,(gam-1.)*( u1(jj,ind,5)-0.5*dens &
      *(velx**2.+vely**2.+velz**2.) ) )
    a=sqrt(gam*pres/dens)
    H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a
  else 
    dens=max(min_dens,u1(ind,jj,1))
    velx=u1(ind,jj,3)/u1(ind,jj,1)
    vely=u1(ind,jj,2)/u1(ind,jj,1)
    velz=u1(ind,jj,4)/u1(ind,jj,1)
    pres=max(min_pres,(gam-1.)*( u1(ind,jj,5)-0.5*dens &
      *(velx**2.+vely**2.+velz**2.) ) )
    a=sqrt(gam*pres/dens)
    H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a
  end if

  !Right Eigenvectors
  Qx(jj,1,1)=1._8
  Qx(jj,2,1)=velx-a
  Qx(jj,3,1)=vely
  Qx(jj,4,1)=velz
  Qx(jj,5,1)=H-velx*a

  Qx(jj,1,2)=1._8
  Qx(jj,2,2)=velx
  Qx(jj,3,2)=vely
  Qx(jj,4,2)=velz
  Qx(jj,5,2)=0.5*(velx**2.+vely**2.+velz**2.)

  Qx(jj,1,3)=0._8
  Qx(jj,2,3)=0._8
  Qx(jj,3,3)=1._8
  Qx(jj,4,3)=0._8
  Qx(jj,5,3)=vely

  Qx(jj,1,4)=0._8
  Qx(jj,2,4)=0._8
  Qx(jj,3,4)=0._8
  Qx(jj,4,4)=1._8
  Qx(jj,5,4)=velz

  Qx(jj,1,5)=1._8
  Qx(jj,2,5)=velx+a
  Qx(jj,3,5)=vely
  Qx(jj,4,5)=velz
  Qx(jj,5,5)=H+velx*a

  !Left Eigenvectors
  Qinv(jj,1,1)=((gam-1.)/(2.*a*a))*(H &
              +a*(velx-a)/(gam-1.))
  Qinv(jj,1,2)=-((gam-1.)/(2.*a*a))*(velx+a/(gam-1.))
  Qinv(jj,1,3)=-((gam-1.)/(2.*a*a))*vely
  Qinv(jj,1,4)=-((gam-1.)/(2.*a*a))*velz
  Qinv(jj,1,5)=((gam-1.)/(2.*a*a))
  
  Qinv(jj,2,1)=((gam-1.)/(2*a*a))*(-2.*H+4.*a*a/(gam-1.))
  Qinv(jj,2,2)=((gam-1.)/(a*a))*velx
  Qinv(jj,2,3)=((gam-1.)/(a*a))*vely
  Qinv(jj,2,4)=((gam-1.)/(a*a))*velz
  Qinv(jj,2,5)=-((gam-1.)/(a*a))

  Qinv(jj,3,1)=-vely
  Qinv(jj,3,2)=0._8
  Qinv(jj,3,3)=1._8
  Qinv(jj,3,4)=0._8
  Qinv(jj,3,5)=0._8

  Qinv(jj,4,1)=-velz
  Qinv(jj,4,2)=0._8
  Qinv(jj,4,3)=0._8
  Qinv(jj,4,4)=1._8
  Qinv(jj,4,5)=0._8

  Qinv(jj,5,1)=((gam-1.)/(2.*a*a))*(H &
              -a*(velx+a)/(gam-1.))
  Qinv(jj,5,2)=((gam-1.)/(2.*a*a))*(-velx+a/(gam-1.))
  Qinv(jj,5,3)=-((gam-1.)/(2.*a*a))*vely
  Qinv(jj,5,4)=-((gam-1.)/(2.*a*a))*velz
  Qinv(jj,5,5)=((gam-1.)/(2.*a*a))
     
  do kk=1,5
    if(dir==1)then
      v(jj,kk)=Qinv(jj,kk,1)*u1(jj,ind,1)+Qinv(jj,kk,2)*u1(jj,ind,2) &
              +Qinv(jj,kk,3)*u1(jj,ind,3)+Qinv(jj,kk,4)*u1(jj,ind,4) &
              +Qinv(jj,kk,5)*u1(jj,ind,5)
    else
      v(jj,kk)=Qinv(jj,kk,1)*u1(ind,jj,1)+Qinv(jj,kk,2)*u1(ind,jj,2) &
              +Qinv(jj,kk,3)*u1(ind,jj,3)+Qinv(jj,kk,4)*u1(ind,jj,4) &
              +Qinv(jj,kk,5)*u1(ind,jj,5)  
    end if  
  end do   
end do

if(dir==1)then
 dr=dx
 dt=dtx
else
 dr=dy
 dt=dty
end if

!------------------------------------------------------------------
!Slope Limiter: S_j
!------------------------------------------------------------------  
do jj=0,n+1  
 do kk=1,5 
   S(jj,kk)=limiter(v(jj,kk)-v(jj-1,kk),v(jj+1,kk)-v(jj,kk))/dr	
 end do
end do

!------------------------------------------------------------------
!(dq/dx)_j
!------------------------------------------------------------------ 
do jj=0,n+1  
 do kk=1,5
   dqdx(jj,kk)=Qx(jj,kk,1)*S(jj,1)+Qx(jj,kk,2)*S(jj,2)+Qx(jj,kk,3)*S(jj,3) &
            +Qx(jj,kk,4)*S(jj,4)+Qx(jj,kk,5)*S(jj,5)
 end do
end do

!------------------------------------------------------------------
!(dq/dt)_j
!------------------------------------------------------------------ 
do jj=0,n+1  
  if(dir==1)then
    dens=max(min_dens,u1(jj,ind,1))
    velx=u1(jj,ind,2)/u1(jj,ind,1)
    vely=u1(jj,ind,3)/u1(jj,ind,1)
    velz=u1(jj,ind,4)/u1(jj,ind,1)
    pres=max(min_pres,(gam-1.)*( u1(jj,ind,5)-0.5*dens &
      *(velx**2.+vely**2.+velz**2.) ) )
    a=sqrt(gam*pres/dens)
    H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a
  else 
    dens=max(min_dens,u1(ind,jj,1))
    velx=u1(ind,jj,3)/u1(ind,jj,1)
    vely=u1(ind,jj,2)/u1(ind,jj,1)
    velz=u1(ind,jj,4)/u1(ind,jj,1)
    pres=max(min_pres,(gam-1.)*( u1(ind,jj,5)-0.5*dens &
      *(velx**2.+vely**2.+velz**2.) ) )
    a=sqrt(gam*pres/dens)
    H=0.5*(velx**2.+vely**2.+velz**2.)+gamul*a*a
  end if

   
  dqdt(jj,1)=-dqdx(jj,2)
  dqdt(jj,2)=-((gam-1.)*H-velx**2.-a**2.)*dqdx(jj,1) &
            -(3.-gam)*velx*dqdx(jj,2)+(gam-1.)*vely*dqdx(jj,3) &
            +(gam-1.)*velz*dqdx(jj,4)-(gam-1.)*dqdx(jj,5)
  dqdt(jj,3)=velx*vely*dqdx(jj,1)-vely*dqdx(jj,2)-velx*dqdx(jj,3)
  dqdt(jj,4)=velx*velz*dqdx(jj,1)-velz*dqdx(jj,2)-velx*dqdx(jj,4)
  dqdt(jj,5)=-0.5*velx*((gam-3.)*H-a**2.)*dqdx(jj,1) &
            -(H-(gam-1.)*velx**2.)*dqdx(jj,2) &
            +(gam-1.)*velx*vely*dqdx(jj,3) &
            +(gam-1.)*velx*velz*dqdx(jj,4)-gam*velx*dqdx(jj,5)

end do
 
!------------------------------------------------------------------
!^q_j
!------------------------------------------------------------------
do jj=0,n+1  
 if(dir==1)then
  do kk=1,5
   qh(jj,kk)=u1(jj,ind,kk)+0.5*(dt*dqdt(jj,kk)+dr*dqdx(jj,kk))
  end do
 else
  qh(jj,1)=u1(ind,jj,1)+0.5*(dt*dqdt(jj,1)+dr*dqdx(jj,1))
  qh(jj,2)=u1(ind,jj,3)+0.5*(dt*dqdt(jj,2)+dr*dqdx(jj,2))
  qh(jj,3)=u1(ind,jj,2)+0.5*(dt*dqdt(jj,3)+dr*dqdx(jj,3))
  qh(jj,4)=u1(ind,jj,4)+0.5*(dt*dqdt(jj,4)+dr*dqdx(jj,4))
  qh(jj,5)=u1(ind,jj,5)+0.5*(dt*dqdt(jj,5)+dr*dqdx(jj,5))
 end if 
end do
 
!---------------------------------------
!Solve Local Riemann Problem: w_j+1/2(0)
!---------------------------------------
do jj=0,n
 !set left and right states for local Riemann problem at j+1/2 cell interface
 rhoL=max(min_dens,qh(jj,1))
 uL=qh(jj,2)/qh(jj,1)
 vL=qh(jj,3)/qh(jj,1)
 wL=qh(jj,4)/qh(jj,1)
 pL=max(min_pres,(gam-1.)*( qh(jj,5)-0.5*rhoL &
    *(uL**2.+vL**2.+wL**2.) ) )
 aL=sqrt(gam*pL/rhoL)
 HL=0.5*(uL**2.+vL**2.+wL**2.)+gamul*aL*aL

 rhoR=max(min_dens,qh(jj+1,1))
 uR=qh(jj+1,2)/qh(jj+1,1)
 vR=qh(jj+1,3)/qh(jj+1,1)
 wR=qh(jj+1,4)/qh(jj+1,1)
 pR=max(min_pres,(gam-1.)*( qh(jj+1,5)-0.5*rhoR &
    *(uR**2.+vR**2.+wR**2.) ) )
 aR=sqrt(gam*pR/rhoR)
 HR=0.5*(uR**2.+vR**2.+wR**2.)+gamul*aR*aR

 call RoeSolver(jj)

 if(debug==1)then
  print*,'Max wave speed=',smax
 end if

end do


!-----------------------------------------------
!f_j
!-----------------------------------------------
do jj=0,nx+1
 u(1)=qh(jj,1)
 u(2)=qh(jj,2)/qh(jj,1)
 u(3)=qh(jj,3)/qh(jj,1)
 u(4)=qh(jj,4)/qh(jj,1)
 u(5)=(gam-1.)*(qh(jj,5)-0.5*u(1)*(u(2)**2.+u(3)**2.+u(4)**2.))
 call flux(u)
   do kk=1,5
     f(jj,kk)=fl(kk)
   end do    
end do

!------------------------------------------------
!^g_i+1/2
!------------------------------------------------
do jj=0,n 
 do kk=1,5
   g(jj,kk)=0.5*(f(jj,kk)+f(jj+1,kk))-0.5*&
        ( K1(jj,kk)*( Q(lambda(jj,1))*alpha(jj,1) ) &
         +K2(jj,kk)*( Q(lambda(jj,2))*alpha(jj,2) ) &
	 +K3(jj,kk)*( Q(lambda(jj,3))*alpha(jj,3) ) &
         +K4(jj,kk)*( Q(lambda(jj,4))*alpha(jj,4) ) &  
         +K5(jj,kk)*( Q(lambda(jj,5))*alpha(jj,5) ) ) 

 end do
end do
 
go to 1000

do jj=0,n
!check for transonic rarefaction waves
if(lambda1L(jj)<0. .and. lambda1R(jj)>0.)then
  do kk=1,5
    g(jj,kk)=f(jj,kk)+(lambda1L(jj)*(lambda1R(jj)&
       -lambda(jj,1))/(lambda1R(jj)-lambda1L(jj)))*alpha(jj,1)*K1(jj,kk) 
  end do
else if(lambda5L(jj)<0. .and. lambda5R(jj)>0.)then
  do kk=1,5
    g(jj,kk)=f(jj+1,kk)-(lambda5R(jj)*(lambda(jj,5)&
       -lambda5L(jj))/(lambda5R(jj)-lambda5L(jj)))*alpha(jj,5)*K5(jj,kk)
  end do
else
  do kk=1,5
    g(jj,kk)=0.5*(f(jj,kk)+f(jj+1,kk))&
               -0.5*( alpha(jj,1)*abs(lambda(jj,1))*K1(jj,kk)&
               +alpha(jj,2)*abs(lambda(jj,2))*K2(jj,kk)&
     	       +alpha(jj,3)*abs(lambda(jj,3))*K3(jj,kk)&
     	       +alpha(jj,4)*abs(lambda(jj,4))*K4(jj,kk)&
     	       +alpha(jj,5)*abs(lambda(jj,5))*K5(jj,kk) )   
  end do
end if

end do
1000 continue



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

end module RoeSolver2D_mod
