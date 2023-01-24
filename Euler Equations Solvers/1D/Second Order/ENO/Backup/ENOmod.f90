module ENO_mod
implicit none

integer,parameter::nt=100
integer,parameter::nx=400
real*8,parameter::del0=0.2
real*8,parameter::cour=0.9 !CFL number 
real*8::dt,dx
real*8::xmin,xmax,x,pressure,Mo
integer::i,k,t
real*8::q1(-1:nx+2,3),q2(-1:nx+2,3),qh(-1:nx+2,3) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::v(-1:nx+2,3) !characteristic variables
real*8::Q(-1:nx+2,3,3),Qinv(-1:nx+2,3,3) !characeteristic matrix
real*8::S(-1:nx+2,3)
real*8::dqdx(-1:nx+2,3),dqdt(-1:nx+2,3)
real*8::gt(-1:nx+2,3),fc(-1:nx+2,3),gamma(-1:nx+2,3),sigma(-1:nx+2,3)
real*8::f(-1:nx+2,3),g(-1:nx+2,3),fl(3)
real*8::u(3)!u=[rho,u,p]
real*8::rhoL,uL,pL,aL,hL   !left state
real*8::rhoR,uR,pR,aR,hR   !right state
real*8::rhoRL,uRL,pRL,aRL,hRL !Roe-averaged quantities
real*8::drho,du,dp,gamil 
real*8::lambda(-1:nx+2,3),delv(-1:nx+2,3) 
real*8::R1(-1:nx+2,3),R2(-1:nx+2,3),R3(-1:nx+2,3)

real*8 ::temp1,temp2,temp3,xc,xt,theta,w,pres,a,vel,dens
real*8::dpdx,dudx,dpdt,dudt

real*8,parameter::premin=1.d-10
real*8,parameter::densmin=1.d-10
real*8,parameter ::gam=5./3. !gam=cp/cv
real*8,parameter ::R=287. !gas constant in SI units
real*8::del!delta parameter for Harten's first order method

contains

subroutine init()

open(unit=10,file='input.txt')
open(unit=11,file='output_euler_ENO.txt')
open(unit=12,file='q1_qh.txt')

!initialization
xmin=0.0
xmax=1.0
gamil=1./(gam-1.)
read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR
dx=(xmax-xmin)/(nx*1.)
dt=0.


do i=-1,nx+2
  x=xmin+(i-0.5)*dx
  if(i<nx/2) then
   q1(i,1)=rhoL
   q1(i,2)=rhoL*uL
   q1(i,3)=0.5*rhoL*uL*uL+gamil*pL
   pressure=pL
  else
   q1(i,1)=rhoR
   q1(i,2)=rhoR*uR
   q1(i,3)=0.5*rhoR*uR*uR+gamil*pR
   pressure=pR
  end if
  if(i>=1 .and. i<=nx) then  
    write(11,*) x,q1(i,1),q1(i,2)/q1(i,1),pressure
  end if    
end do

print*,'Left State=',rhoL,uL,pL
print*,'Right State=',rhoR,uR,pR


end subroutine init


subroutine computeFlux()
 
 integer::i,j,k
 
 do i=-1,nx+2
  !-----------------------------------------------------------------
  !Characteristic Variables: Q_i, Qinv_i v_i
  !-----------------------------------------------------------------
   pres=max((gam-1.)*(q1(i,3)-0.5*q1(i,2)*q1(i,2)/q1(i,1)),premin)
   a=sqrt(gam*pres/q1(i,1))
   vel=q1(i,2)/q1(i,1)
   dens=max(q1(i,1),densmin)
   
   Q(i,1,1)=1
   Q(i,1,2)=0.5*dens/a
   Q(i,1,3)=-0.5*dens/a
   Q(i,2,1)=vel
   Q(i,2,2)=(0.5*dens/a)*(vel+a)
   Q(i,2,3)=-(0.5*dens/a)*(vel-a)
   Q(i,3,1)=0.5*vel*vel
   Q(i,3,2)=(0.5*dens/a)*(0.5*vel*vel+(a*a/(gam-1.))+a*vel)
   Q(i,3,3)=-(0.5*dens/a)*(0.5*vel*vel+(a*a/(gam-1.))-a*vel)
   
   Qinv(i,1,1)=((gam-1.)/(a*a))*(-0.5*vel*vel+a*a/(gam-1.))
   Qinv(i,1,2)=(gam-1.)*vel/(a*a)
   Qinv(i,1,3)=-(gam-1.)/(a*a)
   Qinv(i,2,1)=((gam-1.)/(dens*a))*(0.5*vel*vel-a*vel/(gam-1.))
   Qinv(i,2,2)=((gam-1.)/(dens*a))*(-vel+a/(gam-1.))
   Qinv(i,2,3)=(gam-1.)/(dens*a)
   Qinv(i,3,1)=((gam-1.)/(dens*a))*(-0.5*vel*vel+a*vel/(gam-1.))
   Qinv(i,3,2)=((gam-1.)/(dens*a))*(vel+a/(gam-1.))
   Qinv(i,3,3)=-(gam-1.)/(dens*a)
   
   do j=1,3
     v(i,j)=Qinv(i,j,1)*q1(i,1)+Qinv(i,j,2)*q1(i,2)+Qinv(i,j,3)*q1(i,3)
   end do   
 end do
  
 do i=0,nx+1  
  !------------------------------------------------------------------
  !Slope Limiter: S_i
  !------------------------------------------------------------------
  do j=1,3 
   !if(abs(v(i+1,j)-v(i,j))<abs(v(i,j)-v(i-1,j)))then
     S(i,j)=minmod(v(i+1,j)-v(i,j),v(i,j)-v(i-1,j))/dx	
   !else	
   !  S(i,j)=(v(i,j)-v(i-1,j))/dx
   !end if
  end do
      
 end do
 
 do i=0,nx+1  
  !------------------------------------------------------------------
  !(dq/dx)_i
  !------------------------------------------------------------------
  do j=1,3
    dqdx(i,j)=Q(i,j,1)*S(i,1)+Q(i,j,2)*S(i,2)+Q(i,j,3)*S(i,3)
  end do
 end do
 
 do i=0,nx+1  
  !------------------------------------------------------------------
  !(dq/dt)_i
  !------------------------------------------------------------------
  pres=max((gam-1.)*(q1(i,3)-0.5*q1(i,2)*q1(i,2)/q1(i,1)),premin)
  a=sqrt(gam*pres/q1(i,1))
  vel=q1(i,2)/q1(i,1)
  dens=max(q1(i,1),densmin)
   
  dudx=(1./dens)*(dqdx(i,2)-vel*dqdx(i,1))
  dpdx=(gam-1.)*(dqdx(i,3)-0.5*(q1(i,2)*dudx+vel*dqdx(i,2)))
  
  dqdt(i,1)=-dqdx(i,2)
  dqdt(i,2)=-vel*(dqdx(i,2)+dens*dudx)-dpdx
           
  dpdt=-vel*dpdx-gam*pres*dudx
  dudt=(1./dens)*(dqdt(i,2)-vel*dqdt(i,1))
  
  dqdt(i,3)=(1./(gam-1.))*dpdt+0.5*(q1(i,2)*dudt+vel*dqdt(i,2))
 end do
 
 do i=0,nx+1  
  !------------------------------------------------------------------
  !^q_i
  !------------------------------------------------------------------
  do j=1,3
    qh(i,j)=q1(i,j)+0.5*(dt*dqdt(i,j)+dx*dqdx(i,j))
  end do 
 end do
 
 
 do i=0,nx
  !-----------------------------------------------------------------
  !R_i+1/2, lambda_i+1/2, delv_i+1/2
  !obtain by solving linearized Riemann problem on right cell edge
  !-----------------------------------------------------------------
   rhoL=max(qh(i,1),densmin)
   uL=qh(i,2)/qh(i,1)
   pL=max((gam-1.)*(qh(i,3)-0.5*qh(i,2)*qh(i,2)/qh(i,1)),premin)  
   rhoR=max(qh(i+1,1),densmin)
   uR=qh(i+1,2)/qh(i+1,1)
   pR=max((gam-1.)*(qh(i+1,3)-0.5*qh(i+1,2)*qh(i+1,2)/qh(i+1,1)),premin)
   
   !sound speeds for left and right states
   aR = sqrt((gam*pR)/rhoR) 
   aL = sqrt((gam*pL)/rhoL)
 
   !specific enthalpy for left and right states
   hR=gamil*aR*aR+0.5*uR*uR
   hL=gamil*aL*aL+0.5*uL*uL
   
   !Compute Roe averages
   rhoRL=sqrt(rhoR*rhoL)
   uRL=(sqrt(rhoR)*uR+sqrt(rhoL)*uL)/(sqrt(rhoR)+sqrt(rhoL))
   hRL=(sqrt(rhoR)*hR+sqrt(rhoL)*hL)/(sqrt(rhoR)+sqrt(rhoL))
   aRL=sqrt(max((gam-1.)*(hRL-0.5*uRL*uRL),0.))
   drho=rhoR-rhoL
   du=uR-uL
   dp=pR-pL

   !Compute eigenvalues of Roe-Averaged Jacobian Matrix(i.e.characteristic wave speeds)
   !lambda_i+1/2
   lambda(i,1)=uRL !entropy wave/contact discontinuity propagation velocity
   lambda(i,2)=uRL+aRL
   lambda(i,3)=uRL-aRL

   !Compute wave strengths for characteristic variables
   !delv_i+1/2
   delv(i,1)=drho-dp/(aRL*aRL)
   delv(i,2)=du+dp/(aRL*rhoRL)
   delv(i,3)=du-dp/(aRL*rhoRL)
  
   !Costruct right-eigenvectors of Roe-averaged Jacobian
   !R1_i+1/2
   R1(i,1)=1.0_8
   R1(i,2)=uRL
   R1(i,3)=0.5*uRL*uRL
   !R2_i+1/2
   R2(i,1)=rhoRL/(2.0*aRL)
   R2(i,2)=(rhoRL/(2.0*aRL))*(uRL+aRL)
   R2(i,3)=(rhoRL/(2.0*aRL))*(hRL+uRL*aRL)
   !R3_i+1/2   
   R3(i,1)=-1.0*(rhoRL/(2.0*aRL))
   R3(i,2)=-1.0*(rhoRL/(2.0*aRL))*(uRL-aRL)
   R3(i,3)=-1.0*(rhoRL/(2.0*aRL))*(hRL-uRL*aRL)
  
 end do
   
 dt=0.  
 do i=1,nx-1
   !-----------------------------------------------
   !dt
   !-----------------------------------------------
   dt=max(dt,abs(lambda(i,1)),abs(lambda(i,2)),abs(lambda(i,3)))
 end do
 
 dt=(dx*cour)/dt
 print*,'dt=',dt
   
 do i=0,nx+1
  !-----------------------------------------------
  !f_i
  !-----------------------------------------------
   u(1)=qh(i,1)
   u(2)=qh(i,2)/qh(i,1)
   u(3)=(gam-1.)*(qh(i,3)-0.5*qh(i,2)*qh(i,2)/qh(i,1))
   call flux(u)
   do j=1,3
     f(i,j)=fl(j)
   end do
 end do


 do i=0,nx 
  !------------------------------------------------
  !^g_i+1/2
  !------------------------------------------------
  do j=1,3
	g(i,j)=0.5*(f(i,j)+f(i+1,j)&
	      -R1(i,j)*(psi(lambda(i,1))*delv(i,1))&
		  -R2(i,j)*(psi(lambda(i,2))*delv(i,2))&
		  -R3(i,j)*(psi(lambda(i,3))*delv(i,3)))

  end do
 
 end do
 
 
end subroutine computeFlux

 
subroutine flux(u)
  real*8::u(3)
  fl(1)=u(1)*u(2)
  fl(2)=u(1)*u(2)*u(2)+u(3)
  fl(3)=u(1)*u(2)*(gam*gamil*(u(3)/u(1))+0.5*u(2)*u(2))
end
 
function psi(x) result(qx)

  real*8,intent(in)::x
  real*8::qx
 
  if(abs(x)<del0)then
    qx=(x*x+del0*del0)/(2.*del0)
  else 
    qx=abs(x)
  end if

end function psi

subroutine bound()

!enforce open boundary conditions
  do k=1,3
    q2(0,k)=q2(1,k)	
    q2(-1,k)=q2(1,k)
    q2(nx+1,k)=q2(nx,k)
	q2(nx+2,k)=q2(nx,k)
  end do

end subroutine bound
 


function minmod(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z
  
  z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod
 
 
end module ENO_mod
