module HartenFC_mod
implicit none

integer,parameter::shockSwitch=1 !0:off, 1:on
integer,parameter::nt=100
integer,parameter::nx=200
real*8,parameter::del0=0.1

integer ::i,j,k
real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour

real*8::u1(-1:nx+2,3),u2(-1:nx+2,3) !fluid state vector q=[rho,rho*u,rho*e_T]
real*8::gt(-1:nx+2,3),fc(-1:nx+2,3),gamma(-1:nx+2,3),sigma(-1:nx+2,3)
real*8::f(-1:nx+2,3),g(-1:nx+2,3),fl(3)
real*8::u(3)!u=[rho,u,p]
real*8::rhoL,uL,pL,aL,hL   !left state
real*8::rhoR,uR,pR,aR,hR   !right state
real*8::rhoRL,uRL,pRL,aRL,hRL !Roe-averaged quantities
real*8::dens,vel,pres
real*8::drho,du,dp,smax
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu
real*8::lambda(-1:nx+2,3),delv(-1:nx+2,3) 
real*8::R1(-1:nx+2,3),R2(-1:nx+2,3),R3(-1:nx+2,3)

real*8 ::temp1,temp2,temp3,xc,xt,theta,w

real*8,parameter ::R=287. !gas constant in SI units
real*8::del!delta parameter for Harten's first order method

contains

 subroutine init()
!------------------------------------------------------
!initialization
!------------------------------------------------------

open(unit=10,file='input.txt')
open(unit=11,file='output_euler_harten_FC.txt')
open(unit=12,file='dt.txt')

!Discountinuity in initial condition at the middle of [xmin,xmax]

xmin=0.0
xmax=1.0
cour=0.9
dx=(xmax-xmin)/nx
dt=0.

gam=1.4
gamil=1./(gam+1.)
gamul=1./(gam-1.)
gamel=(gam-1.)/(gam+1.)
gamee=(gam-1.)/(2.*gam)
gamuu=(gam+1.)/(2.*gam)

smax=0.

read(unit=10,fmt=*) rhoL,uL,pL 
read(unit=10,fmt=*) rhoR,uR,pR

do j=-1,nx+2
  if(j<nx/2) then
    u1(j,1)=rhoL
    u1(j,2)=rhoL*uL
    u1(j,3)=0.5*rhoL*uL*uL+pL*gamul
  else
    u1(j,1)=rhoR
    u1(j,2)=rhoR*uR
    u1(j,3)=0.5*rhoR*uR*uR+pR*gamul
  end if
end do

print*,'Left State (rho,u,p) =',rhoL,uL,pL
print*,'Right State (rho,u,p) =',rhoR,uR,pR

end subroutine init



 subroutine computeFlux()
 
 
 do j=-1,nx+1
  !-----------------------------------------------------------------
  !R_j+1/2, lambda_j+1/2, delv_j+1/2
  !obtain by solving linearized Riemann problem on right cell edge
  !-----------------------------------------------------------------
   rhoL=u1(j,1)
   uL=u1(j,2)/u1(j,1)
   pL=max((gam-1.)*(u1(i,3)-0.5*u1(i,2)*u1(i,2)/u1(i,1)),0.)
   !pL=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))  
   rhoR=u1(j+1,1)
   uR=u1(j+1,2)/u1(j+1,1)
   pR=max((gam-1.)*(u1(j+1,3)-0.5*u1(j+1,2)*u1(j+1,2)/u1(j+1,1)),0.)
   
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
   !lambda_j+1/2
   lambda(j,1)=uRL !entropy wave/contact discontinuity propagation velocity
   lambda(j,2)=uRL+aRL
   lambda(j,3)=uRL-aRL

   !Compute wave strengths for characteristic variables
   !delv_j+1/2
   delv(j,1)=drho-dp/(aRL*aRL)
   delv(j,2)=du+dp/(aRL*rhoRL)
   delv(j,3)=du-dp/(aRL*rhoRL)
  
   !Costruct right-eigenvectors of Roe-averaged Jacobian
   !R1_j+1/2
   R1(j,1)=1.0_8
   R1(j,2)=uRL
   R1(j,3)=0.5*uRL*uRL
   !R2_j+1/2
   R2(j,1)=rhoRL/(2.0*aRL)
   R2(j,2)=(rhoRL/(2.0*aRL))*(uRL+aRL)
   R2(j,3)=(rhoRL/(2.0*aRL))*(hRL+uRL*aRL)
   !R3_j+1/2   
   R3(j,1)=-1.0*(rhoRL/(2.0*aRL))
   R3(j,2)=-1.0*(rhoRL/(2.0*aRL))*(uRL-aRL)
   R3(j,3)=-1.0*(rhoRL/(2.0*aRL))*(hRL-uRL*aRL)
  
 end do
   
 
 do j=-1,nx+1
   !-----------------------------------------------
   !Max Wave Speed
   !-----------------------------------------------
   smax=max(smax,abs(lambda(j,1)),abs(lambda(j,2)),abs(lambda(j,3)))
 end do

   
 do j=-1,nx+1 
  !------------------------------------------------
  !~g_j+1/2
  !------------------------------------------------
   do k=1,3
     gt(j,k)=0.5*(Q(lambda(j,k))-(dt/dx)*(lambda(j,k)*lambda(j,k)))*delv(j,k) 
     if(shockSwitch==1)then
       sigma(j,k)=0.5*(1.-Q(lambda(j,k)))
     end if	 
   end do  
 end do
   
 do j=0,nx+1 
  !-----------------------------------------------
  !f(c)_j
  !-----------------------------------------------
   do k=1,3
     fc(j,k)=minmod(gt(j,k),gt(j-1,k))
   end do  
   
   if(shockSwitch==1)then
    do k=1,3
     if(abs(delv(j,k))+abs(delv(j-1,k)) .ne. 0._8)then
      theta=(abs(delv(j,k)-delv(j-1,k)))/(abs(delv(j,k))+abs(delv(j-1,k)))
     else
      theta=0.
     end if 
     w=minmod(delv(j,k)*sigma(j,k),delv(j-1,k)*sigma(j-1,k))
     
     !Modify the entropy wave only (will result in steeper contact discontinuity)
     if(k==2)then	 
       fc(j,k)=fc(j,k)+w*theta
     end if

    end do
   end if
   
  !-----------------------------------------------
  !f_j
  !-----------------------------------------------
   u(1)=u1(j,1)
   u(2)=u1(j,2)/u1(j,1)
   u(3)=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
   call flux(u)
   do k=1,3
     f(j,k)=fl(k)
   end do
 end do

 do j=0,nx  
  !------------------------------------------------
  !gamma_j+1/2
  !------------------------------------------------
  do k=1,3
    if(delv(j,k)>=1.d-30)then
      gamma(j,k)=(fc(j+1,k)-fc(j,k))/delv(j,k)
	else
      gamma(j,k)=0.
    end if	  
  end do
  
  !------------------------------------------------
  !^g_j+1/2
  !------------------------------------------------
  do k=1,3
    g(j,k)=0.5*(f(j,k)+f(j+1,k)&
	      +R1(j,k)*(fc(j,1)+fc(j+1,1)-Q(lambda(j,1)+gamma(j,1))*delv(j,1))&
		  +R2(j,k)*(fc(j,2)+fc(j+1,2)-Q(lambda(j,2)+gamma(j,2))*delv(j,2))&
		  +R3(j,k)*(fc(j,3)+fc(j+1,3)-Q(lambda(j,3)+gamma(j,3))*delv(j,3)) )

  end do
 
 end do
 
 
 end

 
  subroutine flux(u)
    real*8::u(3)
    fl(1)=u(1)*u(2)
    fl(2)=u(1)*u(2)*u(2)+u(3)
    fl(3)=u(1)*u(2)*(gam*gamil*(u(3)/u(1))+0.5*u(2)*u(2))
  end
 
  function Q(x) result(qx)
  real*8,intent(in)::x
  real*8::qx
 
  if(abs(x)<del0)then
    qx=(x*x+del0*del0)/(2.*del0)
  else if(abs(x)>=del0)then
    qx=abs(x)
  end if
 

 end
 
 function minmod(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z
  
  z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
 end
 
subroutine bound()

do k=1,3
    u2(0,k)=u2(1,k)	
    u2(-1,k)=u2(1,k)
    u2(nx+1,k)=u2(nx,k)
    u2(nx+2,k)=u2(nx,k)
end do
 
end subroutine bound

end module HartenFC_mod
