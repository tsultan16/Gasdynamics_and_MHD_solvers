module RiemannHarten_mod
implicit none

real*8::u(3),uleft(3),uright(3) !fluid state vector u=[rho,u,p]
real*8::f(3) !flux vector
real*8::a,h!sound speed,specific enthalpy
real*8::rhoL,uL,pL,aL,hL,sL,ML   !left state
real*8::rhoR,uR,pR,aR,hR,sR,MR   !right state
real*8::rhoRL,uRL,pRL,aRL,hRL !Roe-averaged quantities
real*8::drho,du,dp,gamil 
real*8::lambda(3),dv(3) 
real*8::r1(3),r2(3),r3(3)
real*8::fL(3),fR(3)

real*8 ::temp1,temp2,temp3,xc,xt


real*8,parameter ::gam=5./3. !gam=cp/cv
real*8,parameter ::R=287. !gas constant in SI units
real*8::del!delta parameter for Harten's first order method

contains
 subroutine computeFlux()
 
  !sound speeds for left and right states
  aR = sqrt((gam*pR)/rhoR) 
  aL = sqrt((gam*pL)/rhoL)
   
  uleft=(/rhoL,uL,pL/)
  uright=(/rhoR,uR,pR/)

  !sound speeds for left and right states
  aR = sqrt((gam*pR)/rhoR) 
  aL = sqrt((gam*pL)/rhoL)

  !specific enthalpy for left and right states
  hR=gamil*aR*aR+0.5*uR*uR
  hL=gamil*aL*aL+0.5*uL*uL

  drho=rhoR-rhoL
  du=uR-uL
  dp=pR-pL

  !Compute Roe averages
  rhoRL=sqrt(rhoR*rhoL)
  uRL=(sqrt(rhoR)*uR+sqrt(rhoL)*uL)/(sqrt(rhoR)+sqrt(rhoL))
  hRL=(sqrt(rhoR)*hR+sqrt(rhoL)*hL)/(sqrt(rhoR)+sqrt(rhoL))
  aRL=sqrt((gam-1.)*(hRL-0.5*uRL*uRL))


  !Compute eigenvalues of Roe-Averaged Jacobian Matrix(i.e.characteristic wave speeds)
  lambda(1)=uRL
  lambda(2)=uRL+aRL
  lambda(3)=uRL-aRL

  !Compute wave strengths for characteristic variables
  dv(1)=drho-dp/(aRL*aRL)
  dv(2)=du+dp/(aRL*rhoRL)
  dv(3)=du-dp/(aRL*rhoRL)

  !Costruct right-eigenvectors of Roe-averaged Jacobian
  r1=(/1.0_8,uRL,0.5*uRL*uRL/)
  r2=(rhoRL/(2.0*aRL))*(/1.0_8,uRL+aRL,hRL+uRL*aRL/)
  r3=-1.0*(rhoRL/(2.0*aRL))*(/1.0_8,uRL-aRL,hRL-uRL*aRL/)

  
  !compute flux: f(u(0))
  call flux(uleft)
  fL=f
  call flux(uright)
  fR=f
  f=0.5*(fL+fR)-0.5*(psi(lambda(1))*dv(1)*r1&
    +psi(lambda(2))*dv(2)*r2+psi(lambda(3))*dv(3)*r3)

 end  
 
 subroutine flux(u)
  real*8::u(3)
  f(1)=u(1)*u(2)
  f(2)=u(1)*u(2)*u(2)+u(3)
  f(3)=u(1)*u(2)*(gam*gamil*u(3)/u(1)+0.5*u(2)*u(2))
 end 
 
 function psi(x) result(fx)
  real*8,intent(in)::x
  real*8::fx
  if(abs(x)<del)then
    fx=(x*x+del*del)/(2.*del)
  else if(abs(x)>=del)then
    fx=abs(x)
  end if  
 end
 
end module RiemannHarten_mod
