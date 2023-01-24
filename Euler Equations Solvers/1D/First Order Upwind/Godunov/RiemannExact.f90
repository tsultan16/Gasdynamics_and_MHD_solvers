module RiemannExact_mod
implicit none

real*8 ::rho,u,p,a,s,M !density,velocity,pressure,sound speed,entropy and mach number
real*8 ::f(3) !flux vector 
real*8 ::rhoL,uL,pL,aL,sL,ML   !left state
real*8 ::rhoR,uR,pR,aR,sR,MR   !right state

real*8 ::rho2,u2,p2,a2 !intermediate states
real*8 ::rho3,u3,p3,a3

real*8 ::Sh !shock speed
real*8 ::temp1,temp2,temp3

real*8,parameter ::gam=5./3. !gam=cp/cv
real*8,parameter ::R=287. !gas constant in SI units

contains
 subroutine computeFlux()
 
  !sound speeds for left and right states
  aR = sqrt((gam*pR)/rhoR) 
  aL = sqrt((gam*pL)/rhoL)
  !entropy for left and right states
  sR=(R/(gam-1))*log(pR/(rhoR**gam))
  sL=(R/(gam-1))*log(pL/(rhoL**gam))
  !Mach number for left and right states
  MR=uR/aR
  ML=uL/aL

  call findroot()
  
  !calculate state variables in region 2
  temp1=(gam+1)/(gam-1)
  temp2=p2/pR
  temp3=((gam+1)/(2*gam))*(temp2-1)+1 

  a2=aR*sqrt(temp2*((temp1+temp2)/(1+temp1*temp2)))
  u2=uR+(aR/gam)*(temp2-1)/sqrt(temp3)
  Sh=uR+aR*sqrt(temp3) 
  rho2=(gam*p2)/(a2**2)

  !calculate state variables in region 3
  u3=u2
  p3=p2
  a3=aL*(p3/pL)**((gam-1)/(2*gam))
  rho3=(gam*p3)/(a3**2)

  !compute flux
  call solution()
  f(1)=rho*u !mass flux flux
  f(2)=(rho*u*u)+p !mometum flux
  f(3)=u*(0.5*rho*u*u+p*(gam/(gam-1.))) !total energy flux

  !print*,'left state=',rhoL,uL,pL
  !print*,'right state=',rhoR,uR,pR  
  !print*,'f=',f
  
 end 
 
 
 subroutine findroot() !secant method root solver
   !use pL and pR as initial guesses for p2
   integer ::Niter = 30
   integer ::i, iterCount
 
   real*8 ::p0,p1 !initial guesses
   
   iterCount=0
   
   p0=pL*(1.+0.5)
   p1=pR*(1.+0.1)   
   
   do i=1,Niter
     if(abs((p1-p0)/p0)<0.0000001) then  !stop interations after convergernce criterion satisfied
	   exit
	 else if(abs((p1-p0)/p0)>0.0000001 .and. i<Niter-1) then
	   Niter=Niter+10 !add more iterations if solution does not converge
	   iterCount=iterCount+1
	 else if(abs((p1-p0)/p0)>0.0000001 .and. i<Niter-1 .and. iterCount >5) then
	   print*,'findRoot() unable to converge to the root...'
	   exit
	 end if
     p2=p1-(pf(p1)*(p1-p0)/(pf(p1)-pf(p0)))
	 p0=p1
	 p1=p2
	 
	 !print*,'Iteration# = ',i,', Root = ',p2
   end do
 
 end 

 subroutine solution()  
   real*8 ::b1,b2 !expansion region boundaries
   real*8 ::bc    !location of contact discontinuity
   real*8 ::bs    !location of schock discontinuity
   real*8 ::temp1,xt
   
   xt=0.0  !xt = x/t  
   temp1=(xt)+((gam-1)/2)*uL+aL;
   
   !speeds of expansion fan head and tail(need to be careful here)
   b1=uL-aL
   b2=u3-a3
   !speed of contact
   bc=aR  
   !shock speed
   bs=Sh 
	 
   if(xt<b1) then !Region 4
	 u=uL
	 a=aL
     rho=rhoL
	 p=pL
	 s=sL
	 M=ML
   else if(xt>b1 .and. xt<b2) then !Expansion Region
     u=(2./(gam+1))*temp1
	 a=u-xt
	 p=pL*(a/aL)**((2*gam)/(gam-1))
	 rho=(gam*p)/(a**2)
	 s=(R/(gam-1))*log(p/(rho**gam))
	 M=u/a
   else if(xt>b2 .and. xt<bc) then !Region 3
     u=u3
     a=a3
     rho=rho3	 
   	 p=p3
	 s=(R/(gam-1))*log(p3/(rho3**gam))
	 M=u3/a3
   else if(xt>bc .and. xt<bs) then !Region 2
     u=u2
     a=a2
     rho=rho2	 
   	 p=p2
	 s=(R/(gam-1))*log(p2/(rho2**gam))
	 M=u2/a2
   else if(xt>bs) then !Region 1
     u=uR
	 a=aR
	 rho=rhoR
	 p=pR
	 s=sR
	 M=MR
   end if
 
 end subroutine solution
 
 function pf(x) result(px)
 real*8,intent(in) ::x
 real*8 ::px,temp1,temp2
    
  temp1=uL-uR-(aR/gam)*(((x/pR)-1.)/sqrt((((gam+1)*(x/pR)-1.)/(2*gam))+1)) 
  temp2=1+((gam-1)/(2*aL))*temp1
  px=x-pL*(temp2)**((2*gam)/(gam-1))

 end function pf



end module RiemannExact_mod