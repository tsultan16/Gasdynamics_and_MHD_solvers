!Passive Variable Conservation Equation
!Pure Advcetion: q_t +(v q)_x=0
!Spatial domain: x in [xmin,xmax].
!Cell-interfaces values are approximated using
!the 5th order WENO scheme.
!TVD 3rd Order and 4th Order Runge-Kutta time-integrators are available for use.
!Reference:Jiang & Wu, J Comp. Phys. 150 (1999)

module weno_mod
use global_vars_mod
implicit none

contains
!*********************************************************************
!Subroutine and Function Definitions
!*********************************************************************

subroutine weno_flux(q,flux)
!I/O variables
real::q(-3:nx+4,-3:nx+4),flux(0:nx,0:ny)
!Local variables
real::F(-3:nx+4,-3:ny+4,2),F_plus(-3:nx+4,-3:ny+4,2),F_minus(-3:nx+4,-3:ny+4,2)
real::alphax,alphay
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i,j


!Compute cell-center physical flux and Lax-Freidrichs split flux: F_i,j,F^(+)_i,j,F^(-)_i,j
do i=-3,nx+4
  do j=-3,ny+4
    F(i,j,1)=q(i,j)*v(i,j,1)
    F(i,j,2)=q(i,j)*v(i,j,2)
    
    !local max wave speed
    alphax=max(abs(v(i-1,j,1)),abs(v(i,j,1)),abs(v(i+1,j,1)))
    alphay=max(abs(v(i,j-1,2)),abs(v(i,j,2)),abs(v(i,j+1,2))) 
 
    F_plus(i,j,1)=0.5*(F(i,j,1)+alphax*q(i,j))
    F_minus(i,j,1)=0.5*(F(i,j,1)-alphax*q(i,j))
    F_plus(i,j,2)=0.5*(F(i,j,2)+alphay*q(i,j))
    F_minus(i,j,2)=0.5*(F(i,j,2)-alphay*q(i,j))
  end do
end do


!Compute WENO approximation of cell-interface fluxes
do i=0,nx
  do j=0,ny
    

    !WENO approximation of cell interface flux

    !*********
    !F_i+1/2,j
    !*********
    !Phi function arguments
    a1=F_plus(i-1,j,1)-F_plus(i-2,j,1)
    b1=F_plus(i,j,1)-F_plus(i-1,j,1)
    c1=F_plus(i+1,j,1)-F_plus(i,j,1)
    d1=F_plus(i+2,j,1)-F_plus(i+1,j,1)

    a2=F_minus(i+3,j,1)-F_minus(i+2,j,1)
    b2=F_minus(i+2,j,1)-F_minus(i+1,j,1)
    c2=F_minus(i+1,j,1)-F_minus(i,j,1)
    d2=F_minus(i,j,1)-F_minus(i-1,j,1) 

    flux(i,j,1)=(1./12.)*(-F(i-1,j,1)+7.*F(i,j,1)+7.*F(i+1,j,1)-F(i+2,j,1))-phi(a1,b1,c1,d1)&
         +phi(a2,b2,c2,d2)

    !*********
    !F_i,j+1/2
    !*********
    !Phi function arguments
    a1=F_plus(i,j-1,2)-F_plus(i,j-2,2)
    b1=F_plus(i,j,2)-F_plus(i,j-1,2)
    c1=F_plus(i,j+1,2)-F_plus(i,j,2)
    d1=F_plus(i,j+2,2)-F_plus(i,j+1,2)

    a2=F_minus(i,j+3,2)-F_minus(i,j+2,2)
    b2=F_minus(i,j+2,2)-F_minus(i,j+1,2)
    c2=F_minus(i,j+1,2)-F_minus(i,j,2)
    d2=F_minus(i,j,2)-F_minus(i,j-1,2) 

    flux(i,j,2)=(1./12.)*(-F(i,j-1,2)+7.*F(i,j,2)+7.*F(i,j+1,2)-F(i,j+2,2))-phi(a1,b1,c1,d1)&
         +phi(a2,b2,c2,d2)


  end do
end do

end subroutine weno_flux


function phi(a,b,c,d) result(phix)
!I/0 variables
real,intent(in)::a,b,c,d
real::phix
!Local Variables
real::w0,w2,alpha0,alpha1,alpha2,alphasum,IS0,IS1,IS2

IS0=13.*(a-b)**2+3.*(a-3.*b)**2
IS1=13.*(b-c)**2+3.*(b+c)**2
IS2=13.*(c-d)**2+3.*(3.*c-d)**2

alpha0=1./((eps+IS0)**2)
alpha1=6./((eps+IS1)**2)
alpha2=3./((eps+IS2)**2)

alphasum=alpha0+alpha1+alpha2

w0=alpha0/alphasum
w2=alpha2/alphasum

phix=(1./3.)*w0*(a-2.*b+c)+(1./6.)*(w2-0.5)*(b-2.*c+d)

end function phi

function v(i,j,k) result(vx)
!I/0 variables
integer,intent(in)::i,j,k
real::vx
!Local variables

!Uniform advection to the left
if(k==1)then  !x-component
vx=1.
else if(k==2)then  !y-component
vx=0.
end if

end function v

end module weno_mod
