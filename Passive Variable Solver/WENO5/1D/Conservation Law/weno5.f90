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
real::q(-3:nx+4),flux(0:nx)
!Local variables
real::F(-3:nx+4),F_plus(-3:nx+4),F_minus(-3:nx+4)
real::alpha
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i


!Compute cell-center physical flux and Lax-Freidrichs split flux: F_i,F^(+)_i,F^(-)_i
do i=-3,nx+4
  F(i)=q(i)*v(i)

 !local max wave speed
  alpha=max(abs(v(i-1)),abs(v(i)),abs(v(i+1)))
  
  F_plus(i)=0.5*(F(i)+alpha*q(i))
  F_minus(i)=0.5*(F(i)-alpha*q(i))
end do

do i=0,nx
  !Phi function arguments
  a1=F_plus(i-1)-F_plus(i-2)
  b1=F_plus(i)-F_plus(i-1)
  c1=F_plus(i+1)-F_plus(i)
  d1=F_plus(i+2)-F_plus(i+1)

  a2=F_minus(i+3)-F_minus(i+2)
  b2=F_minus(i+2)-F_minus(i+1)
  c2=F_minus(i+1)-F_minus(i)
  d2=F_minus(i)-F_minus(i-1)

  !WENO approximation of cell interface flux
  !F_i+1/2
  flux(i)=(1./12.)*(-F(i-1)+7.*F(i)+7.*F(i+1)-F(i+2))-phi(a1,b1,c1,d1)&
         +phi(a2,b2,c2,d2)
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

function v(i) result(vx)
!I/0 variables
integer,intent(in)::i
real::vx
!Local variables
real::x

x=xmin+(i-1)*dx

!Uniform advection to the left
!vx=-1.

!go to 101
!Discontinuous velocity field (stationary shock)
if(x<0.5)then
  vx=1.0
else 
  vx=0.25
end if
!101 continue


end function v

end module weno_mod
