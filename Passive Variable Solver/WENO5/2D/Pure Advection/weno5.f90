!Passive Variable Transport Equation
!Pure Advcetion: p_t +v p_x=0
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

subroutine weno_interpolate(q,qint_plus,qint_minus)
!Local Variables
real::q(-3:nx+4,-3:ny+4),qint_plus(0:nx,0:ny,2),qint_minus(0:nx,0:ny,2)
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i,j,k

do i=0,nx
  do j=0,ny

    !----------------------------------------------
    !q_i+1/2,j
    !----------------------------------------------

    !Phi function arguments for up/downwind stencil
    a1=q(i-1,j)-q(i-2,j)
    b1=q(i,j)-q(i-1,j)
    c1=q(i+1,j)-q(i,j)
    d1=q(i+2,j)-q(i+1,j)

    a2=q(i+3,j)-q(i+2,j)
    b2=q(i+2,j)-q(i+1,j)
    c2=q(i+1,j)-q(i,j)
    d2=q(i,j)-q(i-1,j)

    !WENO approximation of cell interface values
    !q_i+1/2,j
    qint_plus(i,j,1)=(1./12.)*(-q(i-1,j)+7.*q(i,j)+7.*q(i+1,j)-q(i+2,j))-phi(a1,b1,c1,d1)
    qint_minus(i,j,1)=(1./12.)*(-q(i-1,j)+7.*q(i,j)+7.*q(i+1,j)-q(i+2,j))+phi(a2,b2,c2,d2)
   
    !----------------------------------------------
    !q_i,j+1/2
    !----------------------------------------------
    !Phi function arguments for up/downwind stencil
    a1=q(i,j-1)-q(i,j-2)
    b1=q(i,j)-q(i,j-1)
    c1=q(i,j+1)-q(i,j)
    d1=q(i,j+2)-q(i,j+1)

    a2=q(i,j+3)-q(i,j+2)
    b2=q(i,j+2)-q(i,j+1)
    c2=q(i,j+1)-q(i,j)
    d2=q(i,j)-q(i,j-1)

    !WENO approximation of cell interface values
    !q_i,j+1/2
    qint_plus(i,j,2)=(1./12.)*(-q(i,j-1)+7.*q(i,j)+7.*q(i,j+1)-q(i,j+2))-phi(a1,b1,c1,d1)
    qint_minus(i,j,2)=(1./12.)*(-q(i,j-1)+7.*q(i,j)+7.*q(i,j+1)-q(i,j+2))+phi(a2,b2,c2,d2)
   
  end do
end do

end subroutine weno_interpolate


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

end module weno_mod
