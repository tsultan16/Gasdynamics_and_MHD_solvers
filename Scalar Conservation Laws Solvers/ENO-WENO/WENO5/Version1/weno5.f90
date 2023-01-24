!Finite-Volume solver for scalar conservation law: u_t+f_x=0 
!where f(u) is the flux function. Spatial domain: x in [xmin,xmax].
!Reconstruction of cell-interface left-right states computed using
!the 5th order WENO scheme.


program scalarHyperbolic
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=200 !number of cells
integer,parameter::nt=300 !total time steps
integer,parameter::flux_option=2 !1:Godunov, 2:Lax-Freidrichs
real,parameter::eps=1.d-7
real,parameter::cfl=0.4 !cfl number

integer::is
real::xmin,xmax,dx,dt,x,t
integer::i,j

real::cr(-1:2,0:2)
real::vbar(-3:nx+3)
real::vL_WENO(-3:nx+3),vR_WENO(-3:nx+3)
real::flux(0:nx) !f_i+1/2 
!*********************************************************************

xmin=0.
xmax=2.
dx=(xmax-xmin)/real(nx)

!time step for linear advection with unit speed
dt=cfl*dx
t=0.
open(unit=10,file='output.txt')

!Set initial state
call initialize()
print*,'Done initialization.'

!Set polynomial interpolation coefficients
cr(-1,0)=11./16.
cr(-1,1)=-7./6.
cr(-1,2)=1./3.
cr(0,0)=1./3.
cr(0,1)=5./6.
cr(0,2)=-1./6.
cr(1,0)=-1./6.
cr(1,1)=5./6.
cr(1,2)=1./3.
cr(2,0)=1./3.
cr(2,1)=-7./6.
cr(2,2)=11./6.


do i=1,nt

  print*,'Time Step,dt=',i,dt
  !Compute cell interface left-right states: v_i+1/2^(+),v_i+1/2^(-) 
  call compute_LR_states_WENO(vL_WENO,vR_WENO)
  
  do j=0,nx
    flux(j)=h(vR_WENO(j),vL_WENO(j+1))
  end do
  

  !Update cell-average value (first order time-integration)
  do j=1,nx
    vbar(j)=vbar(j)-(dt/dx)*(flux(j)-flux(j-1))
  end do  

  call bound()

  t=t+dt

  do j=1,nx
    x=xmin+(j-0.5)*dx
    write(10,*) x,vbar(j),v(x) 
  end do


end do

print*,'Done.'

close(unit=10)


contains
!*********************************************************************
!Subroutine and Function Definitions
!*********************************************************************

subroutine initialize()
integer::i
real::x

!Step function
do i=1,nx
  x=xmin+(i-1)*dx
  if(x>=0. .and. x<=0.1)then
    vbar(i)=1.
  else
    vbar(i)=0.
  end if
end do

call bound()

end subroutine initialize


subroutine compute_LR_states_WENO(vL,vR)

!Local variables
real::Vdiff(-2:nx+4,-2:nx+4)
real::vL(-3:nx+3),vR(-3:nx+3)
integer::i,j,r
real::vright(0:nx+1,0:2),vleft(0:nx+1,0:2)
real::d(3,0:2),dtilde(3,0:2)
real::beta(0:nx+1,3,0:2)
real::alpha(0:2),alphat(0:2),alpha_s,alphat_s
real::omega(0:nx+1,0:2),omegat(0:nx+1,0:2)


!Compute k-th order reconstruction of v_i+1/2 and v_i+1/2
!in each cell using candiate stencils S_r(i)
vleft=0.
vright=0.
do i=0,nx+1
  do r=0,2
    do j=0,2
      !v_i-1/2^(r)
      vleft(i,r)=vleft(i,r)+cr(r-1,j)*vbar(i-r+j)
      !v_i+1/2^(r)
      vright(i,r)=vright(i,r)+cr(r,j)*vbar(i-r+j)
    end do
  end do
end do

!Set optimal weight values
d=0.
!k=1
d(1,0)=1.
!k=2
d(2,0)=2./3.
d(2,1)=1./3.
!k=3
d(3,0)=3./10.
d(3,1)=6./10.
d(3,2)=1./10.

do i=1,3
  do j=0,2
    dtilde(i,j)=d(i,2-j)
  end do
end do

!Set smoothness indicator values
do i=0,nx+1
  !k=2
  beta(i,2,0)=(vbar(i+1)-vbar(i))**2
  beta(i,2,1)=(vbar(i)-vbar(i-1))**2
  !k=3
  beta(i,3,0)=(13./12.)*(vbar(i)-2.*vbar(i+1)+vbar(i+2))**2 &
              +(3./12.)*(3.*vbar(i)-4.*vbar(i+1)+vbar(i+2))**2
  beta(i,3,1)=(13./12.)*(vbar(i-1)-2.*vbar(i)+vbar(i+1))**2 &
              +(3./12.)*(vbar(i-1)-vbar(i+1))**2
  beta(i,3,2)=(13./12.)*(vbar(i-2)-2.*vbar(i-1)+vbar(i))**2 &
              +(3./12.)*(vbar(i-2)-4.*vbar(i-1)+3.*vbar(i))**2 
end do

!Form general weights
do i=0,nx+1
  alpha_s=0.
  alphat_s=0.
  do j=0,2
    alpha(j)=d(3,j)/((eps+beta(i,3,j))**2)
    alphat(j)=dtilde(3,j)*alpha(j)/d(3,j)
    alpha_s=alpha_s+alpha(j)
    alphat_s=alphat_s+alphat(j)
  end do

  do j=0,2
    omega(i,j)=alpha(j)/alpha_s
    omegat(i,j)=alphat(j)/alphat_s
  end do
end do


!Compute weno-reconstructed values at cell edges and store in file
vL=0.
vR=0.
do i=0,nx+1
  do j=0,2
    vL(i)=vL(i)+omega(i,j)*vleft(i,j)
    vR(i)=vR(i)+omegat(i,j)*vright(i,j)
  end do
end do

end subroutine compute_LR_states_WENO

 
function f(v) result(fx)

real,intent(in)::v
real::fx

!linear advection
fx=v

end function f

function h(vL,vR) result(hx)

real,intent(in)::vL,vR
real::hx,alpha

!alpha=max_[vL,vR] (f'(v))
alpha=1.

if(flux_option==1)then
  if(vL<=vR)then
    hx=min(f(vL),f(vR))
  else
    hx=max(f(vL),f(vR))
  end if
else if(flux_option==2)then
  hx=0.5*(f(vL)+f(vR))-0.5*alpha*(vR-vL)
end if

end function h

recursive function Vdiv(i,l) result(vx)
integer,intent(in)::i,l
real::vx

if(l>i+1)then
  vx=(Vdiv(i+1,l)-Vdiv(i,l-1))/((l-i)*dx)   !recursive function call
else if(l==i+1)then
  vx=vbar(i)
end if

end

!Original function
function v(x) result(vx)
real,intent(in)::x
real::vx,xm

xm=mod(t+0.1,xmax)

!step function
if(x>xm-0.1 .and. x<xm)then
 vx=1.
else  
 vx=0.
end if


end function v 


subroutine bound()

integer::j

!Periodic boundary condition
do j=0,3
  vbar(-j)=vbar(nx-j)
  vbar(nx+j)=vbar(1+j)
end do

end subroutine bound

end program scalarHyperbolic
