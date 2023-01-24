!Finite-Volume solver for scalar conservation law: u_t+f_x=0 
!where f(u) is the flux function. Spatial domain: x in [xmin,xmax].
!Reconstruction of cell-interface left-right states computed using
!the 5th order WENO scheme, with (Lax-Friedrichs) flux splitting.
!3rd Order TVD Runge-Kutta time-integrator is used.
!Reference:Jiang & Wu, J Comp. Phys. 150 (1999)

program scalarHyperbolicWENO
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=200 !number of cells
integer,parameter::nt=500 !total time steps
real,parameter::eps=1.d-6
real,parameter::cfl=0.5 !cfl number
real,parameter::RK_option=2 !1:TVD RK3, 2: RK4

integer::is
real::xmin,xmax,dx,dt,x,t
integer::i,j

real::cr(0:2,0:2)
real::v(-3:nx+3)
real::v0(-3:nx+3),v1(-3:nx+3),v2(-3:nx+3)
real::flux0(0:nx),flux1(0:nx),flux2(0:nx),flux3(0:nx)  !f_i+1/2 
!*********************************************************************

xmin=0.
xmax=1.
dx=(xmax-xmin)/real(nx)

!time step for linear advection with unit speed
dt=cfl*dx
t=0.
open(unit=10,file='output.txt')

!Set initial state
call initialize()
print*,'Done initialization.'

!Set polynomial interpolation coefficients
cr(0,0)=1./3.
cr(0,1)=-7./6.
cr(0,2)=11./6.
cr(1,0)=-1./6.
cr(1,1)=5./6.
cr(1,2)=1./3.
cr(2,0)=1./3.
cr(2,1)=5./6.
cr(2,2)=-1./6.
 
do i=1,nt

  print*,'Time Step,dt=',i,dt
  !Compute cell interface left-right states: v_i+1/2^(+),v_i+1/2^(-) 
 
  if(RK_option==1)then
  !Update cell-average value using 3rd order TVD 
  !Runge-Kutta time integration

  !Step 1:
  call compute_flux(flux0,v)
  do j=1,nx
    v0(j)=v(j)-(dt/dx)*(flux0(j)-flux0(j-1))
  end do  

  call bound(v0)

  !Step2:
  call compute_flux(flux1,v0)
  do j=1,nx
    v1(j)=v0(j)-(1./4.)*(dt/dx)*(-3.*(flux0(j)-flux0(j-1))+(flux1(j)-flux1(j-1)))
  end do

  call bound(v1)
  
  !Step3:
  call compute_flux(flux2,v1)
  do j=1,nx 
    v(j)=v1(j)-(1./12.)*(dt/dx)*(-(flux0(j)-flux0(j-1))-(flux1(j)-flux1(j-1))&
          +8.*(flux2(j)-flux2(j-1)))
  end do

  else if(RK_option==2)then
  !Update cell-average value using 4th order 
  !Runge-Kutta time integration

  !Step 1:
  call compute_flux(flux0,v)
  do j=1,nx
    v0(j)=v(j)-(1./2.)*(dt/dx)*(flux0(j)-flux0(j-1))
  end do  

  call bound(v0)

  !Step2:
  call compute_flux(flux1,v0)
  do j=1,nx
    v1(j)=v0(j)-(1./2.)*(dt/dx)*(-(flux0(j)-flux0(j-1))+(flux1(j)-flux1(j-1)))
  end do

  call bound(v1)

  !Step3:
  call compute_flux(flux2,v1)
  do j=1,nx 
    v2(j)=v1(j)-(1./2.)*(dt/dx)*(-(flux1(j)-flux1(j-1))+2.*(flux2(j)-flux2(j-1)))
  end do

  call bound(v2)

  !Step4:
  call compute_flux(flux3,v2)
  do j=1,nx 
    v(j)=v2(j)-(1./6.)*(dt/dx)*((flux0(j)-flux0(j-1))+2.*(flux1(j)-flux1(j-1))&
          -4.*(flux2(j)-flux2(j-1))+(flux3(j)-flux3(j-1)))
  end do


  end if

  call bound(v)

  t=t+dt

  do j=1,nx
    x=xmin+(j-0.5)*dx
    write(10,*) x,v(j),vexact(x) 
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
  if(x>0.1 .and. x<0.2)then
 v(i)=1.
else  
 v(i)=0.
end if
end do

call bound(v)

v0=v

end subroutine initialize

subroutine compute_flux(flux,v)
!Local Variables
real::flux(0:nx),v(-3:nx+3)
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i

do i=0,nx

  a1=fplus(v(i-1))-fplus(v(i-2))
  b1=fplus(v(i))-fplus(v(i-1))
  c1=fplus(v(i+1))-fplus(v(i))
  d1=fplus(v(i+2))-fplus(v(i+1))

  a2=fminus(v(i+3))-fminus(v(i+2))
  b2=fminus(v(i+2))-fminus(v(i+1))
  c2=fminus(v(i+1))-fminus(v(i))
  d2=fminus(v(i))-fminus(v(i-1))

  flux(i)=(1./12.)*(-f(v(i-1))+7.*f(v(i))+7.*f(v(i+1))-f(v(i+2)))&
          -phi(a1,b1,c1,d1)+phi(a2,b2,c2,d2)

end do



end subroutine compute_flux


function phi(a,b,c,d) result(phix)
!Inputs
real,intent(in)::a,b,c,d
!Outputs
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

function f(v) result(fx)

real,intent(in)::v
real::fx

!linear advection
fx=v

end function f

function fplus(v) result(fx)

real,intent(in)::v
real::fx,alpha

!alpha:=max(f'(v)) 
alpha=1 !for linear advection 

!linear advection
fx=0.5*(f(v)+alpha*v)

end function fplus

function fminus(v) result(fx)

real,intent(in)::v
real::fx,alpha

!alpha:=max(f'(v)) 
alpha=1 !for linear advection 

!linear advection
fx=0.5*(f(v)-alpha*v)

end function fminus

function vexact(x) result(vx)
!Inputs
real,intent(in)::x
!Outputs
real::vx

if(x>=0.1+t .and. x<=0.2+t)then
  vx=1.
else
  vx=0.
end if

end function vexact

subroutine bound(v)
integer::j
real::v(-3:nx+3)

!Outflow boundary condition
do j=0,3
  v(-j)=v(1)
  v(nx+j)=v(nx)
end do

end subroutine bound

end program scalarHyperbolicWENO
