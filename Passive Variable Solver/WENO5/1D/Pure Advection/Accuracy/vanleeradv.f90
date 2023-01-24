!Van Leer 2nd Order MUSCL scheme for solving the 1D non-linear advection equation
!u_t+a u_x=0, where a is the wave_speed.


program VanLeerAdvection
implicit none

!***********************************************
!Global Variable Definitions
!***********************************************
integer::nx=25 !number of cells
integer::nt !total time steps
integer,parameter::RK_option=1 !1:1st order 2:2nd order
real,parameter::cfl=0.4 !cfl number (set <=0.5 for stability)
real::x,t,dx,dt,xmin,xmax
integer::i,j
real,parameter::tfinal=0.8
real::error
real,parameter::b=1.2
integer,parameter::niter=12

real,allocatable::u(:),flux(:),u1(:),flux1(:)

!***********************************************
!Initialization
!***********************************************
xmin=0.
xmax=2.0

open(unit=10,file='error_vl.txt')


!**********************************************
!Advect Passive Variable
!**********************************************

do i=1,niter

t=0.
nt=0
dx=(xmax-xmin)/nx
error=0.


allocate(u(-3:nx+3),flux(0:nx),u1(3:nx+3),flux1(0:nx))

print*,'i,nx,dx=',i,nx,dx
call init(u)

do while(t<tfinal)

  nt=nt+1

  !compute time-step size
  dt=cfl*dx/abs(a(0))
  do j=1,nx
    dt=min(dt,cfl*dx/abs(a(j)))
  end do

  !print*,'Time Step,dt =',nt,dt

  if(RK_option==1)then
  
  call compute_flux(u,flux)
  do j=1,nx
    u(j)=u(j)-(dt/dx)*a(j)*(flux(j)-flux(j-1))
  end do

  end if


  if(RK_option==2)then
  !2nd order Runge-Kutta time integration
 
  !*******
  !Step 1:
  !*******
  call compute_flux(u,flux)
  do j=1,nx
    u1(j)=u(j)-(dt/dx)*a(j)*(flux(j)-flux(j-1))
  end do
  call bound(u1)

  !*******
  !Step 2:
  !*******

  call compute_flux(u1,flux1)
  do j=1,nx
    u(j)=0.5*(u(j)+u1(j))-0.5*(dt/dx)*a(j)*(flux1(j)-flux1(j-1))
  end do
 
  end if

  call bound(u)
  t=t+dt
  
end do !end of do while loop

!Compute L1 error
do j=1,nx
  error=error+abs(u(j)-uexact(j))**2
end do
error=sqrt(error/real(nx))

print*,'L1 error=',error
write(10,*) log10(dx),log10(error)
deallocate(u,flux,u1,flux1)

nx=nx*b

end do

close(unit=10)

print*,'Done.'

contains


subroutine init(u)

  !Local Variables
  real::u(-2:nx+3)
  integer::i

  do i=1,nx
    x=xmin+(i-1)*dx
   !Square wave
   !if(x>=0.05 .and. x<=0.6)then
   !   u(i)=1.0
   ! else
   !   u(i)=-1.0
   ! end if
   u(i)=exp(-((x-0.5)/0.15)**2)
  end do

  call bound(u)
  u1=u

end subroutine init

subroutine compute_flux(u,flux)
  !Local variables
  real::u(-3:nx+3),flux(0:nx),S(0:nx+1),a2
  integer::j


  !Compute Limited Slopes
  do j=0,nx+1
    S(j)=minmod(u(j+2)-u(j+1),u(j)-u(j-1))/dx
  end do

  do j=0,nx
    !a_j+1/2    
    a2=0.5*(a(j)+a(j+1)) 
    
    if(a2>=0.)then
      flux(j)=u(j)+0.5*( 1.-min(1.0, (dt/dx)*abs(a2)))*S(j)*dx
    else
      flux(j)=u(j+1)-0.5*( 1.-min(1.0, (dt/dx)*abs(a2)))*S(j+1)*dx
    end if
  end do

end subroutine compute_flux

function limiter(x,y) result(z)
  !Inputs
  real,intent(in)::x,y  
  !Outputs
  real::z
  !Other
  real::b

  b=1. !set b=1 for minmod slope limiter

  if(y>0.0)then
    z=max(0.,min(b*x,y),min(x,b*y))
  else 
    z=min(0.,max(b*x,y),max(x,b*y))
  end if

end function limiter

function minmod(x,y) result(z)
  !Inputs
  real,intent(in)::x,y
  !Outputs
  real::z
  
  z=sign(1.0,x)*max(0.,min(abs(x),sign(1.0,x)*y))
  
end function minmod

subroutine bound(p)
integer::j
real::p(-3:nx+3)

do j=0,3
  !Outflow boundary condition
  p(-j)=p(1)
  p(nx+j)=p(nx)
  !Periodic boundary condition
  !p(-j)=p(nx-j)
  !if(j>0)then
  ! p(nx+j)=p(1+j-1)
  !end if
end do

end subroutine bound

function a(xi) result(ax)
  !Inputs
  integer,intent(in)::xi
  !Outputs
  real::ax,x

  !x=xmin+(xi-0.5)*dx

  !constant wave speed
  ax=1.0

end function a

function uexact(i) result(ux)
!Inputs
integer,intent(in)::i
!Outputs
real::ux

x=xmin+(i-1)*dx

!if(x>=0.05+t .and. x<0.6+t)then
!  ux=1.
!else
!  ux=-1.
!end if

ux=exp(-((x-0.5-t)/0.15)**2)

end function uexact


end program VanLeerAdvection
