program passiveWenoDriver
use global_vars_mod
use weno_mod
implicit none

integer::i,j,k,ts
real::q(-3:nx+4)

!Initialization
open(unit=10,file='output_adv.txt')

xmin=0.
xmax=1.
dx=(xmax-xmin)/real(nx)
t=0.

call initialize(q)
!file output
call fileOut(q)

print*,'Done initialization.'


!Simulation loop (w/ fixed total # of time steps)
do ts=1,nt

  print*,'Time Step:',ts
  !compute dt
  print*,'Computing dt...'
  dt=cfl*dx/abs(v(0))
  do i=1,nx
    dt=min(dt,cfl*dx/abs(v(i)))
  end do

  print*,'dt=',dt

  !Runge-Kutta time integration: q^n --> q^n+1
  if(RKOption==1)then
    call RK3TVD_integrator(q,dt)
  else if(RKOption==2)then
    call RK4_integrator(q,dt) 
  end if 

  t=t+dt

  !file output
  call fileOut(q)

end do

close(unit=10)

print*,'Simulation completed.'


contains

subroutine RK3TVD_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4),dt
!Subroutine local variables
real::q1(-3:nx+4),q2(-3:nx+4)
real::flux0(0:nx),flux1(0:nx),flux2(0:nx)
real::L0,L1,L2

integer::i,j

  !Step 1:
  call weno_flux(q,flux0)
  do i=1,nx   
    L0=-(flux0(i)-flux0(i-1))/dx
    q1(i)=q(i)+dt*L0
  end do  

  call bound(q1)
  
  !Step2:
  call weno_flux(q1,flux1)
  do i=1,nx   
    L0=-(flux0(i)-flux0(i-1))/dx
    L1=-(flux1(i)-flux1(i-1))/dx

    q2(i)=q1(i)+(dt/4.)*(-3.*L0+L1)
   end do

  call bound(q2)

  !Step3:
  call weno_flux(q2,flux2)
  do i=1,nx
    L0=-(flux0(i)-flux0(i-1))/dx
    L1=-(flux1(i)-flux1(i-1))/dx
    L2=-(flux2(i)-flux2(i-1))/dx

    q(i)=q2(i)+(dt/12.)*(-L0-L1+8.*L2)
  end do

  call bound(q)

end subroutine RK3TVD_integrator


subroutine RK4_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4),dt
!Subroutine local variables
real::q1(-3:nx+4),q2(-3:nx+4),q3(-3:nx+4)
real::flux0(0:nx),flux1(0:nx),flux2(0:nx),flux3(0:nx)
real::L0,L1,L2,L3


integer::i,j

  !Step 1:
  call weno_flux(q,flux0)
  do i=1,nx
    L0=-(flux0(i)-flux0(i-1))/dx

    q1(i)=q(i)+(dt/2.)*L0
  end do  

  call bound(q1)
  
  !Step2:
  call weno_flux(q1,flux1)
  do i=1,nx   
    L0=-(flux0(i)-flux0(i-1))/dx
    L1=-(flux1(i)-flux1(i-1))/dx

    q2(i)=q1(i)+(dt/2.)*(-L0+L1)
   end do

  call bound(q2)

  !Step3:
  call weno_flux(q2,flux2)
  do i=1,nx
    L1=-(flux1(i)-flux1(i-1))/dx
    L2=-(flux2(i)-flux2(i-1))/dx

    q3(i)=q2(i)+(dt/2.)*(-L1+2.*L2)
  end do

  call bound(q3)

  !Step4:
  call weno_flux(q3,flux3)
  do i=1,nx
    L0=-(flux0(i)-flux0(i-1))/dx
    L1=-(flux1(i)-flux1(i-1))/dx
    L2=-(flux2(i)-flux2(i-1))/dx
    L3=-(flux3(i)-flux3(i-1))/dx

    q(i)=q3(i)+(dt/6.)*(L0+2.*L1-4.*L2+L3)
  end do

  call bound(q)

end subroutine RK4_integrator

subroutine fileOut(q)
!Subroutine I/O variables
real::q(-3:nx+4)
!Subroutine local variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  write(10,*) x,q(i),qexact(i) 
end do

end subroutine fileOut

subroutine initialize(q)
!I/O variables
real::q(-3:nx+4)
!Local Variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  q(i)=qexact(i)
end do

call bound(q)

end subroutine initialize


subroutine bound(q)
!I/O variables
real::q(-3:nx+4)
!Local variables

!Outflow boundary condition
q(0)=q(1)
q(-1)=q(1)
q(-2)=q(1)
q(-3)=q(1)
q(nx+1)=q(nx)
q(nx+2)=q(nx)
q(nx+3)=q(nx)
q(nx+4)=q(nx)

end subroutine bound

function qexact(i) result(qx)
!I/0 variables
integer,intent(in)::i
real::qx

x=xmin+(i-1)*dx

!Square wave
if(x>=0.1+v(i)*t .and. x<0.3+v(i)*t)then
  qx=1.
else
  qx=0.
end if

end function qexact

end program passiveWenoDriver
