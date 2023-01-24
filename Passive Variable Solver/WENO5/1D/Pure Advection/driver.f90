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
real::qint0_plus(0:nx),qint0_minus(0:nx)
real::qint1_plus(0:nx),qint1_minus(0:nx)
real::qint2_plus(0:nx),qint2_minus(0:nx)
real::L0,L1,L2
real::vavg_plus,vavg_minus
real::temp1,temp2,temp3,temp4


integer::i,j

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
   
    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )
    q1(i)=q(i)+dt*L0
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx   
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )

    L1=-(v(i)/dx)*( (temp1*qint1_plus(i)+temp2*qint1_minus(i))&
          -(temp3*qint1_plus(i-1)+temp4*qint1_minus(i-1)) )

    q2(i)=q1(i)+(dt/4.)*(-3.*L0+L1)
   end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )

    L1=-(v(i)/dx)*( (temp1*qint1_plus(i)+temp2*qint1_minus(i))&
          -(temp3*qint1_plus(i-1)+temp4*qint1_minus(i-1)) )

    L2=-(v(i)/dx)*( (temp1*qint2_plus(i)+temp2*qint2_minus(i))&
          -(temp3*qint2_plus(i-1)+temp4*qint2_minus(i-1)) )

    q(i)=q2(i)+(dt/12.)*(-L0-L1+8.*L2)
  end do

  call bound(q)

end subroutine RK3TVD_integrator


subroutine RK4_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4),dt
!Subroutine local variables
real::q1(-3:nx+4),q2(-3:nx+4),q3(-3:nx+4)
real::qint0_plus(0:nx),qint0_minus(0:nx)
real::qint1_plus(0:nx),qint1_minus(0:nx)
real::qint2_plus(0:nx),qint2_minus(0:nx)
real::qint3_plus(0:nx),qint3_minus(0:nx)
real::L0,L1,L2,L3	
real::vavg_plus,vavg_minus
real::temp1,temp2,temp3,temp4


integer::i,j

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
   
    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )
    q1(i)=q(i)+(dt/2.)*L0
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx   
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )

    L1=-(v(i)/dx)*( (temp1*qint1_plus(i)+temp2*qint1_minus(i))&
          -(temp3*qint1_plus(i-1)+temp4*qint1_minus(i-1)) )

    q2(i)=q1(i)+(dt/2.)*(-L0+L1)
   end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    L1=-(v(i)/dx)*( (temp1*qint1_plus(i)+temp2*qint1_minus(i))&
          -(temp3*qint1_plus(i-1)+temp4*qint1_minus(i-1)) )

    L2=-(v(i)/dx)*( (temp1*qint2_plus(i)+temp2*qint2_minus(i))&
          -(temp3*qint2_plus(i-1)+temp4*qint2_minus(i-1)) )

    q3(i)=q2(i)+(dt/2.)*(-L1+2.*L2)
  end do

  call bound(q3)

  !Step4:
  call weno_interpolate(q3,qint3_plus,qint3_minus)
  do i=1,nx
    vavg_plus=0.5*(v(i)+v(i+1))
    vavg_minus=0.5*(v(i-1)+v(i))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    L0=-(v(i)/dx)*( (temp1*qint0_plus(i)+temp2*qint0_minus(i))&
          -(temp3*qint0_plus(i-1)+temp4*qint0_minus(i-1)) )
 
    L1=-(v(i)/dx)*( (temp1*qint1_plus(i)+temp2*qint1_minus(i))&
          -(temp3*qint1_plus(i-1)+temp4*qint1_minus(i-1)) )

    L2=-(v(i)/dx)*( (temp1*qint2_plus(i)+temp2*qint2_minus(i))&
          -(temp3*qint2_plus(i-1)+temp4*qint2_minus(i-1)) )

    L3=-(v(i)/dx)*( (temp1*qint3_plus(i)+temp2*qint3_minus(i))&
          -(temp3*qint3_plus(i-1)+temp4*qint3_minus(i-1)) )

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

function v(i) result(vx)
!I/0 variables
integer,intent(in)::i
real::vx
!Local variables
real::x

x=xmin+(i-1)*dx

!Uniform advection to the left
!vx=1.

!Discontinuous velocity field (stationary shock)
if(x<0.5)then
  vx=1.0
else 
  vx=0.25
end if

end function v

function qexact(i) result(qx)
!I/0 variables
integer,intent(in)::i
real::qx

x=xmin+(i-1)*dx

!Square wave
if(x>=0.1+v(i)*t .and. x<0.3+v(i)*t)then
  qx=1.
else
  qx=-1.
end if

end function qexact

end program passiveWenoDriver
