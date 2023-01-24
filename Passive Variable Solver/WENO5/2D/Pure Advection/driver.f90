program passiveWenoDriver
use global_vars_mod
use weno_mod
implicit none

integer::i,j,k,ts
real::q(-3:nx+4,-3:ny+4)
character(len=6)::uniti
real::l2err

!Initialization
open(unit=10,file='Output/output_adv.txt')
open(unit=11,file='Output/output_adv1.txt')

xmin=0.
xmax=1.
ymin=0.
ymax=1.0
dx=(xmax-xmin)/real(nx)
dy=(ymax-ymin)/real(ny)
t=0.


call initialize(q)
!file output
call fileOut_xslice(q)
call fileOut_diagslice(q)
call fileOut_2d(q,0)

print*,'Done initialization.'


!Simulation loop (w/ fixed total # of time steps)
do ts=1,nt

  print*,'Time Step:',ts
  !compute dt
  dt=cfl*dx/abs(v(0,0,1))
  do i=1,nx
    do j=1,ny
      dt=min(dt,cfl*dx/abs(v(i,j,1)),cfl*dx/abs(v(i,j,2)))
    end do
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
  call fileOut_xslice(q)
  call fileOut_diagslice(q)
  call fileOut_2d(q,ts)

  l2err=0.
  do i=1,nx
    do j=1,ny
     l2err=l2err+dx*dy*(q(i,j)-qexact(i,j))**2
    end do
  end do

  print*,'L2 error=',l2err


end do

close(unit=10)
close(unit=11)

print*,'Simulation completed.'


contains

subroutine RK3TVD_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4),dt
!Subroutine local variables
real::q1(-3:nx+4,-3:ny+4),q2(-3:nx+4,-3:ny+4)
real::qint0_plus(0:nx,0:ny,2),qint0_minus(0:nx,0:ny,2)
real::qint1_plus(0:nx,0:ny,2),qint1_minus(0:nx,0:ny,2)
real::qint2_plus(0:nx,0:ny,2),qint2_minus(0:nx,0:ny,2)
real::L0,L1,L2
real::vavg_plus(2),vavg_minus(2)
real::temp1(2),temp2(2),temp3(2),temp4(2)
integer::i,j,k

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
   do j=1,ny
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
      q1(i,j)=q(i,j)+dt*L0
    end do
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx 
    do j=1,ny  
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
   
      L1=-(v(i,j,1)/dx)*( (temp1(1)*qint1_plus(i,j,1)+temp2(1)*qint1_minus(i,j,1))&
         -(temp3(1)*qint1_plus(i-1,j,1)+temp4(1)*qint1_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint1_plus(i,j,2)+temp2(2)*qint1_minus(i,j,2))&
         -(temp3(2)*qint1_plus(i,j-1,2)+temp4(2)*qint1_minus(i,j-1,2)) )
   
      q2(i,j)=q1(i,j)+(dt/4.)*(-3.*L0+L1)
     end do
   end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    do j=1,ny
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
   
      L1=-(v(i,j,1)/dx)*( (temp1(1)*qint1_plus(i,j,1)+temp2(1)*qint1_minus(i,j,1))&
         -(temp3(1)*qint1_plus(i-1,j,1)+temp4(1)*qint1_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint1_plus(i,j,2)+temp2(2)*qint1_minus(i,j,2))&
         -(temp3(2)*qint1_plus(i,j-1,2)+temp4(2)*qint1_minus(i,j-1,2)) )
    
      L2=-(v(i,j,1)/dx)*( (temp1(1)*qint2_plus(i,j,1)+temp2(1)*qint2_minus(i,j,1))&
         -(temp3(1)*qint2_plus(i-1,j,1)+temp4(1)*qint2_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint2_plus(i,j,2)+temp2(2)*qint2_minus(i,j,2))&
         -(temp3(2)*qint2_plus(i,j-1,2)+temp4(2)*qint2_minus(i,j-1,2)) )
   
      q(i,j)=q2(i,j)+(dt/12.)*(-L0-L1+8.*L2)
    end do 
  end do

  call bound(q)

end subroutine RK3TVD_integrator


subroutine RK4_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4),dt
!Subroutine local variables
real::q1(-3:nx+4,-3:ny+4),q2(-3:nx+4,-3:ny+4),q3(-3:nx+4,-3:ny+4)
real::qint0_plus(0:nx,0:ny,2),qint0_minus(0:nx,0:ny,2)
real::qint1_plus(0:nx,0:ny,2),qint1_minus(0:nx,0:ny,2)
real::qint2_plus(0:nx,0:ny,2),qint2_minus(0:nx,0:ny,2)
real::qint3_plus(0:nx,0:ny,2),qint3_minus(0:nx,0:ny,2)
real::L0,L1,L2,L3
real::vavg_plus(2),vavg_minus(2)
real::temp1(2),temp2(2),temp3(2),temp4(2)
integer::i,j,k

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
   do j=1,ny
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
   
      q1(i,j)=q(i,j)+(dt/2.)*L0

    end do
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx  
    do j=1,ny 
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
   
      L1=-(v(i,j,1)/dx)*( (temp1(1)*qint1_plus(i,j,1)+temp2(1)*qint1_minus(i,j,1))&
         -(temp3(1)*qint1_plus(i-1,j,1)+temp4(1)*qint1_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint1_plus(i,j,2)+temp2(2)*qint1_minus(i,j,2))&
         -(temp3(2)*qint1_plus(i,j-1,2)+temp4(2)*qint1_minus(i,j-1,2)) )
      
      q2(i,j)=q1(i,j)+(dt/2.)*(-L0+L1)
      
    end do
  end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    do j=1,ny
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do
 
      L1=-(v(i,j,1)/dx)*( (temp1(1)*qint1_plus(i,j,1)+temp2(1)*qint1_minus(i,j,1))&
         -(temp3(1)*qint1_plus(i-1,j,1)+temp4(1)*qint1_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint1_plus(i,j,2)+temp2(2)*qint1_minus(i,j,2))&
         -(temp3(2)*qint1_plus(i,j-1,2)+temp4(2)*qint1_minus(i,j-1,2)) )
    
      L2=-(v(i,j,1)/dx)*( (temp1(1)*qint2_plus(i,j,1)+temp2(1)*qint2_minus(i,j,1))&
         -(temp3(1)*qint2_plus(i-1,j,1)+temp4(1)*qint2_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint2_plus(i,j,2)+temp2(2)*qint2_minus(i,j,2))&
         -(temp3(2)*qint2_plus(i,j-1,2)+temp4(2)*qint2_minus(i,j-1,2)) )

    q3(i,j)=q2(i,j)+(dt/2.)*(-L1+2.*L2)

    end do
  end do

  call bound(q3)

  !Step4:
  call weno_interpolate(q3,qint3_plus,qint3_minus)
  do i=1,nx
    do j=1,ny
      vavg_plus(1)=0.5*(v(i,j,1)+v(i+1,j,1))  !vx_i+1/2,j
      vavg_minus(1)=0.5*(v(i-1,j,1)+v(i,j,1)) !vx_i-1/2,j
      vavg_plus(2)=0.5*(v(i,j,2)+v(i,j+1,2))  !vy_i,j+1/2
      vavg_minus(2)=0.5*(v(i,j-1,2)+v(i,j,2)) !vy_i,j-1/2
      do k=1,2    
        temp1(k)=0.5*(1+sign(1.0,vavg_plus(k)))
        temp2(k)=0.5*(1-sign(1.0,vavg_plus(k)))
        temp3(k)=0.5*(1+sign(1.0,vavg_minus(k)))
        temp4(k)=0.5*(1-sign(1.0,vavg_minus(k)))
      end do

      L0=-(v(i,j,1)/dx)*( (temp1(1)*qint0_plus(i,j,1)+temp2(1)*qint0_minus(i,j,1))&
         -(temp3(1)*qint0_plus(i-1,j,1)+temp4(1)*qint0_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint0_plus(i,j,2)+temp2(2)*qint0_minus(i,j,2))&
         -(temp3(2)*qint0_plus(i,j-1,2)+temp4(2)*qint0_minus(i,j-1,2)) )
   
      L1=-(v(i,j,1)/dx)*( (temp1(1)*qint1_plus(i,j,1)+temp2(1)*qint1_minus(i,j,1))&
         -(temp3(1)*qint1_plus(i-1,j,1)+temp4(1)*qint1_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint1_plus(i,j,2)+temp2(2)*qint1_minus(i,j,2))&
         -(temp3(2)*qint1_plus(i,j-1,2)+temp4(2)*qint1_minus(i,j-1,2)) )
    
      L2=-(v(i,j,1)/dx)*( (temp1(1)*qint2_plus(i,j,1)+temp2(1)*qint2_minus(i,j,1))&
         -(temp3(1)*qint2_plus(i-1,j,1)+temp4(1)*qint2_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint2_plus(i,j,2)+temp2(2)*qint2_minus(i,j,2))&
         -(temp3(2)*qint2_plus(i,j-1,2)+temp4(2)*qint2_minus(i,j-1,2)) )

      L3=-(v(i,j,1)/dx)*( (temp1(1)*qint3_plus(i,j,1)+temp2(1)*qint3_minus(i,j,1))&
         -(temp3(1)*qint3_plus(i-1,j,1)+temp4(1)*qint3_minus(i-1,j,1)) ) &
         -(v(i,j,2)/dy)*( (temp1(2)*qint3_plus(i,j,2)+temp2(2)*qint3_minus(i,j,2))&
         -(temp3(2)*qint3_plus(i,j-1,2)+temp4(2)*qint3_minus(i,j-1,2)) )

      q(i,j)=q3(i,j)+(dt/6.)*(L0+2.*L1-4.*L2+L3)

    end do
  end do

  call bound(q)

end subroutine RK4_integrator

subroutine initialize(q)
!I/O variables
real::q(-3:nx+4,-3:ny+4)
!Local Variables
integer::i,j
real::x,y

do i=1,nx
  do j=1,ny
    x=xmin+(i-1)*dx
    y=ymin+(j-1)*dy
    q(i,j)=qexact(i,j)
  end do
end do

call bound(q)

end subroutine initialize


subroutine bound(q)
!I/O variables
real::q(-3:nx+4,-3:ny+4)
!Local variables
integer::i

!Outflow boundary condition
do j=-3,ny+4
  q(0,j)=q(1,j)
  q(-1,j)=q(1,j)
  q(-2,j)=q(1,j)
  q(-3,j)=q(1,j)
  q(nx+1,j)=q(nx,j)
  q(nx+2,j)=q(nx,j)
  q(nx+3,j)=q(nx,j)
  q(nx+4,j)=q(nx,j)
end do

do j=-3,nx+4
  q(j,0)=q(j,1)
  q(j,-1)=q(j,1)
  q(j,-2)=q(j,1)
  q(j,-3)=q(j,1)
  q(j,ny+1)=q(j,ny)
  q(j,ny+2)=q(j,ny)
  q(j,ny+3)=q(j,ny)
  q(j,ny+4)=q(j,ny)
end do


end subroutine bound

function v(i,j,k) result(vx)
!I/0 variables
integer,intent(in)::i,j,k
real::vx
!Local variables

!Uniform advection to the left
if(k==1)then  !x-component
vx=1.
else if(k==2)then  !y-component
vx=1.
end if

end function v

function qexact(i,j) result(qx)
!I/0 variables
integer,intent(in)::i,j
real::qx
!Local variables
real::x,y

x=xmin+(i-1)*dx
y=ymin+(j-1)*dy

go to 100
!Square wave
if(x>=0.35+v(i,j,1)*t .and. x<0.55+v(i,j,1)*t .and. &
 y>=0.35+v(i,j,2)*t .and. y<0.55+v(i,j,2)*t)then
  qx=1.
else
  qx=-1.
end if
100 continue

!Circle wave
if((x-0.3-v(i,j,1)*t)**2+(y-0.3-v(i,j,2)*t)**2<=0.1**2)then
  qx=1.
else
  qx=-1.
end if


!Gaussian
!px=exp(-((x-0.5-t)/0.15)**2)

end function qexact

subroutine fileOut_xslice(q)
!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4)
!Subroutine local variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  write(10,*) x,q(i,ny/2),qexact(i,ny/2) 
end do

end subroutine fileOut_xslice

subroutine fileOut_diagslice(q)
!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4)
!Subroutine local variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  if(nx .eq. ny)then
    write(11,*) x,q(i,i),qexact(i,i)
  end if 
end do

end subroutine fileOut_diagslice

subroutine fileOut_2d(q,ts)
!I/O variables
real::q(-3:nx+4,-3:ny+4)
integer::ts
!Local variables
integer::i,j
real::x,y
character(len=40)::filename
character(len=6)::uniti


if(ts<10)then
  write(uniti,'(I1.1)') ts
else if(ts>=10 .and. ts<100)then
  write(uniti,'(I2.2)') ts
else if(ts>=100 .and. ts<1000)then
  write(uniti,'(I3.3)') ts
else if(ts>=1000 .and. ts<10000)then
  write(uniti,'(I4.3)') ts
else if(ts>=10000 .and. ts<100000)then
  write(uniti,'(I5.3)') ts
end if
  
filename=trim('Output/t=')//trim(uniti)//trim('.txt')
!print*,'filename=',filename

open(unit=12,file=filename)

do j=1,ny
  y=ymin+(j-1)*dy
  do i=1,nx
    x=xmin+(i-1)*dx
    write(12,*) x,y,q(i,j),qexact(i,j)
  end do
end do

close(12)

end subroutine fileOut_2d

end program passiveWenoDriver
