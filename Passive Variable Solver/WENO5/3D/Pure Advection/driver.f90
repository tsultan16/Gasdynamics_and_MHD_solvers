program passiveWenoDriver
use global_vars_mod
use weno_mod
implicit none

integer::i,j,k,l,ts
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
character(len=6)::uniti
real::l2err

!Initialization
open(unit=10,file='Output/output_adv.txt')
open(unit=11,file='Output/output_adv1.txt')

xmin=0.
xmax=1.
ymin=0.
ymax=1.0
zmin=0.
zmax=1.0
dx=(xmax-xmin)/real(nx)
dy=(ymax-ymin)/real(ny)
dz=(zmax-zmin)/real(nz)
t=0.


call initialize(q)
!file output
!call fileOut_xslice(q)
!call fileOut_diagslice(q)
call fileOut_2dslice_xy(q,0)
call fileOut_2dslice_yz(q,0)
call fileOut_2dslice_xz(q,0)


print*,'Done initialization.'


!Simulation loop (w/ fixed total # of time steps)
do ts=1,nt

  print*,'Time Step:',ts
  !compute dt
  dt=cfl*dx/abs(v(0,0,0,1))
  do i=1,nx
    do j=1,ny
      do k=1,nz
        dt=min(dt,cfl*dx/abs(v(i,j,k,1)),cfl*dy/abs(v(i,j,k,2)),cfl*dz/abs(v(i,j,k,3)))
      end do
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
  call fileOut_2dslice_xy(q,ts)
  call fileOut_2dslice_yz(q,ts)
  call fileOut_2dslice_xz(q,ts)

  l2err=0.
  do i=1,nx
    do j=1,ny
      do k=1,nz
       l2err=l2err+dx*dy*dz*(q(i,j,k)-qexact(i,j,k))**2
      end do
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
real::q(-3:nx+4,-3:ny+4,-3:nz+4),dt
!Subroutine local variables
real::q1(-3:nx+4,-3:ny+4,-3:nz+4),q2(-3:nx+4,-3:ny+4,-3:nz+4)
real::qint0_plus(0:nx,0:ny,0:nz,3),qint0_minus(0:nx,0:ny,0:nz,3)
real::qint1_plus(0:nx,0:ny,0:nz,3),qint1_minus(0:nx,0:ny,0:nz,3)
real::qint2_plus(0:nx,0:ny,0:nz,3),qint2_minus(0:nx,0:ny,0:nz,3)
real::L0,L1,L2
real::vavg_plus(3),vavg_minus(3)
real::temp1(3),temp2(3),temp3(3),temp4(3)
integer::i,j,k,l

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
    do j=1,ny
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
        q1(i,j,k)=q(i,j,k)+dt*L0
      end do
    end do
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx 
    do j=1,ny  
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
   
        L1=-(v(i,j,k,1)/dx)*( (temp1(1)*qint1_plus(i,j,k,1)+temp2(1)*qint1_minus(i,j,k,1))&
           -(temp3(1)*qint1_plus(i-1,j,k,1)+temp4(1)*qint1_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint1_plus(i,j,k,2)+temp2(2)*qint1_minus(i,j,k,2))&
           -(temp3(2)*qint1_plus(i,j-1,k,2)+temp4(2)*qint1_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint1_plus(i,j,k,3)+temp2(3)*qint1_minus(i,j,k,3))&
           -(temp3(3)*qint1_plus(i,j,k-1,3)+temp4(3)*qint1_minus(i,j,k-1,3)) )
   
        q2(i,j,k)=q1(i,j,k)+(dt/4.)*(-3.*L0+L1)
      end do
    end do
  end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    do j=1,ny
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
   
        L1=-(v(i,j,k,1)/dx)*( (temp1(1)*qint1_plus(i,j,k,1)+temp2(1)*qint1_minus(i,j,k,1))&
           -(temp3(1)*qint1_plus(i-1,j,k,1)+temp4(1)*qint1_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint1_plus(i,j,k,2)+temp2(2)*qint1_minus(i,j,k,2))&
           -(temp3(2)*qint1_plus(i,j-1,k,2)+temp4(2)*qint1_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint1_plus(i,j,k,3)+temp2(3)*qint1_minus(i,j,k,3))&
           -(temp3(3)*qint1_plus(i,j,k-1,3)+temp4(3)*qint1_minus(i,j,k-1,3)) )
   
        L2=-(v(i,j,k,1)/dx)*( (temp1(1)*qint2_plus(i,j,k,1)+temp2(1)*qint2_minus(i,j,k,1))&
           -(temp3(1)*qint2_plus(i-1,j,k,1)+temp4(1)*qint2_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint2_plus(i,j,k,2)+temp2(2)*qint2_minus(i,j,k,2))&
           -(temp3(2)*qint2_plus(i,j-1,k,2)+temp4(2)*qint2_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint2_plus(i,j,k,3)+temp2(3)*qint2_minus(i,j,k,3))&
           -(temp3(3)*qint2_plus(i,j,k-1,3)+temp4(3)*qint2_minus(i,j,k-1,3)) )
   
        q(i,j,k)=q2(i,j,k)+(dt/12.)*(-L0-L1+8.*L2)
      end do
    end do 
  end do

  call bound(q)

end subroutine RK3TVD_integrator


subroutine RK4_integrator(q,dt)

!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4),dt
!Subroutine local variables
real::q1(-3:nx+4,-3:ny+4,-3:nz+4),q2(-3:nx+4,-3:ny+4,-3:nz+4),q3(-3:nx+4,-3:ny+4,-3:nz+4)
real::qint0_plus(0:nx,0:ny,0:nz,3),qint0_minus(0:nx,0:ny,0:nz,3)
real::qint1_plus(0:nx,0:ny,0:nz,3),qint1_minus(0:nx,0:ny,0:nz,3)
real::qint2_plus(0:nx,0:ny,0:nz,3),qint2_minus(0:nx,0:ny,0:nz,3)
real::qint3_plus(0:nx,0:ny,0:nz,3),qint3_minus(0:nx,0:ny,0:nz,3)
real::L0,L1,L2,L3
real::vavg_plus(3),vavg_minus(3)
real::temp1(3),temp2(3),temp3(3),temp4(3)
integer::i,j,k,l

  !Step 1:
  call weno_interpolate(q,qint0_plus,qint0_minus)
  do i=1,nx
    do j=1,ny
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
        q1(i,j,k)=q(i,j,k)+(dt/2.)*L0
      end do
    end do
  end do  

  call bound(q1)
  
  !Step2:
  call weno_interpolate(q1,qint1_plus,qint1_minus)
  do i=1,nx 
    do j=1,ny  
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
   
        L1=-(v(i,j,k,1)/dx)*( (temp1(1)*qint1_plus(i,j,k,1)+temp2(1)*qint1_minus(i,j,k,1))&
           -(temp3(1)*qint1_plus(i-1,j,k,1)+temp4(1)*qint1_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint1_plus(i,j,k,2)+temp2(2)*qint1_minus(i,j,k,2))&
           -(temp3(2)*qint1_plus(i,j-1,k,2)+temp4(2)*qint1_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint1_plus(i,j,k,3)+temp2(3)*qint1_minus(i,j,k,3))&
           -(temp3(3)*qint1_plus(i,j,k-1,3)+temp4(3)*qint1_minus(i,j,k-1,3)) )
   
        q2(i,j,k)=q1(i,j,k)+(dt/2.)*(-L0+L1)
      end do
    end do
  end do

  call bound(q2)

  !Step3:
  call weno_interpolate(q2,qint2_plus,qint2_minus)
  do i=1,nx
    do j=1,ny
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
   
        L1=-(v(i,j,k,1)/dx)*( (temp1(1)*qint1_plus(i,j,k,1)+temp2(1)*qint1_minus(i,j,k,1))&
           -(temp3(1)*qint1_plus(i-1,j,k,1)+temp4(1)*qint1_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint1_plus(i,j,k,2)+temp2(2)*qint1_minus(i,j,k,2))&
           -(temp3(2)*qint1_plus(i,j-1,k,2)+temp4(2)*qint1_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint1_plus(i,j,k,3)+temp2(3)*qint1_minus(i,j,k,3))&
           -(temp3(3)*qint1_plus(i,j,k-1,3)+temp4(3)*qint1_minus(i,j,k-1,3)) )
   
        L2=-(v(i,j,k,1)/dx)*( (temp1(1)*qint2_plus(i,j,k,1)+temp2(1)*qint2_minus(i,j,k,1))&
           -(temp3(1)*qint2_plus(i-1,j,k,1)+temp4(1)*qint2_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint2_plus(i,j,k,2)+temp2(2)*qint2_minus(i,j,k,2))&
           -(temp3(2)*qint2_plus(i,j-1,k,2)+temp4(2)*qint2_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint2_plus(i,j,k,3)+temp2(3)*qint2_minus(i,j,k,3))&
           -(temp3(3)*qint2_plus(i,j,k-1,3)+temp4(3)*qint2_minus(i,j,k-1,3)) )
   
        q3(i,j,k)=q2(i,j,k)+(dt/2.)*(-L1+2.*L2)
      end do
    end do 
  end do

  call bound(q3)

  !Step4:
  call weno_interpolate(q3,qint3_plus,qint3_minus)
  do i=1,nx
    do j=1,ny
      do k=1,nz
        vavg_plus(1)=0.5*(v(i,j,k,1)+v(i+1,j,k,1))  !vx_i+1/2,j,k
        vavg_minus(1)=0.5*(v(i-1,j,k,1)+v(i,j,k,1)) !vx_i-1/2,j,k
        vavg_plus(2)=0.5*(v(i,j,k,2)+v(i,j+1,k,2))  !vy_i,j+1/2,k
        vavg_minus(2)=0.5*(v(i,j-1,k,2)+v(i,j,k,2)) !vy_i,j-1/2,k
        vavg_plus(3)=0.5*(v(i,j,k,2)+v(i,j,k+1,2))  !vy_i,j,k+1/2
        vavg_minus(3)=0.5*(v(i,j,k-1,2)+v(i,j,k,2)) !vy_i,j,k-1/2

        do l=1,3    
          temp1(l)=0.5*(1+sign(1.0,vavg_plus(l)))
          temp2(l)=0.5*(1-sign(1.0,vavg_plus(l)))
          temp3(l)=0.5*(1+sign(1.0,vavg_minus(l)))
          temp4(l)=0.5*(1-sign(1.0,vavg_minus(l)))
        end do

        L0=-(v(i,j,k,1)/dx)*( (temp1(1)*qint0_plus(i,j,k,1)+temp2(1)*qint0_minus(i,j,k,1))&
           -(temp3(1)*qint0_plus(i-1,j,k,1)+temp4(1)*qint0_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint0_plus(i,j,k,2)+temp2(2)*qint0_minus(i,j,k,2))&
           -(temp3(2)*qint0_plus(i,j-1,k,2)+temp4(2)*qint0_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint0_plus(i,j,k,3)+temp2(3)*qint0_minus(i,j,k,3))&
           -(temp3(3)*qint0_plus(i,j,k-1,3)+temp4(3)*qint0_minus(i,j,k-1,3)) )
   
        L1=-(v(i,j,k,1)/dx)*( (temp1(1)*qint1_plus(i,j,k,1)+temp2(1)*qint1_minus(i,j,k,1))&
           -(temp3(1)*qint1_plus(i-1,j,k,1)+temp4(1)*qint1_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint1_plus(i,j,k,2)+temp2(2)*qint1_minus(i,j,k,2))&
           -(temp3(2)*qint1_plus(i,j-1,k,2)+temp4(2)*qint1_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint1_plus(i,j,k,3)+temp2(3)*qint1_minus(i,j,k,3))&
           -(temp3(3)*qint1_plus(i,j,k-1,3)+temp4(3)*qint1_minus(i,j,k-1,3)) )
   
        L2=-(v(i,j,k,1)/dx)*( (temp1(1)*qint2_plus(i,j,k,1)+temp2(1)*qint2_minus(i,j,k,1))&
           -(temp3(1)*qint2_plus(i-1,j,k,1)+temp4(1)*qint2_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint2_plus(i,j,k,2)+temp2(2)*qint2_minus(i,j,k,2))&
           -(temp3(2)*qint2_plus(i,j-1,k,2)+temp4(2)*qint2_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint2_plus(i,j,k,3)+temp2(3)*qint2_minus(i,j,k,3))&
           -(temp3(3)*qint2_plus(i,j,k-1,3)+temp4(3)*qint2_minus(i,j,k-1,3)) )
   
        L3=-(v(i,j,k,1)/dx)*( (temp1(1)*qint3_plus(i,j,k,1)+temp2(1)*qint3_minus(i,j,k,1))&
           -(temp3(1)*qint3_plus(i-1,j,k,1)+temp4(1)*qint3_minus(i-1,j,k,1)) ) &
           -(v(i,j,k,2)/dy)*( (temp1(2)*qint3_plus(i,j,k,2)+temp2(2)*qint3_minus(i,j,k,2))&
           -(temp3(2)*qint3_plus(i,j-1,k,2)+temp4(2)*qint3_minus(i,j-1,k,2)) ) &
           -(v(i,j,k,3)/dz)*( (temp1(3)*qint3_plus(i,j,k,3)+temp2(3)*qint3_minus(i,j,k,3))&
           -(temp3(3)*qint3_plus(i,j,k-1,3)+temp4(3)*qint3_minus(i,j,k-1,3)) )
   
        q(i,j,k)=q3(i,j,k)+(dt/6.)*(L0+2.*L1-4.*L2+L3)
      end do
    end do 
  end do

  call bound(q)

end subroutine RK4_integrator

subroutine initialize(q)
!I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
!Local Variables
integer::i,j,k
real::x,y,z

do i=1,nx
  do j=1,ny
    do k=1,nz
      q(i,j,k)=qexact(i,j,k)
    end do
  end do
end do

call bound(q)

end subroutine initialize


subroutine bound(q)
!I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
!Local variables
integer::j,k

!Outflow boundary condition
do j=-3,ny+4
  do k=-3,nz+4
    q(0,j,k)=q(1,j,k)
    q(-1,j,k)=q(1,j,k)
    q(-2,j,k)=q(1,j,k)
    q(-3,j,k)=q(1,j,k)
    q(nx+1,j,k)=q(nx,j,k)
    q(nx+2,j,k)=q(nx,j,k)
    q(nx+3,j,k)=q(nx,j,k)
    q(nx+4,j,k)=q(nx,j,k)
  end do
end do

do j=-3,nx+4
  do k=-3,nz+4
    q(j,0,k)=q(j,1,k)
    q(j,-1,k)=q(j,1,k)
    q(j,-2,k)=q(j,1,k)
    q(j,-3,k)=q(j,1,k)
    q(j,ny+1,k)=q(j,ny,k)
    q(j,ny+2,k)=q(j,ny,k)
    q(j,ny+3,k)=q(j,ny,k)
    q(j,ny+4,k)=q(j,ny,k)
  end do
end do

do j=-3,nx+4
  do k=-3,ny+4
    q(j,k,0)=q(j,k,1)
    q(j,k,-1)=q(j,k,1)
    q(j,k,-2)=q(j,k,1)
    q(j,k,-3)=q(j,k,1)
    q(j,k,nz+1)=q(j,k,nz)
    q(j,k,nz+2)=q(j,k,nz)
    q(j,k,nz+3)=q(j,k,nz)
    q(j,k,nz+4)=q(j,k,nz)
  end do
end do



end subroutine bound

function v(i,j,k,l) result(vx)
!I/0 variables
integer,intent(in)::i,j,k,l
real::vx
!Local variables

!Uniform advection to the left
if(l==1)then  !vx
vx=1.
else if(l==2)then  !vy
vx=1.
else if(l==3)then  !vz
vx=0.
end if

end function v

function qexact(i,j,k) result(qx)
!I/0 variables
integer,intent(in)::i,j,k
real::qx
!Local variables
real::x,y,z

x=xmin+(i-1)*dx
y=ymin+(j-1)*dy
z=zmin+(k-1)*dz

go to 100
!Square wave
if(x>=0.35+v(i,j,k,1)*t .and. x<0.55+v(i,j,k,1)*t .and. &
 y>=0.35+v(i,j,k,2)*t .and. y<0.55+v(i,j,k,2)*t .and. &
 z>=0.35+v(i,j,k,3)*t .and. z<0.55+v(i,j,k,3)*t )then
  qx=1.
else
  qx=-1.
end if
100 continue

!Spherical wave
if((x-0.5-v(i,j,k,1)*t)**2+(y-0.5-v(i,j,k,2)*t)**2+(z-0.5-v(i,j,k,3)*t)**2<=0.1**2)then
  qx=1.
else
  qx=-1.
end if

end function qexact

subroutine fileOut_xslice(q)
!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
!Subroutine local variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  write(10,*) x,q(i,ny/2,nz/2),qexact(i,ny/2,nz/2) 
end do

end subroutine fileOut_xslice

subroutine fileOut_diagslice(q)
!Subroutine I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
!Subroutine local variables
integer::i
real::x

do i=1,nx
  x=xmin+(i-1)*dx
  if(nx .eq. ny)then
    write(11,*) x,q(i,i,nz/2),qexact(i,i,nz/2)
  end if 
end do

end subroutine fileOut_diagslice

subroutine fileOut_2dslice_xy(q,ts)
!I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
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
  
filename=trim('Output/xy_t=')//trim(uniti)//trim('.txt')
!print*,'filename=',filename

open(unit=12,file=filename)

do j=1,ny
  y=ymin+(j-1)*dy
  do i=1,nx
    x=xmin+(i-1)*dx
    write(12,*) x,y,q(i,j,nz/2),qexact(i,j,nz/2)
  end do
end do

close(12)

end subroutine fileOut_2dslice_xy

subroutine fileOut_2dslice_xz(q,ts)
!I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
integer::ts
!Local variables
integer::i,k
real::x,z
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
  
filename=trim('Output/xz_t=')//trim(uniti)//trim('.txt')
!print*,'filename=',filename

open(unit=12,file=filename)

do k=1,nz
  z=zmin+(k-1)*dz
  do i=1,nx
    x=xmin+(i-1)*dx
    write(12,*) x,z,q(i,ny/2,k),qexact(i,ny/2,k)
  end do
end do

close(13)

end subroutine fileOut_2dslice_xz

subroutine fileOut_2dslice_yz(q,ts)
!I/O variables
real::q(-3:nx+4,-3:ny+4,-3:nz+4)
integer::ts
!Local variables
integer::j,k
real::y,z
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
  
filename=trim('Output/yz_t=')//trim(uniti)//trim('.txt')
!print*,'filename=',filename

open(unit=12,file=filename)

do k=1,nz
  z=zmin+(k-1)*dz
  do j=1,ny
    y=ymin+(j-1)*dx
    write(12,*) y,z,q(nx/2,j,k),qexact(nx/2,j,k)
  end do
end do

close(14)

end subroutine fileOut_2dslice_yz

end program passiveWenoDriver
