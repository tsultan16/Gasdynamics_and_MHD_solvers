!Central advection equation (2nd order accurate in time). u_t+f_x=u_t+u_x=0 
!with periodic boundaries

program finiteVolCT
implicit none

integer,parameter::nt=20 !time-steps
integer,parameter::nx=128 !number of cells

integer::i,j
real*8::xmin,xmax,dt,dx
real*8::uFS1(-1:nx+2),uBS1(-1:nx+2),uCS1(-1:nx+2)
real*8::uFS2(-1:nx+2),uBS2(-1:nx+2),uCS2(-1:nx+2)  

!initialization

open(unit=10,file='outputCTFS.txt')
open(unit=11,file='outputCTBS.txt')
open(unit=12,file='outputCTCS.txt')
  
  
xmin=-1.
xmax=1.
dx=(xmax-xmin)/(nx*1.)
dt=0.8*dx
do i=1-2,nx+2
 uFS1(i)=u0(xmin+i*dx)
 uBS1(i)=uFS1(i)
 uCS1(i)=uFS1(i)
 if(i>=1 .and. i<=nx)then 
   write(10,*) xmin+i*dx,uFS1(i)
   write(11,*) xmin+i*dx,uBS1(i)
   write(12,*) xmin+i*dx,uCS1(i)
 end if   
end do
!generate initial condition at second time step using FTBS
do i=1-2,nx+2
 uFS2(i)=uBS1(i)-(dt/dx)*(fu(uBS1(i))-fu(uBS1(i-1)))
 uBS2(i)=uBS1(i)-(dt/dx)*(fu(uBS1(i))-fu(uBS1(i-1)))
 uCS2(i)=uBS1(i)-(dt/dx)*(fu(uBS1(i))-fu(uBS1(i-1)))
 if(i>=1 .and. i<=nx)then 
   write(10,*) xmin+i*dx,uFS2(i)
   write(11,*) xmin+i*dx,uBS2(i)
   write(12,*) xmin+i*dx,uCS2(i)
 end if   
end do



!obtain solution 
do i=3,nt
  do j=1-2,nx+2
    if(mod(i,2)==0) then
      uFS2(j)=uFS2(j)-2*(dt/dx)*(fu(uFS2(j+1))-fu(uFS2(j)))
      uBS2(j)=uBS2(j)-2*(dt/dx)*(fu(uBS2(j))-fu(uBS2(j-1)))
      uCS2(j)=uCS2(j)-(dt/dx)*(fu(uCS2(j+1))-fu(uCS2(j-1)))	  
	  !write to file
      if(j>=1 .and. j<=nx)then 
	    write(10,*) xmin+j*dx,uFS2(j)
	    write(11,*) xmin+j*dx,uBS2(j)
	    write(12,*) xmin+j*dx,uCS2(j)
      end if	
	  !enforce periodic boundary condition
      uFS2(-1)=uFS2(nx-2)
      uFS2(0)=uFS2(nx-1)
      uFS2(nx+1)=uFS2(1)
      uFS2(nx+2)=uFS2(2)
	  uBS2(-1)=uBS2(nx-2)
      uBS2(0)=uBS2(nx-1)
      uBS2(nx+1)=uBS2(1)
      uBS2(nx+2)=uBS2(2)
	  uCS2(-1)=uCS2(nx-2)
      uCS2(0)=uCS2(nx-1)
      uCS2(nx+1)=uCS2(1)
      uCS2(nx+2)=uCS2(2)
	else
      uFS1(j)=uFS1(j)-2*(dt/dx)*(fu(uFS1(j+1))-fu(uFS1(j)))
      uBS1(j)=uBS1(j)-2*(dt/dx)*(fu(uBS1(j))-fu(uBS1(j-1)))
      uCS1(j)=uCS1(j)-(dt/dx)*(fu(uCS1(j+1))-fu(uCS1(j-1)))
	  !write to file
      if(j>=1 .and. j<=nx)then 
	    write(10,*) xmin+j*dx,uFS1(j)
	    write(11,*) xmin+j*dx,uBS1(j)
	    write(12,*) xmin+j*dx,uCS1(j)
      end if
      !enforce periodic boundary condition
      uFS1(-1)=uFS1(nx-2)
      uFS1(0)=uFS1(nx-1)
      uFS1(nx+1)=uFS1(1)
      uFS1(nx+2)=uFS1(2)
	  uBS1(-1)=uBS1(nx-2)
      uBS1(0)=uBS1(nx-1)
      uBS1(nx+1)=uBS1(1)
      uBS1(nx+2)=uBS1(2)
	  uCS1(-1)=uCS1(nx-2)
      uCS1(0)=uCS1(nx-1)
      uCS1(nx+1)=uCS1(1)
      uCS1(nx+2)=uCS1(2)	  
    end if	    
  end do
end do

print*,'Done...'

contains

function u0(x) result(ux)
  real*8,intent(in)::x
  real*8::ux
  !square wave initial profile
  if(x<=1./3. .and. x>=-1./3.) then
    ux=1.
  else
    ux=0.
  end if
  
end function u0

function fu(u) result(fux)
  real*8,intent(in)::u
  real*8::fux
  fux=u
end function fu
  
end program finiteVolCT
