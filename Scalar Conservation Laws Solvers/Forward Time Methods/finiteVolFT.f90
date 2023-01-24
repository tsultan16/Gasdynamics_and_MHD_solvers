!first order accurate forward time approximation for solving linear
!advection equation. u_t+f_x=u_t+u_x=0 
!with periodic boundaries

program finiteVolFT
implicit none

integer,parameter::nt=100 !time-steps
integer,parameter::nx=256 !number of cells

integer::i,j
real*8::xmin,xmax,dt,dx
real*8::uFS(-1:nx+2),uBS(-1:nx+2),uCS(-1:nx+2)
  

!initialization

open(unit=10,file='outputFTFS.txt')
open(unit=11,file='outputFTBS.txt')
open(unit=12,file='outputFTCS.txt')
  
  
xmin=-1.
xmax=1.
dx=(xmax-xmin)/(nx*1.)
dt=0.8*dx
do i=1-2,nx+2
 uFS(i)=u0(xmin+i*dx)
 uBS(i)=uFS(i)
 uCS(i)=uFS(i)
 if(i>=1 .and. i<=nx)then 
   write(10,*) xmin+i*dx,uFS(i)
   write(11,*) xmin+i*dx,uBS(i)
   write(12,*) xmin+i*dx,uCS(i)
 end if   
end do

!obtain solution 
do i=2,nt
  do j=1-2,nx+2
    !forward space
	  !if(j>1 .and. j<nx) then
        uFS(j)=uFS(j)-(dt/dx)*(fu(uFS(j+1))-fu(uFS(j)))
      !else if(j==nx)then
		!u(j)=u(j)-(dt/dx)*(fu(u(1))-fu(u(j)))
	  !end if	
	!else if(option==2)then !backward space
	  !if(j>1 .and. j<nx) then
	    uBS(j)=uBS(j)-(dt/dx)*(fu(uBS(j))-fu(uBS(j-1)))
      !else if(j==1)then
		!u(j)=u(j)-(dt/dx)*(fu(u(j))-fu(u(nx)))
	  !end if	
	!else if(option==3)then !central space
	  !if(j>1 .and. j<nx) then
        uCS(j)=uCS(j)-(dt/dx)*(fu(uCS(j+1))-fu(uCS(j-1)))
      !else if(j==1)then
	    !u(j)=u(j)-(dt/dx)*(fu(u(j+1))-fu(u(nx)))
	  !else if(j==nx)then
	    !u(j)=u(j)-(dt/dx)*(fu(u(1))-fu(u(j-1)))
	  !end if	  
    !end if
	!enforce periodic boundary condition
    uFS(-1)=uFS(nx-2)
    uFS(0)=uFS(nx-1)
    uFS(nx+1)=uFS(1)
    uFS(nx+2)=uFS(2)
	uBS(-1)=uBS(nx-2)
    uBS(0)=uBS(nx-1)
    uBS(nx+1)=uBS(1)
    uBS(nx+2)=uBS(2)
	uCS(-1)=uCS(nx-2)
    uCS(0)=uCS(nx-1)
    uCS(nx+1)=uCS(1)
    uCS(nx+2)=uCS(2)
	
	!write to file
    if(j>=1 .and. j<=nx)then 
	  write(10,*) xmin+j*dx,uFS(j)
	  write(11,*) xmin+j*dx,uBS(j)
	  write(12,*) xmin+j*dx,uCS(j)
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
  
end program finiteVolFT
