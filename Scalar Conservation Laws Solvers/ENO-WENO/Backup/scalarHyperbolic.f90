!Finite-Volume solver for scalar conservation law: u_t+f_x=0 
!where f(u) is the flux function. Spatial domain: x in [xmin,xmax].
!Initial and bouncdary conditions specified.

program scalarHyperbolic
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=50 !number of cells
integer,parameter::k=2 !candidate stencil size (set k<=3)
integer,parameter::nt=1 !total time steps
integer,parameter::r_option=1 !1:ENO, 2:WENO
integer,parameter::flux_option=1 !1:Godunov, 2:Lax-Freidrichs
real,parameter::eps=1.d-7
real,parameter::cfl=0.4 !cfl number

integer::is
real::xmin,xmax,dx,dt,x,t
integer::i,j

real::vbar(-k:nx+k)
real::vL_ENO(-k:nx+k),vR_ENO(-k:nx+k)
real::vL_WENO(-k:nx+k),vR_WENO(-k:nx+k)
real::flux(0:nx) !f_i+1/2 
!*********************************************************************

xmin=-1.
xmax=1.
dx=(xmax-xmin)/real(nx)

!time step for linear advection with unit speed
dt=cfl*dx

open(unit=10,file='output.txt')

!Set initial state
call initialize()

do i=1,nt

  !Compute cell interface left-right states: v_i+1/2^(+),v_i+1/2^(-) 
  if(r_option==1)then
    call compute_LR_states_ENO(vL_ENO,vR_ENO)
  end if

  if(r_option==2)then
    call compute_LR_states_WENO(vL_WENO,vR_WENO)
  end if

  !print*,'j, vL, vR'
  !do j=0,nx
  !  print*,j,vR_WENO(j),vL_WENO(j+1)
  !end do

  !Compute cell-interface numerical flux: f^_j+1/2
  if(r_option==1)then
    do j=0,nx
     flux(j)=h(vR_ENO(j),vL_ENO(j+1))
    end do
  end if

  if(r_option==2)then
    do j=0,nx
     flux(j)=h(vR_WENO(j),vL_WENO(j+1))
    end do
  end if

  !print*,'j     f_j+1/2'
  !do j=0,nx
  ! print*,j,flux(j)
  !end do
  

  do j=1,nx
    x=xmin+(j-0.5)*dx
    write(10,*) x,vbar(j),vR_ENO(j),v(x) 
  end do

  !Update cell-average value (first order time-integration)
  do j=1,nx
    vbar(j)=vbar(j)-(dt/dx)*(flux(j)-flux(j-1))
  end do  

  call bound()

  t=t+dt

!  do j=1,nx
!    x=xmin+(j-0.5)*dx
!    write(10,*) x,vbar(j),vR_ENO(j),v(x) 
!  end do


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
  if(x>=0. .and. x<=0.25)then
    vbar(i)=1.
  else
    vbar(i)=0.
  end if
end do

call bound()

end subroutine initialize


subroutine compute_LR_states_WENO(vL,vR)

!Local variables
real::Vdiff(1-k:nx+k+1,1-k:nx+k+1)
real::vL(-k:nx+k),vR(-k:nx+k)
integer::i,j
real::vright(0:nx+1,0:k-1),vleft(0:nx+1,0:k-1)
real::d(3,0:3-1),dtilde(3,0:3-1)
real::beta(0:nx+1,3,0:2)
real::alpha(0:k-1),alphat(0:k-1),alpha_s,alphat_s
real::omega(0:nx+1,0:k-1),omegat(0:nx+1,0:k-1)

!Compute divided differences
do j=1,k+1
  do i=1-k,nx+1
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

!Compute k-th order reconstruction of v_i+1-2 and v_i+1/2
!in each cell using candiate stencils S_r(i)
do i=0,nx+1
  do j=0,k-1
    vleft(i,j)=p(xmin+(i-1)*dx,j,Vdiff)
    vright(i,j)=p(xmin+i*dx,j,Vdiff)
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

do i=1,k
  do j=0,k-1
    dtilde(i,j)=d(i,k-1-j)
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
  do j=0,k-1
    alpha(j)=d(k,j)/((eps+beta(i,k,j))**2)
    alphat(j)=dtilde(k,j)*alpha(j)/d(k,j)
    alpha_s=alpha_s+alpha(j)
    alphat_s=alphat_s+alphat(j)
  end do

  do j=0,k-1
    omega(i,j)=alpha(j)/alpha_s
    omegat(i,j)=alphat(j)/alphat_s
  end do
end do



!Compute weno-reconstructed values at cell edges and store in file
vL=0.
vR=0.
do i=0,nx+1
  do j=0,k-1
    vL(i)=vL(i)+omega(i,j)*vleft(i,j)
    vR(i)=vR(i)+omegat(i,j)*vright(i,j)
  end do
end do

end subroutine compute_LR_states_WENO


subroutine compute_LR_states_ENO(vL,vR)

!Local Variables
real::vL(-k:nx+k),vR(-k:nx+k)
real::Vdiff(1-k:nx+k+1,1-k:nx+k+1)
integer::i,j
integer::r(0:nx+1)  !eno stencil left-most cell index


!Compute divided differences
do j=1,k+1
  do i=1-k,nx+1
   Vdiff(i,i+j)=Vdiv(i,i+j)  !Vdiff(i,i+j):=V[x_i-1/2,...,x_i+j-1/2]
  end do
end do

!***********************ENO Stencil Selection***************************
!This step involves choosing the "smoothest" 
!stencil from a selection of stencils. Degree of 
!smoothness is measured by the magnitude of the divided differences.
!For a k-point stencils, there will be (k-1) possible stencils to 
!choose from for a given cell. Among these,the one that has the smallest  
!kth order divided difference will be chosen by the ENO condition. In most
!cases, it will turn out that the chosen stencil is also the one that
!does not contain a discontinuous point (i.e. the ENO scheme is 
!designed to pick the stencil that deliberately avoids discontinuities
!to acheive a higher degree of accuracy). 
!************************************************************************ 
do i=0,nx+1
  !Start with two-point stencil
  j=1
  is=i
  do j=2,k
    !apply ENO condition to the two (j+1)-point 
    !stencils formed by adding a point on either 
    !side of the preivious lower oder stencil
    if(abs(Vdiff(is-1,is+j-1))<abs(Vdiff(is,is+j)))then
      is=is-1 
    end if
  end do 
  r(i)=i-is
end do

!Compute k-th order reconstruction of v_i+1-2 and v_i+1/2
!in each cell using the ENO stencil S_r(i)
do i=0,nx+1
    vL(i)=p(xmin+(i-1)*dx,r(i),Vdiff)
    vR(i)=p(xmin+i*dx,r(i),Vdiff)
end do


end subroutine compute_LR_states_ENO
 
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
real::vx


!step function
if(x>t .and. x<t+0.25)then
 vx=1.
else  
 vx=0.
end if


end function v 

!kth order polynomial interpolation function
function p(x,r,Vdiff) result(px)

real,intent(in)::x
integer,intent(in)::r
real,intent(in)::Vdiff(1-k:nx+k+1,1-k:nx+k+1)
real::px,temp1,temp2
integer::i,j,m,l

i=floor((x-xmin)/dx)+1

px=0.
do j=1,k
  temp2=0.
  do m=0,j-1
   temp1=1.
   do l=0,j-1
     if(l .ne. m)then 
       temp1=temp1*(x-(xmin+(i-1-r+l)*dx))
     end if
   end do
   temp2=temp2+temp1
  end do
  px=px+Vdiff(i-r,i-r+j)*temp2
end do

end function p

subroutine bound()

integer::j

!Periodic boundary condition
do j=0,k
  vbar(-j)=vbar(nx-j)
  vbar(nx+j)=vbar(1+j)
end do

end subroutine bound

end program scalarHyperbolic
