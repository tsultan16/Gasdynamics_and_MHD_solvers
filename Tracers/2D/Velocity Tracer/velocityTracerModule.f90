module velocity_tracer_mod_2D
use RoeSolver2D_mod

implicit none

!derived data type (i.e. struct ala c) for tracers
type tr
  real*8 ::current_pos(3),current_vel(3),mass,charge!,exact_pos(3),exact_vel(2);
  integer ::cell_occupied(3)!,exact_cell_occupied(3) !cell occupied by the tracer
  integer ::tracer_id !assigned by order of "creation"
  !integer ::tracer_role !0: passive, 1: active(useless for now)   
  type(tr),pointer ::next   !pointer for linked list capability
end type tr
  

!integer,parameter::option1=1 !1:for fixed # of tracers, 2: variable # of tracers
integer,parameter::N=1 !fixed # of tracers 

integer,parameter::init_option=2! 1:match fluid density distribution 2:uniformly distribute tracers in cells [minCell,maxCell]
integer,parameter::rand_displace=1 !offsets initial tracer placement from cell center by small random amount
integer,parameter::density_option=1 !1:regular 2: CIC
integer,parameter::trboundaryType=3 !1:periodic, 2:outflow, 3:outflow in y & periodic in x
integer,parameter::interpolation_option=1 ;!nearest_cell , 2:linear interpolation

integer ::N_cell !number of tracers in a given cell
real*8 ::tot_fluid_mass
real*8 ::m_cell !fluid mass in given cell

!create tracers
type(tr)::Tracer(N)

real*8::cell_count(-1:nx+2,-1:ny+2)

contains

!-------------------------------------------------------------
subroutine velocityTracerInit()
!-------------------------------------------------------------
!Velocity Tracer Initialization Routine

!Input Variables
!integer::nx
!type(tr)::Tracer(N)
!real*8::dx,xmin,xmax,q1(-1:nx+2,3)

!Local Variables
integer::minCellx,maxCellx,minCelly,maxCelly,i,j,k,l
real*8::rand_num
integer::tr_counter
real*8 ::tot_fluid_mass
real*8 ::m_cell !fluid mass in given cell

open(unit=14,file='Output/velocityTracerDensity.txt')
open(unit=15,file='Output/velocityTracerData.txt')

minCellx=100
maxCellx=100
minCelly=105
maxCelly=105

 
tr_counter=0
k=1

!----------------------------------------------------------------
!Match initial tracer distribution to initial fluid distribution
!----------------------------------------------------------------
if(init_option==1)then
!compute total fluid mass
do i=1,nx
 do l=1,ny
  tot_fluid_mass = tot_fluid_mass+u1(i,l,1)*dx*dy 
 end do
end do

!Distribute tracers according to fluid mass density distribution
do i=1,nx
 do l=1,ny    
  m_cell = u1(i,l,1)*dx*dy !initial fluid mass in ith cell
  rand_num=rand(0)
  !calculate no. of tracers to be placed in the ith cell (equal probablities for rounding Ncell up or down)
  if(rand_num<0.5) then
    N_cell = floor((m_cell/tot_fluid_mass)*N)
  else
    N_cell = ceiling((m_cell/tot_fluid_mass)*N)
  end if 
  tr_counter = tr_counter+N_cell
  if (tr_counter <= N .and. k<=N) then
   do j=1,N_cell	
     if(rand_displace==1) then
       !tracers displaced from cell center by a random amount between(-dx/2,dx/2), this seems to give better agreement between final tracer and fluid density distributions	  
       Tracer(k)%current_pos(1)=xmin+((2.*i-1.)*dx)/2.0 +(2.*rand(0)-1.)*dx/2.
       Tracer(k)%current_pos(2)=ymin+((2.*l-1.)*dy)/2.0 +(2.*rand(0)-1.)*dy/2.
     else
       !Tracers placed at the cell center
       Tracer(k)%current_pos(1)=xmin+(i-1.+(real(j)/real(N_cell)))*dx
       Tracer(k)%current_pos(2)=ymin+(l-1.+(real(j)/real(N_cell)))*dy
     end if  
     Tracer(k)%tracer_id =k	  
     Tracer(k)%cell_occupied(1)=i !1+ floor(Tracer(k)%current_pos(1)/dx)	
     Tracer(k)%cell_occupied(2)=l

     !Tracer(k)%exact_pos=Tracer(k)%current_pos
     !Tracer(k)%exact_cell_occupied=Tracer(k)%cell_occupied
	  
	  !print*,'Tracer#',k,' Pos=',Tracer(k)%current_pos
     write(15,*) k,Tracer(k)%cell_occupied(1),Tracer(k)%cell_occupied(2)
     k=k+1
   end do
  end if	
  !if sum_i(N_cell) exceeds N, then subtract off the excess amount from current N_cell value
  if (tr_counter > N .and. k<=N) then
    N_cell=N_cell-(tr_counter-N)    
    do j=1,N_cell	
      if(rand_displace==1) then
	!tracers displaced from cell center by a small random amount, this seems to give better agreement between final tracer and fluid density distributions	  
        Tracer(k)%current_pos(1)= xmin+(real(2*i-1)*dx)/2.0 +(2.*rand(0)-1.)*dx/2.
        Tracer(k)%current_pos(2)= ymin+(real(2*l-1)*dy)/2.0 +(2.*rand(0)-1.)*dy/2.
      else
	!Tracers placed at the cell center
        Tracer(k)%current_pos(1)= xmin+(real(2*i-1)*dx)/2.0
        Tracer(k)%current_pos(2)= ymin+(real(2*l-1)*dy)/2.0
      end if  
      Tracer(k)%tracer_id =k
      Tracer(k)%cell_occupied(1)=i !1+ floor(Tracer(k)%current_pos(1)/dx)
      Tracer(k)%cell_occupied(2)=l
	 
      !Tracer(k)%exact_pos=Tracer(k)%current_pos
      !Tracer(k)%exact_cell_occupied=Tracer(k)%cell_occupied  
	  !print*,'Tracer#',k,' Pos=',Tracer(k)%current_pos
      write(15,*) k,Tracer(k)%cell_occupied(1),Tracer(k)%cell_occupied(2)
      k=k+1
    end do
  end if
 end do    
end do

!now distribute "leftover" tracers uniformaly across cells
do while(tr_counter < N) 
 do i = 1,nx
  do l=1,ny
    if (k<=N) then 
     if(rand_displace==1) then
     !tracers displaced from cell center by a small random amount, this seems to give better agreement between final tracer and fluid density distributions	  
      Tracer(k)%current_pos(1)= xmin+(real(2.*i-1)*dx)/2.0 +(2.*rand(0)-1.)*dx/2.
      Tracer(k)%current_pos(2)= ymin+(real(2*l-1)*dy)/2.0 +(2.*rand(0)-1.)*dy/2.
     else
      !Tracers placed at the cell center
      Tracer(k)%current_pos= xmin+(real(2.*i-1.)*dx)/2.0
      Tracer(k)%current_pos(2)= ymin+(real(2*l-1)*dy)/2.0 
     end if  
	Tracer(k)%tracer_id =k
	Tracer(k)%cell_occupied(1)=i !1+ floor(Tracer(k)%current_pos(1)/dx)	
        Tracer(k)%cell_occupied(2)=l 
        !Tracer(k)%exact_pos=Tracer(k)%current_pos
	!Tracer(k)%exact_cell_occupied=Tracer(k)%cell_occupied	
		
    !print*,'Tracer#',k,' Pos=',Tracer(k)%current_pos		
	write(15,*) k,Tracer(k)%cell_occupied(1),Tracer(k)%cell_occupied(2)
		  
	k=k+1
	tr_counter=tr_counter+1
	!print*,'Cell #',i,', N_cell= ',N_cell, ', tr_counter= ',tr_counter,', k= ',k 
	
   end if
  end do
 end do
end do



end if

!----------------------------------------------------------------
!Uniformly distribute tracers in cells [minCell,maxCell]
!----------------------------------------------------------------
if(init_option==2)then
do while(tr_counter < N)
 do i=1,nx
  do l=1,ny
   if(i>=minCellx .and. i<=maxCellx .and. l>=minCelly .and. l<=maxCelly) then
    if (k<=N) then 
     if(rand_displace==1) then
       !tracers displaced from cell center by a small random amount, this seems to give better agreement between final tracer and fluid density distributions	  
       Tracer(k)%current_pos(1)= xmin+real(i-0.5)*dx+(rand(0)-0.5)*dx
       Tracer(k)%current_pos(2)= ymin+real(l-0.5)*dy+(rand(0)-0.5)*dy
     else
       !Tracers placed at the cell center
       Tracer(k)%current_pos(1)= xmin+real(i-0.5)*dx
       Tracer(k)%current_pos(2)= ymin+real(l-0.5)*dy 
     end if  
	  
     Tracer(k)%cell_occupied(1)=i !1+ floor(Tracer(k)%current_pos(1)/dx)
     Tracer(k)%cell_occupied(2)=l	 
     Tracer(k)%tracer_id =k
     !Tracer(k)%exact_pos=Tracer(k)%current_pos
     !Tracer(k)%exact_cell_occupied=Tracer(k)%cell_occupied
	
     !print*,'Tracer#',k,' Pos=',Tracer(k)%current_pos
     write(15,*) k,Tracer(k)%cell_occupied(1),Tracer(k)%cell_occupied(2)
		
     k=k+1
     tr_counter=tr_counter+1	
     
    end if
   end if
  end do
 end do
end do

end if

call tracerDensity()

print*,'Completed tracer initialization...' 

end subroutine velocityTracerInit

!-------------------------------------------------------------
subroutine tracerDensity()
!-------------------------------------------------------------

!Input Variables
!integer::nx
!type(tr)::Tracer(N)
!real*8::dx,xmin

!Local Variables
integer ::i,j,cell_index
!real*8::cell_count(-1:nx+2,-1:ny+2)!,exact_cell_count(-1:nx+2,-1:ny+2)
real*8::x,y

cell_count=0
!exact_cell_count=0

do i=1,N
  cell_count(Tracer(i)%cell_occupied(1),Tracer(i)%cell_occupied(2))= &
  cell_count(Tracer(i)%cell_occupied(1),Tracer(i)%cell_occupied(2))+1.
  !exact_cell_count(Tracer(i)%exact_cell_occupied)= &
  !      exact_cell_count(Tracer(i)%exact_cell_occupied)+1.
end do

!write density to file
!do i=1,nx
! do j=1,ny
!  x=xmin+(i-0.5)*dx
!  y=ymin+(j-0.5)*dy
!  write(14,*) x,y,cell_count(i,j) !each row contains the cell position followed by the fraction of tracers occupying that cell
! end do
!end do
     
end subroutine tracerDensity

!-------------------------------------------------------------
subroutine tracerAdvect(dt) 
!-------------------------------------------------------------
!Velocity Tracer Advcetion Routine

!Input Variables
!type(tr)::Tracer(N)
!integer::nx
!integer::dir
real*8::dt!,dx,xmin,xmax,q1(-1:nx+2,3),time

!Local Variables
integer::i,j

do j=1,N
  !calculate interpolated velocity at tracer location
  call velInterpolation(Tracer(j)%current_pos, Tracer(j)%current_vel)
  !update tracer position
  call timeIntegration(Tracer(j)%current_pos,Tracer(j)%current_vel,dt)
 
  !update Cell occupied
  Tracer(j)%cell_occupied(1)=1+floor(Tracer(j)%current_pos(1)/dx)
  Tracer(j)%cell_occupied(2)=1+floor(Tracer(j)%current_pos(2)/dy)
  !Tracer(j)%exact_cell_occupied=1+floor(Tracer(j)%exact_pos/dx)
	
  write(15,*) Tracer(j)%tracer_id,Tracer(j)%cell_occupied,&
	      Tracer(j)%exact_cell_occupied,&
              Tracer(j)%current_pos,Tracer(j)%exact_pos			

end do
call tracerDensity()
  
end subroutine tracerAdvect

!-----------------------------------------------------------
subroutine velInterpolation(tpos,tvel) 
!----------------------------------------------------------- 
!Assign velocity of occupying cell during current time step 

!Input Variables
!integer ::nx
real*8 ::tpos(3),tvel(3)!,q1(-1:nx+2,3),dx


!Local Variables
integer::tcell(2)
real*8::s,v(3) 

tcell(1)=1+floor(tpos(1)/dx)
tcell(2)=1+floor(tpos(2)/dy)

!if(interpolation_option==1)then
!assign fluid velocity of occupying cell
tvel(1)=u1(tcell(1),tcell(2),2)/u1(tcell(1),tcell(2),1) 
tvel(2)=u1(tcell(1),tcell(2),3)/u1(tcell(1),tcell(2),1) 


!if(interpolation_option==2)then
!linear interpolation w/ slope limiter
! if(dir==1)then
!  v(1)=q1(tcell-1,2)/q1(tcell-1,1) 
!  v(2)=q1(tcell,2)/q1(tcell,1) 
!  v(3)=q1(tcell+1,2)/q1(tcell+1,1) 

  !compute slope of interpolation function
!  s=minmod(v(3)-v(2),v(2)-v(1))/dx
  !evalutate interpolation function at tracer position
!  tvel=v(2)+s*(tpos-(tcell-0.5)*dx)
! else
 
! end if
!end if

end subroutine velInterpolation

!-------------------------------------------------------------
subroutine timeIntegration(tpos,tvel,dt)
!-------------------------------------------------------------
! Runge-Kutta Second Order (a.k.a Modified Euler Approximation)

!Input Variables
!integer ::nx
real*8 ::tpos(3),tvel(3),dt!,q1(-1:nx+2,3),xmin,xmax
integer::dir

!Local Variables
real*8::pos1(3),vel1(3)
	     
if(trboundaryType==1) then !periodic in both x and y
 !estimate position after half time-step (predictor step) 
 pos1 = tpos+0.5*dt*tvel;
 !estimate velocity at this location at the initial time
 call velInterpolation(pos1,vel1)
 !estimate position after full time-step (corrector step)
 tpos=tpos+dt*vel1;  
 if(tpos(1)>xmax) then
   tpos(1)=mod(tpos(1),xmax) !enforce periodic boundary conditions
 end if
 if(tpos(1)<xmin) then
   tpos(1)= xmax-mod(abs(tpos(1)),xmax) !enforce periodic boundary conditions
 end if	
 if(tpos(2)>ymax) then
   tpos(2)=mod(tpos(2),ymax) !enforce periodic boundary conditions
 end if
 if(tpos(2)<ymin) then
   tpos(2)= ymax-mod(abs(tpos(2)),ymax) !enforce periodic boundary conditions
 end if	


else if(trboundaryType==2)then ! outflow in both x and y
  if(tpos(1)>xmax .or. tpos(1)<xmin .or. tpos(2)>ymax .or. tpos(2)<ymin)then 
   !dont update tracer position if it's outside the grid
  else
   !estimate position after half time-step (predictor step) 
   pos1 = tpos+0.5*dt*tvel;
   !estimate velocity at this location at the initial time
   call velInterpolation(pos1,vel1)
   !estimate position after full time-step (corrector step)
   tpos = tpos+dt*vel1;  
  end if
else if(trboundaryType==3)then !periodic in x and outflow in y
 if(tpos(2)>ymax .or. tpos(2)<ymin)then
 
 else
  !estimate position after half time-step (predictor step) 
  pos1 = tpos+0.5*dt*tvel;
  !estimate velocity at this location at the initial time
  call velInterpolation(pos1,vel1)
  !estimate position after full time-step (corrector step)
  tpos=tpos+dt*vel1; 
  if(tpos(1)>xmax) then
   tpos(1)=mod(tpos(1),xmax) !enforce periodic boundary conditions
  end if
  if(tpos(1)<xmin) then
   tpos(1)= xmax-mod(abs(tpos(1)),xmax) !enforce periodic boundary conditions
  end if	
 end if
end if

end subroutine timeIntegration

function minmod(x,y) result(z)
  real*8,intent(in)::x,y
  real*8::z
  
  z=sign(1._8,x)*max(0.,min(abs(x),sign(1._8,x)*y))
  
end function minmod


end module velocity_tracer_mod_2D
