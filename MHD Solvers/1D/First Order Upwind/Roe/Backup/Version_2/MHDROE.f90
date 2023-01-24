!Roe's Upwind Method for 1d Ideal MHD Equations 
!(1st Order in space) 

program MHD_Roe
use vars_mod
use init_mod
use RoeSolverMHD_mod
implicit none

!------------------------------------------------------------------
real*8::flux(-1:nx+2,7),u2(-1:nx+2,8)
real*8::dens,velx,vely,velz,pres,magx,magy,magz
integer::jj,kk,ii
real::start,finish
!------------------------------------------------------------------

call cpu_time(start)
call init(u2)
call bound (u2)

!Compute Solution
do ii=1,nt
  print*,'TIME STEP=',ii
  !compute numerical fluxes
  call computeFlux(u2,flux)
  
  print*,'dt=',dt
  
  do jj=1,nx 
    !update cell averages of conerved variables
    do kk=1,7   
     u2(jj,kk)=u2(jj,kk)+(dt/dx)*(flux(jj-1,kk)-flux(jj,kk))
    end do
    x=xmin+(jj-0.5)*dx 
    dens=u2(jj,1)  
    velx=u2(jj,2)/u2(jj,1)
    vely=u2(jj,3)/u2(jj,1)
    velz=u2(jj,4)/u2(jj,1)
    magy=u2(jj,5)
    magz=u2(jj,6)
    magx=u2(jj,8)
    pres=(gam-1.)*(u2(jj,7)-0.5_8*(u2(jj,2)*u2(jj,2)+u2(jj,3)*u2(jj,3) &
          +u2(jj,4)*u2(jj,4))/u2(jj,1)-0.5_8*(magx*magx+magy*magy+magz*magz))
     
    write(11,*) x,dens,velx,vely,velz,pres
    write(12,*) x,magy,magz


  end do

  !call protection(u2)
  call bound(u2)
  !u1=u2

end do 

close(unit=10)
close(unit=11)
close(unit=12)

print*,'Done.'

call cpu_time(finish)

print*,'Time Elapsed (seconds)=',finish-start


contains

subroutine protection(u1)

!Input variables
real*8::u1(-1:nx+2,8)

!Local Variables

do jj=1,nx
  if(u2(jj,1)<tol)then
    u2(jj,1)=premin
  end if
    dens=u1(jj,1)  
    velx=u1(jj,2)/u1(jj,1)
    vely=u1(jj,3)/u1(jj,1)
    velz=u1(jj,4)/u1(jj,1)
    magy=u1(jj,5)
    magz=u1(jj,6)
    magx=u1(jj,8)
    pres=(gam-1.)*(u1(jj,7)-0.5_8*(u1(jj,2)*u1(jj,2)+u1(jj,3)*u1(jj,3) &
          +u1(jj,4)*u1(jj,4))/u1(jj,1)-0.5_8*(magx*magx+magy*magy+magz*magz))
   if(pres<tol)then
     pres=premin
     u1(jj,7)=0.5_8*dens*(velx*velx+vely*vely+velz*velz) &
           +gamul*pres+0.5_8*(magx*magx+magy*magy+magz*magz)
   end if    

    
end do

end subroutine


end program MHD_Roe
