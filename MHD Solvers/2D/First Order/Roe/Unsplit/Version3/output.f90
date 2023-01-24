module output_mod
use vars_mod

contains

subroutine horcut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
integer,parameter::ycell=ny/2
integer::jj
real*8::dens,velx,vely,velz,pres,magx,magy,magz

do jj=1,nx
  x=xmin+(jj-0.5)*dx
  dens=u1(jj,ycell,1)
  velx=u1(jj,ycell,2)/u1(jj,ycell,1)
  vely=u1(jj,ycell,3)/u1(jj,ycell,1)
  velz=u1(jj,ycell,4)/u1(jj,ycell,1)
  magx=u1(jj,ycell,5)
  magy=u1(jj,ycell,6)
  magz=u1(jj,ycell,7)
  pres=(gam-1.)*(u1(jj,ycell,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(12,*) x,dens,velx,vely,velz,pres
  write(13,*) x,magx,magy,magz
end do

end subroutine horcut

subroutine vercut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
integer,parameter::xcell=nx/2
integer::jj
real*8::dens,velx,vely,velz,pres,magx,magy,magz

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(xcell,jj,1)
  velx=u1(xcell,jj,2)/u1(xcell,jj,1)
  vely=u1(xcell,jj,3)/u1(xcell,jj,1)
  velz=u1(xcell,jj,4)/u1(xcell,jj,1)
  magx=u1(xcell,jj,5)
  magy=u1(xcell,jj,6)
  magz=u1(xcell,jj,7)
  pres=(gam-1.)*(u1(xcell,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(14,*) y,dens,velx,vely,velz,pres
  write(15,*) y,magx,magy,magz
end do

end subroutine vercut

subroutine diagcut(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

integer::jj
real*8::vperp,vpar,bpar,bperp
real*8::dens,velx,vely,velz,pres,magx,magy,magz

do jj=1,nx
  y=ymin+(jj-0.5)*dy
  dens=u1(jj,jj,1)
  velx=u1(jj,jj,2)/u1(jj,jj,1)
  vely=u1(jj,jj,3)/u1(jj,jj,1)
  velz=u1(jj,jj,4)/u1(jj,jj,1)
  vpar=(vely+velx)/sqrt(2.)
  vperp=(vely-velx)/sqrt(2.)
  magx=u1(jj,jj,5)
  magy=u1(jj,jj,6)
  magz=u1(jj,jj,7)
  bpar=(magy+magx)/sqrt(2.)
  bperp=(magy-magx)/sqrt(2.)
  pres=(gam-1.)*(u1(jj,jj,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))   
  
  write(16,*) y,dens,vpar,vperp,velz,pres
  write(17,*) y,bpar,bperp,magz
end do

end subroutine diagcut

subroutine fileOutput(iunit,u1)

!Input variables
integer::iunit
real*8::u1(-1:nx+2,-1:ny+2,8)
!local variables
integer::jj,kk
real*8::dens,velx,vely,velz,pres,magx,magy,magz

do kk=1,ny
  y=ymin+(kk-0.5)*dy
  do jj=1,nx
   x=xmin+(jj-0.5)*dx
   dens=u1(jj,kk,1)
   velx=u1(jj,kk,2)/u1(jj,kk,1)
   vely=u1(jj,kk,3)/u1(jj,kk,1)
   velz=u1(jj,kk,4)/u1(jj,kk,1)
   magx=u1(jj,kk,5)
   magy=u1(jj,kk,6)
   magz=u1(jj,kk,7)
   pres=(gam-1.)*(u1(jj,kk,8)-0.5*dens*(velx**2.+vely**2.+velz**2.)&
        -0.5*(magx**2.+magy**2.+magz**2.))    

   write(iunit,*) x,y,dens,velx,vely,magx,magy,pres
  end do
end do

end subroutine fileOutput

end module output_mod
