module bound_mod
use vars_mod

contains
subroutine bound(u1)

!Input Variables
real*8::u1(-1:nx+2,-1:ny+2,8)

!Local Variables
integer::jj,kk,ll

if (boundaryType==1)then
do jj=-1,nx+2
  do ll=1,8
    u1(jj,0,ll)=u1(jj,1,ll)	
    u1(jj,-1,ll)=u1(jj,1,ll)
    u1(jj,ny+1,ll)=u1(jj,ny,ll)
    u1(jj,ny+2,ll)=u1(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u1(0,kk,ll)=u1(1,kk,ll)	
    u1(-1,kk,ll)=u1(1,kk,ll)
    u1(nx+1,kk,ll)=u1(nx,kk,ll)
    u1(nx+2,kk,ll)=u1(nx,kk,ll)
  end do
end do

else if(boundaryType==2)then
do jj=-1,nx+2
  do ll=1,8
    u1(jj,0,ll)=u1(jj,1,ll)	
    u1(jj,-1,ll)=u1(jj,1,ll)
    u1(jj,ny+1,ll)=u1(jj,ny,ll)
    u1(jj,ny+2,ll)=u1(jj,ny,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u1(0,kk,ll)=u1(nx,kk,ll)	
    u1(-1,kk,ll)=u1(nx-1,kk,ll)
    u1(nx+1,kk,ll)=u1(1,kk,ll)
    u1(nx+2,kk,ll)=u1(2,kk,ll)
  end do
end do

else if(boundaryType==3)then

do jj=-1,nx+2
  do ll=1,8
    u1(jj,0,ll)=u1(jj,1,ll)	
    u1(jj,-1,ll)=u1(jj,1,ll)
    u1(jj,ny+1,ll)=u1(jj,ny,ll)
    u1(jj,ny+2,ll)=u1(jj,ny,ll)      
  end do
  u1(jj,0,3)=-u1(jj,1,3)	
  u1(jj,-1,3)=-u1(jj,1,3)
  u1(jj,ny+1,3)=-u1(jj,ny,3)
  u1(jj,ny+2,3)=-u1(jj,ny,3)
  u1(jj,0,6)=-u1(jj,1,6)	
  u1(jj,-1,6)=-u1(jj,1,6)
  u1(jj,ny+1,6)=-u1(jj,ny,6)
  u1(jj,ny+2,6)=-u1(jj,ny,6)

end do

do kk=-1,ny+2
  do ll=1,8
    u1(0,kk,ll)=u1(1,kk,ll)	
    u1(-1,kk,ll)=u1(1,kk,ll)
    u1(nx+1,kk,ll)=u1(nx,kk,ll)
    u1(nx+2,kk,ll)=u1(nx,kk,ll)
  end do
    u1(0,kk,2)=-u1(1,kk,2)	
    u1(-1,kk,2)=-u1(1,kk,2)
    u1(nx+1,kk,2)=-u1(nx,kk,2)
    u1(nx+2,kk,2)=-u1(nx,kk,2)
    u1(0,kk,5)=-u1(1,kk,5)	
    u1(-1,kk,5)=-u1(1,kk,5)
    u1(nx+1,kk,5)=-u1(nx,kk,5)
    u1(nx+2,kk,5)=-u1(nx,kk,5)
end do

else if(boundaryType==4)then
do jj=-1,nx+2
  do ll=1,8
    u1(jj,0,ll)=u1(jj,ny,ll)	
    u1(jj,-1,ll)=u1(jj,ny-1,ll)
    u1(jj,ny+1,ll)=u1(jj,1,ll)
    u1(jj,ny+2,ll)=u1(jj,2,ll)      
  end do
end do

do kk=-1,ny+2
  do ll=1,8
    u1(0,kk,ll)=u1(nx,kk,ll)	
    u1(-1,kk,ll)=u1(nx-1,kk,ll)
    u1(nx+1,kk,ll)=u1(1,kk,ll)
    u1(nx+2,kk,ll)=u1(2,kk,ll)
  end do
end do


end if

end subroutine bound


end module bound_mod
