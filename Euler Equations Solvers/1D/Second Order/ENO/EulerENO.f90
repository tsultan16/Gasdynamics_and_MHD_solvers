!Second Order ENO scheme for Euler's eqautions

program EulerENO
use ENO_mod
implicit none

call init()

!compute solution
do t=1,nt-1
  print*,'% Complete=',(t*1./nt*1.)*100.
  call computeFlux()
  do i=1,nx    
    x=xmin+(i-0.5)*dx
       
    !evolve fluid state
	do k=1,3
	 q2(i,k)=q1(i,k)-(dt/dx)*(g(i,k)-g(i-1,k))
	end do
	pressure=max((gam-1.)*(q2(i,3)-0.5*q2(i,2)*q2(i,2)/q2(i,1)),premin)
	q2(i,1)=max( q2(i,1),densmin) 
	 
    write(11,*) x,q2(i,1),q2(i,2)/q2(i,1),pressure
	write(12,*) x,q1(i,1),q1(i,2),q1(i,3),qh(i,1),qh(i,2),qh(i,3)
    		
  end do

  call bound()
  q1=q2

end do

print*,'Done...'

close(unit=10)
close(unit=11)

 
end program EulerENO
