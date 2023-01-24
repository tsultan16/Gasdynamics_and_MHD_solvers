!Second Order Godunov-type Methods for Euler's eqautions

program Euler
use RoeSolver_mod
implicit none

integer::alt_flag
character(len=6)::uniti
real::start,finish

alt_flag=0

call cpu_time(start)
call init()

!Compute Solution 
do i=1,nt
  print*,'Time Step=',i,', % Complete=',(real(i)/real(nt))*100.
  
  if(alt_flag==0)then
  !-----------------------------------------------
  !x-sweeps
  !-----------------------------------------------
  if(debug==1)print*,'x sweep...'   
  do k=1,ny
    !compute numerical fluxes
    call computeFluxHarten(nx,k,1)

    if(debug==1)then
      print*,'dtx=',dtx
    end if

    do j=1,nx   

      !print*,'Checkpoint 1...xcell=',j,' ycell=',k

      !update cell averages of conserved variables
      do l=1,5    
        u2(j,k,l)=u1(j,k,l)+(dtx/dx)*(g(j-1,l)-g(j,l)) 
      end do
    end do
  end do

  call bound()	
  u1=u2

  !-----------------------------------------------
  !y-sweeps
  !----------------------------------------------- 
  if(debug==1)print*,'y sweep...'   
  do j=1,nx
    !print*,'Checkpoint5...'
    !compute numerical fluxes
    call computeFluxHarten(ny,j,2)
   
    !print*,'Checkpoint6...'    

    if(debug==1)then
     print*,'dty=',dty
    end if
    

    do k=1,ny  

      !print*,'Checkpoint 2...xcell=',j,' ycell=',k

      !update cell averages of conserved variables  
      u2(j,k,1)=u1(j,k,1)+(dty/dy)*(g(k-1,1)-g(k,1))
      u2(j,k,2)=u1(j,k,2)+(dty/dy)*(g(k-1,3)-g(k,3))
      u2(j,k,3)=u1(j,k,3)+(dty/dy)*(g(k-1,2)-g(k,2))
      u2(j,k,4)=u1(j,k,4)+(dty/dy)*(g(k-1,4)-g(k,4))
      u2(j,k,5)=u1(j,k,5)+(dty/dy)*(g(k-1,5)-g(k,5))      
    end do
  end do

  call bound()	
  u1=u2

  alt_flag=1

  else if(alt_flag==1)then

  !-----------------------------------------------
  !y-sweeps
  !-----------------------------------------------  
  
  if(debug==1)print*,'y sweep...'   
  do j=1,nx
    !compute numerical fluxes
    call computeFluxHarten(ny,j,2)
  
    if(debug==1)then
      print*,'dty=',dty
    end if    

    do k=1,ny   

      !print*,'Checkpoint 3...xcell=',j,' ycell=',k

      !update cell averages of conserved variables  
      u2(j,k,1)=u1(j,k,1)+(dty/dy)*(g(k-1,1)-g(k,1))
      u2(j,k,2)=u1(j,k,2)+(dty/dy)*(g(k-1,3)-g(k,3))
      u2(j,k,3)=u1(j,k,3)+(dty/dy)*(g(k-1,2)-g(k,2))
      u2(j,k,4)=u1(j,k,4)+(dty/dy)*(g(k-1,4)-g(k,4))
      u2(j,k,5)=u1(j,k,5)+(dty/dy)*(g(k-1,5)-g(k,5))        
    end do
  end do

  call bound()	
  u1=u2

  !-----------------------------------------------
  !x-sweeps
  !----------------------------------------------- 
  if(debug==1)print*,'x sweep...'     
  do k=1,ny
    !compute numerical fluxes
    call computeFluxHarten(nx,k,1)


    if(debug==1)then
      print*,'dtx=',dtx
    end if
 
    do j=1,nx   
      !print*,'Checkpoint 4...xcell=',j,' ycell=',k

      !update cell averages of conserved variables
      do l=1,5    
        u2(j,k,l)=u1(j,k,l)+(dtx/dx)*(g(j-1,l)-g(j,l)) 
      end do
    end do
  end do

  call bound()	
  u1=u2

  alt_flag=0

  end if  

  if(mod(i,1)==0)then

  if(i<10)then
    write(uniti,'(I1.1)') i
  else if(i>=10 .and. i<100)then
    write(uniti,'(I2.2)') i
  else if(i>=100 .and. i<1000)then
    write(uniti,'(I3.3)') i
  else if(i>=1000 .and. i<10000)then
    write(uniti,'(I4.3)') i
  end if
  filename=trim('Output/t=')//trim(uniti)//trim('.txt')
  open(unit=i+1000,file=filename) 
  call fileOutput(i+1000)
  close(unit=i+1000)

  end if

  call horcut()
  call vercut()
  

end do 

close(unit=10)
close(unit=12)
close(unit=13)

print*,'Done.'
call cpu_time(finish)
print*,'Time Elapsed (seconds)=',finish-start

end program Euler
