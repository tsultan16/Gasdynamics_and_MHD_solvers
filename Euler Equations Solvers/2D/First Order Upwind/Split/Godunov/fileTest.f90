program fileTest
implicit none

character*12 :: filename
integer :: iunit

iunit=1

write(filename,'("TEST",I3,".txt")')iunit
print*,'filename=',filename
open(unit=iunit,file=filename,status='unknown')  
write(iunit,*) iunit,iunit,iunit
close(unit=iunit)


end program fileTest
