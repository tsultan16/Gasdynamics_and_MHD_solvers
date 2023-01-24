module vars_mod
implicit none

integer,parameter::inputOption=1 !1:manual, 2:file
integer,parameter::initOption=5
integer,parameter::debug=0

integer,parameter::nt=300
integer,parameter::nx=200
real*8,parameter::pi=3.14159265359
!------------------------------------------------------------------------------3

!-------------------------------------------------------------------------------

real*8,parameter::pi=3.14152
real*8::tol=1.d-18
real*8::premin=1.d-15,densmin=1.d-15

real*8 ::dt,dx !time step and spatial resolution
real*8::xmin,xmax,x,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu





end module vars_mod
