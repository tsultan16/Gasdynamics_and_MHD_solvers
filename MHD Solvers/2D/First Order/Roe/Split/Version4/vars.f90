module vars_mod
implicit none

integer,parameter::inputOption=1 !1:manual, 2:file
integer,parameter::initOption=0
integer,parameter::debug=0
integer,parameter::outOption=1

integer,parameter::nt=100
integer,parameter::nx=80
integer,parameter::ny=80

integer,parameter::tSkip=1

real*8,parameter::pi=3.14159265359
!-------------------------------------------------------------------------------
integer,parameter::discontinuityDir=4!1:x, 2:y, 3:oblique, 4:partial grid circle 
integer,parameter::boundaryType=1
!BC: 1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic
integer,parameter::perturbationType=2!1:Sinusoidal 2:Random
integer,parameter::KH_test=0 !0:off 1:on
!-------------------------------------------------------------------------------


real*8::tol=1.d-18
real*8::premin=1.d-15,densmin=1.d-15

real*8 ::dt,dx,dy !time step and spatial resolution
real*8::xmin,xmax,ymin,ymax,x,y,cour
real*8 ::gam,gamil,gamul,gamel,gamee,gamuu

character(len=20) :: filename

end module vars_mod
