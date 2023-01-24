module globalVariables_mod
implicit none


integer,parameter::nx=200
integer,parameter::nt=300

real,parameter::cour=0.75
real,parameter::tol=1.d-30
real,parameter::min_pres=1.d-5
real,parameter::min_dens=1.d-5
real,parameter::pi=4.*atan(1.d0)
real,parameter::eps=1.d-6


real::dt,dx,dy !time step and spatial resolution
real::xmin,xmax,ymin,ymax,x,y
real::gam,gamil,gamul,gamel,gamee,gamuu,gam1,gam2,gam3


integer,parameter::initOption=3 !Ryu-Jones (1995) Shock-Tube Tests
integer,parameter::RKOption=2 !1: RK3 2: RK4
integer,parameter::debug=0 !0:off 1:on



!RJ Tests:
! Test 2 Fail, Test 5 fail with RK4 only, 


end module globalVariables_mod
