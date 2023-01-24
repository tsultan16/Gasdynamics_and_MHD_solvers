module globalVariables_mod
implicit none


integer,parameter::nx=300
integer,parameter::nt=100

real,parameter::cour=0.8
real,parameter::tol=1.d-20
real,parameter::min_pres=1.d-5
real,parameter::min_dens=1.d-5
real,parameter::pi=3.14159265359
real,parameter::eps=1.d-6


real::dt,dx,dy !time step and spatial resolution
real::xmin,xmax,ymin,ymax,x,y
real::gam,gamil,gamul,gamel,gamee,gamuu


integer,parameter::initOption=0 !Ryu-Jones (1995) Shock-Tube Tests



end module globalVariables_mod
