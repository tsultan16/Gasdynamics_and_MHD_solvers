module global_vars_mod
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=50 !number of cells in x-direction
integer,parameter::ny=50 !number of cells in y-direction
integer,parameter::nz=50 !number of cells in z-dircetion
integer,parameter::nt=100 !fixed number of time-steps
real,parameter::cfl=0.2 !cfl number
real,parameter::pi=4.*atan(1.d0)
real,parameter::tfinal=0.5
integer::is
real::xmin,xmax,dx
real::ymin,ymax,dy
real::zmin,zmax,dz
real::dt,t
real,parameter::eps=1.d-6

integer,parameter::RKoption=1 !1:RK3 2:RK4


end module global_vars_mod



