module global_vars_mod
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=200 !number of cells in x-direction
integer,parameter::nt=200 !fixed number of time-steps
real,parameter::cfl=0.8 !cfl number
real,parameter::pi=4.*atan(1.d0)
real,parameter::tfinal=0.5
integer::is
real::xmin,xmax,dx,dt,x,t
real,parameter::eps=1.d-6

integer,parameter::RKoption=2 !1:RK3 2:RK4


end module global_vars_mod



