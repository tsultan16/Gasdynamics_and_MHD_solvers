***Version#3***

This code contains an implementation of Roe's first order upwind method for solving the equation of Ideal MHD in 2D. Strang-type dimensional splitting is used.
Div(B) constraint is enforced using the "Constrained Transport technique". This involves evolving area averages of the magnetic field normal component defined at cell interfaces, using the electric field at the cell face boundaries as the advective fluxes. Note that the cell inteface area average magnetic field is evolved independently from the volume averages of the fluid variables. Computations of the advective fluxes for these require volume-averages of the magnetic field. These volume-averages are obtained via linear interpolation of the cell-interface magnetic field. 


