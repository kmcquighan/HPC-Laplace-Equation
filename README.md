# HPC-Laplace-Equation
This provides sample code for how to use to ScaLaPACK high performance linear algebra libraries. The program solves the 2D Laplace equation d_xx + d_yy  = -2*sin(x)*sin(y) with periodic boundary conditions. 

Note: in order to use the included Makefile you will need to change the Makefile.opts to have the correct path to the ScaLaPACK libraries on your machine (or on whatever remote machine you're using). You will also need to have the compiler mpicc.
