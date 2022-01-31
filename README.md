# aemt
Computes a thermal average effective mass tensor at a set of temperatures and chemical potentials
-Be sure to read 'documentation.txt' for how to use the code, and 'compile.sh' for how I compiled it on a linux box (RHEL and ubuntu).  

Given eigenvalues of the FULL 1st Brillouin zone, this code preforms the integrals in equation (10) of
 Hautier et al, How Does Chemistry Influence Electron Effective Mass in Oxides? A
 High-Throughput Computational Analysis, Chemistry of Materials, 2014.
for a set of chemical potentials and temperatures.

Also included are some python scripts for preparing input files plotting some output from the code - either effective masses or energy 
dispersions.

This code is designed to be completly post-process - it will, in theroy, work for any set of k points - they need not be in a grid, or 
even form tetrahedra.  If you need better accuracy, there are other effective mass programs that can preform these integrals more 
accuratly, usually at the cost of needing a specialized set of k points.  

Technical stuff that you should maybe be aware of:
The derivatives are taken by fitting quadratics to a center point and it's two nearest neighbors.
The integrals are evaluated on a rectangular grid.
 !!! IMPORTANT !!! Interpolations from eigenvalues to integration grid points are done using Shepard's method
                   and smoothed with a gaussian profile cross-convlution.  
                   Appropriate sigma values for the smoothing are, ex, (# of grid points / k-point spacing).  
Eigenvalues need to be relativly densly packed in the BZ to get good interpolations.
Nearest neighbor searches are done with the very good, open source program 'kdtree' by John Tsiombikas (https://github.com/jtsiomb/kdtree).  
