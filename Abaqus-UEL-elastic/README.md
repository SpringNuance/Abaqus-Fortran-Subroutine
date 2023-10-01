# UEL-elastic
To define an element for FEA under implicit analysis for linear isotropic materials. UEL provides the freedom to modify the base/shape functions and add additional degree of freedom and their spatial derivatives, which are otherwise not accessible in Abaqus.

The UEL is written for 2D analysis using linear triangular elements.

While changing/modifying the input file, the number of elements (nelem) and state variables (nsdv) have to be modified in the subroutine.
