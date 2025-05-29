# DtN-TDG
This repository implements the DtN-TDG method for scattering by periodic structures. This repository contains the scripts used to carry out the numerical tests of the method described in the paper

_Trefftz Discontinuous Galerkin methods for scattering by periodic structures (2025)_, by Andrea Moiola and Armando Maria Monforte

# Reproducibility Instructions
This directory contains all the code necessary to reproduce the numerical results presented in the paper. All the codes have been tested on MATLAB R2024b release. The MATLAB Partial Differential Equation Toolbox™ is needed for the correct working of the code.

# Contents and Structure

Experiment scripts:
-
The following files contain the main code to run all the experiments in Section 5 of the paper, allowing to reproduce the figures and tables therein:
* The files `pConvTwoFlatNoabs.m`, `pConvTwoFlatAbs.m`, `pConvThreeFlatEps2.m`, `pConvThreeFlatEps10.m`, `pConvCornerSingularities.m`, `pConvMultipleMaterials.m`, `pConvImpenetrableObst.m` and `pConvImpenetrableObstVariableEps.m` are used to derive the _p_-convergence plots in all the numerical experiments of the paper;
* The files `ThetaConvEps2e.m` and `ThetaConvEps10.m` generate the theta-dependent error plots;
* The file `MConvCornerSingularities.m`  generates the _M_-convergence plot.

Other scripts:
-
* The file `SolveSingleProblem.m` is used to run the DtN-TDG method and solve and the scattering problem on a specified domain, plotting only the numerical solution, without convergence or error plots.

Directories:
-
All the auxiliary files are int the `src` directory.
* All the domain configurations used in the paper are implemented in the `GenerateMesh.m` file, which generates the mesh and defines parameters such as the height, the periodicity length and the relative permittivity;
* The `GenerateMeshSol.m` file implements the domain configuration and returns the analytical solution in the implemented cases; 
* The files `MatrixDtNTDG.m` and `rhsDtNTDG.m` implement the linear system as described in Section 4.3, both with and withouth impenetrable obstacles inclusions;
* The files `MeshEdges.m` and `MeshEdgesDir.m` identify the element edges in internal and boundary ones. Theya are compatible with any mesh configuration, as long as a marker for boundary vertices is provided;
* The files `phi_int.m`, `phi_int_bound`, `DtNInt.m` and `DtNInt_m.m` compute all the boundary integrals without numerical quadrature.

The `quadtriangle` directory implements the Duffy quadrature rule on triangular elements and is needed to compute numerical errors, and is taken from https://www.mathworks.com/matlabcentral/fileexchange/72131-quadtriangle.

Change the parameters and new configurations
-
You can easily change the problem parameters such as the wavenumber, the incident angle, the mesh width and the domain configuration in the `SolveSingleProblem.m` file and test the convergence changing the parameters in any of the _p_-convergence scripts. It's easy to add custom configurations with `GenerateMesh.m`, it is only needed to implement a file to generate the mesh and a file to mark the mesh elements depending on the value of the relative permittivity. I suggest to use `MeshDouble.m` and `EpsValDouble.m` as templates.
