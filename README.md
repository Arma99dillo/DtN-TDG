# DtN_TDG
This repository implements the DtN-TDG method for scattering by periodic structures. This repository contains the scripts used to carry out the numerical tests of the method described in the paper

_Trefftz Discontinuous Galerkin methods for scattering by periodic structures (2025)_, by Andrea Moiola and Armando Maria Monforte

# Reproducibility Instructions
This directory contains all the code necessary to reproduce the numerical results presented in the paper. All the codes have been tested on MATLAB R2024b release. The MATLAB pde toolbox is needed for the correct working of the code.

# Contents and Structure

Experiment scripts:
-
The following files contain the main code to run all the experiments in Section 5 of the paper, allowing to reproduce the figures and tables therein:
* The files `ConvergenceTest.m` and `ConvergenceRelative.m` are used to derive the _p_-convergence plots in all the numerical experiments of the paper;
* The file `ThetaConvergence.m` and `MConvergence.m` generate the error plots in Figure 10 and Figure 12 respectively.

Other scripts:
-
* The file `main.m` is used to run the DtN-TDG method and solve and the scattering problem on a specified domain, plotting only the numerical solution, without convergence or error plots; 
* All the configurations used in the paper are implemented in the `GenerateMesh.m` file, which generates the mesh and defines parameters such as the height, the periodicity length and the relative permittivity;
* The `GenerateMeshSol.m` file implements the domain configuration and returns the analytical solution in the implemented cases; 
* The files `MatrixDtNTDG.m` and `rhsDtNTDG.m` implement the linear system as described in Section 4.3, both with and withouth impenetrable obstacles inclusions;

Directories:
-
All the auxiliary files are int the `src` directory. If you want to add your own files it's suggested to do it here.
* The files `MeshEdges.m` and `MeshEdgesDir.m` identify the element edges in internal and boundary ones. Theya are compatible with any mesh configuration, as long as a marker for boundary vertices is provided;
* The files `phi_int.m`, `phi_int_bound`, `DtNInt.m` and `DtNInt_m.m` compute all the boundary integrals without numerical quadrature.

The `quadtriangle` directory implements the Duffy quadrature rule on triangular elements and is needed to compute numerical errors.

Change the parameters and new configurations
-
You can easily change the problem parameters such as the wavenumber, the incident angle, the mesh width and the domain configuration in the `main.m` file and test the convergence in `ConvergenceRelative.m`. In `GenerateMesh.m` you can add custom configuration, you only need to implement a file to generate the mesh and a file to mark the mesh elements depending on the value of the relative permittivity; I suggest to use `MeshDouble.m` and `EpsValDouble.m` as templates.
