# ppf_critical_point_finder_riccati
Code to identify the critical stability (Re, alpha, c) point for plane Poiseuille flow correct to 28 digits.

# Riccati Finder for the Orr-Sommerfeld Equation

This repository contains a Fortran 90 implementation of a fully analytic Riccati shooting method for the Orr-Sommerfeld equation. It is designed to find the critical parameters of plane Poiseuille flow with high precision.

## Features

* **Single-Pass Integration:** Propagates the Riccati matrix $R$ along with its parameter sensitivities ($S^{Re}$, $S^{\alpha}$, $S^c$) simultaneously in one integration pass using a general explicit Runge-Kutta tableau.
* **Analytic Newton Solve:** The inner Newton solve for the wave speed $c$ utilizes the analytic Jacobian $dF/dc$, eliminating the need for finite differences in $c$.
* **Integrated Sensitivities:** Once $c$ converges, the parameter derivatives $dc_i/dRe$ and $dc_i/d\alpha$ are extracted directly from the same integration pass. 
* **Minimal Finite Differences:** The only remaining finite difference in the pipeline is a small step $h_\alpha$ used to compute the $\alpha$-curvature terms ($d^2c_i/d\alpha^2$, $d^2c_i/d\alpha dRe$).
* **Automated LaTeX Output:** Automatically generates a formatted LaTeX table (`convergence_table_riccati.tex`) demonstrating the spectral convergence of the critical point values.

## Prerequisites & Dependencies

To compile and run this code, you will need:
* A modern Fortran compiler (e.g., `gfortran`, `ifort`).
* **The `rkcoeffs` directory:** **(Important)** This program reads Runge-Kutta coefficients at runtime from external text files. You must have a folder named `rkcoeffs` in the same directory as the executable containing the required coefficient files (e.g., `feagin14Aij.txt` and `feagin14B.txt`). If you change to a different tableau, you must change the rkstages=35 to match the number of stages for the new tableau, which is the number of lines in the tableau files ...Aij.txt or ...B.txt.

## Compilation

You can compile the source code using `gfortran` with the following command:

gfortran critical_point_finder_riccati4.f90 -O3 -o riccati_finder
