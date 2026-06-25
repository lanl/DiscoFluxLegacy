# DiscoFlux Legacy

O# (O5149)</br>

This limited-scope basic research code models dislocation transport through the coupling of sets of dislocation density conservation equations to finite deformation crystal plasticity.
In particular, dislocation densities (i.e. densities of line defects in crystals) are tracked per slip system and dislocation character angle.
As such, its main purpose is for testing advancements to the model and framework.
In contrast to most other codes, we include also the relativistic regime of dislocation drag due to phonon wind and account for the dislocation limiting gliding velocities.
The code makes use of the finite element method and explicit time integration.

 © 2026. Triad National Security, LLC. All rights reserved.

## Author

Darby J. Luscher, Jason R. Mayeur, Hashem Mourad, Abigail Hunter, Mark A. Kenamond, Daniel N. Blaschke

## Requirements

A Fortran 2008 compiler.

## Installation

The latest version of this program can be built with either

`make`

or (using the Fortran Package Manager) via

`fpm build --profile release`

## References


1. D. J. Luscher, J.R. Mayeur, H. M. Mourad, A. Hunter, M. A. Kenamond, “Coupling continuum dislocation transport with crystal plasticity for application to shock loading conditions”, Int. J. Plast. 76 (2016) 111
2. J. R. Mayeur, H. M. Mourad, D. J Luscher, A. Hunter, M. A. Kenamond, “Numerical implementation of a crystal plasticity model with dislocation transport for high strain rate applications“, Modelling Simul. Mater. Sci. Eng. 24 (2016) 045013
3. D. N. Blaschke, D. J. Luscher, “Dislocation drag and its influence on elastic precursor decay”, Int. J. Plast. 144 (2021) 103030
