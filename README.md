# DiscoFlux Legacy

O# (O5149)</br>

This limited-scope basic research code models dislocation transport through the coupling of sets of dislocation density conservation equations to finite deformation crystal plasticity.
In particular, dislocation densities (i.e. densities of line defects in crystals) are tracked per slip system and dislocation character angle.
As such, its main purpose is for testing advancements to the model and framework.
In contrast to most other codes, we include also the relativistic regime of dislocation drag due to phonon wind and account for the dislocation limiting gliding velocities.
The code makes use of the finite element method and explicit time integration.

This standalone Fortran code was initially developed for the following papers:

1. D. J. Luscher, J.R. Mayeur, H. M. Mourad, A. Hunter, M. A. Kenamond, “Coupling continuum dislocation transport with crystal plasticity for application to shock loading conditions”, 
[Int. J. Plast. 76 (2016) 111](https://doi.org/10.1016/j.ijplas.2015.07.007)
2. J. R. Mayeur, H. M. Mourad, D. J Luscher, A. Hunter, M. A. Kenamond, “Numerical implementation of a crystal plasticity model with dislocation transport for high strain rate applications“, 
[Modelling Simul. Mater. Sci. Eng. 24 (2016) 045013](https://doi.org/10.1088/0965-0393/24/4/045013)
3. D. N. Blaschke, D. J. Luscher, “Dislocation drag and its influence on elastic precursor decay”, 
[Int. J. Plast. 144 (2021) 103030](https://doi.org/10.1016/j.ijplas.2021.103030) ([arxiv.org/abs/2101.10497](https://arxiv.org/abs/2101.10497))

Users may also be interested in the [DiscoFluxM](https://github.com/lanl/DiscoFluxM) code (written in C++ and implemented within MOOSE framework).

 © 2026. Triad National Security, LLC. All rights reserved.

## Authors

Darby J. Luscher, Jason R. Mayeur, Hashem Mourad, Abigail Hunter, Mark A. Kenamond, Daniel N. Blaschke

## Requirements

* A Fortran 2008 compiler
* the [Fortran standard library](https://stdlib.fortran-lang.org/)
* the [Fortran package manager (fpm)](https://fpm.fortran-lang.org/)
* [Doxygen](https://www.doxygen.nl/) or [Ford](https://forddocs.readthedocs.io/en/stable/)>=7 to build the code documentation

### ... and for the postprocessing python-scripts:

* Python >=3.9,</br>
* [numpy](https://numpy.org/doc/stable/user/) >=1.19,</br>
* [scipy](https://docs.scipy.org/doc/scipy/reference/) >=1.9,</br>
* [matplotlib](https://matplotlib.org/) >=3.3</br>
* [PyDislocDyn](https://github.com/dblaschke-LANL/PyDislocDyn) >=1.2.9 (optional, used for additional plots)

## Installation

The latest version of this program can be built using the Fortran Package Manager via

`fpm build --profile release`

## How to use

The executable, as built by fpm, is called `impactshear` and takes a 'jobname' as its argument.
The associated input file is expected to be located in the current folder and must be named 'input_parameters.jobname.dat'.
Several examples are included in the examples folder.
Results are dumped to text files ending in '.F90txt'.
A postprocessing script is included and can be used to generate a number of plots from the .F90txt files via the following syntax:</br>
`python simple_postprocess.py jobname`
(where 'jobname' is the same argument previously passed to the fortran code `impactshear`)

## License

This program is Open-Source under the BSD-3 License; see [LICENSE](LICENSE).
