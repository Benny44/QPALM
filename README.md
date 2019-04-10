# QPALM

A proximal augmented Lagrangian method for QPs using semismooth Newton direction and exact line search.

## Installation

* To install the mex interface of QPALM, add QPALM and its subfolders to the matlab path. Then run qpalm_make.m
* To install a C-callable library, compile suitesparse, see [here](https://github.com/jluttine/suitesparse). Then run `make` in the QPALM directory.
* To use the Matlab version of QPALM, compile the CHOLMOD mex functions (suitesparse/CHOLMOD/MATLAB/cholmod_make.m), and run QPALM/matlab/mex/PWAlinesearch_setup.m.

## Code Example

Basic demos are available for the different ways to call the solver.
* For the mex interface of QPALM, check out examples/qpalm_mex_demo.m.
* For the C-version of QPALM, check out examples/qpalm_demo.c.
* For the matlab version of QPALM, check out examples/qpalm_matlab_demo.m.

## Tests

In progress

## Benchmarks (random QPs)

![](randomQP.png)

## Contributors

* **Ben Hermans** - *Main developer*
* **Panagiotis Patrinos** - *Codeveloper*
* **Andreas Themelis** - *Theoretical contributions*
## License

TBA
