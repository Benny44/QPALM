# QPALM

[![Build Status](https://travis-ci.com/Benny44/QPALM.svg?branch=master)](https://travis-ci.com/Benny44/QPALM)
[![Coverage Status](https://coveralls.io/repos/github/Benny44/QPALM/badge.svg?branch=master)](https://coveralls.io/github/Benny44/QPALM?branch=master)

A proximal augmented Lagrangian method for (possibly **nonconvex**) QPs using semismooth Newton direction and exact line search.

## Installation

* To install the mex interface of QPALM, add QPALM and its subfolders to the matlab path. Then run qpalm_make.m
* To install a C-callable library, compile suitesparse, see [here](https://github.com/jluttine/suitesparse). Then run compile the QPALM directory. For this you need to link to the BLAS and LAPACK libraries on your computer. You can use BLAS= and BLAS_PATH= options in the make command. For example:
```
make BLAS="-lblas -llapack" BLAS_PATH=path/to/blas
```
* To use the Matlab version of QPALM, compile the CHOLMOD mex functions (suitesparse/CHOLMOD/MATLAB/cholmod_make.m), and run QPALM/matlab/mex/PWAlinesearch_setup.m.

## Code Example

Basic demos are available for the different ways to call the solver.
* For the mex interface of QPALM, check out examples/qpalm_mex_demo.m and examples/qpalm_mex_nonconvex_demo.m.
* For the C-version of QPALM, check out examples/qpalm_demo.c.
* For the matlab version of QPALM, check out examples/qpalm_matlab_demo.m.

## Documentation

You can now find the the documentation [online](https://benny44.github.io/QPALM/).

## Tests

To run the automated tests, do
```
make test BLAS=... BLAS_PATH=...
```

## Benchmarks (random QPs)

![](randomQP.png)

## Contributors

* **Ben Hermans** - *Main developer*
* **Panagiotis Patrinos** - *Codeveloper*
* **Andreas Themelis** - *Theoretical contributions*

## License

TBA
