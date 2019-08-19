# QPALM

[![Coverage Status](https://coveralls.io/repos/github/Benny44/QPALM/badge.svg?branch=master)](https://coveralls.io/github/Benny44/QPALM?branch=master)

Platform | CI Status
---------|:---------
Linux    | [![Linux Build Status](https://travis-ci.org/Benny44/QPALM.svg?env=BADGE=linux&branch=master)](https://travis-ci.com/Benny44/QPALM)
OSX      | [![OSX Build Status](https://travis-ci.org/Benny44/QPALM.svg?env=BADGE=osx&branch=master)](https://travis-ci.com/Benny44/QPALM)
Windows  | [![Build status](https://ci.appveyor.com/api/projects/status/4ep1s4q4arr0l6nb/branch/master?svg=true)](https://ci.appveyor.com/project/Benny44/qpalm/branch/master)

A proximal augmented Lagrangian method for (possibly **nonconvex**) QPs using semismooth Newton direction and exact line search.

## Installation

* To install the mex interface of QPALM, add QPALM and its subfolders to the matlab path. Then run qpalm_make.m
* To install a C-callable library, on Linux/MacOS run buildRelease.sh. Windows support coming soon. For example:
```
chmod 755 buildRelease.sh
./buildRelease.sh
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

To build the debug version and run the automated tests, do
```
chmod 755 buildTest.sh
./buildTest.sh
```

## Benchmarks (random QPs)

![](randomQP.png)

## Contributors

* **Ben Hermans** - *Main developer*
* **Panagiotis Patrinos** - *Codeveloper*
* **Andreas Themelis** - *Theoretical contributions*

## License

QPALM is licensed under GPL v3.0
