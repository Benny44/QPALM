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
* To install a C-callable library, check [Bintray](https://bintray.com/benny44/generic/QPALM) for the binaries.
* To use the Matlab version of QPALM, compile the CHOLMOD mex functions (suitesparse/CHOLMOD/MATLAB/cholmod_make.m), and run QPALM/matlab/mex/PWAlinesearch_setup.m.

## Code Example

Basic demos are available for the different ways to call the solver.
* For the mex interface of QPALM, check out examples/qpalm_mex_demo.m and examples/qpalm_mex_nonconvex_demo.m.
* For the C-version of QPALM, check out examples/qpalm_demo.c.
* For the matlab version of QPALM, check out examples/qpalm_matlab_demo.m.

## Documentation

You can now find the the documentation [online](https://benny44.github.io/QPALM/).

## Tests

The QPALM library is tested and has very high coverage. To build the debug version and run the automated tests yourself, in <span>buildTest.sh</span> change 
```
export SUITESPARSE_ROOT_LIB=path/to/suitesparse_libs
export SUITESPARSE_ROOT_INCLUDE=path/to/suitespare_include
```
and then do
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

## Citing

If you use QPALM in your research, please cite the following paper
```
@inproceedings{hermans2019qpalm,
	author		= {Hermans, B. and Pipeleers, G. and Patrinos, P.},
	booktitle	= {58th IEEE Conference on Decision and Control},
	title		= {{QPALM}: {A} {N}ewton-type {P}roximal {A}ugmented {L}agrangian {M}ethod for {Q}uadratic {P}rograms},
	year		= {2019},
	volume		= {},
	number		= {},
	pages		= {},
	doi			= {},
	issn		= {},
	month		= {Dec.},
}
```

## License

QPALM is licensed under GPL v3.0
