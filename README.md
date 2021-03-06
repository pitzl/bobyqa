# Bounded Optimization by Quadratic Optimization

based on the C version by Éric Thiébaut, 2015-2017
https://github.com/emmt/Algorithms.git

based on the FORTRAN version by Mike Powell, 2009

This provides a C implementation of Mike Powell's BOBYQA algorithm for
minimizing a function of many variables.  The method is *derivatives free*
(only the function values are needed) and accounts for bound constraints on
the variables.  The algorithm is described in:

>  M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
>  Without Derivatives."  Technical report, Department of Applied Mathematics
>  and Theoretical Physics, University of Cambridge (2009).

The present code is based on the original FORTRAN version written by Mike
Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
been converted to C by É. Thiébaut.

BOBYQA builds a quadratic model of the objective function from much less than
`(N+1)(N+2)/2` values of the function (with `N` the number of variables).  The
recommended number of points for building the quadratic model is `2*N+1`.  For
smooth objective functions, BOBYQA is expected to be more efficient than COBYLA
(which exploits a more simple linear model but implements arbitrary inequality
constraints).

In addition to being usable from C code, this version of BOBYQA has a few
improvements over the FORTRAN version:

* any objective function can be used (the function is an argument of the
  method) so different objective functions can co-exist in the same code
  (*e.g.* hierarchical optimization with BOBYQA is possible);

* a return code indicates the success of the method or the reason of the
  failure.

## this version

by Daniel Pitzl, DESY, 2018

The simple Makefile uses the g++ compiler.  
(no extern "C" {} in bobyqa.h)

## make

The code consists in two files [`bobyqa.cc`](./bobyqa.cc) and [`bobyqa.h`](./bobyqa.h)
and is documented in the header [`bobyqa.h`](./bobyqa.h).

To build the library, edit the `Makefile` to suit your preferences and do:
```
make
```

A number of macros can defined for the compiler to adapt the type of variables
used in the code (see [`bobyqa.h`](./bobyqa.h) for details).

## usage

include "bobyqa.h" in your project  
add -L/your-path-to-bobyqa -lbobyqa to your Makefile  
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/your-path-to-bobyqa

## example

The test directory contains an example of least-squares fitting
of a random-generated histogram to a Student's t distribution.

It uses ROOT from CERN for histogramming and visualization.

cd test  
edit the Makefile: path to your ROOT installation  
make  
./student

initial: ![initial values](test/initial.png)  
fitted: ![fitted values](test/final.png)  

```
Hesse Eigenvalues  3.70176e+06  1.30337e+06  571758  441020  744.702
Hesse condition number  4970.8
Hesse Det  9.06005e+26

parameters: fitted values and uncertainties

par 0:  0.000619711 +- 0.00132251
par 1:  0.999723 +- 0.00159345
par 2:  3.03673 +- 0.0366171
par 3:  0.800951 +- 0.00132169
par 4:  0.202071 +- 0.00101357

parameter correlation matrix:

 0             1             2             3             4
 0    1    0.00346015    0.00287422     -0.001493    0.00159949
 1                  1       0.47594     -0.113242     0.0969192
 2                                1     -0.688977      0.766164
 3                                              1     -0.727242
 4                                                            1
```
