# C version of BOBYQA

Éric Thiébaut 2015-2017 

from https://github.com/emmt/Algorithms.git

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

Copyright (c) 2009, Mike Powell (FORTRAN version).

Copyright (c) 2015, Éric Thiébaut (C version).

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

Daniel Pitzl, 2018

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

The test directory constains an example of least-squares fitting of a
random-generated histogram to a Student's t distribution.
It uses ROOT from CERN for histogramming and visualization
cd test
edit the Makefile to set the path to your ROOT installation
make
student
