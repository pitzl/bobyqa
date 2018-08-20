/*
 * bobyqa.h -
 *
 * Definitions for Mike Powell's BOBYQA algorithm for minimizing a function of
 * many variables.  The method is "derivatives free" (only the function values
 * are needed) and accounts for bound constraints on the variables.  The
 * algorithm is described in:
 *
 *   M.J.D. Powell, "The BOBYQA Algorithm for Bound Constrained Optimization
 *   Without Derivatives."  Technical report, Department of Applied Mathematics
 *   and Theoretical Physics, University of Cambridge (2009).
 *
 * The present code is based on the original FORTRAN version written by Mike
 * Powell who kindly provides his code on demand (at mjdp@cam.ac.uk) and has
 * been converted to C by É. Thiébaut.
 *
 * Copyright (c) 2009, Mike Powell (FORTRAN version).
 * Copyright (c) 2015, Éric Thiébaut (C version).
 *
 * Read the accompanying `LICENSE` file for details.
 */

#ifndef _BOBYQA_H
#define  _BOBYQA_H 1

#ifndef LOGICAL
# define LOGICAL int
#endif

#ifndef INTEGER
# define INTEGER long
#endif

#ifndef REAL
# define REAL      double
#endif

/* Prototype of the objective function assumed by the BOBYQA routine.  The
   returned value is the function value at X the current variables, N is the
   number of variables and DATA is anything needed by the function (unused by
   BOBYQA itself). */

  typedef REAL bobyqa_objfun( const INTEGER n, const REAL* x, void* data );

/* BOBYQA seeks the least value of a function of many variables, by applying a
   trust region method that forms quadratic models by interpolation.  There is
   usually some freedom in the interpolation conditions, which is taken up by
   minimizing the Frobenius norm of the change to the second derivative of the
   model, beginning with the zero matrix.  The values of the variables are
   constrained by upper and lower bounds.  The arguments of the subroutine are
   as follows.

   N must be set to the number of variables and must be at least two.  NPT is
   the number of interpolation conditions.  Its value must be in the interval
   [N+2,(N+1)(N+2)/2].  Choices that exceed 2*N+1 are not recommended.

   OBJFUN is provided by the user to compute the objective function value at
   the values of the variables X(1),X(2),...,X(N), which are generated
   automatically by BOBYQA in a way that satisfies the bounds given in XL and
   XU.  DATA is anything needed by the function and which is passed as is to
   OBJFUN by BOBYQA.

   Initial values of the variables must be set in X(1),X(2),...,X(N).  They
   will be changed to the values that give the least calculated F.  For
   I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper bounds,
   respectively, on X(I).  The construction of quadratic models requires XL(I)
   to be strictly less than XU(I) for each I.  Further, the contribution to a
   model from changes to the I-th variable is damaged severely by rounding
   errors if XU(I)-XL(I) is too small.

   RHOBEG and RHOEND must be set to the initial and final values of a trust
   region radius, so both must be positive with RHOEND no greater than RHOBEG.
   Typically, RHOBEG should be about one tenth of the greatest expected change
   to a variable, while RHOEND should indicate the accuracy that is required in
   the final values of the variables.  An error return occurs if any of the
   differences XU(I)-XL(I), I=1,...,N, is less than 2*RHOBEG.

   The value of IPRINT should be set to 0, 1, 2 or 3, which controls the amount
   of printing.  Specifically, there is no output if IPRINT=0 and there is
   output only at the return if IPRINT=1.  Otherwise, each new value of RHO is
   printed, with the best vector of variables so far and the corresponding
   value of the objective function.  Further, each new value of F with its
   variables are output if IPRINT=3.

   MAXFUN must be set to an upper bound on the number of calls of OBJFUN.

   The array W will be used for working space.  Its length must be at least
   (NPT+5)*(NPT+N)+3*N*(N+5)/2.  Upon successful return, the first element of W
   will be set to the function value at the solution. */

  extern int bobyqa( const INTEGER n, const INTEGER npt,
                     bobyqa_objfun* objfun, void* data,
                     REAL* x, const REAL* xl, const REAL* xu,
                     const REAL rhobeg, const REAL rhoend,
                     const INTEGER iprint, const INTEGER maxfun, REAL* w );

// Possible values returned by BOBYQA

#define BOBYQA_SUCCESS                 (0) // algorithm converged
#define BOBYQA_BAD_NPT                (-1) // NPT is not in the required interval
#define BOBYQA_TOO_CLOSE              (-2) // insufficient space between the bounds
#define BOBYQA_ROUNDING_ERRORS        (-3) // too much cancellation in a denominator
#define BOBYQA_TOO_MANY_EVALUATIONS   (-4) // maximum number of function evaluations exceeded
#define BOBYQA_STEP_FAILED            (-5) // a trust region step has failed to reduce Q

#endif // _BOBYQA_H
