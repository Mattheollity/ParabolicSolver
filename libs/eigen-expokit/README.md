# eigen-expokit

Port of `expokit` (https://www.maths.uq.edu.au/expokit/) to
header-file-only C++ routines using
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).

# Expokit

This package provides Eigen implementations of some routines contained
in [EXPOKIT](http://www.maths.uq.edu.au/expokit). Those routines allow
an efficient calculation of the action of matrix exponentials on vectors
for large sparse matrices. For more details about the methods see
*R.B. Sidje, ACM Trans. Math. Softw., 24(1):130-156, 1998*
(or [its preprint](http://www.maths.uq.edu.au/expokit/paper.pdf)).

# Routines

| | | |
|---|---|---|
| [`padm(A, p)`](./eexpokit/padm.hpp) | `exp(A)` | Pade, scaling and squaring |
| [`expv(t, A, v)`](./eexpokit/expv.hpp) | `exp(t*A)*v` | Krylov subspace projection |
| [`expmv(t, A, v)`](./eexpokit/expmv.hpp) | `exp(t*A)*v` | scaled truncated Taylor series |
| [`chbv(t, H, x)`](./eexpokit/chbv.hpp) | `exp(t*H)*x` | Chebyshev, H symmetric negative definite |
| [`mexpv(t, A, v)`](./eexpokit/mexpv.hpp) | `exp(t*A)*v` | Markov chains, keeps `w` on the simplex |
| [`phiv(t, A, u, v)`](./eexpokit/phiv.hpp) | `exp(t*A)*v + t*phi(t*A)*u` | Krylov, solves `w' = Aw + u` |

All of them take dense or sparse Eigen types, and scalar types other
than the builtin ones: the tests run them over `float`, `double`, `long
double`, both complexes, and [cppduals](https://gitlab.com/tesch1/cppduals)
dual numbers, which differentiate through the exponential.
[example.hpp](./eexpokit/example.hpp) shows their use.

# Against Eigen's own exp()

Eigen computes `exp(A)` by Pade, and to apply it to a vector it has to
form the whole matrix.  `expv()` and `expmv()` never do: the operator
enters only through products with a vector, so they cost what the
sparsity costs.  2D Laplacian, `exp(A)*v`, clang 21 -O3 on an M2:

| n | `A.exp()*v` | `expv` | `expmv` | `padm(A,13)` |
|-----:|------------:|-------:|--------:|-------------:|
| 64 | 0.13 ms | 0.006 ms | 0.39 ms | 0.3 ms |
| 256 | 6.1 ms | 0.14 ms | 0.45 ms | 13.9 ms |
| 1024 | 342 ms | 0.51 ms | 0.56 ms | 817 ms |
| 2025 | 2474 ms | 0.97 ms | 0.85 ms | 6226 ms |

Median of repeated runs after a warmup; the sub-millisecond entries move
by a factor of two between runs.  Every column lands within 3e-15 of the
others, so that is a straight 2500x at n=2025 for the same answer.  Use
`padm()` when you want the matrix itself.

`expmv` does not get cheaper with n because below about n=1000 it is not
measuring the matrix at all: it carries a fixed cost of roughly 0.4 ms,
independent of the problem, from computing its own scaling thresholds.
It bisects a Taylor tail bound for each degree it considers, about
169000 dependent divisions, in place of the table of precomputed
constants the published algorithm uses.  So the `expmv` column is that
setup until the matrix is large enough to outweigh it, and calling
`expmv` in a loop over small problems pays it every time.

`padm` is slower than Eigen's own `exp()` because of how the Pade
numerator and denominator are evaluated, not the degree: both end up at
degree 13 for these matrices.  `padm` uses plain Horner in A squared,
which costs 15 matrix products at p=13, while Eigen groups the terms
around A2, A4 and A6 and needs 6.  Eigen also picks the lowest degree
that will do, where `padm` uses the p it is given, so `padm(A,6)` is
about a third quicker than `padm(A,13)` here.

Accuracy is not uniform, and the tests record what to expect:

| | relative error |
|---|---|
| `padm(A,13)`, `expv`, `expmv` | 1e-16 .. 3e-15 |
| `padm(A,6)`, the default | 1e-8 .. 3e-6 |
| `chbv`, H symmetric negative definite | ~7e-12 |
| `chbv`, anything else | up to 100% |
| `expv`, strongly non-normal A | up to 2e-3, and `err` will not say so |

The last line is the one to watch: on a bidiagonal matrix with a hump of
1e19, `expv` returns 2e-3 while reporting 1e-7, and `expmv` returns
3e-15.  Read `hump`, and prefer `expmv` when it is large.

Those are relative errors, `norm(computed - reference) / norm(reference)`,
worst case over the matrices in the tests: diagonals and rotation
generators, where the exponential is known in closed form, and otherwise
Eigen's own `exp()` as the reference.  Each figure is measured, not
derived from the methods, and the tests assert the same quantity with 10
to 100 times the measured error as the tolerance, so a real regression
trips them and noise does not.  `eexpokit/test_accuracy.cpp` has the
per-routine numbers, and the scalar type barely moves them: what matters
is the Pade degree, and whether the matrix is in the regime a routine
documents.

The original `expokit` is not bundled.  `-DEXPOKIT_BUILD_FORTRAN=ON`
fetches it from https://www.maths.uq.edu.au/expokit/ and builds the
reference Fortran library, its samples, and a test against its
`c1024.crs` matrix.  It is non-commercial-use-only; see NOTICE.

# Using it

The library is header-only, so it needs no build.  Either fetch it:

```cmake
include(FetchContent)
FetchContent_Declare(eigen-expokit
  GIT_REPOSITORY https://gitlab.com/tesch1/eigen-expokit.git
  GIT_TAG        v0.1.0)
FetchContent_MakeAvailable(eigen-expokit)

target_link_libraries(your_target PRIVATE eigen-expokit::eexpokit)
```

or install it once and find it:

```sh
git clone https://gitlab.com/tesch1/eigen-expokit.git && cd eigen-expokit
cmake -Bbuild -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build --target install
```

```cmake
find_package(eigen-expokit 0.1 REQUIRED)
target_link_libraries(your_target PRIVATE eigen-expokit::eexpokit)
```

Either way you then write

```cpp
#include "eexpokit/expv.hpp"
Eigen::VectorXd w = eexpokit::expv(t, A, v);   // exp(t*A)*v
```

Eigen is the consumer's to provide, so any 3.3 or 5.x you already have
works; the package config does not go looking for one.  As a
subdirectory the library target is all you get: `EXPOKIT_BUILD_TESTING`
and `EXPOKIT_INSTALL` default off unless this is the top-level project,
so nothing is downloaded.

# Building it

Only needed to run the tests and examples, which fetch Eigen 5.0.1 and
googletest:

```sh
cmake -Bbuild . && cmake --build build && ctest --test-dir build
./build/example_double
```

`-DEXPOKIT_EIGEN_VERSION=3.3.7` tests against the older Eigen,
`-DEXPOKIT_BUILD_FORTRAN=ON` adds the reference Fortran and the tests
against it, and `-DEXPOKIT_BUILD_DUAL_TESTS=OFF` drops the cppduals
download.  Eigen 5 requires C++14; googletest requires C++17, so the
tests build as C++17 while the library itself still only asks for 14.


# License

MPL-2.0, but see NOTICE: these routines are ports of `expokit`, whose
original notice reads:

```
                              NOTICE

Permission to use, copy, modify, and distribute EXPOKIT and its
supporting documentation for non-commercial purposes, is hereby
granted without fee, provided that this permission message and
copyright notice appear in all copies. Approval must be sought for
commercial purposes as testimony of its usage in applications.

Neither the Institution (University of Queensland) nor the Author
make any representations about the suitability of this software for
any purpose.  This software is provided ``as is'' without express or
implied warranty.

The work resulting from EXPOKIT has been published in ACM-Transactions 
on Mathematical Software, 24(1):130-156, 1998.

The bibtex record of the citation:

@ARTICLE{EXPOKIT,
        AUTHOR  = {Sidje, R. B.},
        TITLE   = {{\sc Expokit.} {A} Software Package for
		  Computing Matrix Exponentials},
        JOURNAL = {ACM Trans. Math. Softw.},
        VOLUME  = {24},
        NUMBER  = {1},
        PAGES   = {130-156}
        YEAR    = {1998}
}

Certain elements of the current software may include inadequacies
that may be corrected at any time, as they are discovered. The Web 
always contains the latest updates.

                            Original Author:
                 Roger B. Sidje <rbs@maths.uq.edu.au>
           Department of Mathematics, University of Queensland 
     Brisbane, QLD-4072, Australia, (c) 1996-2006 All Rights Reserved
```
