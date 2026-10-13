# Changelog

## v0.1.0

- first packaged release: install and export `eigen-expokit::eexpokit`
  with a package config, so downstreams can `find_package` it.
- `padm()`, `chbv()`, `expv()`, `mexpv()`, `phiv()` and `loadcrs()` from
  EXPOKIT, plus `expmv()` after Al-Mohy and Higham, for dense and sparse
  Eigen types and for user-defined scalars.
- `mexpv()` was a stub that returned its input; it now runs the Markov
  projection and keeps the result on the simplex.  `phiv()` is new.
- the routines that act on a vector share one shape: `(t, A, v, tol, m)`
  returning the vector, with `expv()`'s and `mexpv()`'s error estimate and
  hump moved to an optional trailing `expv_stats`.  `chbv()` gained a `t`.
- `tol` no longer decides the arithmetic type, so a double tolerance can be
  passed against a float matrix.
- `expv()` used one Krylov vector fewer than it had built when the Arnoldi
  process broke down early, which was wrong for small systems; it also
  ignored the caller's `m` and left `t` out of the `n == 1` case.
- the reference Fortran EXPOKIT is fetched on demand behind
  `EXPOKIT_BUILD_FORTRAN`, not bundled.
- MPL-2.0, matching Eigen.
