// eigen-expokit: exp(t*A)*v by scaled truncated Taylor series
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Algorithm after Al-Mohy & Higham, SIAM J. Sci. Comput. 33(2):488-511, 2011.
//
#pragma once
#include "eexpokit/expokit.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <Eigen/Core>

namespace eexpokit {

namespace internal {

/// sum_{k>m} x^k/k!, the tail dropped by a degree-m Taylor sum.
template <class Real>
Real taylor_tail(Real x, int m)
{
  Real term(1), sum(0);
  // callers keep x small, so the terms are long dead by the last k
  for (int k = 1; k <= m + 64; k++) {
    term = term * x / Real(k);
    if (k > m)
      sum = sum + term;
  }
  return sum;
}

/// Largest x <= \a cap whose degree-m Taylor tail is within \a tol, bisected.
template <class Real>
Real taylor_theta(int m, Real tol, Real cap)
{
  if (taylor_tail(cap, m) <= tol)
    return cap;
  Real lo(0), hi = cap;
  for (int i = 0; i < 50; i++) {
    Real mid = (lo + hi) / Real(2);
    if (taylor_tail(mid, m) <= tol)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

}

/** \brief w = exp(t*A)*v, by a truncated Taylor series applied s times to v.
  * \a A_ enters only through matrix-vector products.
  *
  * The scaling s and degree m minimize the s*m products needed to pull the
  * Taylor tail at ||t*A||/s under \a tol; each sum stops as soon as its terms
  * stop contributing.  \a v_ may have several columns.
  *
  * Unlike expv(), the accuracy holds up for non-normal \a A_, at the cost of
  * more products once ||t*A|| is large.
  */
template <class tType, class DerivedA, class Derivedv,
          class Scalar = typename Eigen::ScalarBinaryOpTraits<typename DerivedA::Scalar,
                                                              typename Derivedv::Scalar>::ReturnType,
          class Real = typename Eigen::NumTraits<Scalar>::Real>
typename Derivedv::PlainObject
expmv(tType t,
      const Eigen::EigenBase<DerivedA> & A_,
      const Eigen::EigenBase<Derivedv> & v_,
      typename internal::nodeduce<Real>::type tol = Eigen::NumTraits<Real>::epsilon())
{
  const DerivedA & A = A_.derived();
  const Derivedv & v = v_.derived();
  typedef typename Derivedv::PlainObject MatrixV;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> DenseA;
  // not using namespace std, so that ADL finds user scalar overloads
  using std::abs;
  using std::ceil;
  using std::max;
  using std::min;

  // the longest sum worth running, the largest ||t*A||/s the sums are
  // allowed to see, and the most scaling steps to sit through.  Terms of a
  // sum reach exp(mu_max)*||w||, so mu_max also caps the cancellation.
  const int m_max = 40;
  const Real mu_max(4);
  const long s_max = 1L << 20;

  if (A.rows() != A.cols() || A.cols() != v.rows() || A.rows() == 0)
    throw std::invalid_argument("Misshapen arguments to expmv().");

  DenseA w = v;
  // either norm bounds the tail in its own vector norm, so take the smaller
  Real bnorm = Real(abs(t)) * min(norm(A, 1), norm(A, Eigen::Infinity));
  if (!(bnorm > Real(0)))
    return internal::unwrap<MatrixV>(w);

  // Pick the (s,m) of least cost.  The exponents are integers, so keep them
  // out of Scalar arithmetic, as padm() does.
  int m = 0;
  long s = 1, best = std::numeric_limits<long>::max();
  for (int mm = 1; mm <= m_max; mm++) {
    Real theta = internal::taylor_theta(mm, tol, mu_max);
    if (!(theta > Real(0)))
      continue;
    Real ss = ceil(bnorm / theta);
    if (!(ss <= Real(s_max)))
      continue;
    long si = max(1L, (long)ss);
    if (si * mm < best) {
      best = si * mm;
      s = si;
      m = mm;
    }
  }
  if (m == 0)
    throw std::invalid_argument("expmv(): ||t*A|| too large for the requested tolerance.");

  eigen_assert(s >= 1 && m >= 1 && m <= m_max);

  DenseA term, next;
  for (long i = 0; i < s; i++) {
    term = w;
    Real prev(0);
    for (int k = 1; k <= m; k++) {
      // the scalar is applied after the product, not folded into it: Eigen
      // 3.3 cannot fold one into a product whose lhs is a nullary expression
      next = A * term;
      next *= Scalar(t) / Real(s * k);
      term.swap(next);
      w += term;
      // two small terms in a row, in case one alone was an accident
      Real cur = term.norm();
      if (cur + prev <= tol * w.norm())
        break;
      prev = cur;
    }
  }

  return internal::unwrap<MatrixV>(w);
}

}
