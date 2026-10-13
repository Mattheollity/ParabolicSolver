// eigen-expokit: exp(t*A)*v + t*phi(t*A)*u by Krylov subspace projection
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Algorithm from Sidje, ACM Trans. Math. Softw. 24(1):130-156, 1998.
//
#pragma once
#include "eexpokit/expokit.hpp"
#include "eexpokit/expv.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace eexpokit {

namespace internal {

/// [[A,u],[0,0]], dense.
template <class Scalar, class DerivedA, class Derivedu,
          typename std::enable_if<!IsSparse<DerivedA>::value,int>::type = 0>
Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>
augment(const Eigen::EigenBase<DerivedA> & A_, const Derivedu & u)
{
  const DerivedA & A = A_.derived();
  Eigen::Index n = A.rows();
  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Ah(n+1,n+1);
  Ah.setZero();
  Ah.topLeftCorner(n,n) = A.template cast<Scalar>();
  Ah.col(n).head(n) = u;
  return Ah;
}

/// [[A,u],[0,0]], sparse: nnz(A)+n triplets, never an n-by-n dense pass.  The
/// storage order is A's, since expv() spends its time in products with this.
template <class Scalar, class DerivedA, class Derivedu,
          typename std::enable_if<IsSparse<DerivedA>::value,int>::type = 0>
Eigen::SparseMatrix<Scalar, DerivedA::IsRowMajor ? Eigen::RowMajor : Eigen::ColMajor>
augment(const Eigen::EigenBase<DerivedA> & A_, const Derivedu & u)
{
  typedef typename DerivedA::PlainObject SparseA;
  typedef Eigen::SparseMatrix<Scalar, DerivedA::IsRowMajor ? Eigen::RowMajor
                                                           : Eigen::ColMajor> SparseAh;
  typedef Eigen::Triplet<Scalar> Trip;
  const SparseA & A = A_.derived();
  Eigen::Index n = A.rows();
  std::vector<Trip> trips;
  trips.reserve((size_t)A.nonZeros() + (size_t)n);
  for (Eigen::Index k = 0; k < A.outerSize(); k++)
    for (typename SparseA::InnerIterator it(A,k); it; ++it)
      trips.push_back(Trip(it.row(), it.col(), Scalar(it.value())));
  for (Eigen::Index i = 0; i < n; i++)
    if (u(i) != Scalar(0))
      trips.push_back(Trip(i, n, u(i)));
  SparseAh Ah(n+1,n+1);
  Ah.setFromTriplets(trips.begin(), trips.end());
  return Ah;
}

}

/** \brief w = exp(t*A)*v + t*phi(t*A)*u, phi(z) = (exp(z)-1)/z, ie the solution
  * of <tt>w' = Aw + u</tt>, <tt>w(0) = v</tt>.  This is expv() on the augmented
  * <tt>[[A,u],[0,0]]</tt> and <tt>[v;1]</tt>, so it costs one copy of \a A_ plus
  * n, and \a stats describes that augmented problem.
  */
template <class tType, class DerivedA, class Derivedu, class Derivedv,
          class Scalar = typename Eigen::ScalarBinaryOpTraits<
            typename DerivedA::Scalar,
            typename Eigen::ScalarBinaryOpTraits<typename Derivedu::Scalar,
                                                 typename Derivedv::Scalar>::ReturnType>::ReturnType,
          class Real = typename Eigen::NumTraits<Scalar>::Real>
typename Derivedv::PlainObject
phiv(tType t,
     const Eigen::EigenBase<DerivedA> & A_,
     const Eigen::EigenBase<Derivedu> & u_,
     const Eigen::EigenBase<Derivedv> & v_,
     typename internal::nodeduce<Real>::type tol
       = 1.L/(1ull << (std::numeric_limits<Real>::digits / 2)),
     Eigen::Index m = 30,
     expv_stats<Real> * stats = nullptr)
{
  const DerivedA & A = A_.derived();
  typedef typename Derivedv::PlainObject MatrixV;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> DenseV;
  // not using namespace std, so that ADL finds user scalar overloads
  using std::max;

  Eigen::Index n = A.rows();
  if (A.rows() != A.cols() || A.cols() != v_.rows()
      || u_.rows() != v_.rows() || u_.cols() != 1 || v_.cols() != 1 || n == 0) {
    throw std::invalid_argument("Misshapen arguments to phiv().");
  }

  DenseV u = u_.derived().template cast<Scalar>();
  DenseV v = v_.derived().template cast<Scalar>();

  // w is linear in (u,v) together, so both are normalized before the augmented
  // problem is built and the answer is scaled back.  That keeps the u column
  // from adding stiffness that is not in A, keeps ||[v;1]|| >= 1 so expv() can
  // never mistake the operand for zero, and makes w proportional to (u,v) the
  // way the exact answer is.
  Real scale = max(u.norm(), v.norm());
  if (!(scale > Real(0)))
    scale = Real(1);

  DenseV vh(n+1);
  vh.head(n) = v / Scalar(scale);
  vh(n) = Scalar(1);
  DenseV uh = u / Scalar(scale);

  DenseV wh = expv(t, internal::augment<Scalar>(A_, uh), vh, tol, m, stats);
  DenseV w = Scalar(scale) * wh.head(n);
  return internal::unwrap<MatrixV>(w);
}

}
