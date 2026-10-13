// eigen-expokit: exp(A) by irreducible Pade approximation, scaling and squaring
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Algorithm from Sidje, ACM Trans. Math. Softw. 24(1):130-156, 1998.
//
#pragma once
#include "eexpokit/expokit.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>

// TODO, FIX this A:
//A=[-0.444450578393624, -0.0452058962756795;
//   0.107939911590861, 0.257741849523849]

namespace eexpokit {

/** \brief exp(A), by the irreducible (p,p)-degree Pade approximation with
  * scaling and squaring.  \a p defaults to 6.
  */
template <class DerivedA>
typename DerivedA::PlainObject
padm(const Eigen::EigenBase<DerivedA> & A_, int p = 6)
{
  const DerivedA & Au = A_.derived();
  typedef typename DerivedA::PlainObject::RealScalar Real;
  typedef typename DerivedA::Scalar Scalar;
  typedef typename DerivedA::PlainObject MatrixA;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorA;
  // not using namespace std, so that ADL finds user scalar overloads
  using std::log;
  using std::max;
  using std::pow;
  using std::trunc;

  eigen_assert(Au.rows() == Au.cols());
  int n = Au.rows();

  // Pade coefficients (0-based instead of 0-based as in the literature)
  VectorA c(p+1,1);
  c(0) = 1;
  for (int k = 0; k < p; k++)
    c(k+1) = c(k)*Real((p+1.-(k+1.))/((k+1.)*(2.*p+1.-(k+1.))));

  // Scaling.  The exponent is an integer, so keep it out of Scalar
  // arithmetic: mixing it in breaks scalar types whose Real is not a
  // builtin (eg forward-mode AD types, which would also pick up a
  // spurious derivative from the norm).
  MatrixA A; // scaled
  int s = 0;
  Real nrm = norm(Au, Eigen::Infinity);
  if (nrm > Real(0.5)) {
    s = (int)max(Real(0), trunc(log(nrm)/log(Real(2))) + Real(2));
    A = Real(pow(2., -s)) * Au;
  }
  else
    A = Au;

  // Horner evaluation of the irreducible fraction (see ref. above)

  MatrixA I_(n,n);
  I_.setIdentity();
  MatrixA A2 = A*A;
  MatrixA Q = c(p)*I_;
  MatrixA P = c(p-1)*I_;
  MatrixA E;
  int odd = 1;
  for (int k = p-1; k >= 0; k--) {
    if (odd == 1)
      Q = Q*A2 + c(k)*I_;
    else
      P = P*A2 + c(k)*I_;
    odd = 1-odd;
  }
  if (odd == 1) {
    Q = Q*A;
    Q = Q - P;
    E = -(I_ + 2*Ldiv(Q,P));
  }
  else {
    P = P*A;
    Q = Q - P;
    E = I_ + 2*Ldiv(Q,P);
  }

  // Squaring
  for (int k = 0; k < s; k++)
    E = E*E;

  return E;
}

}
