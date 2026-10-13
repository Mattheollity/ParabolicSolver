// eigen-expokit: exp(t*A)*v for Markov chains
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
#include <iostream>
#include <limits>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace eexpokit {

/** \brief As expv(), constraining \c w to a probability vector: every step
  * clips the negative components and rescales the sum back to 1.  Assumes
  * <tt>Q e = 0</tt> and that \a A_ is Q transposed; use expv() otherwise.
  *
  * \warning Nothing here checks that assumption.  Where it does not hold the
  * clipping changes the answer and only \c stats->err notices, jumping to
  * O(1); reach for expv() there.
  */
template <class tType, class DerivedA, class Derivedv,
          class Scalar = typename Eigen::ScalarBinaryOpTraits<typename DerivedA::Scalar,
                                                              typename Derivedv::Scalar>::ReturnType,
          class Real = typename Eigen::NumTraits<Scalar>::Real>
typename Derivedv::PlainObject
mexpv(tType t,
      const Eigen::EigenBase<DerivedA> & A_, const Eigen::EigenBase<Derivedv> & v_,
      typename internal::nodeduce<Real>::type tol
        = 1.L/(1ull << (std::numeric_limits<Real>::digits / 2)),
      Eigen::Index m = 30,
      expv_stats<Real> * stats = nullptr)
{
  const DerivedA & A = A_.derived();
  const Derivedv & v = v_.derived();
  typedef typename Derivedv::PlainObject MatrixV;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> DenseV;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> DenseA;
  // not using namespace std, so that ADL finds user scalar overloads
  using std::abs;
  using std::ceil;
  using std::exp;
  using std::floor;
  using std::log10;
  using std::max;
  using std::min;
  using std::pow;
  using std::sqrt;

  Eigen::Index n = A.rows();
  m = max<Eigen::Index>(1, min<Eigen::Index>(n,m));
  if (stats)
    *stats = expv_stats<Real>();

  // check shape.  One column, unlike expv(): the constraint is on a single
  // probability vector.
  if (A.rows() != A.cols() || A.cols() != v.rows() || v.cols() != 1 || n == 0) {
    throw std::invalid_argument("Misshapen arguments to mexpv().");
  }

  if (n == 1) {
    // Q e = 0 makes a 1x1 chain the identity, so v is already the answer
    if (stats)
      stats->hump = 1;
    return v;
  }

  Real anorm = norm(A, Eigen::Infinity);
  int mxrej = 10;   Real btol  = 1.0e-7;
  Real gamma = 0.9; Real delta = 1.2;
  Eigen::Index mb = m;    Real t_out   = abs(t);
  int nstep = 0;    Real t_new   = 0;
  Real t_now = 0;   Real s_error = 0;
  Real rndoff= anorm * std::numeric_limits<Real>::epsilon();
  Real normv = v.norm();

  // The original has no early exit here, but its first step size divides by
  // both of these, so guard the degenerate case.  Only at exactly zero: a
  // relative test against tol throws away small operands that are perfectly
  // representable, and both answers below are exact.
  if (anorm == Real(0) || normv == Real(0)) {
    if (stats)
      stats->hump = 1;
    if (normv == Real(0))
      return 0*v;
    return v;
  }

  //C---  obtain the very first stepsize ...
  int k1 = 2; Real xm = 1./m; Real beta = normv;
  Real fact = pow((m+1.)/exp(1.), (m+1))*sqrt(2*Real(EIGEN_PI)*(m+1));
  t_new = (1/anorm)*pow((fact*tol)/(4*beta*anorm), xm);
  Real s = pow(10., floor(log10(t_new))-1);
  t_new = ceil(t_new/s)*s;
  Real sgn = sign(t);
  nstep = 0;
  Real avnorm(0);

  DenseA w = v;
  Real hump = normv;
  Eigen::Index mx;
  Real err_loc = -1;
  //C---  step-by-step integration ...
  while (t_now < t_out) {
    nstep = nstep + 1;
    Real t_step = min( t_out-t_now,t_new );
    DenseA V(n,m+1);
    V.setZero();
    DenseA H(m+2,m+2);
    H.setZero();

    V.col(0) = (1./beta)*w;
    //C---  Arnoldi loop ...
    for (int j = 0; j < m; j++) {
      DenseV p = A*V.col(j);
      for (int i = 0; i <= j; i++) {
        H(i,j) = (V.col(i).adjoint()*p)(0,0);
        p = p - H(i,j) * V.col(i);
      }
      s = p.norm();
      if (s < btol) {
        k1 = 0;
        mb = j + 1;
        t_step = t_out-t_now;
        break;
      }
      H(j+1,j) = s;
      V.col(j+1) = (1./s)*p;
    }
    if (k1 != 0) {
      H(m+1,m) = 1;
      avnorm = (A*V.col(m)).norm();
    }
    //C---  set 1 for the 2-corrected scheme ...
    int ireject = 0;
    DenseA F;
    //C---  loop while ireject<mxreject until the tolerance is reached ...
    while (ireject <= mxrej) {
      //C---  compute w = beta*V*exp(t_step*H)*e1 ..
      mx = mb + k1;
      F = (sgn*t_step*H.block(0,0,mx,mx)).exp();
      //C---  error estimate ...
      if (k1 == 0) {
        err_loc = btol;
        break;
      }
      else {
        Real phi1 = abs( beta*F(m+0,0) );
        Real phi2 = abs( beta*F(m+1,0) * avnorm );
        if (phi1 > 10*phi2) {
          err_loc = phi2;
          xm = 1./m;
        }
        else if (phi1 > phi2) {
          err_loc = (phi1*phi2)/(phi1-phi2);
          xm = 1./m;
        }
        else {
          err_loc = phi1;
          xm = 1./(m-1.);
        }
      }
      if (err_loc <= delta * t_step*tol) {
        break;
      }
      else {
        //C---  reject the step-size if the error is not acceptable ...
        t_step = gamma * t_step * pow(t_step*tol/err_loc,xm);
        s = pow(10., floor(log10(t_step)) - 1.);
        t_step = ceil(t_step/s) * s;
        if (ireject == mxrej) {
          std::cout << "The requested tolerance is too high: " << tol << "\n";
          std::cout << "err_loc=" << err_loc << " delta=" << delta << "\n";
          if (stats)
            stats->err = std::numeric_limits<Real>::infinity();
          return internal::unwrap<MatrixV>(w);
        }
        ireject = ireject + 1;
      }
    }

    //C---  now update w = beta*V*exp(t_step*H)*e1 ...
    mx = mb + max( 0,k1-1 );
    w = V.block(0,0,V.rows(),mx)*(beta*F.block(0,0,mx,1));

    //C---  Markov model constraints ...
    for (Eigen::Index i = 0; i < n; i++)
      if (Eigen::numext::real(w(i,0)) < Real(0))
        w(i,0) = Scalar(0);
    Real psum = w.cwiseAbs().sum();
    if (!(psum > Real(0))) {
      // nothing left to renormalize: A_ was not a transition rate matrix
      if (stats)
        stats->err = std::numeric_limits<Real>::infinity();
      return internal::unwrap<MatrixV>(w);
    }
    w /= Scalar(psum);
    // how far this step drifted off the simplex, which is an error the plain
    // expv() estimate knows nothing about
    Real roundoff = abs(Real(1) - psum) / Real(n);

    // beta normalizes the next Arnoldi start, so it is taken after the clipping
    beta = w.norm();
    hump = max(hump,beta);

    //C---  update the time covered ...
    t_now = t_now + t_step;
    //C---  suggested value for the next stepsize ...
    t_new = gamma * t_step * pow(t_step*tol/err_loc,xm);
    s = pow(10., floor(log10(t_new))-1);
    t_new = ceil(t_new/s) * s;

    err_loc = max(err_loc,rndoff);
    err_loc = max(err_loc,roundoff);
    s_error = s_error + err_loc;
  }
  if (stats) {
    stats->err = s_error;
    stats->hump = hump / normv;
  }
  return internal::unwrap<MatrixV>(w);
}

}
