// eigen-expokit: exp(H)*x by uniform rational Chebyshev approximation
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Algorithm from Sidje, ACM Trans. Math. Softw. 24(1):130-156, 1998.
//
#pragma once
#include "eexpokit/expokit.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>

#include <type_traits>

namespace eexpokit {

/** \brief exp(H)*x, by the type (14,14) uniform rational Chebyshev
  * approximation.  About 14 digits when \a H_ is symmetric negative definite,
  * possibly poor otherwise.
  */
template <class DerivedH, class Derivedx>
typename Derivedx::PlainObject
chbv(const Eigen::EigenBase<DerivedH> & H_, const Eigen::EigenBase<Derivedx> & x_)
{
  const DerivedH & H = H_.derived();
  const Derivedx & x = x_.derived();
  typedef typename Derivedx::PlainObject::RealScalar Real;
  typedef std::complex<Real> Complex;
  typedef typename std::conditional<IsSparse<DerivedH>::value,
                                    Eigen::SparseMatrix<Complex>,
                                    Eigen::Matrix<Complex,DerivedH::RowsAtCompileTime,
                                                  DerivedH::ColsAtCompileTime>>::type MatrixCH;
  typedef typename std::conditional<IsSparse<Derivedx>::value,
                                    Eigen::SparseMatrix<Complex>,
                                    Eigen::Matrix<Complex, Derivedx::RowsAtCompileTime,
                                                  Derivedx::ColsAtCompileTime>>::type MatrixCx;
  Eigen::SparseMatrix<Complex> SparseC;

  eigen_assert(H.rows() == H.cols());
  eigen_assert(H.rows() == x.rows());


  Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic> alpha(8,1);
  Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic> theta(8,1);

  // Coefficients and poles of the partial fraction expansion ...
  std::complex<double> i(0,1);
  alpha(0) =  0.183216998528140087E-11;
  alpha(1) =  0.557503973136501826E+02 - i*0.204295038779771857E+03;
  alpha(2) = -0.938666838877006739E+02 + i*0.912874896775456363E+02;
  alpha(3) =  0.469965415550370835E+02 - i*0.116167609985818103E+02;
  alpha(4) = -0.961424200626061065E+01 - i*0.264195613880262669E+01;
  alpha(5) =  0.752722063978321642E+00 + i*0.670367365566377770E+00;
  alpha(6) = -0.188781253158648576E-01 - i*0.343696176445802414E-01;
  alpha(7) =  0.143086431411801849E-03 + i*0.287221133228814096E-03;

  theta(1) = -0.562314417475317895E+01 + i*0.119406921611247440E+01;
  theta(2) = -0.508934679728216110E+01 + i*0.358882439228376881E+01;
  theta(3) = -0.399337136365302569E+01 + i*0.600483209099604664E+01;
  theta(4) = -0.226978543095856366E+01 + i*0.846173881758693369E+01;
  theta(5) =  0.208756929753827868E+00 + i*0.109912615662209418E+02;
  theta(6) =  0.370327340957595652E+01 + i*0.136563731924991884E+02;
  theta(7) =  0.889777151877331107E+01 + i*0.166309842834712071E+02;

  int p = 7;
  theta = -theta;
  alpha = -alpha;

  Eigen::Index n = H.rows();
  MatrixCH I_(n,n);
  I_.setIdentity();
  MatrixCx y(n,1);
  y = alpha(0) * x;

  bool both_real =
    !Eigen::NumTraits<typename DerivedH::Scalar>::IsComplex &&
    !Eigen::NumTraits<typename Derivedx::Scalar>::IsComplex;

  if (both_real) {
    // both are real
    for (int i = 1; i <= p; i++) {
      //y = y + (H-theta(i)*I_) \ (alpha(i)*x);
      y = y + Ldiv(H.template cast<Complex>() - theta(i)*I_, alpha(i)*x);
    }
    y = y.real().template cast<Complex>();
  }
  else {
    // H or x or both are complex
    theta.conservativeResize(16,1);
    alpha.conservativeResize(16,1);
    theta.block(8,0,7,1) = theta.block(1,0,7,1).conjugate();
    alpha.block(8,0,7,1) = alpha.block(1,0,7,1).conjugate();
    alpha.array() *= Real(0.5);
    for (int i = 1; i <= 2*p; i++) {
      //y = y + (H-theta(i)*I_) \ (alpha(i)*x);
      y = y + Ldiv(H.template cast<Complex>() - theta(i)*I_, alpha(i)*x);
    }
  }
  return internal::maybe_real<typename Derivedx::PlainObject>(y);
}

/** \brief exp(t*H)*x, chbv() with the time scaling folded in.  The scaling
  * costs a fraction of a percent of the solves it feeds.
  */
template <class tType, class DerivedH, class Derivedx>
typename Derivedx::PlainObject
chbv(tType t, const Eigen::EigenBase<DerivedH> & H_, const Eigen::EigenBase<Derivedx> & x_)
{
  typedef typename DerivedH::Scalar Scalar;
  return chbv((Scalar(t) * H_.derived()).eval(), x_);
}

}
