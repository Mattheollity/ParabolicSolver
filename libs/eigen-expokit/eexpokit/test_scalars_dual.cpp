// eigen-expokit: padm/chbv with a forward-mode AD scalar (cppduals)
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// duals must be included before anything that pulls in Eigen/Core, or libc++
// resolves hypot() inside std::complex<dual<T>> before cppduals declares it.
//
#include <duals/dual_eigen>
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "eexpokit/padm.hpp"
#include "eexpokit/chbv.hpp"
#include "eexpokit/expv.hpp"

namespace {

typedef duals::dual<double> Dual;
typedef Eigen::Matrix<Dual,Eigen::Dynamic,Eigen::Dynamic> MatrixXdu;
typedef Eigen::Matrix<Dual,Eigen::Dynamic,1> VectorXdu;

/// diag(a_i), each entry carrying dpart 1, so the seeded direction is dA = I
MatrixXdu seeded_diag(int n) {
  MatrixXdu A = MatrixXdu::Zero(n,n);
  for (int i = 0; i < n; i++) A(i,i) = Dual(-(i+1.)/3., 1.);
  return A;
}

/// 1-D Laplacian with dpart 1 on the diagonal: same dA = I direction
MatrixXdu seeded_lap1d(int n) {
  MatrixXdu A = MatrixXdu::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = Dual(-2., 1.);
    if (i)     A(i,i-1) = Dual(1.);
    if (i<n-1) A(i,i+1) = Dual(1.);
  }
  return A;
}

/// d/ds exp(A + sI) at s=0 is exp(A), so every dpart must equal its own rpart
template <class M> void dpart_equals_rpart(const M & E, double tol) {
  for (Eigen::Index j = 0; j < E.cols(); j++)
    for (Eigen::Index i = 0; i < E.rows(); i++)
      EXPECT_NEAR(duals::dpart(E(i,j)), duals::rpart(E(i,j)), tol) << "(" << i << "," << j << ")";
}

} // namespace

// this is the check the padm scaling-exponent regression failed
TEST(dual, padm_derivative_of_diagonal) {
  MatrixXdu A = seeded_diag(5);
  dpart_equals_rpart(MatrixXdu(eexpokit::padm(A)), 1e-5);    // measured 4.1e-7
  dpart_equals_rpart(MatrixXdu(eexpokit::padm(A,13)), 1e-13); // measured 1.1e-16
}

TEST(dual, padm_rpart_is_the_undifferentiated_answer) {
  MatrixXdu E = eexpokit::padm(seeded_diag(5), 13);
  for (int i = 0; i < 5; i++)
    EXPECT_NEAR(duals::rpart(E(i,i)), std::exp(-(i+1.)/3.), 1e-14);
}

// closed form off the diagonal: d/dth [[cos,-sin],[sin,cos]]
TEST(dual, padm_derivative_of_rotation) {
  Dual th(1.0, 1.0);
  MatrixXdu A(2,2);
  A << Dual(0), -th, th, Dual(0);
  MatrixXdu E = eexpokit::padm(A, 13);
  const double c = std::cos(1.0), s = std::sin(1.0);
  EXPECT_NEAR(duals::rpart(E(0,0)),  c, 1e-14);
  EXPECT_NEAR(duals::rpart(E(0,1)), -s, 1e-14);
  EXPECT_NEAR(duals::dpart(E(0,0)), -s, 1e-14);
  EXPECT_NEAR(duals::dpart(E(0,1)), -c, 1e-14);
  EXPECT_NEAR(duals::dpart(E(1,0)),  c, 1e-14);
  EXPECT_NEAR(duals::dpart(E(1,1)), -s, 1e-14);
}

TEST(dual, chbv_derivative) {
  MatrixXdu A = seeded_lap1d(8);
  VectorXdu v = VectorXdu::Ones(8);
  dpart_equals_rpart(VectorXdu(eexpokit::chbv(A,v)), 1e-10);
}

TEST(dual, chbv_agrees_with_padm) {
  MatrixXdu A = seeded_lap1d(8);
  VectorXdu v = VectorXdu::Ones(8);
  VectorXdu y = eexpokit::chbv(A,v), ref = eexpokit::padm(A,13)*v;
  for (int i = 0; i < 8; i++) {
    EXPECT_NEAR(duals::rpart(y(i)), duals::rpart(ref(i)), 1e-10);
    EXPECT_NEAR(duals::dpart(y(i)), duals::dpart(ref(i)), 1e-10);
  }
}

// rparts are right, dparts blow up (1e4 at n=10, NaN at n=40): Eigen's own
// MatrixBase::exp() misbehaves on the dual Hessenberg matrix expv builds.
TEST(dual, DISABLED_expv_derivative) {
  const int n = 12;
  MatrixXdu A = seeded_lap1d(n);
  VectorXdu v = VectorXdu::Ones(n), ref = eexpokit::padm(A,13)*v;
  VectorXdu w = eexpokit::expv(Dual(1.), A, v);
  for (int i = 0; i < n; i++) {
    EXPECT_NEAR(duals::rpart(w(i)), duals::rpart(ref(i)), 1e-12);
    EXPECT_NEAR(duals::dpart(w(i)), duals::dpart(ref(i)), 1e-12);
  }
}
