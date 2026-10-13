// eigen-expokit: expv tests
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/expv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

// diagonal A, so exp(t*A)*v is known exactly
void diagonal(int n, double t) {
  MatrixXd A = MatrixXd::Zero(n,n);
  VectorXd exact(n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -(i + 1.0);
    exact(i) = std::exp(t * -(i + 1.0));
  }
  auto r = eexpokit::expv(t, A, VectorXd::Ones(n));
  EXPECT_LT((r - exact).norm(), 1e-12) << "n=" << n << " t=" << t;
}

} // namespace

// the Krylov space spans everything once n is small, so these are exact.
// Getting one vector too few here used to make n=2 return exp(-1.5) twice.
TEST(expv, diagonal_small) {
  for (int n : {1, 2, 3, 4, 5, 8, 17})
    diagonal(n, 1.0);
}

TEST(expv, diagonal_large) {
  diagonal(40, 1.0);
  diagonal(70, 1.0);
}

TEST(expv, respects_t) {
  for (double t : {0.25, 0.5, 2.0})
    for (int n : {1, 3, 40})
      diagonal(n, t);
}

// m is the Krylov dimension the caller asked for, and used to be ignored
TEST(expv, honours_krylov_dimension) {
  const int n = 60;
  MatrixXd A = MatrixXd::Random(n,n) * 0.1;
  VectorXd v = VectorXd::Ones(n);
  VectorXd ref = A.exp() * v;
  for (int m : {5, 15, 30, 50})
    EXPECT_LT((eexpokit::expv(1.0, A, v, 1e-10, m) - ref).norm(), 1e-8) << "m=" << m;
}

TEST(expv, matches_expm) {
  for (int n : {2, 6, 25}) {
    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd v = VectorXd::Random(n);
    VectorXd ref = A.exp() * v;
    EXPECT_LT((eexpokit::expv(1.0, A, v) - ref).norm(), 1e-10 * ref.norm()) << "n=" << n;
  }
}

TEST(expv, misshapen_throws) {
  EXPECT_THROW(eexpokit::expv(1.0, MatrixXd::Zero(3,4), VectorXd::Ones(4)), std::invalid_argument);
  EXPECT_THROW(eexpokit::expv(1.0, MatrixXd::Zero(3,3), VectorXd::Ones(2)), std::invalid_argument);
}

// the original has no small-norm early exit; ours guards only exact zero, so a
// small but representable operand still gets integrated rather than truncated
TEST(expv, small_operands_are_not_truncated) {
  const int n = 6;
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) A(i,i) = -(i + 1.0);
  VectorXd tiny = VectorXd::Ones(n) * 1e-20;
  VectorXd w = eexpokit::expv(1.0, A, tiny);
  EXPECT_GT(w.norm(), 0);
  for (int i = 0; i < n; i++)
    EXPECT_NEAR(w(i) / 1e-20, std::exp(-(i + 1.0)), 1e-12);

  // an A whose norm is under any sane tol still exponentiates to about I
  MatrixXd wee = MatrixXd::Identity(n,n) * 1e-30;
  VectorXd v = VectorXd::Ones(n);
  EXPECT_LT((eexpokit::expv(1.0, wee, v) - v).norm(), 1e-12);
}

TEST(expv, exact_zero_is_exact) {
  const int n = 4;
  VectorXd v = VectorXd::Random(n);
  EXPECT_EQ((eexpokit::expv(1.0, MatrixXd::Zero(n,n), v) - v).norm(), 0);
  EXPECT_EQ(eexpokit::expv(1.0, MatrixXd::Identity(n,n), VectorXd::Zero(n)).norm(), 0);
}
