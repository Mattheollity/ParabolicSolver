// eigen-expokit: expmv tests
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#include "gtest/gtest.h"
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/expmv.hpp"
#include "eexpokit/expv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXcd;
using Eigen::VectorXcd;

namespace {

// diagonal A, so exp(t*A)*v is known exactly
void diagonal(int n, double t) {
  MatrixXd A = MatrixXd::Zero(n,n);
  VectorXd exact(n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -(i + 1.0);
    exact(i) = std::exp(t * -(i + 1.0));
  }
  VectorXd w = eexpokit::expmv(t, A, VectorXd::Ones(n));
  // relative: at t<0 these entries reach e^40 and beyond
  EXPECT_LT((w - exact).norm(), 1e-13 * exact.norm()) << "n=" << n << " t=" << t;
}

// random A, rescaled to the requested inf-norm
MatrixXd scaled_random(int n, double anorm) {
  MatrixXd A = MatrixXd::Random(n,n);
  return A * (anorm / eexpokit::norm(A, Eigen::Infinity));
}

} // namespace

TEST(expmv, diagonal_small) {
  for (int n : {1, 2, 3, 4, 5, 8, 17})
    diagonal(n, 1.0);
}

TEST(expmv, diagonal_large) {
  diagonal(40, 1.0);
  diagonal(70, 1.0);
}

TEST(expmv, respects_t) {
  for (double t : {0.25, 0.5, 2.0, -1.0, -3.5})
    for (int n : {1, 2, 3, 40})
      diagonal(n, t);
}

// nothing is applied at t=0, so v comes back untouched
TEST(expmv, zero_t_returns_v) {
  VectorXd v = VectorXd::Random(9);
  VectorXd w = eexpokit::expmv(0.0, scaled_random(9, 10.0), v);
  EXPECT_EQ((w - v).norm(), 0.0);
}

TEST(expmv, matches_expm) {
  for (int n : {2, 6, 25, 60}) {
    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd v = VectorXd::Random(n);
    VectorXd ref = A.exp() * v;
    EXPECT_LT((eexpokit::expmv(1.0, A, v) - ref).norm(), 1e-12 * ref.norm()) << "n=" << n;
  }
}

// the scaling is the whole point: a wrong s shows up as accuracy collapse
// at large norm, not at small
TEST(expmv, norm_regimes) {
  for (double anorm : {0.1, 1.0, 10.0, 100.0}) {
    for (int n : {2, 6, 25, 60}) {
      MatrixXd A = scaled_random(n, anorm);
      VectorXd v = VectorXd::Random(n);
      VectorXd ref = A.exp() * v;
      EXPECT_LT((eexpokit::expmv(1.0, A, v) - ref).norm(), 1e-12 * ref.norm())
        << "|A|inf=" << anorm << " n=" << n;
    }
  }
}

TEST(expmv, sparse_matches_dense) {
  const int n = 50;
  Eigen::SparseMatrix<double> S(n,n);
  std::vector<Eigen::Triplet<double> > T;
  for (int i = 0; i < n; i++) {
    T.push_back(Eigen::Triplet<double>(i,i,-2.0));
    if (i+1 < n) {
      T.push_back(Eigen::Triplet<double>(i,i+1,1.0));
      T.push_back(Eigen::Triplet<double>(i+1,i,1.0));
    }
  }
  S.setFromTriplets(T.begin(), T.end());
  MatrixXd D(S);
  VectorXd v = VectorXd::Random(n);
  for (double t : {1.0, -2.5}) {
    VectorXd ws = eexpokit::expmv(t, S, v);
    VectorXd ref = (t*D).exp() * v;
    EXPECT_LT((ws - eexpokit::expmv(t, D, v)).norm(), 1e-12 * ref.norm()) << "t=" << t;
    EXPECT_LT((ws - ref).norm(), 1e-12 * ref.norm()) << "t=" << t;
  }
}

TEST(expmv, agrees_with_expv) {
  for (double anorm : {1.0, 10.0}) {
    for (int n : {3, 20, 45}) {
      MatrixXd A = scaled_random(n, anorm);
      VectorXd v = VectorXd::Random(n);
      VectorXd w = eexpokit::expmv(1.0, A, v);
      EXPECT_LT((w - eexpokit::expv(1.0, A, v, 1e-14)).norm(), 1e-9 * w.norm())
        << "|A|inf=" << anorm << " n=" << n;
    }
  }
}

TEST(expmv, complex_scalars) {
  for (int n : {1, 3, 12}) {
    MatrixXcd A = MatrixXcd::Random(n,n) * 3.0;
    VectorXcd v = VectorXcd::Random(n);
    VectorXcd ref = A.exp() * v;
    EXPECT_LT((eexpokit::expmv(1.0, A, v) - ref).norm(), 1e-12 * ref.norm()) << "n=" << n;
  }
}

// v may be a block of vectors, each getting its own exp(t*A)*v
TEST(expmv, several_columns) {
  const int n = 15;
  MatrixXd A = scaled_random(n, 8.0);
  MatrixXd V = MatrixXd::Random(n,3);
  MatrixXd W = eexpokit::expmv(1.0, A, V);
  for (int j = 0; j < V.cols(); j++)
    EXPECT_LT((W.col(j) - eexpokit::expmv(1.0, A, VectorXd(V.col(j)))).norm(), 1e-13) << "j=" << j;
}

TEST(expmv, misshapen_throws) {
  EXPECT_THROW(eexpokit::expmv(1.0, MatrixXd::Zero(3,4), VectorXd::Ones(4)), std::invalid_argument);
  EXPECT_THROW(eexpokit::expmv(1.0, MatrixXd::Zero(3,3), VectorXd::Ones(2)), std::invalid_argument);
}
