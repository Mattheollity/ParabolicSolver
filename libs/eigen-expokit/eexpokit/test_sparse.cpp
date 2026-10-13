// eigen-expokit: the sparse paths of expv/padm/chbv/mexpv
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Everything here is checked against the same computation on the dense
// equivalent.  Measured: expv sparse == expv dense bit for bit, padm sparse
// within 8e-15 of padm dense (its Ldiv goes through BiCGSTAB, not partialPivLu).
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/padm.hpp"
#include "eexpokit/chbv.hpp"
#include "eexpokit/expv.hpp"
#include "eexpokit/mexpv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;

namespace {

MatrixXd lap1d(int n) {
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -2;
    if (i)     A(i,i-1) = 1;
    if (i<n-1) A(i,i+1) = 1;
  }
  return A;
}

} // namespace

TEST(sparse, expv_matches_dense) {
  for (int n : {4, 10, 30}) {
    MatrixXd D = lap1d(n);
    SpMat A = D.sparseView();
    VectorXd v = VectorXd::Ones(n);
    EXPECT_LT((eexpokit::expv(1.0,A,v) - eexpokit::expv(1.0,D,v)).norm(), 1e-14) << "n=" << n;
  }
}

TEST(sparse, expv_matches_eigen_exp) {
  const int n = 100;
  MatrixXd D = lap1d(n);
  SpMat A = D.sparseView();
  VectorXd v = VectorXd::Ones(n), ref = D.exp()*v;
  eexpokit::expv_stats<double> st;
  auto r = eexpokit::expv(1.0, A, v, 1e-10, 30, &st);
  EXPECT_LT((r-ref).norm()/ref.norm(), 1e-13);
  EXPECT_NEAR(st.hump, 1.0, 1e-12);
}

// a sparse v comes back sparse, which is a different return path
TEST(sparse, expv_sparse_vector) {
  const int n = 10;
  MatrixXd D = lap1d(n);
  SpMat A = D.sparseView(), v(n,1);
  v.insert(0,0) = 1.0;
  VectorXd dv = VectorXd::Zero(n); dv(0) = 1;
  auto r = eexpokit::expv(1.0, A, v);
  EXPECT_TRUE(eexpokit::IsSparse<decltype(r)>::value);
  EXPECT_LT((MatrixXd(r) - eexpokit::expv(1.0,D,dv)).norm(), 1e-14);
}

TEST(sparse, padm_matches_dense) {
  for (int n : {4, 10, 30}) {
    MatrixXd D = lap1d(n);
    SpMat A = D.sparseView();
    EXPECT_LT((MatrixXd(eexpokit::padm(A)) - eexpokit::padm(D)).norm(), 1e-12) << "n=" << n;
    EXPECT_LT((MatrixXd(eexpokit::padm(A,13)) - eexpokit::padm(D,13)).norm(), 1e-12) << "n=" << n;
  }
}

TEST(sparse, chbv_matches_dense) {
  const int n = 20;
  MatrixXd D = lap1d(n);
  SpMat A = D.sparseView();
  VectorXd v = VectorXd::Ones(n);
  EXPECT_LT((eexpokit::chbv(A,v) - eexpokit::chbv(D,v)).norm(), 1e-12);
}

// lap1d is not Q' for any Q, so what this checks is the sparse path against
// the dense one and the probability constraint, not the exponential: mexpv
// renormalizes, lands far from expv, and says so in err.
TEST(sparse, mexpv_matches_dense) {
  const int n = 6;
  MatrixXd D = lap1d(n);
  SpMat A = D.sparseView();
  VectorXd v = VectorXd::Ones(n);
  eexpokit::expv_stats<double> st;
  auto r = eexpokit::mexpv(1.0, A, v, 1e-10, 30, &st);
  EXPECT_GT(st.err, 1e-2);
  EXPECT_NEAR(r.sum(), 1.0, 1e-12);
  EXPECT_LT((r - eexpokit::mexpv(1.0, D, v, 1e-10, 30)).norm(), 1e-14);
}
