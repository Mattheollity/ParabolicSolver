// eigen-expokit: mexpv tests
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// The chain here is a birth-death random walk, so its stationary distribution
// is (lambda/mu)^i normalized and Q e = 0 holds exactly.  Measured on it, over
// both supported Eigen versions: |sum(w)-1| <= 5.3e-16, and w within 2.2e-15
// of exp(t*A)*v from Eigen's own .exp() for t <= 10.  Tolerances carry ~50x.
//
#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/expv.hpp"
#include "eexpokit/mexpv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef std::complex<double> cdouble;
typedef Eigen::Matrix<cdouble,Eigen::Dynamic,Eigen::Dynamic> Matc;
typedef Eigen::Matrix<cdouble,Eigen::Dynamic,1> Vecc;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matf;
typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vecf;

namespace {

const double TOL = 1e-13, SUMTOL = 1e-12;
const double UP = 1.0, DOWN = 2.0;

/// transition rate matrix of a birth-death walk: rows sum to zero.  mexpv
/// wants its transpose.
MatrixXd rateQ(int n) {
  MatrixXd Q = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    if (i+1 < n) Q(i,i+1) = UP;
    if (i)       Q(i,i-1) = DOWN;
    Q(i,i) = -Q.row(i).sum();
  }
  return Q;
}

/// the walk's stationary distribution, up to normalization (UP/DOWN)^i
VectorXd stationary(int n) {
  VectorXd pi(n);
  for (int i = 0; i < n; i++)
    pi(i) = std::pow(UP/DOWN, i);
  return pi / pi.sum();
}

MatrixXd lap1d(int n) {
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -2;
    if (i)     A(i,i-1) = 1;
    if (i<n-1) A(i,i+1) = 1;
  }
  return A;
}

VectorXd e1(int n) {
  VectorXd v = VectorXd::Zero(n);
  v(0) = 1;
  return v;
}

double relv(const VectorXd & a, const VectorXd & b) { return (a-b).norm() / b.norm(); }

/// non-negative and summing to 1, which is the whole point of mexpv
void expect_probability(const VectorXd & w, const char * what) {
  EXPECT_GE(w.minCoeff(), 0.0) << what;
  EXPECT_NEAR(w.sum(), 1.0, SUMTOL) << what;
}

} // namespace

TEST(mexpv, is_a_probability_vector) {
  for (int n : {2, 6, 20, 50})
    for (double t : {0.5, 1.0, 10.0, 100.0}) {
      std::string what = "n=" + std::to_string(n) + " t=" + std::to_string(t);
      expect_probability(eexpokit::mexpv(t, MatrixXd(rateQ(n).transpose()), e1(n),
                                         1e-12, 30), what.c_str());
    }
}

TEST(mexpv, matches_dense_expm) {
  for (int n : {2, 6, 20, 50}) {
    MatrixXd A = rateQ(n).transpose();
    for (double t : {0.5, 1.0, 10.0}) {
      VectorXd ref = MatrixXd((t*A).exp()) * e1(n);
      EXPECT_LT(relv(eexpokit::mexpv(t,A,e1(n),1e-12,30), ref), TOL)
        << "n=" << n << " t=" << t;
    }
  }
}

// where the chain is genuine, expv() lands on a probability vector by itself
// and the two must agree
TEST(mexpv, agrees_with_expv) {
  for (int n : {2, 6, 20, 50}) {
    MatrixXd A = rateQ(n).transpose();
    for (double t : {0.5, 1.0, 10.0})
      EXPECT_LT(relv(eexpokit::mexpv(t,A,e1(n),1e-12,30),
                     eexpokit::expv(t,A,e1(n),1e-12,30)), TOL) << "n=" << n << " t=" << t;
  }
}

// long enough and the walk forgets where it started.  The renormalization
// pins w to the simplex, so this beats the dense .exp(), which drifts off it
// by 3e-14 at t=200.
TEST(mexpv, reaches_the_stationary_distribution) {
  for (int n : {4, 6, 10}) {
    MatrixXd A = rateQ(n).transpose();
    VectorXd w = eexpokit::mexpv(200.0, A, e1(n), 1e-12, 30);
    expect_probability(w, "stationary");
    EXPECT_LT(relv(w, stationary(n)), SUMTOL) << "n=" << n;
  }
}

TEST(mexpv, uniform_start_default_tolerance) {
  for (int n : {6, 20}) {
    MatrixXd A = rateQ(n).transpose();
    VectorXd v = VectorXd::Ones(n)/n;
    VectorXd w = eexpokit::mexpv(1.0, A, v);
    expect_probability(w, "uniform");
    EXPECT_LT(relv(w, VectorXd(MatrixXd(A.exp())*v)), TOL) << "n=" << n;
  }
}

// t<0 runs the chain backwards.  While exp(t*A)*v is still non-negative that
// is just the exponential, and mexpv reproduces it.
TEST(mexpv, respects_negative_t) {
  const int n = 6;
  MatrixXd A = rateQ(n).transpose();
  VectorXd v = VectorXd::Ones(n)/n;
  for (double t : {-0.1, -0.25, -0.5}) {
    VectorXd w = eexpokit::mexpv(t, A, v, 1e-12, 30);
    expect_probability(w, "negative t");
    EXPECT_LT(relv(w, VectorXd(MatrixXd((t*A).exp())*v)), TOL) << "t=" << t;
  }
}

// past t=-0.6 or so the reverse-time solution leaves the positive cone, and
// clipping it back is no longer the exponential of anything.  expv() keeps
// sum(w)=1 either way, so the sum is not what tells the two apart.
TEST(mexpv, reverse_time_leaves_the_cone) {
  const int n = 6;
  MatrixXd A = rateQ(n).transpose();
  VectorXd v = VectorXd::Ones(n)/n;
  VectorXd w = eexpokit::mexpv(-1.0, A, v, 1e-12, 30);
  VectorXd e = eexpokit::expv(-1.0, A, v, 1e-12, 30);
  expect_probability(w, "reverse time");
  EXPECT_LT(e.minCoeff(), -0.1);
  EXPECT_GT(relv(w, e), 0.1);
}

// A 1-D Laplacian is not Q' for any Q: its columns sum to -2 at the ends, so
// exp(t*A)*v leaves the simplex and the renormalization pulls it somewhere
// else entirely.  The step's drift off the simplex goes into err, which is the
// only warning the caller gets.  Documenting the divergence, not endorsing it.
TEST(mexpv, breaks_where_its_assumption_does) {
  for (int n : {6, 20}) {
    MatrixXd A = lap1d(n);
    eexpokit::expv_stats<double> st;
    VectorXd w = eexpokit::mexpv(1.0, A, e1(n), 1e-12, 30, &st);
    expect_probability(w, "lap1d");
    VectorXd e = eexpokit::expv(1.0, A, e1(n), 1e-12, 30);
    EXPECT_LT(e.sum(), 0.9) << "n=" << n;
    EXPECT_GT(relv(w, e), 0.1) << "n=" << n;
    EXPECT_GT(st.err, 1e-2) << "n=" << n;
    // where the assumption does hold, err stays at the Krylov floor
    eexpokit::expv_stats<double> ok;
    eexpokit::mexpv(1.0, MatrixXd(rateQ(n).transpose()), e1(n), 1e-12, 30, &ok);
    EXPECT_LT(ok.err, 1e-6) << "n=" << n;
  }
}

TEST(mexpv, stats_are_optional) {
  MatrixXd A = rateQ(10).transpose();
  VectorXd v = e1(10);
  eexpokit::expv_stats<double> st;
  VectorXd plain = eexpokit::mexpv(1.0, A, v);
  VectorXd with = eexpokit::mexpv(1.0, A, v, 1e-10, 30, &st);
  EXPECT_EQ((plain - with).norm(), 0);
  EXPECT_GT(st.err, 0);
  EXPECT_NEAR(st.hump, 1.0, 1e-9);
}

TEST(mexpv, misshapen_throws) {
  EXPECT_THROW(eexpokit::mexpv(1.0, MatrixXd::Zero(3,4), VectorXd::Ones(4)),
               std::invalid_argument);
  EXPECT_THROW(eexpokit::mexpv(1.0, MatrixXd::Zero(3,3), VectorXd::Ones(2)),
               std::invalid_argument);
  // the constraint is on one probability vector, so a block of them is out
  EXPECT_THROW(eexpokit::mexpv(1.0, MatrixXd::Zero(3,3), MatrixXd::Ones(3,2)),
               std::invalid_argument);
}

TEST(mexpv, sparse_matches_dense) {
  for (int n : {6, 20, 50}) {
    MatrixXd D = rateQ(n).transpose();
    SpMat A = D.sparseView();
    VectorXd w = eexpokit::mexpv(1.0, A, e1(n), 1e-12, 30);
    expect_probability(w, "sparse");
    EXPECT_LT(relv(w, eexpokit::mexpv(1.0,D,e1(n),1e-12,30)), 1e-14) << "n=" << n;
  }
}

TEST(mexpv, float_scalar) {
  const int n = 20;
  MatrixXd D = rateQ(n).transpose();
  Matf A = D.cast<float>();
  Vecf v = e1(n).cast<float>();
  // a double tol against a float matrix must compile, and must not drive Real
  Vecf w = eexpokit::mexpv(1.0, A, v, 1e-6, 30);
  EXPECT_TRUE(w.allFinite());
  EXPECT_GE(w.minCoeff(), 0.0f);
  EXPECT_NEAR((double)w.sum(), 1.0, 1e-6);
  EXPECT_LT(relv(w.cast<double>(), eexpokit::mexpv(1.0,D,e1(n),1e-12,30)), 1e-5);
}

// the clipping tests the real part, so a chain with real rates behaves as it
// does in double
TEST(mexpv, complex_scalar) {
  const int n = 20;
  MatrixXd D = rateQ(n).transpose();
  Matc A = D.cast<cdouble>();
  Vecc v = e1(n).cast<cdouble>();
  Vecc w = eexpokit::mexpv(1.0, A, v, 1e-12, 30);
  EXPECT_NEAR(std::abs(w.sum() - cdouble(1,0)), 0.0, SUMTOL);
  VectorXd ref = eexpokit::mexpv(1.0, D, e1(n), 1e-12, 30);
  EXPECT_LT((w - ref.cast<cdouble>()).norm()/ref.norm(), TOL);
}
