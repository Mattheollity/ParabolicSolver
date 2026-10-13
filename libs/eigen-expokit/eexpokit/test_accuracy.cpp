// eigen-expokit: measured accuracy of padm/chbv/expv
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Worst relative error measured over the matrices used below (Apple clang 21,
// arm64, Eigen 5.0.1 and 3.3.7; float figures in parentheses):
//
//   padm p=6    1.4e-8 .. 2.9e-6     (2.8e-7 .. 5.6e-6)
//   padm p=13   1.3e-16 .. 2.2e-15   (0 .. 8.3e-7)
//   chbv        6.9e-12 .. 7.5e-12   (2.5e-6 .. 3.3e-6)
//   expv        3.0e-16 .. 5.0e-16   (9.0e-8 .. 1.1e-7)
//
// padm p=6 is bounded by the Pade truncation and not by the scalar type: 1e-9
// on the negative definite Laplacian, 2.9e-6 on a rotation generator with
// theta=7.  The tolerances asserted here are those figures with ~10x margin.
//
// chbv is only accurate on symmetric negative definite input; see
// chbv_accuracy.off_regime for what it does elsewhere.
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/padm.hpp"
#include "eexpokit/chbv.hpp"
#include "eexpokit/expv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

/// the 1-D Laplacian: symmetric negative definite, the regime chbv wants
MatrixXd lap1d(int n) {
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -2;
    if (i)     A(i,i-1) = 1;
    if (i<n-1) A(i,i+1) = 1;
  }
  return A;
}

MatrixXd rotgen(double th) { MatrixXd A(2,2); A << 0, -th, th, 0; return A; }
MatrixXd rotexp(double th) {
  MatrixXd A(2,2);
  A << std::cos(th), -std::sin(th), std::sin(th), std::cos(th);
  return A;
}

double rel(const MatrixXd & a, const MatrixXd & b)
{ return (a-b).cwiseAbs().maxCoeff() / b.cwiseAbs().maxCoeff(); }

} // namespace

TEST(padm_accuracy, diagonal_closed_form) {
  const int n = 5;
  MatrixXd A = MatrixXd::Zero(n,n), X = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) { A(i,i) = -(i+1.)/3.; X(i,i) = std::exp(A(i,i)); }
  EXPECT_LT(rel(eexpokit::padm(A), X), 1e-6);      // measured 1.4e-7
  EXPECT_LT(rel(eexpokit::padm(A,13), X), 1e-13);  // measured 1.6e-16
}

// the worst padm p=6 case here: measured 2.9e-6 at theta=7, 1.4e-8 at theta=1
TEST(padm_accuracy, rotation_closed_form) {
  for (double th : {0.5, 1.0, 3.0, 7.0}) {
    EXPECT_LT(rel(eexpokit::padm(rotgen(th)), rotexp(th)), 1e-4) << "theta=" << th;
    EXPECT_LT(rel(eexpokit::padm(rotgen(th),13), rotexp(th)), 1e-13) << "theta=" << th;
  }
}

TEST(padm_accuracy, matches_eigen_exp) {
  for (int n : {8, 20}) {
    MatrixXd A = lap1d(n), X = A.exp();
    EXPECT_LT(rel(eexpokit::padm(A), X), 1e-8) << "n=" << n;      // measured 1.0e-9
    EXPECT_LT(rel(eexpokit::padm(A,13), X), 1e-13) << "n=" << n;  // measured 1.8e-15
  }
}

// p is the only accuracy knob padm has, so it had better do something
TEST(padm_accuracy, p13_beats_p6) {
  MatrixXd A = rotgen(7.0), X = rotexp(7.0);
  EXPECT_LT(rel(eexpokit::padm(A,13), X), rel(eexpokit::padm(A), X) / 1e6);
}

TEST(chbv_accuracy, symmetric_negative_definite) {
  for (int n : {8, 20}) {
    MatrixXd A = lap1d(n);
    VectorXd v = VectorXd::Ones(n), ref = A.exp()*v;
    EXPECT_LT((eexpokit::chbv(A,v)-ref).norm()/ref.norm(), 1e-10) << "n=" << n; // measured 7.1e-12
  }
}

// documented caveat: away from the negative real axis chbv is not an
// approximation at all.  100% relative error here, against padm's 1.4e-15.
TEST(chbv_accuracy, off_regime) {
  MatrixXd A = rotgen(20.0);
  VectorXd v = VectorXd::Ones(2), ref = MatrixXd(A.exp())*v;
  EXPECT_GT((eexpokit::chbv(A,v)-ref).norm()/ref.norm(), 0.1);
  EXPECT_LT((eexpokit::padm(A,13)*v-ref).norm()/ref.norm(), 1e-13);

  MatrixXd P(2,2); P << 2, 0, 0, 3;   // positive spectrum: measured 2.4e-4
  VectorXd w = VectorXd::Ones(2), rp = MatrixXd(P.exp())*w;
  EXPECT_GT((eexpokit::chbv(P,w)-rp).norm()/rp.norm(), 1e-6);
}

TEST(expv_accuracy, matches_eigen_exp) {
  for (int n : {8, 20}) {
    MatrixXd A = lap1d(n);
    VectorXd v = VectorXd::Ones(n), ref = A.exp()*v;
    EXPECT_LT((eexpokit::expv(1.0,A,v)-ref).norm()/ref.norm(), 1e-13) << "n=" << n; // measured 3e-16
  }
}

// hump is max||exp(sA)|| over [0,t] normalised by ||v||: 1 for a contraction
TEST(expv_accuracy, hump_of_a_contraction) {
  MatrixXd A = lap1d(8);
  eexpokit::expv_stats<double> st;
  eexpokit::expv(1.0, A, VectorXd::Ones(8), 1e-10, 30, &st);
  EXPECT_NEAR(st.hump, 1.0, 1e-12);
  EXPECT_LT(st.err, 1e-6);
}
