// eigen-expokit: padm/chbv/expv with long double
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// One scalar type per translation unit, like the examples: each type costs
// ~0.9GB of compiler memory, mostly the EigenSolver that norm() instantiates.
// Tolerances are the errors measured in test_accuracy.cpp with ~10x margin.
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/padm.hpp"
#include "eexpokit/chbv.hpp"
#include "eexpokit/expv.hpp"

typedef long double Scalar;
typedef Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> Mat;
typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> Vec;

namespace {

const double TOL_P6 = 1e-5, TOL = 1e-10;

/// 1-D Laplacian: symmetric negative definite, the regime chbv wants
Mat lap1d(int n) {
  Mat A = Mat::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = Scalar(-2);
    if (i)     A(i,i-1) = Scalar(1);
    if (i<n-1) A(i,i+1) = Scalar(1);
  }
  return A;
}

double rel(const Mat & a, const Mat & b)
{ return (double)(a-b).cwiseAbs().maxCoeff() / (double)b.cwiseAbs().maxCoeff(); }
double relv(const Vec & a, const Vec & b)
{ return (double)(a-b).norm() / (double)b.norm(); }

} // namespace

TEST(scalars_ldouble, padm_diagonal) {
  const int n = 5;
  Mat A = Mat::Zero(n,n), X = Mat::Zero(n,n);
  for (int i = 0; i < n; i++) { A(i,i) = Scalar(-(i+1.)/3.); X(i,i) = std::exp(A(i,i)); }
  EXPECT_LT(rel(Mat(eexpokit::padm(A)), X), TOL_P6);
  EXPECT_LT(rel(Mat(eexpokit::padm(A,13)), X), TOL);
}

TEST(scalars_ldouble, padm_matches_eigen_exp) {
  Mat A = lap1d(8), X = A.exp();
  EXPECT_LT(rel(Mat(eexpokit::padm(A)), X), TOL_P6);
  EXPECT_LT(rel(Mat(eexpokit::padm(A,13)), X), TOL);
}

TEST(scalars_ldouble, chbv_matches_eigen_exp) {
  Mat A = lap1d(8);
  Vec v = Vec::Ones(8), ref = A.exp()*v;
  EXPECT_LT(relv(Vec(eexpokit::chbv(A,v)), ref), TOL);
}

TEST(scalars_ldouble, expv_matches_eigen_exp) {
  Mat A = lap1d(8);
  Vec v = Vec::Ones(8), ref = A.exp()*v;
  eexpokit::expv_stats<long double> st;
  auto r = eexpokit::expv(1.0, A, v, TOL, 30, &st);
  EXPECT_LT(relv(Vec(r), ref), TOL);
  EXPECT_NEAR((double)st.hump, 1.0, 1e-9);
}
