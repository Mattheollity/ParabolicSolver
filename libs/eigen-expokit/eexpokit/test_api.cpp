// eigen-expokit: the shape of the public API
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "eexpokit/chbv.hpp"
#include "eexpokit/expmv.hpp"
#include "eexpokit/expv.hpp"
#include "eexpokit/mexpv.hpp"
#include "eexpokit/padm.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace {

// symmetric negative definite, which is chbv's regime
MatrixXd snd(int n) {
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -2;
    if (i) { A(i,i-1) = 1; A(i-1,i) = 1; }
  }
  return A;
}

} // namespace

TEST(api, chbv_t_equals_scaled_h) {
  MatrixXd H = snd(6);
  VectorXd x = VectorXd::Ones(6);
  for (double t : {0.5, 1.0, 2.0})
    EXPECT_LT((eexpokit::chbv(t, H, x) - eexpokit::chbv(MatrixXd(t*H), x)).norm(), 1e-14) << "t=" << t;
}

// the diagnostics are opt-in, and asking for them must not change the answer
TEST(api, stats_are_optional) {
  MatrixXd A = snd(10);
  VectorXd v = VectorXd::Ones(10);
  eexpokit::expv_stats<double> st;
  VectorXd plain = eexpokit::expv(1.0, A, v);
  VectorXd with = eexpokit::expv(1.0, A, v, 1e-10, 30, &st);
  EXPECT_EQ((plain - with).norm(), 0);
  EXPECT_GT(st.hump, 0);
}

// a double tolerance against a float matrix used to fail to compile: tol drove
// the internal real type instead of the scalar doing it
TEST(api, tolerance_does_not_pick_the_scalar_type) {
  typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matf;
  typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vecf;
  Matf A = snd(6).cast<float>();
  Vecf v = Vecf::Ones(6);
  eexpokit::expv_stats<float> st;
  EXPECT_TRUE(eexpokit::expv(1.0, A, v, 1e-6, 30, &st).allFinite());
  EXPECT_TRUE(eexpokit::expmv(1.0, A, v, 1e-6).allFinite());
  EXPECT_TRUE(eexpokit::chbv(1.0, A, v).allFinite());
}

// all four actions read (t, A, v), and all of them return the vector
TEST(api, actions_agree) {
  MatrixXd A = snd(12);
  VectorXd v = VectorXd::Ones(12);
  VectorXd a = eexpokit::expv(1.0, A, v);
  VectorXd b = eexpokit::expmv(1.0, A, v);
  VectorXd c = eexpokit::chbv(1.0, A, v);
  VectorXd d = eexpokit::padm(A, 13) * v;
  EXPECT_LT((a-d).norm(), 1e-12);
  EXPECT_LT((b-d).norm(), 1e-12);
  EXPECT_LT((c-d).norm(), 1e-10);
}
