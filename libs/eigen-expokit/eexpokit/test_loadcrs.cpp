// eigen-expokit: tests against the CRS matrix shipped with EXPOKIT
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// Only built when EXPOKIT_BUILD_FORTRAN has fetched the EXPOKIT tree.
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "eexpokit/loadcrs.hpp"
#include "eexpokit/expv.hpp"
#include "eexpokit/chbv.hpp"

TEST(loadcrs, c1024) {
  Eigen::SparseMatrix<double> A = eexpokit::loadcrs<double>(EEXPOKIT_CRS_FILE);

  ASSERT_EQ(A.rows(), 1024);
  ASSERT_EQ(A.cols(), 1024);
  EXPECT_EQ(A.nonZeros(), 11264);
  EXPECT_GT(eexpokit::norm(A, Eigen::Infinity), 0);
  EXPECT_TRUE(std::isfinite(eexpokit::norm(A, 1)));

  Eigen::MatrixXd v = Eigen::MatrixXd::Identity(A.rows(), 1);

  eexpokit::expv_stats<double> st;
  auto ret = eexpokit::expv(1., A, v, 1e-10, 30, &st);
  EXPECT_TRUE(ret.allFinite());
  EXPECT_GT(ret.norm(), 0);
  EXPECT_LT(st.err, 1e-10);

  Eigen::MatrixXd y = eexpokit::chbv(A, v);
  EXPECT_TRUE(y.allFinite());
  EXPECT_LT((y - ret).norm(), 1e-8);
}
