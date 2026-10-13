// eigen-expokit: tests
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#include "gtest/gtest.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "eexpokit/expokit.hpp"

#define _EXPECT_TRUE(...)  {typedef __VA_ARGS__ tru; EXPECT_TRUE(tru::value); }
#define _EXPECT_FALSE(...) {typedef __VA_ARGS__ fal; EXPECT_FALSE(fal::value); }

TEST(templates, IsSparse) {
  Eigen::MatrixXd dense;
  Eigen::VectorXd dvec;
  Eigen::SparseMatrix<double> sparse(1,1);

  using eexpokit::IsSparse;

  _EXPECT_TRUE(IsSparse<Eigen::SparseMatrix<double>>);
  _EXPECT_TRUE(IsSparse<decltype(sparse)>);
  _EXPECT_FALSE(IsSparse<decltype(dense)>);
  _EXPECT_FALSE(IsSparse<decltype(dvec)>);

  _EXPECT_TRUE(IsSparse<decltype(sparse * sparse)>);
  _EXPECT_FALSE(IsSparse<decltype(dense * sparse)>);
  _EXPECT_FALSE(IsSparse<decltype(sparse * dense)>);
  _EXPECT_FALSE(IsSparse<decltype(dense * dense)>);

  _EXPECT_TRUE(IsSparse<decltype(dvec * sparse)>);
  _EXPECT_FALSE(IsSparse<decltype(sparse * dvec)>);
  _EXPECT_FALSE(IsSparse<decltype(dvec * dvec)>);
}

TEST(types, sparse_sparse) {
  //Eigen::SparseMatrix<>
}

TEST(norm, various) {
  int N = 5;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(N,N);
  Eigen::SparseMatrix<double> a = A.sparseView();
  double tol = 1e-12;
  using eexpokit::norm;

  EXPECT_NEAR(norm(A,1), norm(A.transpose(), Eigen::Infinity), tol);
  EXPECT_NEAR(norm(A,1), norm(a,1), tol);
  EXPECT_NEAR(norm(A,Eigen::Infinity), norm(a.transpose(),1), tol);
  EXPECT_NEAR(norm(A,Eigen::Infinity), norm(a,Eigen::Infinity), tol);
}

TEST(Ldiv, types) {
  using eexpokit::Ldiv;

  Eigen::SparseMatrix<double> s1(1,1), s2(1,1), rs(1,1);
  Eigen::MatrixXd d1(1,1), d2(1,1), rd(1,1);

  // sparse-sparse
  rs = Ldiv(s1, s2);

  // sparse-dense
  rd = Ldiv(s1, d2);

  // dense-dense
  rd = Ldiv(d1, d2);

  // dense-sparse
  //rd = Ldiv(d1, s1);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  std::cout.precision(20);
  std::cerr.precision(20);
  return RUN_ALL_TESTS();
}
