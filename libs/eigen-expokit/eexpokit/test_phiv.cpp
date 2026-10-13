// eigen-expokit: phiv tests
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
// The reference is exp(t*[[A,u],[0,0]])*[v;1] built with Eigen's own .exp().
// Measured against it over n <= 60, |t| <= 3 and both supported Eigen versions:
// 2.9e-15 relative in double, 3.7e-7 in float with tol=1e-6.  The tolerances
// here carry ~30x on those.
//
#include "gtest/gtest.h"
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include "eexpokit/expv.hpp"
#include "eexpokit/phiv.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef std::complex<double> cdouble;
typedef Eigen::Matrix<cdouble,Eigen::Dynamic,Eigen::Dynamic> Matc;
typedef Eigen::Matrix<cdouble,Eigen::Dynamic,1> Vecc;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic> Matf;
typedef Eigen::Matrix<float,Eigen::Dynamic,1> Vecf;

namespace {

const double TOL = 1e-13;

MatrixXd lap1d(int n) {
  MatrixXd A = MatrixXd::Zero(n,n);
  for (int i = 0; i < n; i++) {
    A(i,i) = -2;
    if (i)     A(i,i-1) = 1;
    if (i<n-1) A(i,i+1) = 1;
  }
  return A;
}

// exp(t*[[A,u],[0,0]])*[v;1], whose first n rows are what phiv approximates
VectorXd phiref(double t, const MatrixXd & A, const VectorXd & u, const VectorXd & v) {
  Eigen::Index n = A.rows();
  MatrixXd Ah = MatrixXd::Zero(n+1,n+1);
  Ah.topLeftCorner(n,n) = A;
  Ah.col(n).head(n) = u;
  VectorXd vh(n+1);
  vh.head(n) = v;
  vh(n) = 1;
  return (MatrixXd((t*Ah).exp())*vh).head(n);
}

double relv(const VectorXd & a, const VectorXd & b) { return (a-b).norm() / b.norm(); }

} // namespace

TEST(phiv, matches_augmented_expm) {
  for (int n : {2, 6, 25, 60}) {
    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd u = VectorXd::Random(n), v = VectorXd::Random(n);
    for (double t : {0.25, 1.0, 2.0, -1.0})
      EXPECT_LT(relv(eexpokit::phiv(t,A,u,v,1e-12,30), phiref(t,A,u,v)), TOL)
        << "n=" << n << " t=" << t;
  }
}

// the default tol is sqrt(eps), and the Krylov error lands well below it
TEST(phiv, default_tolerance) {
  for (int n : {6, 25, 60}) {
    MatrixXd A = MatrixXd::Random(n,n);
    VectorXd u = VectorXd::Random(n), v = VectorXd::Random(n);
    EXPECT_LT(relv(eexpokit::phiv(1.0,A,u,v), phiref(1.0,A,u,v)), TOL) << "n=" << n;
  }
}

TEST(phiv, zero_u_is_expv) {
  for (int n : {1, 2, 3, 5, 8, 17, 40, 70}) {
    MatrixXd A = lap1d(n);
    VectorXd v = VectorXd::Ones(n);
    VectorXd w = eexpokit::phiv(1.0, A, VectorXd::Zero(n), v, 1e-10, 30);
    EXPECT_LT(relv(w, eexpokit::expv(1.0,A,v,1e-10,30)), TOL) << "n=" << n;
  }
}

// A=0 leaves the augmented operator nilpotent, so this one is exact
TEST(phiv, zero_A_is_v_plus_tu) {
  for (int n : {1, 3, 10}) {
    VectorXd u = VectorXd::Random(n), v = VectorXd::Random(n);
    for (double t : {0.5, 3.0, -2.0})
      EXPECT_LT((eexpokit::phiv(t, MatrixXd::Zero(n,n), u, v) - (v + t*u)).norm(), 1e-14)
        << "n=" << n << " t=" << t;
  }
}

TEST(phiv, diagonal_closed_form) {
  for (int n : {1, 4, 9, 30}) {
    MatrixXd A = MatrixXd::Zero(n,n);
    VectorXd u(n), v(n), exact(n);
    for (int i = 0; i < n; i++) {
      double a = -(i+1.0)/2;
      A(i,i) = a;
      u(i) = 1.0 + i;
      v(i) = 1.0/(i+1.0);
      exact(i) = std::exp(a)*v(i) + u(i)*(std::exp(a)-1)/a;
    }
    EXPECT_LT(relv(eexpokit::phiv(1.0,A,u,v,1e-12,30), exact), TOL) << "n=" << n;
  }
}

TEST(phiv, honours_krylov_dimension) {
  const int n = 60;
  MatrixXd A = MatrixXd::Random(n,n) * 0.1;
  VectorXd u = VectorXd::Ones(n), v = VectorXd::Ones(n);
  VectorXd ref = phiref(1.0,A,u,v);
  for (int m : {5, 15, 30, 50})
    EXPECT_LT(relv(eexpokit::phiv(1.0,A,u,v,1e-10,m), ref), 1e-8) << "m=" << m;
}

// normalizing the operands before building the augmented problem is what keeps
// the answer proportional to (u,v).  The dense reference cannot be used out
// here: at c=1e30 exponentiating the augmented matrix overflows, and phiv is
// the one still right.
TEST(phiv, scales_with_its_operands) {
  const int n = 20;
  MatrixXd A = MatrixXd::Random(n,n);
  VectorXd u = VectorXd::Random(n), v = VectorXd::Random(n);
  VectorXd a = eexpokit::phiv(1.0,A,u,v,1e-12,30);
  for (double c : {1e-30, 1e-6, 1e6, 1e30}) {
    VectorXd b = eexpokit::phiv(1.0, A, VectorXd(c*u), VectorXd(c*v), 1e-12, 30);
    EXPECT_LT(relv(VectorXd(c*a), b), TOL) << "c=" << c;
  }
}

// ||u|| and ||v|| far apart, plus the degenerate ends of that.  Carrying ||v||
// in the last slot of the operand instead of normalizing would have left expv()
// looking at an operand of norm 1e-40 here, and returning zero.
TEST(phiv, lopsided_operands) {
  const int n = 12;
  MatrixXd A = lap1d(n);
  VectorXd o = VectorXd::Ones(n), z = VectorXd::Zero(n);
  for (double c : {1e-40, 1e-20, 1.0, 1e20}) {
    VectorXd v = c*o;
    EXPECT_LT(relv(eexpokit::phiv(1.0,A,o,v,1e-12,30), phiref(1.0,A,o,v)), TOL)
      << "c=" << c;
  }
  EXPECT_LT(relv(eexpokit::phiv(1.0,A,o,z,1e-12,30), phiref(1.0,A,o,z)), TOL);
  EXPECT_LT(relv(eexpokit::phiv(1.0,A,VectorXd(1e-20*o),o,1e-12,30),
                 phiref(1.0,A,VectorXd(1e-20*o),o)), TOL);
  EXPECT_EQ(eexpokit::phiv(1.0,A,z,z).norm(), 0);
}

TEST(phiv, stats_are_optional) {
  MatrixXd A = lap1d(10);
  VectorXd u = VectorXd::Ones(10), v = VectorXd::Ones(10);
  eexpokit::expv_stats<double> st;
  VectorXd plain = eexpokit::phiv(1.0, A, u, v);
  VectorXd with = eexpokit::phiv(1.0, A, u, v, 1e-10, 30, &st);
  EXPECT_EQ((plain - with).norm(), 0);
  EXPECT_GT(st.hump, 0);
}

TEST(phiv, misshapen_throws) {
  VectorXd v3 = VectorXd::Ones(3), v4 = VectorXd::Ones(4);
  EXPECT_THROW(eexpokit::phiv(1.0, MatrixXd::Zero(3,4), v4, v4), std::invalid_argument);
  EXPECT_THROW(eexpokit::phiv(1.0, MatrixXd::Zero(3,3), v3, v4), std::invalid_argument);
  EXPECT_THROW(eexpokit::phiv(1.0, MatrixXd::Zero(3,3), v4, v3), std::invalid_argument);
  EXPECT_THROW(eexpokit::phiv(1.0, MatrixXd::Zero(0,0), v3, v3), std::invalid_argument);
}

// the augmented operator is assembled from triplets for a sparse A, so this
// path never densifies it.  Measured: same answer as dense to 7.4e-17.
TEST(phiv, sparse_matches_dense) {
  typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpRow;
  for (int n : {4, 10, 30, 100}) {
    MatrixXd D = lap1d(n);
    SpMat A = D.sparseView();
    SpRow R = D.sparseView();
    VectorXd u = VectorXd::Ones(n)*0.5, v = VectorXd::Ones(n);
    VectorXd ref = eexpokit::phiv(1.0,D,u,v,1e-12,30);
    EXPECT_LT(relv(eexpokit::phiv(1.0,A,u,v,1e-12,30), ref), 1e-14) << "n=" << n;
    EXPECT_LT(relv(eexpokit::phiv(1.0,R,u,v,1e-12,30), ref), 1e-14) << "n=" << n;
  }
}

// a sparse v comes back sparse, as it does out of expv()
TEST(phiv, sparse_vector) {
  const int n = 10;
  MatrixXd D = lap1d(n);
  SpMat A = D.sparseView(), v(n,1);
  v.insert(0,0) = 1.0;
  VectorXd dv = VectorXd::Zero(n);
  dv(0) = 1;
  VectorXd u = VectorXd::Ones(n)*0.25;
  auto r = eexpokit::phiv(1.0, A, u, v);
  EXPECT_TRUE(eexpokit::IsSparse<decltype(r)>::value);
  EXPECT_LT(relv(VectorXd(MatrixXd(r)), eexpokit::phiv(1.0,D,u,dv)), 1e-14);
}

TEST(phiv, float_scalar) {
  const int n = 8;
  Matf A = lap1d(n).cast<float>();
  Vecf u = Vecf::Ones(n)*0.5f, v = Vecf::Ones(n);
  // a double tol against a float matrix must compile, and must not drive Real
  Vecf w = eexpokit::phiv(1.0, A, u, v, 1e-6, 30);
  EXPECT_TRUE(w.allFinite());
  VectorXd ref = phiref(1.0, lap1d(n), VectorXd::Ones(n)*0.5, VectorXd::Ones(n));
  EXPECT_LT(relv(w.cast<double>(), ref), 1e-5);
}

TEST(phiv, complex_scalar) {
  const int n = 8;
  Matc A = lap1d(n).cast<cdouble>();
  A(0,0) = cdouble(-2, 0.5);
  Vecc u = Vecc::Ones(n)*cdouble(0.5,0), v = Vecc::Ones(n);
  Matc Ah = Matc::Zero(n+1,n+1);
  Ah.topLeftCorner(n,n) = A;
  Ah.col(n).head(n) = u;
  Vecc vh(n+1);
  vh.head(n) = v;
  vh(n) = cdouble(1,0);
  Vecc ref = (Matc(Ah.exp())*vh).head(n);
  Vecc w = eexpokit::phiv(1.0, A, u, v, 1e-12, 30);
  EXPECT_LT((w-ref).norm()/ref.norm(), TOL);
}
