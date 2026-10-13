// eigen-expokit: shared helpers: norm(), Ldiv(), expv_stats
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#pragma once
#include <iostream>
#include <type_traits>
#include <utility>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>

namespace eexpokit {

/// argument for the norm() function to request a fro-norm.
const static int Frobenius = -2;

/// optional diagnostics from expv(), mexpv() and phiv().
template <class Real>
struct expv_stats {
  /// error estimate.  Not a bound, see expv().
  Real err = 0;
  /// max||exp(sA)|| over the integration, 1 when well conditioned.
  Real hump = 0;
};

/// sign(um) function
template <typename T> int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}

/// Check if a matrix is sparse or not
template <typename T, typename S = decltype(std::declval<T>().makeCompressed())> static std::true_type sparseTest(int);
template <typename T> static std::false_type sparseTest(...);
template <typename T, typename U = void> struct IsSparse : decltype(sparseTest<typename T::PlainObject>(0)) {};

/// Sparse matrices norm_1()
template <typename Derived>
typename Derived::RealScalar
norm_1(const Eigen::SparseMatrixBase<Derived> & A)
{
  typedef typename Derived::RealScalar Real;
  typedef typename Eigen::Matrix<Real,1,Eigen::Dynamic> RowVectorX;
  return (RowVectorX::Ones(A.rows()) * A.cwiseAbs()).maxCoeff();
}

/// Sparse matrices norm()
template <typename Derived>
typename Derived::RealScalar
norm_inf(const Eigen::SparseMatrixBase<Derived> & A)
{
  typedef typename Derived::RealScalar Real;
  typedef typename Eigen::Matrix<Real,Eigen::Dynamic,1> VectorX;
  return (A.cwiseAbs() * VectorX::Ones(A.cols())).maxCoeff();
}

/** \brief L-norm of a sparse matrix or vector, \a L being 1, 2,
  * Eigen::Infinity or Frobenius.  Entrywise for vectors, induced for
  * matrices; \c L==2 on a sparse matrix throws std::invalid_argument.
  */
template <typename Derived>
typename Derived::RealScalar
norm(const Eigen::SparseMatrixBase<Derived> & A, int L = 2)
{
  typedef typename Derived::RealScalar Real;
  // not using namespace std, so that ADL finds user scalar overloads
  using std::abs;
  using std::pow;

  // vector norm
  if (A.rows() == 1 || A.cols() == 1 || L == Frobenius) {
    if (L == Frobenius || L == 2)
      return A.norm();
    Real v(0);
    for (Eigen::Index k = 0; k < A.outerSize(); ++k)
      for (typename Derived::InnerIterator it(A.derived(),k); it; ++it)
        //for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(A,k); it; ++it)
        v += pow(abs(it.value()), L);
    return pow(v, 1./L);
  }

  // matrix norm
  if (L == 2) {
    // no.
    throw std::invalid_argument("l2 matrix norm on sparse matrices not available.");
  }
  else if (L == 1) {
    // L1
    return norm_1(A);
  }
  else if (L == Eigen::Infinity) {
    // inf
    return norm_inf(A);
  }

  std::cout << "Unimplemented Sparse Norm, L=" << L << "\n";
  // punt, give them something
  return A.norm();
}

/// Dense matrix norm_inf()
template <typename DerivedA>
typename DerivedA::RealScalar
norm_1(const Eigen::MatrixBase<DerivedA> & A_)
{
  const DerivedA & A = A_.derived();
  return A.cwiseAbs().colwise().sum().maxCoeff();
}

/// Dense matrix norm_inf()
template <typename DerivedA>
typename DerivedA::RealScalar
norm_inf(const Eigen::MatrixBase<DerivedA> & A_)
{
  const DerivedA & A = A_.derived();
  return A.cwiseAbs().rowwise().sum().maxCoeff();
}

/** \brief L-norm of a dense matrix or vector, \a L being 1, 2,
  * Eigen::Infinity or Frobenius.  Entrywise for vectors, induced for
  * matrices, \c L==2 being the largest singular value.
  */
template <typename DerivedA>
typename DerivedA::RealScalar
norm(const Eigen::MatrixBase<DerivedA> & A_, int L = 2)
{
  const DerivedA & A = A_.derived();
  using std::pow;

  // vector norm
  if (A.rows() == 1 || A.cols() == 1 || L == Frobenius) {
    if (L == Frobenius || L == 2)
      return A.norm();
    return pow(A.array().abs().pow(L).sum(),1./L);
  }

  // matrix norm.  This branch instantiates an eigensolver even when the
  // caller only ever asks for L==1 or Infinity, and dominates the compile
  // time of everything that calls norm().
  if (L == 2) {
    typedef typename std::conditional<Eigen::NumTraits<typename DerivedA::Scalar>::IsComplex,
                                      Eigen::ComplexEigenSolver<typename DerivedA::PlainObject>,
                                      Eigen::EigenSolver<typename DerivedA::PlainObject> >::type
      EigenSolverType;
    if (A.rows() == A.cols()) {
      EigenSolverType ed(A, false);
      return ed.eigenvalues().array().abs().maxCoeff();
    }
    else {
      EigenSolverType ed(A.transpose() * A, false);
      return ed.eigenvalues().array().abs().sqrt().maxCoeff();
    }
  }
  else if (L == 1) {
    // L1
    return norm_1(A_);
  }
  else if (L == Eigen::Infinity) {
    // inf
    return norm_inf(A_);
  }

  // punt
  throw std::invalid_argument("TODO.");
  return pow(A.array().abs().pow(L).sum(),1./L);
}

/** \brief Solve <tt>D X = N</tt> for X, sparse \a D_, iteratively (BiCGSTAB).
  * On decomposition failure it reports on std::cout and returns \a D_.
  */
template <class DerivedD, class DerivedN,
          typename std::enable_if<IsSparse<DerivedD>::value,int>::type = 0>
typename DerivedN::PlainObject
Ldiv(const Eigen::EigenBase<DerivedD> & D_, const DerivedN & N)
{
  const DerivedD & D = D_.derived();
  typedef typename Eigen::SparseMatrixBase<DerivedD>::PlainObject MatrixT;

  // Pick one
  //typedef Eigen::SparseLU<MatrixC> SolverType;
  typedef Eigen::BiCGSTAB<MatrixT> SolverType; // *
  //typedef Eigen::LeastSquaresConjugateGradient<MatrixC> SolverType;
  //typedef Eigen::ConjugateGradient<MatrixC> SolverType;
  //typedef Eigen::SimplicialLDLT<MatrixC> SolverType;
  //typedef Eigen::SparseQR<MatrixC, Eigen::COLAMDOrdering<int> > SolverType;

  SolverType solver;
  solver.compute(D);
  if (solver.info() != Eigen::Success) {
    std::cout << "decomposition faild\n";
    return D;
  }

  auto R = solver.solve(N);

  if(solver.info() != Eigen::Success) {
    // solving failed
    std::cout << "solving faild\n";
  }
  return R;
}

/** \brief Solve <tt>D X = N</tt> for X, dense \a D_.
  */
template <typename DerivedD, typename DerivedN,
          typename std::enable_if<!IsSparse<DerivedD>::value,int>::type = 0>
typename DerivedN::PlainObject
Ldiv(const Eigen::EigenBase<DerivedD> & D_, const Eigen::EigenBase<DerivedN> & N_)
{
  const DerivedD & D = D_.derived();
  return D.derived().partialPivLu().solve(N_.derived());
}

namespace internal {
/// keeps its argument out of template deduction
template <class T> struct nodeduce { typedef T type; };

template <class D, class C> inline D unwrap(const Eigen::MatrixBase<C> & c) {
  return c.sparseView().template cast<typename D::Scalar>();
}
template <class D, class C> inline D unwrap(const Eigen::SparseMatrixBase<C> & c) {
  return c.template cast<typename D::Scalar>();;
}

template <class D, class C,
          typename std::enable_if<Eigen::NumTraits<typename D::Scalar>::IsComplex,int>::type = 0>
inline D maybe_real(const C & c) { return c; }
template <class D, class C,
          typename std::enable_if<!Eigen::NumTraits<typename D::Scalar>::IsComplex,int>::type = 0>
inline D maybe_real(const C & c) { return c.real(); }

}

}
