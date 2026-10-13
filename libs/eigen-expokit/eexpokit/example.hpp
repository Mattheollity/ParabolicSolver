// eigen-expokit: example driver, one scalar type per translation unit
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#pragma once
#include <chrono>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>

#include "expv.hpp"
#include "mexpv.hpp"
#include "chbv.hpp"
#include "padm.hpp"

using namespace eexpokit;

class Timer
{
public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_>
      (clock_::now() - beg_).count(); }

private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
};

/// The negative 5-point Laplacian on a \a g x \a g grid: symmetric negative
/// definite, the regime chbv() is accurate in.
template <class T>
Eigen::SparseMatrix<T> laplacian2d(int g)
{
  Eigen::SparseMatrix<T> A(g*g, g*g);
  std::vector<Eigen::Triplet<T>> tr;
  tr.reserve(5*g*g);
  for (int j = 0; j < g; j++)
    for (int i = 0; i < g; i++) {
      int k = j*g + i;
      tr.emplace_back(k, k, T(-4));
      if (i > 0)   tr.emplace_back(k, k-1, T(1));
      if (i < g-1) tr.emplace_back(k, k+1, T(1));
      if (j > 0)   tr.emplace_back(k, k-g, T(1));
      if (j < g-1) tr.emplace_back(k, k+g, T(1));
    }
  A.setFromTriplets(tr.begin(), tr.end());
  return A;
}

template <class T>
void demo()
{
  typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> MatrixX;
  Timer timer;

  std::cout << "Building the matrix ...\n";
  auto A = laplacian2d<T>(32);
  A = A; //  + T(0,2) * A;

  std::cout << "A(" << A.rows() << "," << A.cols() << ")\n";
  std::cout << "Computing w = exp(A)e_1 ...\n";
  auto n = A.rows();
  auto v = MatrixX::Identity(n,1);

  std::cout << "A(1:5) =" << (A.block(0,0,5,1)).transpose() << "\n";

  std::cout << " *** expv ***\n";
  timer.reset();
  expv_stats<typename MatrixX::RealScalar> st;
  auto ret = expv(1,A,v,1e-10,30,&st);
  std::cout << "time   =" << timer.elapsed() << " ";
  std::cout << "err    =" << st.err << "\n";

  std::cout << " *** chbv ***\n";
  timer.reset();
  auto cy = chbv(A,v);
  std::cout << "time   =" << timer.elapsed() << " ";
  //std::cout << "y(1:5) =" << (cy.block(0,0,5,1)).transpose() << "\n";
  std::cout << "err    =" << (cy - ret).norm() << "\n";

  if (A.rows() < 513) {
    std::cout << " *** padm ***\n";
    timer.reset();
    auto E = padm(A);
    auto py = E*v;
    std::cout << "time   =" << timer.elapsed() << "\n";
    MatrixX F;
    timer.reset();
    F = A; F = F.exp();
    std::cout << "time   =" << timer.elapsed() << "\n";
    timer.reset();
    F = A; F = padm(F);
    std::cout << "time   =" << timer.elapsed() << " ";
    //std::cout << "y(1:5) =" << (py.block(0,0,5,1)).transpose() << "\n";
    std::cout << "err    =" << (py - ret).norm() << "\n";
    std::cout << "err2   =" << (E - F).norm() << "\n";
  }
  else
    std::cout << "-- Skipping padm --\n";

  //std::cout << " *** mexpv***\n";
  //ret = mexpv(1,A,v);
}
