// eigen-expokit: read a matrix in Compressed Row Storage format
// https://gitlab.com/tesch1/eigen-expokit
//
// SPDX-License-Identifier: MPL-2.0
// Copyright (C) 2019 Michael Tesch <tesch1@gmail.com>
//
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <stdio.h>

namespace eexpokit {

/** \brief Read a matrix in Compressed Row Storage format.  Returns an empty
  * matrix if \a filename cannot be opened.
  */
template <class Scalar>
Eigen::SparseMatrix<Scalar>
loadcrs(const char * filename)
{
  auto fid = fopen(filename, "r");
  if (!fid) {
    perror("could not open file\n");
    Eigen::SparseMatrix<Scalar> A;
    return A;
  }
  int n;  fscanf(fid,"%u", &n);
  int nz; fscanf(fid,"%u", &nz);

  eigen_assert(nz >= n);

  Eigen::VectorXi ia(nz,1), ja(nz,1);
  Eigen::VectorXd ra(nz,1);
  ia.setZero();
  ja.setZero();
  ra.setZero();

  for (int i = 0; i <= n; i++) {
    fscanf(fid,"%u",&ia(i));
    ia(i)--;
  }

  for (int i = 0; i < nz; i++) {
    int x;
    double d;
    fscanf(fid,"%u",&x);
    ja(i) = x - 1;
    fscanf(fid,"%lf",&d);
    ra(i) = d;
  }
  fclose(fid);

  // expand row indices ...
  for (int i = n-1; i >= 0; i--) {
    for (int k = 0; k < ia(i+1)-ia(i); k++) {
      eigen_assert(nz > 0);
      ia(nz-1) = i;
      nz = nz-1;
    }
  }

  eigen_assert(nz == 1);

  Eigen::SparseMatrix<Scalar> A(n,n);
  for (int i = 0; i < ia.rows(); i++) {
    eigen_assert(ia(i) < n);
    eigen_assert(ja(i) < n);
    eigen_assert(ia(i) >= 0);
    eigen_assert(ja(i) >= 0);
    //if (i < 10 || i > ia.rows()-10)
    //  std::cout << "adding (" << ia(i) << "," << ja(i) << ") = " << ra(i) << "\n";
    A.insert(ia(i), ja(i)) = ra(i);
  }
  return A;
}

}
