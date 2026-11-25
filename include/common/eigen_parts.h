#ifndef EIGEN_PARTS_
#define EIGEN_PARTS_

#include "common.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/SparseLU/SparseLU.h"

typedef Eigen::SparseMatrix<double> spmat;
typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vect;
typedef Eigen::Triplet<double> triplet;
typedef Eigen::SparseLU<spmat> SparseLU;

inline void concate_vector(const vector<vect> & vv, vect & v)
{
    vector<int> csum = {0};
    for(vect el : vv)
    csum.push_back(csum.back()+el.rows());

    v = vect(csum.back());

    for(int i=0;i<vv.size();i++)
    v.block(csum[i], 0, csum[i+1]-csum[i], 1) = vv[i];
}

inline void split_vector(vector<vect> & vv, const vect & v)
{
    vector<int> csum = {0};
    for(vect el : vv)
    csum.push_back(csum.back()+el.rows());

    for(int i=0;i<vv.size();i++)
    vv[i] = v.block(csum[i], 0, csum[i+1]-csum[i], 1);
}

#endif