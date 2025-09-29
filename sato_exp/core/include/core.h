#ifndef CORE_H
#define CORE_H

#include <iostream>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;        // long-type matrix
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;                                      // long-type vector
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld; // long double-type matrix
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;                               // long double-type vector

extern VectorXld B;
extern VectorXld s;
extern MatrixXli basis;
extern MatrixXld mu;
extern MatrixXld R;
extern NTL::ZZ vol;
extern MatrixXld hmu;
extern VectorXld hB;

extern NTL::vec_ZZ B_num, B_den;
extern NTL::mat_ZZ mu_num, mu_den;

/**
 * @brief computes GSA-slope
 *
 * @param n rank of lattice
 * @return long double GSA-slope
 */
long double sl(const long n);

/**
 * @brief computes root of Hermite-factor
 *
 * @param n rank of lattice
 * @return long double root of Hermite-factor
 */
long double rhf(const long n);

/**
 * @brief Compute GSO-informations of lattice basis
 *
 */
void computeGSO();

#endif // !LAT_PY_H
