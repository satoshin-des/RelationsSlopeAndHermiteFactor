#include <iostream>
#include <vector>

#include <eigen3/Eigen/Dense>

typedef Eigen::Matrix<long, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXli;        // long-type matrix
typedef Eigen::Matrix<long, 1, Eigen::Dynamic> VectorXli;                                      // long-type vector
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXld; // long double-type matrix
typedef Eigen::Matrix<long double, 1, Eigen::Dynamic> VectorXld;                               // long double-type vector

/**
 * @brief Updates GSO-information with swapping of the lattice basis vectors \bm{b}_{k-1} and \bm{b}_{k}
 *
 * @param k index
 * @param n rank of lattice
 */
void updateSwapGSO(const long k, const long n);

/**
 * @brief Applies LLL-reduction
 *
 * @param delta reduction parameter
 * @param end end index
 */
void LLL(const double delta, const long end, const long n);

extern "C"
{
    /**
     * @brief Applies BKZ-reduction to lattice basis
     *
     * @param basis_ptr lattice basis matrix
     * @param delta reduction parameter to Lovasz condition
     * @param beta block size
     * @param pruning if make use of pruning or not
     * @param output_sl output GSA-slope or not
     * @param output_rhf output root of Hermite-factor or not
     * @param n rank of lattice
     * @param m null of lattice
     */
    void BKZ(
        long **basis_ptr,
        const double delta,
        const long beta,
        const long max_loops,
        const bool pruning,
        const long n,
        const long m);
}
