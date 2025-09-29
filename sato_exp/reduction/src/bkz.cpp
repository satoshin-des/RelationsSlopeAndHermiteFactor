#include "reduction.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "core.h"
#include "svp.h"

extern "C" void BKZ(
    long **basis_ptr,
    const double delta,
    const long beta,
    const long max_loops,
    const bool pruning,
    const long n,
    const long m)
{
    long z, i, j, num_tour = 0, k = 0, h, d, l, p;
    long double radius;
    long double log_vol = 0.0L;
    FILE *log_sl, *log_hf;
    VectorXli t, coeff_vec, v;
    MatrixXld err_mat = MatrixXld::Zero(n, n);
    NTL::mat_ZZ basis_ntl;
    NTL::vec_RR B_ntl;
    NTL::mat_RR mu_ntl;

    basis_ntl.SetDims(n, m);
    basis = MatrixXli::Zero(n, m);

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis.coeffRef(i, j) = basis_ptr[i][j];
            basis_ntl[i][j] = NTL::to_ZZ(basis_ptr[i][j]);
        }
    }

    log_sl = fopen("sl_log.csv", "w");
    log_hf = fopen("hf_log.csv", "w");
    fprintf(log_sl, "val\n");
    fprintf(log_hf, "val\n");

    computeGSO();

    for (i = 0; i < n; ++i)
    {
        log_vol += log2l(B.coeff(i));
    }
    log_vol *= 0.5L;

    fprintf(log_sl, "%.15Lf\n", sl(n));
    fprintf(log_hf, "%.15Lf\n", -2.0L * (log2l(B.coeff(0)) * 0.5L - log_vol / n) / (n - 1.0L));

    for (z = k = 0; z < n - 2;)
    {
        if (k == n - 1)
        {
            fprintf(log_sl, "%.15Lf\n", sl(n));
            fprintf(log_hf, "%.15Lf\n", -2.0L * (log2l(B.coeff(0)) * 0.5L - log_vol / n) / (n - 1.0L));

            k = 0;
            ++num_tour;
            if (max_loops > -1)
            {
                if (num_tour > max_loops)
                {
                    break;
                }
            }
        }
        ++k;

        if (n < k + beta)
        {
            l = n;
        }
        else
        {
            l = k + beta - 1;
        }
        if (l < n)
        {
            h = l + 1;
        }
        else
        {
            h = n;
        }
        d = l - k + 1;

        radius = delta * B.coeff(k - 1);
        if (enumSV(coeff_vec, radius, mu, B, pruning, k - 1, l))
        {
            v = coeff_vec * basis.block(k - 1, 0, d, m);

            z = 0;

            basis_ntl.SetDims(n + 1, m);
            for (j = 0; j < n + 1; ++j)
            {
                if (j < k - 1)
                {
                    for (p = 0; p < m; ++p)
                    {
                        basis_ntl[j][p] = basis.coeff(j, p);
                    }
                }
                else if (j == k - 1)
                {
                    for (p = 0; p < m; ++p)
                    {
                        basis_ntl[j][p] = v.coeff(p);
                    }
                }
                else
                {
                    for (p = 0; p < m; ++p)
                    {
                        basis_ntl[j][p] = basis.coeff(j - 1, p);
                    }
                }
            }
            NTL::LLL_FP(basis_ntl, 0.99);
            for (j = 0; j < n; ++j)
            {
                for (p = 0; p < m; ++p)
                {
                    basis.coeffRef(j, p) = NTL::to_long(basis_ntl[j + 1][p]);
                }
            }

            computeGSO();
        }
        else
        {
            LLL(delta, h, n);

            ++z;
        }
    }

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < m; ++j)
        {
            basis_ptr[i][j] = basis.coeff(i, j);
        }
    }

    fclose(log_sl);
    fclose(log_hf);
}
