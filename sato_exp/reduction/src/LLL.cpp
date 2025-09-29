#include "reduction.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

#include <eigen3/Eigen/Dense>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>

#include "core.h"

void LLL(const double delta, const long end, const long n)
{
    long j, k, q;

    for (k = 1; k < end;)
    {
        for (j = k - 1; j > -1; --j)
        {
            if (fabsl(mu.coeff(k, j)) > 0.5)
            {
                q = round(mu.coeff(k, j));
                basis.row(k) -= q * basis.row(j);
                mu.row(k).head(j + 1) -= static_cast<long double>(q) * mu.row(j).head(j + 1);
            }
        }

        if ((k > 0) && (B.coeff(k) < (delta - mu.coeff(k, k - 1) * mu.coeff(k, k - 1)) * B.coeff(k - 1)))
        {
            basis.row(k - 1).swap(basis.row(k));
            updateSwapGSO(k, n);

            --k;
        }
        else
        {
            ++k;
        }
    }
}
