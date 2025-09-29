#include "core.h"

#include <iostream>
#include <cmath>

#include <eigen3/Eigen/Dense>

long double sl(const long n)
{
    long double S = 0, T = 0;
    for (long i = 0; i < n; ++i)
    {
        S += (i + 1) * log2l(B.coeff(i));
        T += log2l(B.coeff(i));
    }

    return 6.0L * (S - (n + 1) * T * 0.5L) / (n * (n * n - 1));
}
