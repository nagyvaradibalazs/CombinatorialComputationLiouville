#include "summatory_liouville.h"

#include <cmath>

namespace ComputationLiouville
{
    int64_t CombinatorialLiouville::compute_L(int64_t x)
    {
        int64_t s1, s2, s3, s4;

        double logx = std::log(x);
        double loglogx = std::log(logx);
        double value = std::pow(x, 0.4) * std::pow(loglogx / logx, 0.6);

        int64_t y = static_cast<int64_t>(std::floor(value));

        s1 = compute_trivial_sum_s1(x, y);
        s2 = compute_dependent_variable_sum_s2(x, y);
        s3 = compute_independent_variable_sum_s3(x, y);
        s4 = compute_trivial_sum_s4(x, y);

        return s1 - s2 - s3 + s4;
    }
}