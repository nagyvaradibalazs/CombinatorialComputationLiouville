#ifndef SUMMATORY_LIOUVILLE_H
#define SUMMATORY_LIOUVILLE_H

#include "data_types.h"

#include <cstdint>
#include <algorithm>

using namespace DataTypes;

namespace ComputationLiouville
{
    class CombinatorialLiouville
    {
    public:
        // Computes L(x).
        int64_t compute_L(int64_t x);

    private:
        // Computes the trival case S_1 = sum_{m <= y} mu(m) * floor(sqrt(x / m))
        static int64_t compute_trivial_sum_s1(int64_t x, int64_t y);

        // Computes the trivial case S_4 = sum_{k <= y} lambda(k) * sum_{m <= y} mu(m) * (floor(x / ym) - Delta(x / ym))
        static int64_t compute_trivial_sum_s4(int64_t x, int64_t y);

        // Computes S_2 = sum_{l=mn; y < l < x / y} L(x / l) sum_{m | l; m <= y} mu(m)
        static int64_t compute_dependent_variable_sum_s2(int64_t x, int64_t y);

        // Computes S_3 = sum_{m <= y} mu(m) sum_{k <= y} lambda(k) * floor(x / km)
        static int64_t compute_independent_variable_sum_s3(int64_t x, int64_t y);
    };

    class DependentVariable
    {
    public:
        // Computes sum_{n0 < n <= n1} sum_{m|l, m <= a} mu(m) in segments of choice
        int64_t compute_sum_of_divisor_sums(int64_t n0, int64_t n1, int64_t a, int64_t segment_length);

    private:
        // Computes sum_{m|l, m <= a} mu(m).
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t compute_divisor_sum_with_limit(PrimeFactorization &factorization_of_l, int64_t a);

        // Helper function
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t helper_compute_divisor_sum_with_limit(PrimeFactorization &factorization_of_l, int64_t prim_index, int64_t m1, int64_t m2, int64_t a, int64_t l);
    };

    class IndependentVariable
    {
    };
}

#endif // SUMMATORY_LIOUVILLE_H