// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "summatory_liouville.h"
#include "sieve.h"

namespace ComputationLiouville
{
    int64_t CombinatorialLiouville::compute_dependent_variable_sum_s2(int64_t x, int64_t y)
    {
        // Calculating the exact value of y
        double logx = std::log(x);
        double loglogx = std::log(logx);
        double y_exact = std::pow(x, 0.4) * std::pow(loglogx / logx, 0.6);

        DependentVariable dependent_variable;
        ComputationalSieves::SegmentedSieve sieve;

        int64_t result = 0;

        int64_t l_min = std::floor(y_exact) + 1;
        int64_t l_max = std::floor(x / y_exact);

        int64_t x_over_l_min = std::floor(x / l_max);
        int64_t x_over_l_max = std::floor(x / l_min);

        int64_t temp_L = 0;

        int64_t segment_length = std::floor(std::sqrt(x / y_exact));

        for (int64_t n = 1; n <= x_over_l_max; n += segment_length)
        {
            int64_t end = std::min(n + segment_length, x_over_l_max + 1);
            int64_t d = end - n;

            std::vector<int> lambda_segment = sieve.list_lambda_segmented(n, d);

            for (int64_t i = 0; i < d; i++)
            {
                temp_L += lambda_segment[i];

                if (n + i >= x_over_l_min)
                {
                    int64_t temp = std::floor(x / (n + i));
                    int64_t l_to_sum_min = std::floor(x / (n + i + 1));
                    int64_t l_to_sum_max = std::min(l_max, temp);

                    int64_t temp2 = dependent_variable.compute_sum_of_divisor_sums(l_to_sum_min, l_to_sum_max, y, segment_length);
                    result += temp_L * temp2;
                }
            }
        }

        return result;
    }

    int64_t DependentVariable::compute_sum_of_divisor_sums(int64_t n0, int64_t n1, int64_t a, int64_t segment_length)
    {
        ComputationalSieves::SegmentedSieve sieve;

        int64_t result = 0;

        for (int64_t n = n0 + 1; n <= n1; n += segment_length)
        {
            // We have segments [n, n + segment_length)
            int64_t end = std::min(n + segment_length, n1 + 1);
            int64_t d = end - n;

            std::vector<PrimeFactorization> prime_factorization_segment = sieve.list_prime_factorization_segmented(n, d);

            for (int64_t i = 0; i < d; i++)
            {
                // std::cout << "n=" << prime_factorization_segment[i].get_product_of_factors() << ", ";
                result += compute_divisor_sum_with_limit(prime_factorization_segment[i], a);
            }
        }

        return result;
    }

    int64_t DependentVariable::compute_divisor_sum_with_limit(PrimeFactorization &factorization_of_l, int64_t a)
    {
        int64_t l = 1;
        for (PrimeFactor factor : factorization_of_l.get_prime_factors())
        {
            l *= factor.prime;
        }

        return helper_compute_divisor_sum_with_limit(factorization_of_l, factorization_of_l.get_prime_factors().size() - 1, 1, 1, a, l);
    }

    int64_t DependentVariable::helper_compute_divisor_sum_with_limit(PrimeFactorization &factorization_of_l, int64_t prime_index, int64_t m1, int64_t m2, int64_t a, int64_t l)
    {
        if (m1 > a)
        {
            return 0;
        }
        else if (prime_index < 0)
        {
            return 1;
        }
        else if (m2 * a >= l)
        {
            return 0;
        }

        // Prime factors are in increasing order, so to get the maximum we take the last entry
        Prime p = factorization_of_l.get_prime_factors()[prime_index].prime;
        return helper_compute_divisor_sum_with_limit(factorization_of_l, prime_index - 1, m1, p * m2, a, l) - helper_compute_divisor_sum_with_limit(factorization_of_l, prime_index - 1, m1 * p, m2, a, l);
    }
}