// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "summatory_liouville.h"
#include "sieve.h"

namespace ComputationLiouville
{
    ComputationalSieves::SegmentedSieve sieve;

    int64_t CombinatorialLiouville::compute_trivial_sum_s1(int64_t x, int64_t y)
    {
        int64_t result = 0;

        int64_t segment_length = std::sqrt(y);

        for (int64_t n = 1; n <= y; n += segment_length)
        {
            int64_t end = std::min(n + segment_length, y + 1);
            int64_t d = end - n;

            std::vector<int> mu_segment = sieve.list_mu_segmented(n, d);

            for (int64_t i = 0; i < d; i++)
            {
                int64_t current_number = n + i;
                int mu = mu_segment[i];
                if (mu == 0)
                {
                    continue;
                }

                int64_t root = std::floor(std::sqrt(x / current_number));
                result += mu * root;
            }
        }
        return result;
    }

    int64_t CombinatorialLiouville::compute_trivial_sum_s4(int64_t x, int64_t y)
    {
        int64_t result = 0;

        int64_t segment_length = std::sqrt(y);

        // Calculating the exact value of y
        double logx = std::log(x);
        double loglogx = std::log(logx);
        double y_exact = std::pow(x, 0.4) * std::pow(loglogx / logx, 0.6);

        // Calculating inner sum: sum_{m <= y} mu(m) * floor(x / ym)
        int64_t inner_sum = 0;

        for (int64_t n = 1; n <= y; n += segment_length)
        {
            int64_t end = std::min(n + segment_length, y + 1);
            int64_t d = end - n;

            std::vector<int> mu_segment = sieve.list_mu_segmented(n, d);

            for (int64_t i = 0; i < d; i++)
            {
                int64_t current_number = n + i;
                int mu = mu_segment[i];
                if (mu == 0)
                {
                    continue;
                }

                int64_t floor = std::floor(x / (y_exact * current_number));
                inner_sum += mu * floor;
            }
        }

        for (int64_t n = 1; n <= y; n += segment_length)
        {
            int64_t end = std::min(n + segment_length, y + 1);
            int64_t d = end - n;

            std::vector<int> lambda_segment = sieve.list_lambda_segmented(n, d);

            for (int64_t i = 0; i < d; i++)
            {
                result += lambda_segment[i] * inner_sum;
            }
        }

        return result;
    }
}
