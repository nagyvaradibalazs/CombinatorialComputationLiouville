// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "sieve.h"

namespace ComputationalSieves
{
    int64_t SegmentedSieve::multiple_with_lower_limit(int64_t k, int64_t c)
    {
        if (c % k == 0)
        {
            return c;
        }
        else
        {
            return c + k - (c % k);
        }
    }

    std::vector<Prime> SegmentedSieve::list_primes(int64_t n)
    {
        int64_t m = current_max;

        if (m < n)
        {
            std::vector<bool> extended_prime_mark = std::vector<bool>(n - m, true);

            // We mark multiples of primes we already know as composites
            for (Prime p : primes)
            {
                for (int64_t i = multiple_with_lower_limit(p, m + 1); i <= n; i += p)
                {
                    extended_prime_mark[i - (m + 1)] = false;
                }
            }

            // We then apply the method of sieve of Eratosthenes for remaining primes.
            int64_t segment_length = std::sqrt(n);

            for (int64_t i = m + 1; i <= segment_length; i++)
            {
                // If i is prime we mark its multiples
                if (extended_prime_mark[i - (m + 1)])
                {
                    for (int64_t j = i * i; j <= n; j += i)
                    {
                        extended_prime_mark[j - (m + 1)] = false;
                    }
                }
            }

            for (int64_t i = 0; i < extended_prime_mark.size(); i++)
            {
                if (extended_prime_mark[i])
                {
                    primes.push_back(m + 1 + i);
                }
            }

            current_max = n;
        }
        return primes;
    }

    std::vector<Prime> SegmentedSieve::list_primes_segmented(int64_t n, int64_t d)
    {
        int64_t segment_length = static_cast<int64_t>(std::sqrt(n + d));
        list_primes(segment_length);

        std::vector<bool> prime_mark = std::vector<bool>(d, true);
        std::vector<Prime> segmented_primes_result;

        for (Prime p : primes)
        {
            if (p > segment_length)
            {
                break;
            }

            if (p >= n)
            {
                segmented_primes_result.push_back(p);
            }

            for (int64_t i = multiple_with_lower_limit(p, n); i < n + d; i += p)
            {
                prime_mark[i - n] = false;
            }
        }

        for (int64_t i = 0; i < d; i++)
        {
            if (prime_mark[i])
            {
                segmented_primes_result.push_back(n + i);
            }
        }
        return segmented_primes_result;
    }

    std::vector<int> SegmentedSieve::list_mu_segmented(int64_t n, int64_t d)
    {
        int64_t segment_length = static_cast<int64_t>(std::sqrt(n + d));
        list_primes(segment_length);

        std::vector<int64_t> product_of_prime_factors_mark = std::vector<int64_t>(d, 1);
        std::vector<int> mu_result = std::vector<int>(d, 1);

        for (Prime p : primes)
        {
            if (p > segment_length)
            {
                break;
            }

            // Multiples of p^2 are 0 since are not square-free
            int64_t p_square = p * p;
            for (int64_t i = multiple_with_lower_limit(p_square, n); i < n + d; i += p_square)
            {
                mu_result[i - n] = 0;
            }

            for (int64_t i = multiple_with_lower_limit(p, n); i < n + d; i += p)
            {
                mu_result[i - n] *= -1;
                if (mu_result[i - n] != 0)
                {
                    product_of_prime_factors_mark[i - n] *= p;
                }
            }
        }

        // We account for large prime factors that require one more flip
        for (int64_t i = n; i < n + d; i++)
        {
            if (product_of_prime_factors_mark[i - n] != i)
            {
                mu_result[i - n] *= -1;
            }
        }

        return mu_result;
    }

    std::vector<int> SegmentedSieve::list_lambda_segmented(int64_t n, int64_t d)
    {
        int64_t segment_length = static_cast<int64_t>(std::sqrt(n + d));
        list_primes(segment_length);

        std::vector<int64_t> product_of_prime_factors_mark = std::vector<int64_t>(d, 1);
        std::vector<int> lambda_result = std::vector<int>(d, 1);

        for (Prime p : primes)
        {
            if (p > segment_length)
            {
                break;
            }

            int64_t p_power = p;

            while (p_power < n + d)
            {
                for (int64_t i = multiple_with_lower_limit(p_power, n); i < n + d; i += p_power)
                {
                    lambda_result[i - n] *= -1;
                    product_of_prime_factors_mark[i - n] *= p;
                }

                p_power *= p;
            }
        }

        // We account for large prime factors that require one more flip
        for (int64_t i = n; i < n + d; i++)
        {
            if (product_of_prime_factors_mark[i - n] != i)
            {
                lambda_result[i - n] *= -1;
            }
        }

        return lambda_result;
    }

    std::vector<PrimeFactorization> SegmentedSieve::list_prime_factorization_segmented(int64_t n, int64_t d)
    {
        int64_t segment_length = static_cast<int64_t>(std::sqrt(n + d));
        list_primes(segment_length);

        std::vector<int64_t> product_of_prime_factors_mark = std::vector<int64_t>(d, 1);
        std::vector<PrimeFactorization> factorization_result = std::vector<PrimeFactorization>(d, PrimeFactorization());

        for (Prime p : primes)
        {
            if (p > segment_length)
            {
                break;
            }

            int64_t p_power = p;
            Exponent k = 1;

            while (p_power < n + d)
            {
                for (int64_t i = multiple_with_lower_limit(p_power, n); i < n + d; i += p_power)
                {
                    // If the next exponent does not divide anymore, we have p with largest exponent, we can add
                    if (i % (p_power * p) != 0)
                    {
                        factorization_result[i - n].add_prime_factor(p, k);
                        product_of_prime_factors_mark[i - n] *= p_power;
                    }
                }

                p_power *= p;
                k += 1;
            }
        }

        // We account for large prime factors that require one more addition
        for (int64_t i = n; i < n + d; i++)
        {
            if (product_of_prime_factors_mark[i - n] != i)
            {
                int64_t p = i / product_of_prime_factors_mark[i - n];
                factorization_result[i - n].add_prime_factor(p, 1);
            }
        }

        return factorization_result;
    }
}