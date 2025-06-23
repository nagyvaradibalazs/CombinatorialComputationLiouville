// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef SIEVE_H
#define SIEVE_H

#include "data_types.h"

#include <vector>

using namespace DataTypes;

namespace ComputationalSieves
{
    class SegmentedSieve
    {
    private:
        // Current upper bound
        int64_t current_max;

        // List of all primes up to this bound
        std::vector<Prime> primes;

        // Computes first instance of a multiple of K such that it is at least C
        static int64_t multiple_with_lower_limit(int64_t k, int64_t c);

    public:
        // Constructor
        SegmentedSieve() : current_max(1) {}

        // Returns list of primes up to N
        std::vector<Prime> list_primes(int64_t n);

        // Returns list of primes in [N, N + D)
        std::vector<Prime> list_primes_segmented(int64_t n, int64_t d);

        // Returns list of values of mu(n) for n in [N, N + D)
        std::vector<int> list_mu_segmented(int64_t n, int64_t d);

        // Returns list of values of lambda(n) for n in [N, N + D)
        std::vector<int> list_lambda_segmented(int64_t n, int64_t d);

        // Returns prime factorization of integers in [N, N + D)
        std::vector<PrimeFactorization> list_prime_factorization_segmented(int64_t n, int64_t d);
    };
}

#endif // SIEVE_H