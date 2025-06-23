// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef SUMMATORY_LIOUVILLE_H
#define SUMMATORY_LIOUVILLE_H

#include "data_types.h"

#include <cstdint>
#include <algorithm>
#include <functional>

#define INFTY LLONG_MAX
#define MINFTY LLONG_MIN

using namespace DataTypes;

namespace ComputationLiouville
{
    enum COMP_MODE
    {
        LINEAR_APPR,
        BRUTE_FORCE
    };

    enum DIAG_MODE
    {
        DIAGONAL,
        NON_DIAG
    };

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
    public:
        int64_t dd_sum(int64_t A0, int64_t A1, int64_t B0, int64_t B1, int64_t N, int64_t D, COMP_MODE comp_mode, int64_t a, int64_t b, DIAG_MODE diag_mode);

    private:
        // Returns congurent t a mod q representation
        static int64_t mod(int64_t a, int64_t q);

        // Returns largest integer <= n congurent to a mod q
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t largest_int_congurent(int64_t n, int64_t a, int64_t q);

        // Returns sum of g(n) over n in I with n = r mod q
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t sum_inter(std::vector<int64_t> G, int64_t r, Interval I, int64_t b, int64_t q);

        // Returns tree tables of vals of g: G(n), rho(r), sigma(r)
        // Algorithm by Helfgott & Thompson, 2023.
        static std::tuple<std::vector<int64_t>, std::vector<int64_t>, std::vector<int64_t>> sum_table(std::function<int64_t(int64_t)> g, int64_t q, int64_t b, int64_t a0);

        // Returns the sum of values of g(n) for n = b + kq for integers k with sk < 0
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t ray_sum(std::function<int64_t(int64_t)> g, int64_t q, int64_t b, int s);

        // Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m) (floor(alpha0 + alpha1 d) + floor(alpha2 m)) by separation of variables
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t linear_sum(std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t a, int64_t b, int64_t N, int64_t d0, int64_t m0);

        // Computation of L - L1 for a0(m-m0) + r0 = -1 mod q, with q > 1.
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t special_1(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b);

        // Computing L - L1 for a0(m-m0) + r0 = 0 mod q, with q > 1
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t special_0b(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_beta, int s_delta);

        // Computation of L - L1 for q = 1
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t special_00(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_delta);

        // Computation of L1 - L2 for where a0(m-m0) + r0 = 0 mod q.
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t special_0a(std::vector<int64_t> G, int64_t q, int64_t a, int64_t a_inv, int64_t r0, int64_t b, double Q, int s_beta, int s_delta);

        // Returns sum_{(d,m) in [d0 - a, d0 + a) x [m0 - b, m0 + b)} f(d)g(m)floor(N/dm)
        // Algorithm by Helfgott & Thompson, 2023.
        static int64_t sum_by_linear_approximation(std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t N, int64_t d0, int64_t m0, int64_t a, int64_t b);

        // Computes the sum f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by brute force
        static int64_t brute_force_double_sum(int64_t d0, int64_t d1, int64_t m0, int64_t m1, std::vector<int> f, std::vector<int> g, int64_t N);

        // Computes the sum of f(d)g(m)floor(N/dm) over the rectangle [d0, d1) x [m0, m1) by linear approximation
        static int64_t double_sum(int64_t d0, int64_t d1, int64_t m0, int64_t m1, int64_t a, int64_t b, std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t N);

        // Returns a tuple (p, p', q, s) of integers so that abs(alpha - p/q) <= 1/(qQ) with gcd(p, q) = 1, q <= Q, and pp' = 1 mod q, while s = sgn(alpha - p/q)
        // Algorithm by Helfgott & Thompson, 2023.
        static std::tuple<int128_t, int128_t, int128_t, int> approximation_by_reduces_fractions(Fraction alpha, int64_t Q);
    };
}

#endif // SUMMATORY_LIOUVILLE_H