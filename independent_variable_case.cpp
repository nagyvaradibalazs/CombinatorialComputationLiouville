// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "summatory_liouville.h"
#include "sieve.h"

namespace ComputationLiouville
{
    int64_t CombinatorialLiouville::compute_independent_variable_sum_s3(int64_t x, int64_t y)
    {
        IndependentVariable independent_variable;

        int64_t result = 0;

        int64_t A1 = y + 1, B1 = y + 1;
        int64_t C = 10, D = 8;

        while (A1 >= 2 * std::pow(6 * C * C * C * x, 0.25) && A1 >= std::sqrt(y) + 1 && A1 >= 2 * D)
        {
            int64_t A = A1 - 2 * (A1 / (2 * D));
            double cube_root_term = std::cbrt((double)A / (6 * x));
            int64_t a = A * cube_root_term;

            while (B1 >= 2 * C * std::cbrt(6 * x / A) && B1 >= std::sqrt(y) + 1 && B1 >= 2 * D)
            {
                int64_t B = B1 - 2 * (B1 / (2 * D));
                int64_t b = B * cube_root_term;

                int64_t D = (1 + y / std::max(2 * a, 2 * b)) * std::max(2 * a, 2 * b);
                DIAG_MODE diag_mode = (A == B && A1 == B1) ? DIAGONAL : NON_DIAG;

                result += independent_variable.dd_sum(A, A1, B, B1, x, D, LINEAR_APPR, a, b, diag_mode);
                B1 = B;
            }

            result += independent_variable.dd_sum(A, A1, 1, B1, x, std::sqrt(y) + 1, BRUTE_FORCE, 0, 0, NON_DIAG);
            A1 = A;
            B1 = A;
        }

        result += independent_variable.dd_sum(1, A1, 1, B1, x, std::sqrt(y) + 1, BRUTE_FORCE, 0, 0, DIAGONAL);

        return result;
    }

    int64_t IndependentVariable::mod(int64_t a, int64_t q)
    {
        if (a >= 0 || a % q == 0)
        {
            return a % q;
        }
        else
            return a % q + q;
    }

    int64_t IndependentVariable::largest_int_congurent(int64_t n, int64_t a, int64_t q)
    {
        if (n == MINFTY || n == INFTY)
        {
            return n;
        }
        return n - mod(n - a, q);
    }

    int64_t IndependentVariable::sum_inter(std::vector<int64_t> G, int64_t r, Interval I, int64_t b, int64_t q)
    {
        int64_t I0 = I.start, I1 = I.end;
        if (I0 <= I1)
        {
            if (I0 != MINFTY && I0 != INFTY)
            {
                I0--;
            }
            int64_t r0 = largest_int_congurent(I0, r, q);
            int64_t r1 = largest_int_congurent(std::min(I1, b - 1), r, q);
            if (r0 > r1 || r1 < -b)
            {
                return 0;
            }
            else if (r0 < -b)
            {
                return G[r1 + b];
            }
            return G[r1 + b] - G[r0 + b];
        }
        return 0;
    }

    std::tuple<std::vector<int64_t>, std::vector<int64_t>, std::vector<int64_t>> IndependentVariable::sum_table(std::function<int64_t(int64_t)> g, int64_t q, int64_t b, int64_t a0)
    {
        std::vector<int64_t> G(2 * b), rho(q), sigma(q + 1, 0);
        for (int64_t m = -b; m < -b + q; m++)
        {
            G[m + b] = g(m);
        }
        for (int64_t m = -b + q; m < b; m++)
        {
            G[m + b] = G[m + b - q] + g(m);
        }
        int64_t r = mod(a0 * (b - q), q);
        int64_t a = mod(a0, q);

        for (int64_t m = b - q; m < b; m++)
        {
            rho[r] = G[m + b];
            r += a;
            if (r >= q)
            {
                r -= q;
            }
        }

        for (int64_t r = 1; r < q; r++)
        {
            sigma[r + 1] = sigma[r] + rho[q - r];
        }
        return std::make_tuple(G, rho, sigma);
    }

    int64_t IndependentVariable::ray_sum(std::function<int64_t(int64_t)> g, int64_t q, int64_t b, int s)
    {
        int64_t S = 0;
        if (s < 0)
        {
            for (int64_t m = q; m < b; m += q)
            {
                S += g(m);
            }
        }
        else if (s > 0)
        {
            for (int64_t m = q; m <= b; m += q)
            {
                S += g(-m);
            }
        }
        return S;
    }

    int64_t IndependentVariable::linear_sum(std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t a, int64_t b, int64_t N, int64_t d0, int64_t m0)
    {
        int64_t S_f0 = 0, S_f1 = 0;
        int64_t S_g0 = 0, S_g1 = 0;
        for (int64_t d = -a; d < a; d++)
        {
            S_f0 += f(d);
            S_f1 += f(d) * (int64_t)(Fraction((int128_t)N * (d0 - d), (int128_t)d0 * d0 * m0).floor());
        }
        for (int64_t m = -b; m < b; m++)
        {
            S_g0 += g(m);
            S_g1 += g(m) * (int64_t)(Fraction((int128_t)N * (-m), (int128_t)d0 * m0 * m0).floor());
        }
        return S_f0 * S_g1 + S_f1 * S_g0;
    }

    int64_t IndependentVariable::special_1(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b)
    {
        int64_t r = a_inv * (-1 - r0);

        int64_t gamma1 = d * (-R0 * q - (r0 + 1) + a * m0);
        Interval J(-a * d, gamma1, N * q);
        J.shift(-m0);

        return sum_inter(G, r, Interval(MINFTY, INFTY), b, q) - sum_inter(G, r, J, b, q);
    }

    int64_t IndependentVariable::special_0b(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_beta, int s_delta)
    {
        Interval I;
        if (s_delta > 0)
        {
            I = Interval(MINFTY, std::ceil(-Q) - 1);
        }
        else if (s_delta < 0)
        {
            I = Interval(std::floor(-Q) + 1, INFTY);
        }
        else if (s_beta < 0)
        {
            I = Interval(MINFTY, INFTY);
        }

        int64_t gamma1 = d * (-R0 * q - r0 + a * m0);
        Interval J(-a * d, gamma1, N * q);
        J.shift(-m0);

        J.intersect(I);
        return sum_inter(G, -r0 * a_inv, I, b, q) - sum_inter(G, -r0 * a_inv, J, b, q);
    }

    int64_t IndependentVariable::special_00(std::vector<int64_t> G, int64_t N, int64_t q, int64_t a, int64_t a_inv, int64_t R0, int64_t r0, int64_t m0, int64_t d, int64_t b, double Q, int s_delta)
    {
        Interval I, I_c(MINFTY, INFTY);
        if (s_delta > 0)
        {
            I = Interval(MINFTY, std::ceil(-Q) - 1);
            I_c = Interval(std::ceil(-Q), INFTY);
        }
        else if (s_delta < 0)
        {
            I = Interval(std::floor(-Q) + 1, INFTY);
            I_c = Interval(MINFTY, std::floor(-Q));
        }

        std::vector<Interval> J(2);
        for (int j = 0; j <= 1; j++)
        {
            if (a != 0)
            {
                int64_t gamma1 = d * (-R0 - (r0 + j) + a * m0);
                J[j] = Interval(-a * d, gamma1, N);
                J[j].shift(-m0);
            }
            else
            {
                J[j] = Interval(N / (d * (R0 + r0 + j)) - m0 + 1, INFTY);
            }
        }

        J[0].intersect(I);
        J[1].intersect(I_c);
        return sum_inter(G, 0, Interval(MINFTY, INFTY), b, q) - sum_inter(G, 0, J[0], b, q) - sum_inter(G, 0, J[1], b, q);
    }

    int64_t IndependentVariable::special_0a(std::vector<int64_t> G, int64_t q, int64_t a, int64_t a_inv, int64_t r0, int64_t b, double Q, int s_beta, int s_delta)
    {
        int r = -r0 * a_inv;

        Interval I;
        if (0 < r0 && r0 < q)
        {
            if (s_delta > 0)
            {
                I = Interval(std::ceil(-Q), INFTY);
            }
            else if (s_delta < 0)
            {
                I = Interval(MINFTY, std::floor(-Q));
            }
            else if (s_beta >= 0)
            {
                I = Interval(MINFTY, INFTY);
            }
        }
        else
        {
            if (s_beta < 0)
            {
                if (s_delta < 0)
                {
                    return sum_inter(G, r, Interval(MINFTY, std::floor(-Q)), b, q) + sum_inter(G, r, Interval(1, INFTY), b, q);
                }
                else if (s_delta > 0)
                {
                    return sum_inter(G, r, Interval(MINFTY, -1), b, q) + sum_inter(G, r, Interval(std::ceil(-Q), INFTY), b, q);
                }
            }
            else if (s_beta > 0)
            {
                if (s_delta > 0)
                {
                    I = Interval(std::ceil(-Q), -1);
                }
                else if (s_delta < 0)
                {
                    I = Interval(1, std::floor(-Q));
                }
            }
        }
        return sum_inter(G, r, I, b, q);
    }

    int64_t IndependentVariable::sum_by_linear_approximation(std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t N, int64_t d0, int64_t m0, int64_t a, int64_t b)
    {
        Fraction alpha = Fraction(-N, (int128_t)d0 * m0 * m0);

        int64_t S = linear_sum(f, g, a, b, N, d0, m0);

        std::tuple<int128_t, int128_t, int128_t, int> tuple = approximation_by_reduces_fractions(alpha, 2 * b);
        int64_t a0 = (int64_t)std::get<0>(tuple);
        int64_t a0_inv = (int64_t)std::get<1>(tuple);
        int64_t q = (int64_t)std::get<2>(tuple);
        int sgn_delta = std::get<3>(tuple);
        Fraction delta = alpha - Fraction(a0, q);

        int64_t raySum = ray_sum(g, q, b, sgn_delta);

        std::tuple<std::vector<int64_t>, std::vector<int64_t>, std::vector<int64_t>>
            tuple2 = sum_table(g, q, b, a0);
        std::vector<int64_t> G = std::get<0>(tuple2);
        std::vector<int64_t> rho = std::get<1>(tuple2);
        std::vector<int64_t> sigma = std::get<2>(tuple2);

        for (int64_t d = -a; d < a; d++)
        {
            if (f(d) != 0)
            {
                int64_t d_ = d0 + d;

                Fraction R0 = Fraction(N * (d0 - d), d0 * d0 * m0);
                Fraction R0_frac = R0.fractional_part();
                int64_t r0 = (int64_t)((R0_frac * Fraction(q)).round());
                Fraction beta = R0_frac - Fraction(r0, q);
                int sgn_beta = beta.sign();

                double Q(1);
                if (sgn_delta != 0)
                {
                    Q = (double)((beta / delta).numerical());
                }

                int64_t T1 = 0;
                if (q > 1)
                {
                    T1 += special_1(G, N, q, a0, a0_inv, (int64_t)(R0.floor()), r0, m0, d_, b);

                    T1 += special_0b(G, N, q, a0, a0_inv, (int64_t)(R0.floor()), r0, m0, d_, b, Q, sgn_beta, sgn_delta);
                }
                else
                {
                    T1 += special_00(G, N, q, a0, a0_inv, (int64_t)(R0.floor()), r0, m0, d_, b, Q, sgn_delta);
                }

                int64_t T2 = 0;
                T2 += sigma[r0];
                if (0 < r0 && r0 < q)
                {
                    T2 += raySum;
                }

                T2 += special_0a(G, q, a0, a0_inv, r0, b, Q, sgn_beta, sgn_delta);

                S += (T1 + T2) * f(d);
            }
        }
        return S;
    }

    int64_t IndependentVariable::brute_force_double_sum(int64_t d0, int64_t d1, int64_t m0, int64_t m1, std::vector<int> f, std::vector<int> g, int64_t N)
    {
        int64_t S = 0;
        for (int64_t d = d0; d < d1; d++)
        {
            int64_t T = 0;
            for (int64_t m = m0; m < m1; m++)
            {
                int64_t k = N / (d * m);
                T += k * g[m - m0];
            }
            S += T * f[d - d0];
        }
        return S;
    }

    int64_t IndependentVariable::double_sum(int64_t d0, int64_t d1, int64_t m0, int64_t m1, int64_t a, int64_t b, std::function<int64_t(int64_t)> f, std::function<int64_t(int64_t)> g, int64_t N)
    {
        int64_t S = 0;
        int64_t mm = m0, dm = d0;
        for (int64_t dm = d0; dm < d1; dm += 2 * a)
        {
            int64_t dp = std::min(dm + 2 * a, d1);
            int64_t d_mid = (dm + dp) / 2, d_len = (dp - dm) / 2;
            std::function<int64_t(int64_t)> f_ = [f, d0, dm, d_len](int64_t d)
            {
                return f(d + dm - d0 + d_len);
            };
            for (int64_t mm = m0; mm < m1; mm += 2 * b)
            {
                int64_t mp = std::min(mm + 2 * b, m1);
                int64_t m_mid = (mm + mp) / 2, m_len = (mp - mm) / 2;
                std::function<int64_t(int64_t)> g_ = [g, m0, mm, m_len](int64_t m)
                {
                    return g(m + mm - m0 + m_len);
                };
                S += sum_by_linear_approximation(f_, g_, N, d_mid, m_mid, d_len, m_len);
            }
        }
        return S;
    }

    int64_t IndependentVariable::dd_sum(int64_t A0, int64_t A1, int64_t B0, int64_t B1, int64_t N, int64_t D, COMP_MODE comp_mode, int64_t a, int64_t b, DIAG_MODE diag_mode)
    {
        ComputationalSieves::SegmentedSieve sieve;

        int64_t S = 0;

        std::vector<int> mu_d, mu_m;
        std::function<int64_t(int64_t)> m_d, m_m;
        std::vector<int> lambda_d, lambda_m;
        std::function<int64_t(int64_t)> l_d, l_m;
        for (int64_t d0 = A0; d0 < A1; d0 += D)
        {
            int64_t d1 = std::min(d0 + D, A1);
            lambda_d = sieve.list_lambda_segmented(d0, d1 - d0);
            l_d = [lambda_d](int64_t d)
            {
                return lambda_d[d];
            };
            if (diag_mode == NON_DIAG)
            {
                mu_d = sieve.list_mu_segmented(d0, d1 - d0);
                m_d = [mu_d](int64_t d)
                {
                    return mu_d[d];
                };
            }

            for (int64_t m0 = B0; m0 < B1; m0 += D)
            {
                int64_t m1 = std::min(m0 + D, B1);
                mu_m = sieve.list_mu_segmented(m0, m1 - m0);
                m_m = [mu_m](int64_t m)
                {
                    return mu_m[m];
                };
                if (diag_mode == NON_DIAG)
                {
                    lambda_m = sieve.list_lambda_segmented(m0, m1 - m0);
                    l_m = [lambda_m](int64_t m)
                    {
                        return lambda_m[m];
                    };
                }

                switch (comp_mode)
                {
                case LINEAR_APPR:
                    S += double_sum(d0, d1, m0, m1, a, b, l_d, m_m, N);
                    if (diag_mode == NON_DIAG)
                    {
                        S += double_sum(d0, d1, m0, m1, a, b, m_d, l_m, N);
                    }
                    break;
                case BRUTE_FORCE:

                    S += brute_force_double_sum(d0, d1, m0, m1, lambda_d, mu_m, N);
                    if (diag_mode == NON_DIAG)
                    {
                        S += brute_force_double_sum(m0, m1, d0, d1, lambda_m, mu_d, N);
                    }
                    break;
                }
            }
        }
        return S;
    }

    std::tuple<int128_t, int128_t, int128_t, int> IndependentVariable::approximation_by_reduces_fractions(Fraction alpha, int64_t Q)
    {
        int sgn_alpha = 1;
        if (alpha.is_negative())
        {
            sgn_alpha = -1;
            alpha.negate();
        }

        int128_t p[2]{0, 1};
        int128_t q[2]{1, 0};
        int s = 1;

        while (q[1] <= Q)
        {
            int128_t a = alpha.floor();
            int128_t new_p = a * p[1] + p[0];
            int128_t new_q = a * q[1] + q[0];
            p[0] = p[1];
            p[1] = new_p;
            q[0] = q[1];
            q[1] = new_q;

            if (alpha.is_integral() && q[1] <= Q)
            {
                return std::make_tuple(sgn_alpha * p[1], (-sgn_alpha * s * q[0]) % q[1], q[1], 0);
            }

            alpha = alpha.fractional_part();
            alpha.invert();
            s *= -1;
        }
        return std::make_tuple(sgn_alpha * p[0], (-sgn_alpha * s * q[1]) % q[0], q[0], sgn_alpha * s);
    }
}