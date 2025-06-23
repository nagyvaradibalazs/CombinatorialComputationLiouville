// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>

#include <boost/multiprecision/cpp_int.hpp>

typedef int Exponent;
typedef int64_t Prime;

typedef boost::multiprecision::int128_t int128_t;

namespace DataTypes
{
    // Structure stores a prime and its exponent
    struct PrimeFactor
    {
        Prime prime;
        Exponent exponent;

        PrimeFactor(Prime p, Exponent exp) : prime(p), exponent(exp)
        {
        }

        // Returns factor as an integer
        operator int64_t()
        {
            return (int64_t)std::pow(prime, exponent);
        }
    };

    struct PrimeFactorization
    {
    private:
        // List of all prime factors
        std::vector<PrimeFactor> prime_factors;

        // Product of these prime factors
        int64_t product_of_factors;

    public:
        // Construct
        PrimeFactorization() : product_of_factors(1)
        {
        }

        // Returns prime factors
        const std::vector<PrimeFactor> get_prime_factors() const
        {
            return prime_factors;
        }

        // Returns product of prime factors
        const int64_t get_product_of_factors() const
        {
            return product_of_factors;
        }

        // Adds prime factor - with update
        void add_prime_factor(Prime p, Exponent k)
        {
            prime_factors.push_back(PrimeFactor(p, k));
            product_of_factors *= (int64_t)std::pow(p, k);
        }
    };

    struct Interval
    {
        int64_t start;
        int64_t end;

        // Creates interval [start, end], by default is empty set
        Interval(int64_t start = 1, int64_t end = 0);

        // Creates the integer interval by roots of ax^2 + bx + c by quadratic formula
        Interval(int64_t a, int64_t b, int64_t c);

        // Shifts interval by a.
        void shift(int64_t a);

        // Intersects with another interval
        void intersect(Interval &I);
    };

    struct Fraction
    {
        int128_t num;
        int128_t denom;

        Fraction(int128_t num, int128_t denom = 1) : num(num), denom(denom)
        {
        }

        double numerical() const;

        int128_t floor() const;

        int128_t round() const;

        Fraction fractional_part() const;

        int sign() const;

        bool is_integral() const;

        bool is_negative() const;

        bool is_zero() const;

        void invert();

        void negate();

        Fraction operator*(Fraction f);
        Fraction operator+(Fraction f);
        Fraction operator-(Fraction f);
        Fraction operator/(Fraction f);
    };
}

#endif // DATA_TYPES_H