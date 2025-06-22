#ifndef DATA_TYPES_H
#define DATA_TYPES_H

#include <vector>
#include <cmath>
#include <cstdint>

typedef int Exponent;
typedef int64_t Prime;

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
}

#endif // DATA_TYPES_H