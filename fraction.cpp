// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "data_types.h"

namespace DataTypes
{
    double Fraction::numerical() const
    {
        if (denom != 0)
        {
            return double(num) / double(denom);
        }
        return 1;
    }

    int128_t Fraction::floor() const
    {
        if (num < 0 && num % denom != 0)
        {
            return (num / denom) - 1;
        }
        return num / denom;
    }

    int128_t Fraction::round() const
    {
        Fraction q(2 * num + denom, 2 * denom);
        return q.floor();
    }

    Fraction Fraction::fractional_part() const
    {
        int128_t temp = num % denom;
        if (temp < 0)
        {
            temp += denom;
        }
        return Fraction(temp, denom);
    }

    bool Fraction::is_negative() const
    {
        if (denom > 0)
        {
            return num < 0;
        }
        else
            return num > 0;
    }

    bool Fraction::is_zero() const
    {
        return num == 0;
    }

    int Fraction::sign() const
    {
        if (is_negative())
        {
            return -1;
        }
        else if (is_zero())
        {
            return 0;
        }
        return 1;
    }

    bool Fraction::is_integral() const
    {
        return num % denom == 0;
    }

    void Fraction::invert()
    {
        int128_t num_temp = num;
        num = denom;
        denom = num_temp;
    };

    void Fraction::negate()
    {
        num *= -1;
    }

    Fraction Fraction::operator*(Fraction f)
    {
        return Fraction(num * f.num, denom * f.denom);
    }

    Fraction Fraction::operator+(Fraction f)
    {
        return Fraction(num * f.denom + f.num * denom, denom * f.denom);
    }

    Fraction Fraction::operator-(Fraction f)
    {
        return Fraction(num * f.denom - f.num * denom, denom * f.denom);
    }

    Fraction Fraction::operator/(Fraction f)
    {
        return Fraction(num * f.denom, denom * f.num);
    }
}