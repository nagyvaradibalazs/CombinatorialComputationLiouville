// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "data_types.h"

namespace DataTypes
{
    Interval::Interval(int64_t start, int64_t end) : start(start), end(end) {};

    Interval::Interval(int64_t a, int64_t b, int64_t c)
    {
        int64_t discriminant = b * b - 4 * a * c;
        if (discriminant >= 0)
        {
            int64_t sqrt_term = std::sqrt(discriminant);
            if (a < 0)
            {
                start = std::ceil((double)(-b + sqrt_term) / (2 * a));
                end = std::floor((double)(-b - sqrt_term) / (2 * a));
            }
            else if (a > 0)
            {
                // Check if discriminant is a square or not
                if (sqrt_term * sqrt_term != discriminant)
                {
                    start = std::ceil((double)(-b - sqrt_term) / (2 * a));
                    end = std::floor((double)(-b + sqrt_term) / (2 * a));
                }
                else
                {
                    start = std::floor((double)(-b - sqrt_term) / (2 * a)) + 1;
                    end = std::ceil((double)(-b + sqrt_term) / (2 * a)) - 1;
                }
            }
        }
    };

    void Interval::shift(int64_t a)
    {
        if (start != LLONG_MIN && start != LLONG_MAX)
        {
            start += a;
        }
        if (end != LLONG_MIN && end != LLONG_MAX)
        {
            end += a;
        }
    };

    void Interval::intersect(Interval &I)
    {
        start = std::max(start, I.start);
        end = std::min(end, I.end);
    };
}