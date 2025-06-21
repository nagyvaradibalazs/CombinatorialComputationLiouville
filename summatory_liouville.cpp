#include "summatory_liouville.h"

#include <cstdint>

int64_t compute_L(int64_t x)
{
    int64_t result = 0;

    for (int64_t i = 0; i < x; i++)
    {
        result += i;
    }
    return result;
}