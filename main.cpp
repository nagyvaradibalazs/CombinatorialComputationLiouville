// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "summatory_liouville.h"

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

int main()
{
    ComputationLiouville::CombinatorialLiouville combinatorial_liouville;

    int base, initial_power, max_power;

    std::cout << "Base: ";
    std::cin >> base;

    std::cout << "Initial power: ";
    std::cin >> initial_power;

    std::cout << "Max power: ";
    std::cin >> max_power;

    for (int i = initial_power; i <= max_power; i++)
    {
        int64_t x = std::pow(base, i);

        std::ofstream outfile("result_" + std::to_string(base) + "_" + std::to_string(i) + ".txt");

        if (!outfile.is_open())
        {
            std::cerr << "Failed to open output file.\n";
            return 1;
        }

        std::chrono::steady_clock::time_point start = std::chrono::high_resolution_clock::now();

        int64_t result = combinatorial_liouville.compute_L(x);

        std::chrono::steady_clock::time_point end = std::chrono::high_resolution_clock::now();

        double duration = static_cast<double>(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()) / 1000000;

        outfile << "L(" << x << ") = " << result << " | Time: " << std::fixed << std::setprecision(6) << duration << "s\n";

        std::cout << base << "^" << i << " done" << std::endl;

        outfile.close();
    }

    return 0;
}
