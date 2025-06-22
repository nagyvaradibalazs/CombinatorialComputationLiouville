#include "summatory_liouville.h"

#include <chrono>
#include <fstream>
#include <string>
#include <iostream>

int main()
{
    ComputationLiouville::CombinatorialLiouville combinatorial_liouville;

    int64_t x;

    std::cout << "To compute L(x), please specify x: ";

    std::cin >> x;

    std::chrono::steady_clock::time_point start, end;

    start = std::chrono::high_resolution_clock::now();

    int64_t result = combinatorial_liouville.compute_L(x);

    end = std::chrono::high_resolution_clock::now();

    double duration = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::ofstream outfile("output_" + std::to_string(x) + ".txt");

    if (outfile.is_open())
    {
        outfile << "Computing L(" << x << "):" << std::endl
                << "L(" << x << ") = " << result << std::endl
                << "Computation took " << duration << " time.";
        outfile.close();
    }

    return 0;
}