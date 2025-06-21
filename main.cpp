#include "summatory_liouville.h"

#include <chrono>
#include <fstream>
#include <string>

int main()
{
    int64_t x = 1000000000;

    std::chrono::steady_clock::time_point start, end;

    start = std::chrono::high_resolution_clock::now();

    int64_t result = compute_L(x);

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