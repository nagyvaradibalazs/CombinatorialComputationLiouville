#include <chrono>
#include <fstream>
#include <string>

int main()
{
    int64_t n = 1000;

    std::chrono::steady_clock::time_point start, end;

    start = std::chrono::high_resolution_clock::now();

    std::string s = std::to_string(n);

    end = std::chrono::high_resolution_clock::now();

    double duration = (double)std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

    std::ofstream outfile("output_" + s + ".txt");

    if (outfile.is_open())
    {
        outfile << s << " takes " << duration << " time " << std::endl;
        outfile.close();
    }

    return 0;
}