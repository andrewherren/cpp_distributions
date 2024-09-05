#include <cpp11.hpp>
#include <gig.h>
#include <random>

[[cpp11::register]]
cpp11::writable::doubles sample_gig_cpp(int n, double p, double a, double b, int random_seed = -1) {
    std::mt19937 rng;
    if (random_seed == -1) {
        std::random_device rd;
        rng = std::mt19937(rd());
    } else {
        rng = std::mt19937(random_seed);
    }

    cpp11::writable::doubles output(n);
    for (int i = 0; i < n; i++) {
        output[i] = CustomDists::SampleGIG(p, a, b, rng);
    }
    return output;
}
