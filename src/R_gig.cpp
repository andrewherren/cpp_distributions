#include <cpp11.hpp>
#include <gig.h>

[[cpp11::register]]
cpp11::writable::doubles sample_gig_cpp(int n, double p, double a, double b) {
    cpp11::writable::doubles output(n);
    for (int i = 0; i < n; i++) {
        output[i] = CustomDists::SampleGIG(p, a, b);
    }
    return output;
}
