#include <gig.h>
#include <cassert>

namespace CustomDists {

double SampleGIG(double p, double a, double b) {
  assert((a > 0) || (b > 0));
  assert(!((p == 0) && ((a == 0) || (b == 0))));
  if (a == 0) {
    // TODO: IG(p, b/2) sampler
  } else if (b == 0) {
    // TODO: Gamma(p, a/2) sampler (where a/2 is a rate parameter)
  } else {
    // Transform parameters to the "alternative" parameterization 
    // presented on wikipedia and in Hörmann and Leydold (2014)
    double alpha = std::sqrt(a/b);
    double beta = std::sqrt(a*b);
    // Here GIG(p, a, b) is now GIG(p, alpha, beta)
    // Now, since alpha is a scaling parameter of the distribution, 
    // we sample Z ~ GIG(p, 1, beta) and then scale X = alpha*Z
    double Z = SampleGIGHormanLeydold(p, beta);
    double X = alpha*Z;
    return X;
  }
}

double SampleGIGHormanLeydold(double p, double beta) {
  // TODO: placeholder to run the "sampler" end-to-end in R
  // This function will eventually be a 
  return 1.;
}

} // namespace CustomDists