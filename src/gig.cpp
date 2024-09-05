#include <gig.h>
#include <meta.h>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace CustomDists {

double SampleGIG(double p, double a, double b, std::mt19937& rng) {
  assert((a > 0) || (b > 0));
  assert(!((p == 0) && ((a == 0) || (b == 0))));
  if (a == 0) {
    return SampleInverseGamma(p, b/2., rng);
  } else if (b == 0) {
    return SampleGamma(p, a/2, rng);
  } else {
    // Transform parameters to the "alternative" parameterization 
    // presented on wikipedia and in Hörmann and Leydold (2014)
    double alpha = std::sqrt(b/a);
    double beta = std::sqrt(a*b);
    // Here GIG(p, a, b) is now GIG(p, alpha, beta)
    // Now, since alpha is a scaling parameter of the distribution, 
    // we sample Z ~ GIG(p, 1, beta) and then rescale X = alpha*Z
    double Z = SampleGIGHormanLeydold(p, beta, rng);
    return alpha*Z;
  }
}

double SampleGamma(double shape, double rate, std::mt19937& rng) {
  double gamma_scale = 1./rate;
  std::gamma_distribution<double> gamma_dist = std::gamma_distribution<double>(shape, gamma_scale);
  return gamma_dist(rng);
}

double SampleInverseGamma(double shape, double scale, std::mt19937& rng) {
  double gamma_scale = 1./scale;
  std::gamma_distribution<double> gamma_dist = std::gamma_distribution<double>(shape, gamma_scale);
  return (1/gamma_dist(rng));
}

double GenInvGaussLogQuasiPDF(double x, double p, double beta) {
  return (p - 1)*std::log(x) - beta*(x + 1./x)/2.;
}

double GenInvGaussMode(double p, double beta) {
  if (p < 1) {
    return beta / (std::sqrt((p - 1.)*(p - 1.) + beta*beta) + (1 - p));
  } else {
    return (std::sqrt((1. - p)*(1. - p) + beta*beta) - (1 - p)) / beta;
  }
}

double SampleGIGHormanLeydoldRatioUnifModeShift(double p, double beta, double m, std::mt19937& rng) {
  // Initial computations
  double a2 = -2. * (p + 1) / beta - m;
  double a1 = 2. * m * (p - 1) / beta - 1;
  double p1 = a1 - (a2 * a2) / 3;
  double q1 = 2 * (a2 * a2 * a2) / 27. - a2 * a1 / 3 + m;
  double phi = std::acos(-q1 * std::sqrt(-27. / (p1 * p1 * p1)) / 2);
  double s1 = -std::sqrt(-4. * p1 / 3.);
  double root1 = s1 * std::cos(phi / 3. + pi_constant / 3.) - a2 / 3.;
  double root2 = -s1 * std::cos(phi / 3.) - a2 / 3.;
  double lm = GenInvGaussLogQuasiPDF(m, p, beta);
  double d1 = GenInvGaussLogQuasiPDF(root1, p, beta) - lm;
  double d2 = GenInvGaussLogQuasiPDF(root2, p, beta) - lm;
  // Bounding rectange for the rejection sampler
  double vmin = (root1 - m) * std::exp(0.5 * d1);
  double vmax = (root2 - m) * std::exp(0.5 * d2);
  double umax = 1.; // sqrt(h(m)) = sqrt(g(m) / g(m)) = 1
  // Rejection sampler
  int iter = 0;
  int accept = false;
  double u, v, result;
  std::uniform_real_distribution<double> unif_dist(0., 1.);
  while (!accept) {
    u = unif_dist(rng) * umax;
    v = unif_dist(rng) * (vmax - vmin) + vmin;
    result = v / u + m;
    accept = 2 * std::log(u) <= (GenInvGaussLogQuasiPDF(result, p, beta) - lm);
    iter++;
  }
  return result;
}

double SampleGIGHormanLeydoldRatioUnifNoModeShift(double p, double beta, double m, std::mt19937& rng) {
  // Initial computations
  double xplus = ((1 + p) + std::sqrt((1 + p) * (1 + p) + (beta * beta))) / beta;
  // Bounding rectange for the rejection sampler
  double umax = std::exp(0.5 * GenInvGaussLogQuasiPDF(m, p, beta));
  double vmin = 0.;
  double vmax = xplus * std::exp(0.5 * GenInvGaussLogQuasiPDF(xplus, p, beta));
  // Rejection sampler
  int iter = 0;
  int accept = false;
  double u, v, result;
  std::uniform_real_distribution<double> unif_dist(0., 1.);
  while (!accept) {
    u = unif_dist(rng) * umax;
    v = unif_dist(rng) * (vmax - vmin) + vmin;
    result = v / u;
    accept = 2 * std::log(u) <= (GenInvGaussLogQuasiPDF(result, p, beta));
    iter++;
  }
  return result;
}

double SampleGIGHormanLeydoldNoRatioUnif(double p, double beta, double m, std::mt19937& rng) {
  // Initial setup
  double x0 = beta / (1 - p);
  double xs = std::max(x0, 2. / beta);
  double k1 = std::exp(GenInvGaussLogQuasiPDF(m, p, beta));
  double A1 = k1 * x0;
  double k2, A2;
  if (x0 < 2 / beta) {
    k2 = std::exp(-beta);
    if (p > 0) {
      A2 = k2 * (std::pow(2. / beta, p) - std::pow(x0, p)) / p;
    } else {
      A2 = k2 * std::log(2. / (beta * beta));
    }
  } else {
    k2 = 0.;
    A2 = 0.;
  }
  double k3 = std::pow(xs, (p - 1));
  double A3 = 2 * k3 * std::exp(-xs * beta / 2) / beta;
  double A = A1 + A2 + A3;
  // Rejection sampler
  int iter = 0;
  int accept = false;
  double u, v, result, h;
  std::uniform_real_distribution<double> unif_dist(0., 1.);
  while (!accept) {
    u = unif_dist(rng);
    v = unif_dist(rng) * A;
    if (v <= A1) {
      result = k2 * v / A1;
      h = k1;
    } else if ((v > A1) && (v <= (A1 + A2))) {
      if (p > 0) {
        result = std::pow(std::pow(x0, p) + (v - A1) * p / k2, (1. / p));
      } else {
        result = beta * std::exp((v - A1) * std::exp(beta));
      }
      h = k2 * std::pow(result, p - 1);
    } else {
      double z = std::exp(-xs * beta / 2) - beta * (v - A1 - A2) / (2 * k3);
      result = (-2. / beta) * std::log(z);
      h = k3 * std::exp(-result * beta / 2);
    }
    accept = std::log(u * h) <= (GenInvGaussLogQuasiPDF(result, p, beta));
    iter++;
  }
  return result;
}

double SampleGIGHormanLeydold(double p, double beta, std::mt19937& rng) {
  // TODO: handle the case where p < 0 using the fact that X ~ GIG(p,b) 
  // implies that 1/X ~ GIG(-p,b). For now, we will just require that p >= 0.
  assert(p >= 0);
  // Compute mode
  double m = GenInvGaussMode(p, beta);
  // Figure out which of the three samplers presented in Hörmann and Leydold (2014) to use
  bool mode_shift = false;
  bool ratio_unif = true;
  if ((p >= 1) || (beta > 1)) {
    mode_shift = true;
  } else if (beta >= (std::min(0.5, 2 * std::sqrt(1.-p) / 3))) {
    mode_shift = false;
  } else {
    ratio_unif = false;
  }
  // Run the algorithms
  if (ratio_unif) {
    if (mode_shift) {
      return SampleGIGHormanLeydoldRatioUnifModeShift(p, beta, m, rng);
    } else {
      return SampleGIGHormanLeydoldRatioUnifNoModeShift(p, beta, m, rng);
    }
  } else {
    return SampleGIGHormanLeydoldNoRatioUnif(p, beta, m, rng);
  }
}

} // namespace CustomDists