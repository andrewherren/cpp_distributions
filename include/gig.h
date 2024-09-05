/*!
 * Copyright (c) 2024 customdists authors. All rights reserved.
 * Licensed under the MIT License. See LICENSE file in the project root for license information.
 */
#ifndef CUSTOMDISTS_GIG_H_
#define CUSTOMDISTS_GIG_H_
#include <random>
#include <boost/math/special_functions/bessel.hpp>

namespace CustomDists {

/*!
 * \brief Sample from a Generalized Inverse Gaussian distribution. 
 * Uses the default parameterization described in 
 * https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution
 * 
 * \param p First parameter of the distribution, which reduces to a shape parameter 
 * in the degenerate case where either of `a` or `b` is 0 (in which case this 
 * distribution collapses to inverse gamma and gamma, respectively). Can be any real 
 * number but if `p = 0` then both `a > 0` and `b > 0` is required.
 * \param a Second parameter of the distribution. Must be nonnegative and cannot be zero 
 * if either of `p` or `b` is zero.
 * \param b Third parameter of the distribution. Must be nonnegative and cannot be zero 
 * if either of `p` or `a` is zero.
 * \param rng Reference to a C++ random number generator initialized in outer program scope.
 * \return double 
 */
double SampleGIG(double p, double a, double b, std::mt19937& rng);

/*!
 * \brief Sample from a Generalized Inverse Gaussian distribution in the "hard" case
 * in which both `a` and `b` are nonzero, so the distribution does not simplify to a 
 * Gamma or Inverse Gamma. This relies on the sampler proposed in Hörmann and Leydold
 * (2014) [1] which uses an "alternative parameterization" of the GIG distribution. 
 * Rather than GIG(p, a, b), we define alpha = sqrt(a/b) and beta = sqrt(ab) and 
 * parameterize the distribution as GIG(p, alpha, beta). Since alpha reduces to a scale 
 * parameter under this parameterization, we can sample from GIG(p, 1, beta) and then multiply 
 * by alpha.
 * 
 * [1] Wolfgang Hörmann and Josef Leydold. 2014. Generating generalized inverse Gaussian random variates. 
 * Statistics and Computing 24, 4 (July 2014), 547–557. https://doi.org/10.1007/s11222-013-9387-3
 * 
 * \param p First parameter of the alternative parameterization of GIG used in Hörmann and Leydold (2014). 
 * Can be any real number.
 * \param beta Second parameter of the alternative parameterization of GIG used in Hörmann and Leydold (2014). 
 * \param rng Reference to a C++ random number generator initialized in outer program scope.
 * \return double 
 */
double SampleGIGHormanLeydold(double p, double beta, std::mt19937& rng);

/*!
 * \brief Sample from a Gamma distribution with rate parameterization. 
 * This is a simplification of the GIG when b = 0.
 * 
 * \param shape Shape parameter for the Gamma distribution
 * \param rate Rate parameter for the Gamma distribution
 * \param rng Reference to a C++ random number generator initialized in outer program scope.
 * \return double 
 */
double SampleGamma(double shape, double rate, std::mt19937& rng);

/*!
 * \brief Sample from an Inverse Gamma distribution. 
 * This is a simplification of the GIG when a = 0.
 * 
 * \param shape Shape parameter for the Inverse Gamma distribution
 * \param scale Scale parameter for the Inverse Gamma distribution
 * \param rng Reference to a C++ random number generator initialized in outer program scope.
 * \return double 
 */
double SampleInverseGamma(double shape, double scale, std::mt19937& rng);

} // namespace CustomDists

#endif // CUSTOMDISTS_GIG_H_
