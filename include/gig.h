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
 * \return double 
 */
double SampleGIG(double p, double a, double b);

double SampleGIGHormanLeydold(double p, double beta);

} // namespace CustomDists

#endif // CUSTOMDISTS_GIG_H_
