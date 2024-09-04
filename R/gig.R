#' Take n draws from a Generalized Inverse Gaussian (GIG) 
#' distribution, using the default parameterization in
#' https://en.wikipedia.org/wiki/Generalized_inverse_Gaussian_distribution 
#'
#' @param n Number of samples to draw
#' @param p First parameter of the distribution, which reduces to a shape parameter 
#' in the degenerate case where either of `a` or `b` is 0 (in which case this 
#' distribution collapses to inverse gamma and gamma, respectively). Can be any real 
#' number but if `p = 0` then both `a > 0` and `b > 0` is required.
#' @param a Second parameter of the distribution. Must be nonnegative and cannot be zero 
#' if either of `p` or `b` is zero.
#' @param b Third parameter of the distribution. Must be nonnegative and cannot be zero 
#' if either of `p` or `a` is zero.
#'
#' @return Vector of length `n`
#' @export
#'
#' @examples
#' n <- 100
#' p <- 1.5
#' a <- 3.4
#' b <- 0.6
#' gig_draws <- sample_gig(n,p,a,b)
#' hist(gig_draws)
sample_gig <- function(n, p, a, b) {
    warning("sample_gig has not yet been implemented. The output is constant. Do not rely on this function if this warning is displayed at runtime.")
    sample_gig_cpp(n, p, a, b)
}