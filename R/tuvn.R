#' Univariate Truncated Normal Distribution
#'
#' Random number generator for the truncated normal distribution. To
#'     calculate the density please use dtuvn.
#'
#' @param n number of samples
#' @param mean mean
#' @param sd standard deviation
#' @param lower lower bound
#' @param upper upper bound
#' @return a vector of generated samples
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#'
#' @examples
#' # sample from truncated normal with mean 10, sd 20, and lower = 10, upper = 20
#' rtuvn(n = 1, mean = 10, sd = 20, lower = 10, upper = 20)
#'
#' @export

rtuvn = function(n = 1, mean = 0, sd = 1, lower, upper)
{
  # transform the boundaries
  a <- (lower - mean)/sd
  b <- (upper - mean)/sd

  # generate n samples from TN(0,1;a,b)
  Z <- rep(0,n)
  for (i in 1:n)
  {
    temp <- sample_tuvsn(a,b)
    Z[i] <- temp$x
  }

  # transform the data back
  samp <- sd*Z + mean

  return(samp)
}

#' Univariate Truncated Normal Distribution
#'
#' Calculates density for the truncated normal distribution.
#'
#' @param x value of sample
#' @param mean mean
#' @param sd standard deviation
#' @param lower lower bound
#' @param upper upper bound
#' @return a vector of generated samples
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#'
#' @export

dtuvn = function(x, mean, sd, lower, upper)
{
  value = dnorm(x, mean, sd) /
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  return(ifelse(x < lower | x > upper, 0, value))
}
