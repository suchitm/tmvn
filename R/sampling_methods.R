#' TUVN normal rejection sampling
#'
#' Normal rejection sampling for the truncated standard normal distribution
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncation; default is infinity
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# normal rejection
norm_rej_r = function(a, b = Inf)
{
  acc = 0
  repeat
  {
    x = rnorm(1)
    acc = acc + 1
    if (x >= a & x <= b)
      return(list(x = x, acc = acc))
  }
}

#' TUVN Half-normal rejection sampling
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncation; default is infinity
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# half-normal rejection
halfnorm_rej_r = function(a, b)
{
  acc = 0
  repeat
  {
    x = rnorm(1)
    acc = acc + 1
    if (abs(x) >= a & abs(x) <= b)
      return(list(x = abs(x), acc = acc))
  }
}

#' TUVN uniform rejection sampling
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncation; default is infinity
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# uniform rejection
unif_rej_r = function(a, b)
{
  acc = 0
  repeat
  {
    x = runif(1, a, b)
    u = runif(1)

    # different cases for the ratio
    if (0 >= a & 0 <= b) {rho = exp(-x^2/2)}
    if (a > 0){rho = exp(-(x^2 - a^2) / 2)}
    if (b < 0){rho = exp(-(x^2 - b^2) / 2)}

    acc = acc + 1
    if (u <= rho)
      return(list(x = x, acc = acc))
  }
}

#' TUVN exponential rejection sampling
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncation; default is infinity
#' @param lambda rate parameter for the exponential distribution. Optimal value is
#'     calculated by default. Refer to Lemma 2.1.1.2 in Li and Ghosh (2015).
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export
# exponential rejection
exp_rej_r <- function(a, b = Inf, lambda = NULL)
{
  if (is.null(lambda))
  {
    lambda = (a + sqrt(a^2 + 4))/2
  } else {
    lambda = lambda
  }

  acc = 0
  repeat
  {
    x = rweibull(1, shape=1, scale = 1/lambda) + a
    # x = rexp(1, rate = lambda) + a
    u = runif(1)
    rho = exp(-(x - lambda)^2/2)
    acc = acc + 1

    if (u <= rho & x < b)
      return(list(x = x, acc = acc))
  }
}
