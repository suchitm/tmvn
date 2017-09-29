###############################################################################
# Rejection sampling by case - list of cases given on page 4 of
# Li & Ghosh (2015)
###############################################################################

#' Sample Case 1
#'
#' generate a sample from case 1 where the bounds are \eqn{[a, \infty)}.
#'    The cases are listed on page 4 of Li & Ghosh (2015).
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# Case 1: [a,infty)
sample_case1_r = function(a, b = Inf)
{
  if (a <= 0)
  {
    samp = norm_rej_r(a = a, b = b)
  } else if (a < 0.25696) {
      samp = halfnorm_rej_r(a = a, b = b)
  } else {
      samp = exp_rej_r(a = a, b = b, lambda = NULL)
  }
  return(samp)
}

#' Sample Case 2
#'
#' generate a sample from case 2 where the bounds are \eqn{0 \in [a, b]}.
#'    The cases are listed on page 4 of Li & Ghosh (2015).
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# Case 2: 0 in [a,b], a<0<b
sample_case2_r = function(a, b)
{
  if (b > lower_b_r(a))
  {
    samp = norm_rej_r(a = a, b = b)
  } else {
    samp = unif_rej_r(a = a, b = b)
  }
  return(samp)
}

#' Sample Case 3
#'
#' generate a sample from case 3 where the bounds are \eqn{[a, b], a > 0}.
#'    The cases are listed on page 4 of Li & Ghosh (2015).
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# Case 3: [a,b], a>0
sample_case3_r = function(a, b)
{
  if (a < 0.25696)
  {
    blower1 = lower_b1_r(a)
    if (b <= blower1)
    {
      samp = unif_rej_r(a = a, b = b)
    } else {
      samp = halfnorm_rej_r(a = a, b = b)
    }
  } else {
    blower2 = lower_b2_r(a)
    if (b <= blower2)
    {
      samp = unif_rej_r(a = a, b = b)
    } else {
      samp = exp_rej_r(a = a, b = b)
    }
  }
  return(samp)
}

#' Sample Case 4
#'
#' generate a sample from case 4 where the bounds are \eqn{(-\infty, b]}.
#'    This case is symmetric to case 1 in the function sample_case1.
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# Case 4: (-infty,b] (symmetric to Case 1)
sample_case4_r = function(a, b)
{
  temp = sample_case1_r(a= -b, b = -a)
  samp = list(x = -temp$x, acc = temp$acc)
  return(samp)
}

#' Sample Case 5
#'
#' generate a sample from case 5 where the bounds are \eqn{[a, b], b < 0}.
#'    This case is symmetric to case 3 in the function sample_case3.
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# Case 5: [a,b], b<=0 (symmetric to Case 3)
sample_case5_r <- function(a, b)
{
  temp = sample_case3_r(a = -b, b = -a)
  samp = list(x = -temp$x, acc = temp$acc)
  return(samp)
}

#' Full Rejection sampling steps
#'
#' Generate a sample from a truncated univariate standard normal distribution,
#'     TN(0, 1; a, b).
#'
#' @param a Lower bound of truncation
#' @param b Upper bound of truncattion
#'
#' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

## Final rejection sampling
sample_tuvsn_r <- function(a, b)
{
  if (a == -Inf | b == Inf)
  {
    if (b == Inf)
    {
      samp <- sample_case1_r(a=a,b=b)
    } else {
        samp <- sample_case4_r(a=a,b=b)
    }
  } else {
    if(a >= 0) {
      samp = sample_case3_r(a = a, b = b)
    } else if (b <= 0) {
      samp = sample_case5_r(a = a, b = b)
    } else {
      samp = sample_case2_r(a = a, b = b)
    }
  }
  return(samp)
}

