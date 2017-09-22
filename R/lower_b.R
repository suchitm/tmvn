################################################################################
# functions for lower bounds on b that are listed on page 4 and 5 of
# Li and Ghosh
################################################################################

#' lower bound for b - normal vs. uniform
#'
#' Calculates the lower bound in Lemma 2.1.1.4 in Li & Ghosh (2015). It is the
#'     bound for b case 2.
#'
#' @param a Lower bound of truncated univariate normal
#'
#' @return value of the bound for b
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# lower bound of b for normal vs unif in 0 in [a,b]
lower_b <- function(a)
{
  sqrt(2 * pi) + a
}

#' lower bound for b - half-normal vs uniform
#'
#' Calculates the lower bound in 2.1.1.5 of Li and Ghosh. Called b_1(a) in the
#'     case breakdown
#'
#' @param a Lower bound of truncated univariate normal
#'
#' @return value of the bound for b
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# lower bound of b for half vs unif in [a,b]>=0
lower_b1 = function(a)
{
  sqrt(pi / 2) * exp(a^2 / 2) + a
}

#' bound for b - exponential vs uniform sampling
#'
#' Calculates the bound for b in lemma 2.1.1.6 of Li & Ghosh. Called b_2(a)
#'     in the case breakdown.
#'
#' @param a Lower bound of truncated univariate normal
#'
#' @return value of the bound for b
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#' @export

# lower bound of b for exp vs unif in [a,b]>=0
lower_b2 = function(a)
{
  lambda = (a + sqrt(a^2 + 4)) / 2
  return(a + exp(1/2) / lambda * exp((a^2 - a * sqrt(a^2 + 4)) / 4))
}
