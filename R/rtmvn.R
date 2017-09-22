#' Truncated Multivariate Normal Distribution
#'
#' Random vector generation for the truncated multivariate normal distribution
#'     using a Gibbs sampler.
#'
#' @param n number of samples to be generated
#' @param Mean mean vector
#' @param Sigma covariance matrix
#' @param D matrix of linear constraints
#' @param lower vector of lower bounds
#' @param upper vector of upper bounds
#' @param init vector of initial values for the Gibbs sampler. Must satisfy
#'     the linear constraints.
#'
#' @return a matrix of samples with each column being an idependent sample.
#'
#' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
#'     truncated multivariate normal and student-t distributions subject to
#'     linear inequality constraints. Journal of Statistical Theory and
#'     Practice, 9(4), 712-732.
#'
#' @examples
#' Mean = rep(0,2)
#' rho = 0.5
#' Sigma = matrix(c(10,rho,rho,0.1),2,2)
#' D = matrix(c(1,1,1,-1),2,2)
#' varp = Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2] # var of the sum
#' varm = Sigma[1,1]+Sigma[2,2]-2*Sigma[1,2] # var of the diff
#' sd = c(sqrt(varp),sqrt(varm))
#' lower = -1.5*sd; upper = 1.5*sd; init = rep(0,2)
#' n = 20
#' rtmvn(n,Mean,Sigma,D,lower,upper,init)
#'
#' @export
#'

# Truncated Multivariate Normal Sampler
rtmvn = function(n, Mean, Sigma, D, lower, upper, init)
{
  inits_test = D %*% init
  if((prod(inits_test >= lower & inits_test <= upper)) == 0)
  {
    stop("initial value outside bounds. \n")
  }

  if(is.vector(D) == TRUE)
  {
    Rtilde = t(as.matrix(D))
    lower = as.vector(lower)
    upper = as.vector(upper)
  } else {
    Rtilde = D
  }

  # standardizing the problem
  a = lower - Rtilde %*% Mean
  b = upper - Rtilde %*% Mean
  Sigma_chol = t(chol(Sigma))
  R = Rtilde %*% Sigma_chol

  p = ncol(R)

  z = forwardsolve(l = Sigma_chol, x = init - Mean)

  keep_x = matrix(0, p, n)

  #*****************************************************
  # Gibbs Sampler - seep page 10 of Li & Ghosh (2015)
  #****************************************************
  for (i in 1:n)
  {
    for (j in 1:p)
    {
      rj = as.vector(R[, j])
      Rj = as.matrix(R[, -j])
      zj = as.vector(z[-j])
      a_temp = a - Rj %*% zj
      b_temp = b - Rj %*% zj

      # cases due to equation 6 in paper
      pos = rj > 0
      neg = rj < 0

      if(sum(pos) == 0)
      {
        lower_pos = -Inf
        upper_pos = Inf
      } else {
        lower_pos = max(a_temp[pos] / rj[pos])
        upper_pos = min(b_temp[pos] / rj[pos])
      }

      if(sum(neg) == 0)
      {
        upper_neg = Inf
        lower_neg = -Inf
      } else {
        upper_neg = min(a_temp[neg] / rj[neg])
        lower_neg = max(b_temp[neg]/ rj[neg])
      }

      lower_j = max(lower_pos, lower_neg)
      upper_j = min(upper_pos, upper_neg)

      z[j] <- rtuvn(lower = lower_j, upper = upper_j)
    }

    # for the original results
    x <- Sigma_chol %*% z + Mean
    keep_x[, i] <- x
  }
  return(keep.x)
}
