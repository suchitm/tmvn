#' Truncated Multivariate Student t Distribution
#'
#' This function generates random vectors form the truncates multivariate
#'     Student t distribution. It uses random sampling from rtmvn since the
#'     Student t distribution can be represented as a scale-mixture of normals.
#'
#' @param n number of samples to be generated
#' @param Mean mean vector
#' @param Sigma covariance matrix
#' @param nu degress of freedom for the t-distribution
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
#' lower = -1.5*sd; upper = 1.5*sd; int = rep(0,2)
#' nu = 5
#' n = 20
#' rtmvt(n,Mean,Sigma,nu,D,lower,upper,int)

#'
#' @export
#'

# Truncated Multivariate Student-T Sampler
rtmvt = function(n, Mean, Sigma, nu, D, lower, upper, init)
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

  a = lower - Rtilde %*% Mean
  b = upper - Rtilde %*% Mean
  Sigma_chol = t(chol(Sigma))
  R = Rtilde %*% Sigma_chol
  x = forwardsolve(l = Sigma_chol, x = init - Mean)

  p = ncol(R)
  keep_t = matrix(0, p, n)

  for (i in 1:n)
  {
    u = rchisq(1, df = nu)
    denom = sqrt(u / nu)
    lw = a * denom
    up = b * denom
    z0 = x * denom
    z = rtmvn(n = 1, Mean = rep(0, p), Sigma = diag(1, p), D = R, lower = lw,
               upper = up, init = z0)

    x = z / denom

    w = Sigma_chol %*% x + Mean
    keep_t[, i] = w
  }
  return(keep_t)
}



