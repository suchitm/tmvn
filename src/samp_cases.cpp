#include <Rcpp.h>
using namespace Rcpp;

#include "include/sampling_methods.hpp"
#include "include/lower_b.hpp"

//##############################################################################
// Rejection sampling by case - list of cases given on page 4 of
// Li & Ghosh (2015)
//##############################################################################

//' Sample Case 1
//'
//' generate a sample from case 1 where the bounds are \eqn{[a, \infty)}.
//'    The cases are listed on page 4 of Li & Ghosh (2015).
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return the value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_case1(double a, double b)
{
  double samp;

  if(a <= 0)
    samp = norm_rej(a, b);
  else if (a < 0.25696)
    samp = halfnorm_rej(a, b);
  else
    samp = exp_rej(a, b);

  return(samp);
}

//' Sample Case 2
//'
//' generate a sample from case 2 where the bounds are \eqn{0 \in [a, b]}.
//'    The cases are listed on page 4 of Li & Ghosh (2015).
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_case2(double a, double b)
{
  double samp;
  double this_lower_b = lower_b(a);

  if(b > this_lower_b)
    samp = norm_rej(a, b);
  else
    samp = unif_rej(a, b);

  return(samp);
}

//' Sample Case 3
//'
//' generate a sample from case 3 where the bounds are \eqn{[a, b], a > 0}.
//'    The cases are listed on page 4 of Li & Ghosh (2015).
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_case3(double a, double b)
{
  double samp;

  if(a < 0.25696)
  {
    double blower1 = lower_b1(a);
    if(b <= blower1)
      samp = unif_rej(a, b);
    else
      samp = halfnorm_rej(a, b);
  }
  else
  {
    double blower2 = lower_b2(a);
    if(b <= blower2)
      samp = unif_rej(a, b);
    else
      samp = exp_rej(a, b);
  }
  return(samp);
}

//' Sample Case 4
//'
//' generate a sample from case 4 where the bounds are \eqn{(-\infty, b]}.
//'    This case is symmetric to case 1 in the function sample_case1.
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return value of the samples
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_case4(double a, double b)
{
  double temp = sample_case1(-b, -a);
  return(-temp);
}

//' Sample Case 5
//'
//' generate a sample from case 5 where the bounds are \eqn{[a, b], b < 0}.
//'    This case is symmetric to case 3 in the function sample_case3.
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_case5(double a, double b)
{
  double temp = sample_case3(-b, -a);
  return(-temp);
}

//' Full Rejection sampling steps
//'
//' Generate a sample from a truncated univariate standard normal distribution,
//'     TN(0, 1; a, b).
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncattion
//'
//' @return value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double sample_tuvsn(double a, double b)
{
  double samp;
  if( (a == R_NegInf) || (b == R_PosInf) )
  {
    if(b == R_PosInf)
      samp = sample_case1(a, b);
    else
      samp = sample_case4(a, b);
  }
  else
  {
    if(a >= 0)
      samp = sample_case3(a, b);
    else if (b <= 0)
      samp = sample_case5(a, b);
    else
      samp = sample_case2(a, b);
  }
  return(samp);
}




