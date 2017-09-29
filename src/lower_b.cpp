#include <Rcpp.h>
using namespace Rcpp;

//' lower bound for b - normal vs. uniform
//'
//' Calculates the lower bound in Lemma 2.1.1.4 in Li & Ghosh (2015). It is the
//'     bound for b case 2.
//'
//' @param a Lower bound of truncated univariate normal
//'
//' @return value of the bound for b
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double lower_b(double a)
{
  return( sqrt(2.0 * M_PI) + a );
}

//' lower bound for b - half-normal vs uniform
//'
//' Calculates the lower bound in 2.1.1.5 of Li and Ghosh. Called b_1(a) in the
//'     case breakdown
//'
//' @param a Lower bound of truncated univariate normal
//'
//' @return value of the bound for b
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double lower_b1(double a)
{
  return( sqrt(M_PI / 2.0) * exp(a*a / 2.0) + a );
}

//' bound for b - exponential vs uniform sampling
//'
//' Calculates the bound for b in lemma 2.1.1.6 of Li & Ghosh. Called b_2(a)
//'     in the case breakdown.
//'
//' @param a Lower bound of truncated univariate normal
//'
//' @return value of the bound for b
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double lower_b2(double a)
{
  double lambda = a / 2.0 + sqrt(a*a + 4.0) / 2.0;
  return( a + exp(0.5) / lambda * exp((a*a - a * sqrt(a*a + 4.0)) / 4.0) );
}


