#include <Rcpp.h>
using namespace Rcpp;

//' TUVN normal rejection sampling
//'
//' Normal rejection sampling for the truncated standard normal distribution
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncation
//' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double norm_rej(double a, double b)
{
  // first sample from a random normal
  double x = rnorm(1)[0];

  // keep going if p not in (a, b)
  while( (x < a) || (x > b) )
  {
    x = rnorm(1)[0];
  }
  return(x);
}

//' TUVN Half-normal rejection sampling
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncation
//'
//' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double halfnorm_rej(double a, double b)
{
  double x = rnorm(1)[0];
  double abs_x = std::abs(x);

  // keep going if abs_p not in (a, b)
  while( (abs_x < a) || (abs_x > b) )
  {
    x = rnorm(1)[0];
    abs_x = std::abs(x);
  }
  return(abs_x);
}

//' TUVN uniform rejection sampling
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncation
//'
//' @return \item{x}{value of the sample} \item{acc}{number of proposals for acceptance}
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double unif_rej(double a, double b)
{
  while(true)
  {
    double x = runif(1, a, b)[0];
    double u = runif(1)[0];
    double rho;

    // cases for the ratio
    if( (0 >= a) && (0 <= b)) {rho = exp(-1 * (x*x) / 2.0);}
    if (a > 0) {rho = exp( -1 * (x*x - a*a) / 2.0);}
    if (b < 0) {rho = exp(-1 * (x*x - b*b) / 2.0);}

    // accept step
    if(u <= rho) {return(x);}
  }
}

//' TUVN exponential rejection sampling
//'
//' @param a Lower bound of truncation
//' @param b Upper bound of truncation
//' @return value of the sample
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//' @export
// [[Rcpp::export]]
double exp_rej(double a, double b)
{
  // using the optimal lambda define in the paper
  double lambda = (a + sqrt(a*a + 4.0)) / 2.0;

  // loop to generate the sample
  while(true)
  {
    double x = rweibull(1, 1, 1/lambda)[0] + a;
    double u = runif(1)[0];
    double rho = exp(-1 * (x - lambda) * (x - lambda) / 2.0);

    if(u <= rho) {return(x);}
  }
}







