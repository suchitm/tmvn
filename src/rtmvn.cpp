# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "include/rtuvn.hpp"

//' Truncated Multivariate Normal Distribution
//'
//' Random vector generation for the truncated multivariate normal distribution
//'     using a Gibbs sampler.
//'
//' @param n number of samples to be generated
//' @param p dimension of the distribution to sample
//' @param R  matrix of linear constraints
//' @param a vector of lower bounds
//' @param b vector of upper bounds
//' @param z vector of initial values for the Gibbs sampler.
//'
//' @return a matrix of samples with each column being an idependent sample.
//'
//' @references Li, Y., & Ghosh, S. K. (2015). Efficient sampling methods for
//'     truncated multivariate normal and student-t distributions subject to
//'     linear inequality constraints. Journal of Statistical Theory and
//'     Practice, 9(4), 712-732.
//'
//' @examples
//' Mean = rep(0,2)
//' rho = 0.5
//' Sigma = matrix(c(10,rho,rho,0.1),2,2)
//' D = matrix(c(1,1,1,-1),2,2)
//' varp = Sigma[1,1]+Sigma[2,2]+2*Sigma[1,2] # var of the sum
//' varm = Sigma[1,1]+Sigma[2,2]-2*Sigma[1,2] # var of the diff
//' sd = c(sqrt(varp),sqrt(varm))
//' lower = -1.5*sd; upper = 1.5*sd; init = rep(0,2)
//' n = 20
//' rtmvn(n,Mean,Sigma,D,lower,upper,init)
//'
//' @export
// [[Rcpp::export]]
arma::mat rtmvn_gibbs(int n, int p, arma::vec Mean, arma::mat Sigma_chol,
                   arma::mat R, arma::vec a, arma::vec b, arma::vec z)
{
  arma::mat keep_x(p, n);
  int nrow_R = R.n_rows;
  arma::vec temp;

  double lower_pos, upper_pos, lower_neg, upper_neg;
  int sum_pos, sum_neg;

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      // loop to create indicies for subsetting
      arma::uvec idx(p - 1);
      for(int k = 0; k < p; k++)
      {
        if(k == j)
          continue;
        else if (k > j)
          idx(k - 1) = k;
        else
          idx(k) = k;
      }

      arma::vec rj = R.col(j);

      arma::mat Rj = R.cols(idx);
      arma::vec zj = z.rows(idx);
      arma::vec Rj_zj = Rj * zj;

      arma::vec a_temp = a - Rj_zj;
      arma::vec b_temp = b - Rj_zj;

      arma::vec pos(nrow_R), neg(nrow_R);
      for(int k = 0; k < nrow_R; k++)
      {
        if (rj(k) < 0)
        {
          pos(k) = 0;
          neg(k) = 1;
        }
        else
        {
          pos(k) = 1;
          neg(k) = 0;
        }
      }

      sum_pos = sum(pos);
      sum_neg = sum(neg);
      arma::uvec pos_idx(sum_pos), neg_idx(sum_neg);

      int lp = 0; int ln = 0;
      for(int k = 0; k < nrow_R; k++)
      {
        if(pos(k) == 1)
        {
          pos_idx(lp) = k;
          lp += 1;
        }
        if(neg(k) == 1)
        {
          neg_idx(ln) = k;
          ln += 1;
        }
      }

      if(sum(pos) == 0 )
      {
        lower_pos = R_NegInf;
        upper_pos = R_PosInf;
      }
      else
      {
        temp = a_temp.rows(pos_idx) / rj.rows(pos_idx);
        lower_pos = max(temp);

        temp = b_temp.rows(pos_idx) / rj.rows(pos_idx);
        upper_pos = min(temp);
      }

      if(sum(neg) == 0)
      {
        upper_neg = R_PosInf;
        lower_neg = R_NegInf;
      }
      else
      {
        temp = a_temp.rows(neg_idx) / rj.rows(neg_idx);
        upper_neg = min(temp);

        temp = b_temp.rows(neg_idx) / rj.rows(neg_idx);
        lower_neg = min(temp);
      }

      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);

      z(j) = rtuvn(1, 0, 1, lower_j, upper_j)[0];
    }
    keep_x.col(i) = Sigma_chol * z + Mean;
  }
  return(keep_x);
}
