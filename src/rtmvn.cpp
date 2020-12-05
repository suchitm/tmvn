# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "include/rtuvn.hpp"

//' Gibbs sampler for the Truncated Multivariate Normal Distribution
//'
//' Random vector generation for the truncated multivariate normal distribution
//'     using a Gibbs sampler.
//'
//' @param n number of samples to be generated
//' @param p dimension of the distribution to sample
//' @param R matrix of linear constraints
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
//' @export
// [[Rcpp::export]]
arma::mat rtmvn_gibbs(int n, int p, arma::vec Mean, arma::mat Sigma_chol,
                      arma::mat R, arma::vec a, arma::vec b, arma::vec z)
{
  arma::mat keep_x(p, n);
  arma::vec temp;

  // std::cout << "R: \n" << R << endl;

  double lower_pos, upper_pos, lower_neg, upper_neg;
  int sum_pos, sum_neg;

  for(int i = 0; i < n; i++)
  {
    for(int j = 0; j < p; j++)
    {
      // loop to create indicies for subsetting
      // 1 to p, excluding j
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

      // subsetting the matricies
      arma::vec rj = R.col(j);
      arma::mat Rj = R.cols(idx);
      arma::vec zj = z.rows(idx);
      arma::vec Rj_zj = Rj * zj;

      arma::vec a_temp = a - Rj_zj;
      arma::vec b_temp = b - Rj_zj;

      arma::uvec pos = rj > 0;
      arma::uvec neg = rj < 0;

      // cout << "Pos: \n" << pos << endl;
      // cout << "Rj: \n" << rj << endl;
      // cout << "Neg: \n" << neg << endl;

      sum_pos = sum(pos);
      sum_neg = sum(neg);

      // which are negtive and positive
      arma::uvec neg_idx = find(neg == 1);
      arma::uvec pos_idx = find(pos == 1);

      // cout << "pos_idx: \n" << pos_idx << endl;
      // cout << "neg_idx: \n" << neg_idx << endl;

      if(sum(pos) == 0)
      {
        lower_pos = R_NegInf;
        upper_pos = R_PosInf;
      }
      else
      {
        temp = a_temp.rows(pos_idx) / rj.rows(pos_idx);
        lower_pos = max(temp);

        // // debug
        // std::cout << "-------------------- Lower Pos ---------------" << endl;
        // std::cout << "A_temp: \n" << a_temp.rows(pos_idx) << endl;
        // std::cout << "Rj.rows: \n" << rj.rows(pos_idx) << endl;
        // std::cout << "Temp: \n" << temp << endl;

        temp = b_temp.rows(pos_idx) / rj.rows(pos_idx);
        upper_pos = min(temp);

        // // debug
        // std::cout << "-------------------- Upper Pos ---------------" << endl;
        // std::cout << "A_temp: \n" << b_temp.rows(pos_idx) << endl;
        // std::cout << "Rj.rows: \n" << rj.rows(pos_idx) << endl;
        // std::cout << "Temp: \n" << temp << endl;
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

        // // debug
        // std::cout << "-------------------- Upper Neg ---------------" << endl;
        // std::cout << "A_temp: \n" << a_temp.rows(neg_idx) << endl;
        // std::cout << "Rj.rows: \n" << rj.rows(neg_idx) << endl;
        // std::cout << "Temp: \n" << temp << endl;

        temp = b_temp.rows(neg_idx) / rj.rows(neg_idx);
        lower_neg = max(temp);

        // // debug
        // std::cout << "-------------------- Lower Neg ---------------" << endl;
        // std::cout << "B_temp: \n" << b_temp.rows(neg_idx) << endl;
        // std::cout << "Rj.rows: \n" << rj.rows(neg_idx) << endl;
        // std::cout << "Temp: \n" << temp << endl;
      }

      double lower_j = std::max(lower_pos, lower_neg);
      double upper_j = std::min(upper_pos, upper_neg);

      // cout << "lower_pos =" << lower_pos << "lower_neg" << lower_neg << endl;
      // cout << "upper_pos =" << upper_pos << "upper_neg" << upper_neg << endl;
      // cout << "lower_j = " << lower_j << " upper_j = " << upper_j << endl;
      //
      // cout << "***********************************************" << endl;
      // cout << "iteration: " << i << " variable: " << j << endl;
      // cout << "***********************************************" << endl;


      // cout << "Z: \n " << z << endl;

      // univariate sample
      z(j) = rtuvn(1, 0, 1, lower_j, upper_j)[0];

      // cout << "Z_j: " << z(j) << endl;
    }
    keep_x.col(i) = Sigma_chol * z + Mean;
  }
  return(keep_x);
}
