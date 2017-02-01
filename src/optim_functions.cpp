
// notes ----------------------------------------------------------------------

// questions ------------------------------------------------------------------


#include <Rcpp.h>
using namespace Rcpp;

// function for fygeson and ritov original method -----------------------------


// [[Rcpp::export]]
NumericVector fyg_rit_orig(NumericVector surv_cens,
                           NumericVector beta_res,
                           NumericMatrix x_mat) {
  
  // additional parameters
  NumericVector sum_val(x_mat.ncol());
  
  for (int obs_i = 0; obs_i < surv_cens.size(); ++obs_i){
    for (int obs_j = 0; obs_j < surv_cens.size(); ++obs_j){
     
      sum_val = sum_val + surv_cens[obs_i] * ( x_mat(obs_i, _) - x_mat(obs_j, _) ) * (1 - 1 * (beta_res[obs_i] > beta_res[obs_j]));
      
    }
  }
  
  return sum_val;
}


// function for fygeson and ritov smooth --------------------------------------


// [[Rcpp::export]]
NumericVector fyg_rit_smooth(NumericVector surv_cens,
                             NumericVector beta_res,
                             NumericMatrix x_mat, 
                             double h) {
  
  // additional parameters
  NumericVector sum_val(x_mat.ncol());
  
  for (int obs_i = 0; obs_i < surv_cens.size(); ++obs_i){
    for (int obs_j = 0; obs_j < surv_cens.size(); ++obs_j){
      
      sum_val = sum_val + surv_cens[obs_i] * ( x_mat(obs_i, _) - x_mat(obs_j, _) ) * (1 - R::pnorm( (beta_res[obs_i] - beta_res[obs_j]) / h , 0, 1, 1, 0) );
      
    }
  }
  
  return sum_val;
}


// function to calculate A_n ------------------------------------------------

// [[Rcpp::export]]
NumericVector A_n_calc(NumericVector surv_cens,
                       NumericVector beta_res,
                       NumericMatrix x_mat, 
                       double h) {
  
  // additional parameters
  NumericMatrix A_n(x_mat.ncol(), x_mat.ncol());
  double sum_val;
  
  for (int cov_l = 0; cov_l < x_mat.ncol(); ++cov_l){
    for (int cov_m = 0; cov_m <= cov_l; ++cov_m){
      
      for (int obs_i = 0; obs_i < surv_cens.size(); ++obs_i){
        for (int obs_j = 0; obs_j < surv_cens.size(); ++obs_j){
          
          sum_val = sum_val + surv_cens[obs_i] * ( x_mat(obs_i, cov_l) - x_mat(obs_j, cov_l) ) * ( x_mat(obs_i, cov_m) - x_mat(obs_j, cov_m) ) * R::dnorm( (beta_res[obs_i] - beta_res[obs_j]) / h , 0, 1, 0);

        }
      }
      
      A_n(cov_l, cov_m) = sum_val;
      
      sum_val = 0;
    }
  }
  
  return A_n;
}

// function to calculate V_n ------------------------------------------------



// [[Rcpp::export]]
NumericVector V_n_calc(NumericVector surv_cens,
                       NumericVector beta_res,
                       NumericMatrix x_mat, 
                       double h) {
  
  // additional parameters
  NumericMatrix V_n(x_mat.ncol(), x_mat.ncol());
  double eps;
  double sum_val;
  
  for (int cov_l = 0; cov_l < x_mat.ncol(); ++cov_l){
    for (int cov_m = 0; cov_m <= cov_l; ++cov_m){
      
      for (int obs_i = 0; obs_i < surv_cens.size(); ++obs_i){
        for (int obs_j = 0; obs_j < surv_cens.size(); ++obs_j){
          for (int obs_k = 0; obs_k < surv_cens.size(); ++obs_k){
            
            if (obs_j != obs_k){
              
              eps = (surv_cens[obs_i] * (1 - R::pnorm( (beta_res[obs_i] - beta_res[obs_j]) / h , 0, 1, 1, 0) ) -
                surv_cens[obs_j] * (1 - R::pnorm( (beta_res[obs_j] - beta_res[obs_i]) / h , 0, 1, 1, 0) ) ) * 
                (surv_cens[obs_i] * (1 - R::pnorm( (beta_res[obs_i] - beta_res[obs_k]) / h , 0, 1, 1, 0) ) -
                surv_cens[obs_k] * (1 - R::pnorm( (beta_res[obs_k] - beta_res[obs_i]) / h , 0, 1, 1, 0) ) );
              
              sum_val = sum_val + ( x_mat(obs_i, cov_l) - x_mat(obs_j, cov_l) ) * ( x_mat(obs_i, cov_m) - x_mat(obs_k, cov_m) ) * eps;
              
              eps = 0;
            }
            
          }
        }
      }
      
      V_n(cov_l, cov_m) = sum_val;
      
      sum_val = 0;
    }
  }
  
  return V_n;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
