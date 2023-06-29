# https://stackoverflow.com/questions/57760655/generating-multivariate-gaussian-distribution-in-rcpp
# https://github.com/veradjordjilovic/screenMin/blob/master/SimulationStudyFunctions.R
Rcpp::sourceCpp(code='
  #include <RcppArmadillo.h>
  // [[Rcpp::depends(RcppArmadillo)]]
  
  using namespace Rcpp;
  
  // [[Rcpp::export]]
  arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma){
    int ncols = sigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
  }'
)
