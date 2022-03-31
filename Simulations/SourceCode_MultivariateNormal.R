### Function for simulating from a multivariate normal distribution with big n
### This is done via C++
###

## C++ code to sample from a high-dimensional multivariate normal distribution
rtools <- "C:\\Rtools\\bin"
gcc <- "C:\\Rtools\\gcc-4.6.3\\bin"
path <- strsplit(Sys.getenv("PATH"), ";")[[1]]
new_path <- c(rtools, gcc, path)
new_path <- new_path[!duplicated(tolower(new_path))]
Sys.setenv(PATH = paste(new_path, collapse = ";"))

code <- '
  using namespace Rcpp;
  int n = as<int>(n_);
  arma::vec mu = as<arma::vec>(mu_);
  arma::mat sigma = as<arma::mat>(sigma_);
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return wrap(arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma));
'
rmvnorm.rcpp <- cxxfunction(signature(n_="integer", mu_="numeric", 
                                      sigma_="matrix"), code,
                            plugin="RcppArmadillo", verbose=TRUE)