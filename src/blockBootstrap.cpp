#include <Rcpp.h>
using namespace Rcpp;

// Helper: Compute default block length
// [[Rcpp::export]]
int compute_default_block_length(const Rcpp::NumericMatrix &x) {
  int n = x.nrow();
  Rcpp::NumericVector ac1(x.ncol());
  for (int j = 0; j < x.ncol(); ++j) {
    double mu = mean(x(_, j));
    double num = 0.0, denom = 0.0;
    
    for (int t = 1; t < n; ++t) {
      num += (x(t, j) - mu) * (x(t - 1, j) - mu);
      denom += (x(t, j) - mu) * (x(t, j) - mu);
    }
    
    ac1[j] = (denom != 0.0) ? std::abs(num / denom) : 0.0;
  }
  
  double rho1 = mean(ac1);
  int candidate = std::floor(10.0 / (1.0 - rho1));
  return std::max(5, std::min(candidate, (int) std::sqrt(n)));
}

// [[Rcpp::export]]
Rcpp::List blockBootstrap_cpp(SEXP xSEXP,
                          SEXP n_boot_spec,
                          SEXP block_length_spec,
                          SEXP num_blocks_spec,
                          const int num_boots,
                          const std::string &block_type,
                          SEXP p_spec,
                          const bool overlap) {
  
  // All inputs are validated in tsbs.R before blockBootstrap_cpp is called.
  
  NumericMatrix x(xSEXP);
  
  int n = x.nrow();
  int d = x.ncol();
  
  // ---- Block length ----
  
  // If our Rcpp function takes const Rcpp::NumericVector &block_length_spec as 
  // input, and it receives NULL, it will fail. So if we instead specify 
  // SEXP block_length_spec and then do 
  //
  //  int block_length = Rf_isNull(block_length_spec)
  //  ? compute_default_block_length(x)
  //    : as<int>(block_length_spec);
  //
  // then that allows us to accept the NULL input, compute a valid value if NULL, 
  // and specify integer type at that point for further processing.
  int block_length = Rf_isNull(block_length_spec)
    ? compute_default_block_length(x)
      : as<int>(block_length_spec);

  // Number of blocks
  int num_blocks = (num_blocks_spec.size() > 0 && !Rcpp::NumericVector::is_na(num_blocks_spec[0]))
    ? static_cast<int>(num_blocks_spec[0]) : -1;
  
  // Final output length
  int n_boot;
  if (num_blocks > 0) {
    n_boot = num_blocks * block_length;
  } else if (n_boot_spec.size() > 0 && !Rcpp::NumericVector::is_na(n_boot_spec[0])) {
    n_boot = static_cast<int>(n_boot_spec[0]);
  } else {
    n_boot = n;
  }
  
  Rcpp::List boots(num_boots);
  for (int b = 0; b < num_boots; ++b) {
    Rcpp::NumericMatrix sample(n_boot, d);
    int pos = 0;
    
    while (pos < n_boot) {
      int start_idx;
      int current_block_len;
      
      if (block_type == "moving") {
        start_idx = static_cast<int>(R::runif(0, n - block_length + 1));
        current_block_len = block_length;
        
      } else if (block_type == "stationary") {
        if (p <= 0.0 || p >= 1.0) {
          Rcpp::stop("'p' must be in the range (0, 1) for 'stationary' block_type.");
        }
        current_block_len = static_cast<int>(R::rgeom(p)) + 1;
        if (current_block_len > block_length) {
          current_block_len = block_length;
        }
        if (overlap) {
          start_idx = static_cast<int>(R::runif(0, n));
        } else {
          int max_start = (n / block_length) * block_length;
          int num_starts = max_start / block_length;
          int block_num = static_cast<int>(R::runif(0, num_starts));
          start_idx = block_num * block_length;
        }
        
      } else {
        Rcpp::stop("Unsupported block_type. Use 'moving' or 'stationary'.");
      }
      
      int rows_to_copy = std::min(current_block_len, n_boot - pos);
      for (int i = 0; i < rows_to_copy; ++i) {
        int idx = (start_idx + i) % n;
        for (int col = 0; col < d; ++col) {
          sample(pos + i, col) = x(idx, col);
        }
      }
      pos += rows_to_copy;
    }
    boots[b] = sample;
  }
  
  return boots;
}
