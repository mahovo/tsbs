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

// +++++++++++++ Compute weights for tapered blocks +++++++++++++
// Cosine tapering window function.
Rcpp::NumericVector cosine_weights(int len) {
  Rcpp::NumericVector w(len);
  for (int i = 0; i < len; ++i) {
    w[i] = 0.5 * (1 - std::cos(2 * M_PI * i / (len - 1)));
  }
  return w;
}


// Bartlett (triangular) window. 
// Linearly tapers from 1 in the center to 0 at both ends.
Rcpp::NumericVector bartlett_weights(int block_length) {
  Rcpp::NumericVector w(block_length);
  double N = static_cast<double>(block_length - 1);
  for (int i = 0; i < block_length; ++i) {
    w[i] = 1.0 - std::abs((i - N / 2.0) / (N / 2.0));
  }
  return w;
}


// Tukey window (a.k.a. tapered cosine window).
// Allows tuning of the taper via alpha ∈ [0, 1]:
// alpha = 0 → rectangular (no taper)
// alpha = 1 → Hann (fully tapered) 
Rcpp::NumericVector tukey_weights(int block_length, double alpha = 0.5) {
  Rcpp::NumericVector w(block_length);
  double N = static_cast<double>(block_length - 1);
  
  for (int i = 0; i < block_length; ++i) {
    double n = static_cast<double>(i);
    if (n < alpha * N / 2.0) {
      w[i] = 0.5 * (1.0 + std::cos(M_PI * ((2.0 * n) / (alpha * N) - 1.0)));
    } else if (n <= N * (1.0 - alpha / 2.0)) {
      w[i] = 1.0;
    } else {
      w[i] = 0.5 * (1.0 + std::cos(M_PI * ((2.0 * n) / (alpha * N) - (2.0 / alpha) + 1.0)));
    }
  }
  
  return w;
}



// [[Rcpp::export]]
Rcpp::List blockBootstrap_cpp(
  SEXP xSEXP,
  SEXP n_boot_spec,
  SEXP block_length_spec,
  const std::string &bs_type,
  const std::string &block_type,
  const std::string &taper_type,
  const double &tukey_alpha,
  SEXP num_blocks_spec,
  const int num_boots,
  SEXP p,
  const double stationary_max_percentile,
  const double stationary_max_fraction_of_n
) {
  
  // All inputs are validated in tsbs().R before blockBootstrap_cpp() is called.
  
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
      : Rcpp::as<int>(block_length_spec);
  
  if (block_length == 1 && block_type == "tapered") {
   Rcpp::stop("Can not use block type \"tapered\" when block length is 1. See `?tsbs::tsbs`");
  }

  // Final output length
  int n_boot;
  int num_blocks;
  if (!Rf_isNull(n_boot_spec)) {
    n_boot = Rcpp::as<int>(n_boot_spec);
  } else if (!Rf_isNull(num_blocks_spec)) {
    num_blocks = Rcpp::as<int>(num_blocks_spec);
    n_boot = num_blocks * block_length;
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
      
      if (bs_type == "moving") {
        current_block_len = block_length;
        
        if (block_type == "overlapping") {
          start_idx = static_cast<int>(R::runif(0, n - block_length + 1));
        } else if (block_type == "non-overlapping") {
          // Pick a block number from non-overlapping partitions
          int max_start = (n / block_length) * block_length;
          int num_starts = max_start / block_length;
          int block_num = static_cast<int>(R::runif(0, num_starts));
          start_idx = block_num * block_length;
        } else if (block_type == "tapered") {
          start_idx = static_cast<int>(R::runif(0, n - block_length + 1));
          current_block_len = block_length;
          Rcpp::NumericVector w = cosine_weights(block_length);
          
          int rows_to_copy = std::min(current_block_len, n_boot - pos);
          for (int i = 0; i < rows_to_copy; ++i) {
            int idx = (start_idx + i) % n;
            for (int col = 0; col < d; ++col) {
              sample(pos + i, col) = x(idx, col) * w[i];
            }
          }
          pos += rows_to_copy;
        }
        
      } else if (bs_type == "stationary") {
        double prob = Rcpp::as<double>(p);
        
        // Compute maximum block length
        int max_from_percentile = R::qgeom(stationary_max_percentile, prob, 1, 0) + 1;
        int max_from_n_fraction = static_cast<int>(stationary_max_fraction_of_n * n);
        int max_block_length = std::min(max_from_percentile, max_from_n_fraction);
        max_block_length = std::max(1, std::min(max_block_length, n)); // keep within [1, n]
        
        // Draw random block length
        current_block_len = static_cast<int>(R::rgeom(prob)) + 1;
        current_block_len = std::min(current_block_len, max_block_length);
        
        if (block_type == "overlapping") {
          start_idx = static_cast<int>(R::runif(0, n));
          
        } else if (block_type == "non-overlapping") {
          int max_start = (n / current_block_len) * current_block_len;
          int num_starts = max_start / current_block_len;
          int block_num = static_cast<int>(R::runif(0, num_starts));
          start_idx = block_num * current_block_len;
          
        } else if (block_type == "tapered") {
          start_idx = static_cast<int>(R::runif(0, n));
          
          Rcpp::NumericVector w;
          if (taper_type == "cosine") {
            w = cosine_weights(current_block_len);
          } else if (taper_type == "bartlett") {
            w = bartlett_weights(current_block_len);
          } else if (taper_type == "tukey") {
            w = tukey_weights(current_block_len, tukey_alpha);
          } else {
            Rcpp::stop("Unknown taper_type.");
          }
          
          int rows_to_copy = std::min(current_block_len, n_boot - pos);
          for (int i = 0; i < rows_to_copy; ++i) {
            int idx = (start_idx + i) % n;
            for (int col = 0; col < d; ++col) {
              sample(pos + i, col) = x(idx, col) * w[i];
            }
          }
          pos += rows_to_copy;
          continue;  // Skip generic copy logic below, already done
        }
      } else {
        Rcpp::stop("Unsupported type. Use 'moving' or 'stationary'.");
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
    boots[b] = sample( Range(0, n_boot - 1), _ );
  }
  
  return boots;
}


