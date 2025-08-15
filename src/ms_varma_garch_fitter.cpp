// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

// --- MAIN C++ FITTING FUNCTION ---

//' @title Fit a general MS-VARMA-GARCH model via the EM Algorithm in C++
 //' @description Internal C++ function to orchestrate the EM estimation.
 //'              This version delegates all statistical calculations (likelihood,
 //'              parameter estimation) to helper functions in R to correctly
 //'              interface with TMB-based packages like tsgarch.
 //' @param y A (T x k) matrix of (potentially differenced) time series data.
 //' @param M The number of states.
 //' @param spec A list of model specifications from R.
 //' @param model_type "univariate" or "multivariate".
 //' @param control A list with max_iter and tol.
 //' @return A raw list with estimated parameters and results.
 // [[Rcpp::export]]
 Rcpp::List fit_ms_varma_garch_cpp(
     const arma::mat& y,
     int M,
     Rcpp::List spec,
     std::string model_type,
     Rcpp::List control
 ) {
   // --- 1. SETUP & INITIALIZATION ---
   int T = y.n_rows;
   int k = y.n_cols;
   int max_iter = Rcpp::as<int>(control["max_iter"]);
   double tol = Rcpp::as<double>(control["tol"]);
   
   // Parameters to be estimated
   Rcpp::List model_fits(M);
   arma::mat P(M, M);
   arma::mat smooth_probs(T, M);
   
   // Initialize parameters
   P.fill(1.0 / M);
   for(int j=0; j<M; ++j) { P(j,j) = 0.9; }
   P = P.each_row([](arma::rowvec& r){ r /= arma::sum(r); });
   
   // Initialize model_fits with starting parameters provided in the spec.
   for (int j = 0; j < M; ++j) {
     model_fits[j] = Rcpp::as<Rcpp::List>(spec[j])["start_pars"];
   }
   
   double log_lik_old = -std::numeric_limits<double>::infinity();
   
   // Prepare R helper functions
   Rcpp::Environment global = Rcpp::Environment::global_env();
   Rcpp::Function calculate_loglik_vector_r = global["calculate_loglik_vector_r"];
   Rcpp::Function estimate_arma_weighted_r = global["estimate_arma_weighted_r"];
   Rcpp::Function estimate_garch_weighted_r = global["estimate_garch_weighted_r"];
   
   // --- 2. EM ALGORITHM LOOP ---
   for (int iter = 0; iter < max_iter; ++iter) {
     Rcpp::checkUserInterrupt();
     
     // --- E-STEP: Calculate Likelihoods in R, Filter/Smooth in C++ ---
     arma::mat cond_logliks(T, M, arma::fill::zeros);
     
     for (int j = 0; j < M; ++j) {
       // Call R helper to get the log-likelihood vector for state j
       arma::vec ll_vec = Rcpp::as<arma::vec>(calculate_loglik_vector_r(y, model_fits[j], spec[j]));
       cond_logliks.col(j) = ll_vec;
     }
     arma::mat cond_dens = arma::exp(cond_logliks);
     
     // Hamilton Filter & Kim Smoother (fast matrix algebra in C++)
     arma::mat filt_probs(T, M);
     arma::mat pred_probs(T, M);
     arma::vec lik_contrib(T);
     arma::vec pred_init(M, arma::fill::ones);
     pred_init /= M;
     arma::vec joint_dens = pred_init % cond_dens.row(0).t();
     lik_contrib(0) = arma::sum(joint_dens);
     filt_probs.row(0) = joint_dens.t() / lik_contrib(0);
     for (int t = 1; t < T; ++t) {
       arma::vec pred_t = P.t() * filt_probs.row(t - 1).t();
       pred_probs.row(t) = pred_t.t();
       joint_dens = pred_t % cond_dens.row(t).t();
       lik_contrib(t) = arma::sum(joint_dens);
       filt_probs.row(t) = joint_dens.t() / lik_contrib(t);
     }
     smooth_probs.row(T - 1) = filt_probs.row(T - 1);
     for (int t = T - 2; t >= 0; --t) {
       arma::vec ratio = smooth_probs.row(t + 1).t() / pred_probs.row(t + 1).t();
       smooth_probs.row(t) = filt_probs.row(t) % (P * ratio).t();
     }
     
     // --- M-STEP (2-Step via R calls): Update Parameters ---
     // Update transition matrix P
     arma::mat trans_mat(M, M, arma::fill::zeros);
     for (int t = 0; t < T - 1; ++t) {
       arma::mat num = P % (filt_probs.row(t).t() * (smooth_probs.row(t + 1) / pred_probs.row(t + 1)));
       trans_mat += num / arma::accu(num);
     }
     P = trans_mat.each_col() / arma::sum(trans_mat, 1);
     
     for (int j = 0; j < M; ++j) {
       arma::vec weights_j = smooth_probs.col(j);
       Rcpp::List state_spec_j = Rcpp::as<Rcpp::List>(spec[j]);
       
       // M-Step 1: Update Mean Parameters by calling R helper
       Rcpp::List new_arma_fit = Rcpp::as<Rcpp::List>(estimate_arma_weighted_r(y, weights_j, state_spec_j));
       
       // M-Step 2: Update Variance Parameters by calling R helper
       arma::mat residuals = Rcpp::as<arma::mat>(new_arma_fit["residuals"]);
       Rcpp::List new_garch_fit = Rcpp::as<Rcpp::List>(estimate_garch_weighted_r(residuals, weights_j, state_spec_j));
       
       // Update the full parameter set for the state
       Rcpp::as<Rcpp::List>(model_fits[j])["arma_pars"] = new_arma_fit["coefficients"];
       Rcpp::as<Rcpp::List>(model_fits[j])["garch_pars"] = new_garch_fit["coefficients"];
     }
     
     // --- 3. CHECK CONVERGENCE ---
     double log_lik_new = arma::sum(arma::log(lik_contrib));
     if (std::abs(log_lik_new - log_lik_old) < tol) {
       log_lik_old = log_lik_new;
       break;
     }
     log_lik_old = log_lik_new;
   }
   
   // --- 4. RETURN RESULTS ---
   return Rcpp::List::create(
     Rcpp::Named("model_fits") = model_fits,
     Rcpp::Named("P") = P,
     Rcpp::Named("log_likelihood") = log_lik_old,
     Rcpp::Named("smoothed_probabilities") = smooth_probs,
     Rcpp::Named("convergence") = Rcpp::List::create(Rcpp::Named("final_loglik") = log_lik_old)
   );
 }