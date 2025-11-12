// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <chrono> // For timing
#include <iomanip> // For std::setw and std::setfill

// ---- MAIN C++ FITTING FUNCTION ----

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <chrono> // For timing
#include <iomanip> // For std::setw and std::setfill

// --- MAIN C++ FITTING FUNCTION ---
//' @title Fit a general MS-ARMA-GARCH model via the EM Algorithm in C++
//' @description Internal C++ function to orchestrate the EM estimation.
//'              Delegates all statistical calculations (likelihood,
//'              parameter estimation) to helper functions in R to interface
//'              with the tsgarch and tsmarch packages.
//' @param y A (T x k) matrix of (differenced) time series data.
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
  int max_iter = Rcpp::as<int>(control["max_iter"]);
  double tol = Rcpp::as<double>(control["tol"]);
  
  Rcpp::List model_fits(M);
  arma::mat P(M, M);
  arma::mat smooth_probs(T, M);
  Rcpp::List all_warnings; 
  
  P.fill(1.0 / M);
  for(int j=0; j<M; ++j) { P(j,j) = 0.9; }
  P = P.each_row([](arma::rowvec& r){ r /= arma::sum(r); });
  
  for (int j = 0; j < M; ++j) {
    model_fits[j] = Rcpp::as<Rcpp::List>(spec[j])["start_pars"];
  }
  
  double log_lik_old = -std::numeric_limits<double>::infinity();
  
  // Rcpp::Environment global = Rcpp::Environment::global_env();
  // Rcpp::Function calculate_loglik_vector_r = global["calculate_loglik_vector_r"];
  // Rcpp::Function perform_m_step_parallel_r = global["perform_m_step_parallel_r"];
  Rcpp::Environment pkg_env = Rcpp::Environment::namespace_env("tsbs");
  Rcpp::Function calculate_loglik_vector_r = pkg_env["calculate_loglik_vector_r"];
  Rcpp::Function perform_m_step_parallel_r = pkg_env["perform_m_step_parallel_r"];
  
  // --- 2. EM ALGORITHM LOOP ---
  for (int iter = 0; iter < max_iter; ++iter) {
    auto start = std::chrono::high_resolution_clock::now();
    Rcpp::checkUserInterrupt();
    
    // --- E-STEP ---
    arma::mat cond_logliks(T, M, arma::fill::zeros);
    for (int j = 0; j < M; ++j) {
      arma::vec ll_vec = Rcpp::as<arma::vec>(calculate_loglik_vector_r(y, model_fits[j], spec[j], model_type));
      cond_logliks.col(j) = ll_vec;
    }
    arma::mat cond_dens = arma::exp(cond_logliks);
    
    // Hamilton Filter & Kim Smoother
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
    
    
    // DIAGNOSTIC begin
    double log_lik_before_mstep = arma::sum(arma::log(lik_contrib));
    // DIAGNOSTIC end
    
    
    // --- M-STEP ---
    // Update transition matrix P
    arma::mat trans_mat(M, M, arma::fill::zeros);
    for (int t = 0; t < T - 1; ++t) {
      arma::mat num = P % (filt_probs.row(t).t() * (smooth_probs.row(t + 1) / pred_probs.row(t + 1)));
      trans_mat += num / arma::accu(num);
    }
    P = trans_mat.each_col() / arma::sum(trans_mat, 1);
    
    // --- Call the single parallel R function for the M-step ---
    model_fits = Rcpp::as<Rcpp::List>(perform_m_step_parallel_r(y, smooth_probs, spec, model_type));
    
    // === DIAGNOSTIC: Check if M-step improved likelihood ===
    Rcpp::Rcout << "\n=== M-STEP CHECK (Iteration " << (iter + 1) << ") ===" << std::endl;
    
    // Recompute likelihoods with NEW parameters
    arma::mat cond_logliks_new(T, M, arma::fill::zeros);
    for (int j = 0; j < M; ++j) {
      arma::vec ll_vec = Rcpp::as<arma::vec>(
        calculate_loglik_vector_r(y, model_fits[j], spec[j], model_type)
      );
      cond_logliks_new.col(j) = ll_vec;
    }
    
    // Recompute total LL with new parameters
    arma::mat cond_dens_new = arma::exp(cond_logliks_new);
    arma::vec pred_init_new(M, arma::fill::ones);
    pred_init_new /= M;
    arma::vec lik_contrib_new(T);
    arma::vec joint_dens_new = pred_init_new % cond_dens_new.row(0).t();
    lik_contrib_new(0) = arma::sum(joint_dens_new);
    arma::vec filt_t = joint_dens_new / lik_contrib_new(0);
    
    for (int t = 1; t < T; ++t) {
      arma::vec pred_t = P.t() * filt_t;
      joint_dens_new = pred_t % cond_dens_new.row(t).t();
      lik_contrib_new(t) = arma::sum(joint_dens_new);
      filt_t = joint_dens_new / lik_contrib_new(t);
    }
    
    double log_lik_after_mstep = arma::sum(arma::log(lik_contrib_new));
    
    Rcpp::Rcout << "  LL before M-step: " << log_lik_before_mstep << std::endl;
    Rcpp::Rcout << "  LL after M-step:  " << log_lik_after_mstep << std::endl;
    Rcpp::Rcout << "  Change: " << (log_lik_after_mstep - log_lik_before_mstep) << std::endl;
    
    if (log_lik_after_mstep < log_lik_before_mstep - 1e-6) {
      Rcpp::Rcout << "  *** WARNING: M-step DECREASED LL! ***" << std::endl;
    }
    
    // Print parameter values for diagnosis
    for (int j = 0; j < M; ++j) {
      Rcpp::Rcout << "  State " << (j+1) << " params after M-step:" << std::endl;
      Rcpp::List state_pars = model_fits[j];
      
      // Print GARCH parameters
      if (state_pars.containsElementNamed("garch_pars")) {
        Rcpp::List garch_pars = state_pars["garch_pars"];
        Rcpp::Rcout << "    garch_pars length: " << garch_pars.size() << std::endl;
        
        for (int s = 0; s < garch_pars.size(); ++s) {
          Rcpp::List series_pars = garch_pars[s];
          Rcpp::Rcout << "      Series " << (s+1) << ": ";
          if (series_pars.containsElementNamed("alpha1")) {
            Rcpp::Rcout << "alpha1=" << Rcpp::as<double>(series_pars["alpha1"]) << " ";
          }
          if (series_pars.containsElementNamed("beta1")) {
            Rcpp::Rcout << "beta1=" << Rcpp::as<double>(series_pars["beta1"]) << " ";
          }
          Rcpp::Rcout << std::endl;
        }
      }
      
      // Print DCC parameters
      if (state_pars.containsElementNamed("alpha_1")) {
        Rcpp::Rcout << "    DCC alpha_1: " << Rcpp::as<double>(state_pars["alpha_1"]) << std::endl;
      }
      if (state_pars.containsElementNamed("beta_1")) {
        Rcpp::Rcout << "    DCC beta_1: " << Rcpp::as<double>(state_pars["beta_1"]) << std::endl;
      }
    }
    
    // Update cond_logliks for next E-step
    cond_logliks = cond_logliks_new;
    
    
    // ================ END Check if M-step improved likelihood  ===============
    
    
    // --- 3. CHECK CONVERGENCE & PROVIDE FEEDBACK ---
    
    // === DIAGNOSTIC begin ===

    // Next line replaced with diagnostic code.
    //double log_lik_new = arma::sum(arma::log(lik_contrib));
    double log_lik_new = log_lik_after_mstep;  // Use the post-M-step LL
    
    // === DIAGNOSTIC end ===
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    int total_seconds = static_cast<int>(elapsed.count());
    int hours = total_seconds / 3600;
    int minutes = (total_seconds % 3600) / 60;
    int seconds = total_seconds % 60;
    
    Rcpp::Rcout << "EM Iteration: " << iter + 1 
                << ", Log-Likelihood: " << log_lik_new 
                << ", Duration: "
                << std::setw(2) << std::setfill('0') << hours << ":"
                << std::setw(2) << std::setfill('0') << minutes << ":"
                << std::setw(2) << std::setfill('0') << seconds << std::endl;
    
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
    Rcpp::Named("convergence") = Rcpp::List::create(Rcpp::Named("final_loglik") = log_lik_old),
    Rcpp::Named("warnings") = Rcpp::List() // Warning collection removed for simplicity
  );
}


