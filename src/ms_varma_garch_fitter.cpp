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
    Rcpp::List control,
    Rcpp::Nullable<Rcpp::List> diagnostics = R_NilValue,
    bool verbose = false
) {
  // --- 1. SETUP & INITIALIZATION ---
  int T = y.n_rows;
  int max_iter = Rcpp::as<int>(control["max_iter"]);
  double tol = Rcpp::as<double>(control["tol"]);
  
  Rcpp::List model_fits(M);
  arma::mat P(M, M);
  arma::mat smooth_probs(T, M);
  Rcpp::List all_warnings; 
  
  // Initialize diagnostics if provided
  bool collect_diagnostics = diagnostics.isNotNull();
  Rcpp::List diag_collector;
  if (collect_diagnostics) {
    diag_collector = Rcpp::as<Rcpp::List>(diagnostics);
  }
  
  P.fill(1.0 / M);
  for(int j=0; j<M; ++j) { P(j,j) = 0.9; }
  P = P.each_row([](arma::rowvec& r){ r /= arma::sum(r); });
  
  for (int j = 0; j < M; ++j) {
    model_fits[j] = Rcpp::as<Rcpp::List>(spec[j])["start_pars"];
  }
  
  double log_lik_old = -std::numeric_limits<double>::infinity();
  
  Rcpp::Environment pkg_env = Rcpp::Environment::namespace_env("tsbs");
  Rcpp::Function calculate_loglik_vector_r = pkg_env["calculate_loglik_vector_r"];
  //Rcpp::Function perform_m_step_parallel_r = pkg_env["perform_m_step_parallel_r"];
  Rcpp::Function perform_m_step_r = pkg_env["perform_m_step_r"];
  
  // Extract DCC control parameters
  double dcc_threshold = Rcpp::as<double>(control["dcc_boundary_threshold"]);
  std::string dcc_criterion = Rcpp::as<std::string>(control["dcc_boundary_criterion"]);
  
  // Load diagnostic functions if needed
  // Rcpp::Function add_em_iteration_diagnostic = R_NilValue;
  // Rcpp::Function add_parameter_evolution = R_NilValue;
  // Rcpp::Function add_diagnostic_warning = R_NilValue;
  
  if (collect_diagnostics) {
    // add_em_iteration_diagnostic = pkg_env["add_em_iteration_diagnostic"];
    // add_parameter_evolution = pkg_env["add_parameter_evolution"];
    // add_diagnostic_warning = pkg_env["add_diagnostic_warning"];
    Rcpp::Function add_em_iteration_diagnostic = pkg_env["add_em_iteration_diagnostic"];
    Rcpp::Function add_parameter_evolution = pkg_env["add_parameter_evolution"];
    Rcpp::Function add_diagnostic_warning = pkg_env["add_diagnostic_warning"];
  }
  
  // --- 2. EM ALGORITHM LOOP ---
  for (int iter = 0; iter < max_iter; ++iter) {
    auto start = std::chrono::high_resolution_clock::now();
    Rcpp::checkUserInterrupt();
    
    Rcpp::Rcout << "EM Iteration " << iter + 1 << "... ";
    
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
    
    double log_lik_before_mstep = arma::sum(arma::log(lik_contrib));
    
    // --- M-STEP ---
    // Update transition matrix P
    arma::mat trans_mat(M, M, arma::fill::zeros);
    for (int t = 0; t < T - 1; ++t) {
      arma::mat num = P % (filt_probs.row(t).t() * (smooth_probs.row(t + 1) / pred_probs.row(t + 1)));
      trans_mat += num / arma::accu(num);
    }
    P = trans_mat.each_col() / arma::sum(trans_mat, 1);
    
    // Call M-step with diagnostics if collecting
    /*
    if (collect_diagnostics) {
      model_fits = Rcpp::as<Rcpp::List>(perform_m_step_parallel_r(
        y, 
        smooth_probs, 
        spec, 
        model_type,
        Rcpp::Named("diagnostics") = diag_collector,
        Rcpp::Named("iteration") = iter + 1,
        Rcpp::Named("verbose") = verbose
      ));
    } else {
      model_fits = Rcpp::as<Rcpp::List>(perform_m_step_parallel_r(
      
        y, 
        smooth_probs, 
        spec, 
        model_type,
        Rcpp::Named("verbose") = verbose
      ));
    }
     */
    Rcpp::List m_step_result;
    
    if (collect_diagnostics) {
      m_step_result = Rcpp::as<Rcpp::List>(perform_m_step_r(
        y, 
        smooth_probs, 
        spec, 
        model_type,
        Rcpp::Named("diagnostics") = diag_collector,
        Rcpp::Named("iteration") = iter + 1,
        Rcpp::Named("verbose") = verbose,
        Rcpp::Named("dcc_threshold") = dcc_threshold,
        Rcpp::Named("dcc_criterion") = dcc_criterion
      ));
      
      // Extract fits and diagnostics
      model_fits = Rcpp::as<Rcpp::List>(m_step_result["fits"]);
      diag_collector = Rcpp::as<Rcpp::List>(m_step_result["diagnostics"]);
      
    } else {
      m_step_result = Rcpp::as<Rcpp::List>(perform_m_step_r(
        y, 
        smooth_probs, 
        spec, 
        model_type,
        Rcpp::Named("verbose") = verbose,
        Rcpp::Named("dcc_threshold") = dcc_threshold,
        Rcpp::Named("dcc_criterion") = dcc_criterion
      ));
      
      // Extract just fits (no diagnostics)
      model_fits = Rcpp::as<Rcpp::List>(m_step_result["fits"]);
    }
    
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
    double ll_change = log_lik_after_mstep - log_lik_before_mstep;
    
    // Update cond_logliks for next E-step
    cond_logliks = cond_logliks_new;
    
    // Calculate iteration duration
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    double duration_sec = elapsed.count();
    
    // Output to console
    if (verbose) {
      Rcpp::Rcout << "\n=== EM ITERATION " << (iter + 1) << " ===" << std::endl;
      Rcpp::Rcout << "  LL before M-step: " << log_lik_before_mstep << std::endl;
      Rcpp::Rcout << "  LL after M-step:  " << log_lik_after_mstep << std::endl;
      Rcpp::Rcout << "  Change: " << ll_change << std::endl;
      
      if (ll_change < -1e-6) {
        Rcpp::Rcout << "  *** WARNING: M-step DECREASED LL! ***" << std::endl;
      }
    }
    
    int total_seconds = static_cast<int>(duration_sec);
    int hours = total_seconds / 3600;
    int minutes = (total_seconds % 3600) / 60;
    int seconds = total_seconds % 60;
    
    //Rcpp::Rcout << "EM Iteration: " << iter + 1 
    Rcpp::Rcout << "Log-Likelihood: " << log_lik_after_mstep 
                << " (Duration: "
                << std::setw(2) << std::setfill('0') << hours << ":"
                << std::setw(2) << std::setfill('0') << minutes << ":"
                << std::setw(2) << std::setfill('0') << seconds << ")" 
                << std::endl;
    
    // Collect diagnostics
    if (collect_diagnostics) {
      Rcpp::Function add_em_iteration_diagnostic = pkg_env["add_em_iteration_diagnostic"];
      Rcpp::Function add_parameter_evolution = pkg_env["add_parameter_evolution"];
      Rcpp::Function add_diagnostic_warning = pkg_env["add_diagnostic_warning"];
      
      diag_collector = Rcpp::as<Rcpp::List>(add_em_iteration_diagnostic(
        diag_collector, // diagnostics
        iter + 1, // iteration
        log_lik_before_mstep, // log_lik_before
        log_lik_after_mstep, // log_lik_after
        ll_change, // ll_change
        duration_sec, // duration_sec
        model_fits, // parameters
        false  // converged flag
      ));
      
      // Collect parameter evolution for each state
      for (int j = 0; j < M; ++j) {
        diag_collector = Rcpp::as<Rcpp::List>(add_parameter_evolution(
          diag_collector,
          iter + 1,
          j + 1,
          model_fits[j]
        ));
      }
      
      // Add warning if LL decreased
      if (ll_change < -1e-6) {
        diag_collector = Rcpp::as<Rcpp::List>(add_diagnostic_warning(
          diag_collector,
          iter + 1,
          "ll_decrease",
          "M-step decreased log-likelihood",
          Rcpp::List::create(Rcpp::Named("decrease") = ll_change)
        ));
      }
    }
    
    // Update log_lik_old
    double log_lik_new = log_lik_after_mstep;
    
    // --- 3. CHECK CONVERGENCE ---
    if (std::abs(log_lik_new - log_lik_old) < tol) {
      if (verbose) {
        Rcpp::Rcout << "\n=== CONVERGED at iteration " << (iter + 1) << " ===" << std::endl;
      }
      
      if (collect_diagnostics) {
        // Mark last iteration as converged
        Rcpp::List last_iter = diag_collector["em_iterations"];
        if (last_iter.size() > 0) {
          Rcpp::List final_iter = last_iter[last_iter.size() - 1];
          final_iter["converged"] = true;
          last_iter[last_iter.size() - 1] = final_iter;
          diag_collector["em_iterations"] = last_iter;
        }
      }
      
      log_lik_old = log_lik_new;
      break;
    }
    log_lik_old = log_lik_new;
  }
  
  // --- 4. RETURN RESULTS ---
  Rcpp::List result = Rcpp::List::create(
    Rcpp::Named("model_fits") = model_fits,
    Rcpp::Named("P") = P,
    Rcpp::Named("log_likelihood") = log_lik_old,
    Rcpp::Named("smoothed_probabilities") = smooth_probs,
    Rcpp::Named("convergence") = Rcpp::List::create(Rcpp::Named("final_loglik") = log_lik_old),
    Rcpp::Named("warnings") = Rcpp::List()
  );
  
  if (collect_diagnostics) {
    result["diagnostics"] = diag_collector;
  }
  
  return result;
}


