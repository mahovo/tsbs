// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

const double log2pi = std::log(2.0 * M_PI);

 //' Calculate the log-density of a multivariate normal distribution
 //' @param x vector
 //' @param mean vector 
 //' @param sigma matrix
 //' @return Log-density value
 //' @keywords internal
 double dmvnorm_arma_log(const arma::vec& x, const arma::vec& mean, const arma::mat& sigma) {
   int k = x.n_elem;
   // rooti stands for root inverse.
   // The Cholesky decomposition is often considered a type of matrix square root. 
   // It decomposes a positive-definite matrix into the product of a triangular 
   // matrix and its transpose.
   // E.g Σ=U′U (where U is the Cholesky factor, or "root").
   // The inverse of the Cholesky root factor is a key component for 
   // efficiently calculating both the determinant and the Mahalanobis distance 
   // in the log-density formula.
   arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
   double rootisum = arma::sum(log(rooti.diag()));
   double constants = -(static_cast<double>(k)/2.0) * log2pi;
   arma::vec z = rooti * (x - mean);
   return constants + rootisum - 0.5 * arma::sum(z % z);
 }
 
 //' Fit a 2-State MS-VAR(1) Model via the EM Algorithm in C++
 //' @param y A (T x k) matrix of time series data.
 //' @param max_iter Maximum number of EM iterations.
 //' @param tol Convergence tolerance for the log-likelihood.
 //' @return A list with estimated parameters.
 // [[Rcpp::export]]
 Rcpp::List fit_msvar_cpp(const arma::mat& y, int max_iter = 100, double tol = 1e-6) {
   
   // --- 1. Setup & Initialization ---
   int T = y.n_rows;
   int k = y.n_cols;
   int p = 1; // AR(1) model
   
   // Create a new matrix containing only the rows from p+1 to the end (row T)
   // Autoregressive model of order p needs p previous observations to predict 
   // the current one. This means the first p observations in dataset cannot be 
   // used as targets because they don't have a complete history of p lags. 
   // The first value you can actually predict is the one at time t = p+1 (row p).
   // So y_target is a matrix with T - p rows.
   arma::mat y_target = y.rows(p, T - 1);
   
   // Declare new, empty matrix of predictors named X_lagged with dimensions 
   // matching the y_target matrix.
   // T - p is the number of rows in y_target.
   // 1 + p * k columns: 1 for the intercept, and for each of the p lags, we 
   // need all k variables.
   arma::mat X_lagged(T - p, 1 + p * k);
   X_lagged.col(0).ones();
   for (int i = 0; i < p; ++i) {
     X_lagged.cols(1 + i * k, k + i * k) = y.rows(p - 1 - i, T - 2 - i);
   }
   
   // Effective sample size. The size of the estimation sample.
   int T_eff = y_target.n_rows;
   
   // Initial parameter estimates from a single-regime VAR
   // Note: "arma" stands for "Armadillo"! Not "AutoRegressive Moving Average"!
   arma::mat beta_full = arma::solve(X_lagged, y_target);
   arma::mat resid_full = y_target - X_lagged * beta_full;
   arma::mat sigma_full = arma::cov(resid_full);
   
   // State-dependent parameters (start with small perturbations)
   arma::mat beta1 = beta_full * 0.9;
   arma::mat beta2 = beta_full * 1.1;
   arma::mat sigma1 = sigma_full * 0.9;
   arma::mat sigma2 = sigma_full * 1.1;
   
   // Transition probability matrix P
   arma::mat P = {{0.9, 0.1}, {0.1, 0.9}};
   
   double log_lik_old = -std::numeric_limits<double>::infinity();
   
   // --- 2. EM Algorithm Loop ---
   arma::mat smooth_probs(T_eff, 2);
   for (int i = 0; i < max_iter; ++i) {
     
     // --- E-STEP: Hamilton Filter & Kim Smoother ---
     arma::mat pred_probs(T_eff, 2);
     arma::mat filt_probs(T_eff, 2);
     arma::mat lik_contrib(T_eff, 1);
     
     // t = 1
     arma::vec pred_init = {0.5, 0.5};
     arma::vec dens(2);
     // use exp(log(...)) structure to handle extremely small numerical values
     dens(0) = exp(dmvnorm_arma_log(y_target.row(0).t(), (X_lagged.row(0) * beta1).t(), sigma1));
     dens(1) = exp(dmvnorm_arma_log(y_target.row(0).t(), (X_lagged.row(0) * beta2).t(), sigma2));
     
     arma::vec joint_dens = pred_init % dens;
     lik_contrib(0) = arma::sum(joint_dens);
     filt_probs.row(0) = joint_dens.t() / lik_contrib(0);
     
     // t = 2...T
     for (int t = 1; t < T_eff; ++t) {
       arma::vec pred_t = P.t() * filt_probs.row(t - 1).t();
       pred_probs.row(t) = pred_t.t();
       
       dens(0) = exp(dmvnorm_arma_log(y_target.row(t).t(), (X_lagged.row(t) * beta1).t(), sigma1));
       dens(1) = exp(dmvnorm_arma_log(y_target.row(t).t(), (X_lagged.row(t) * beta2).t(), sigma2));
       
       joint_dens = pred_t % dens;
       lik_contrib(t) = arma::sum(joint_dens);
       filt_probs.row(t) = joint_dens.t() / lik_contrib(t);
     }
     
     // Kim Smoother (Backward Pass)
     smooth_probs.row(T_eff - 1) = filt_probs.row(T_eff - 1);
     for (int t = T_eff - 2; t >= 0; --t) {
       arma::vec ratio = smooth_probs.row(t + 1).t() / pred_probs.row(t + 1).t();
       smooth_probs.row(t) = filt_probs.row(t) % (P * ratio).t();
     }
     
     // --- M-STEP: Update Parameters ---
     // Weights for State 1 and State 2
     arma::vec w1 = smooth_probs.col(0);
     arma::vec w2 = smooth_probs.col(1);
     
     // Update betas using weighted least squares
     arma::mat Xw1 = X_lagged.each_col() % w1;
     beta1 = arma::solve(Xw1.t() * X_lagged, Xw1.t() * y_target);
     
     arma::mat Xw2 = X_lagged.each_col() % w2;
     beta2 = arma::solve(Xw2.t() * X_lagged, Xw2.t() * y_target);
     
     // Update sigmas using weighted covariance
     arma::mat resid1 = y_target - X_lagged * beta1;
     arma::mat resid2 = y_target - X_lagged * beta2;
     sigma1 = (resid1.t() * (resid1.each_col() % w1)) / arma::sum(w1);
     sigma2 = (resid2.t() * (resid2.each_col() % w2)) / arma::sum(w2);
     
     // Update transition matrix P
     arma::mat trans_mat(2, 2, arma::fill::zeros);
     for (int t = 0; t < T_eff - 1; ++t) {
       arma::mat num = P % (filt_probs.row(t).t() * (smooth_probs.row(t + 1) / pred_probs.row(t + 1)));
       trans_mat += num / arma::accu(num);
     }
     P = trans_mat.each_col() / arma::sum(trans_mat, 1);
     
     // --- 3. Check Convergence ---
     // double log_lik_new = arma::sum(arma::log(lik_contrib));
     double log_lik_new = arma::as_scalar(arma::sum(arma::log(lik_contrib)));
     if (std::abs(log_lik_new - log_lik_old) < tol) {
       break;
     }
     log_lik_old = log_lik_new;
   }
   
   return Rcpp::List::create(
     Rcpp::Named("beta1") = beta1,
     Rcpp::Named("beta2") = beta2,
     Rcpp::Named("sigma1") = sigma1,
     Rcpp::Named("sigma2") = sigma2,
     Rcpp::Named("P") = P,
     Rcpp::Named("log_likelihood") = log_lik_old,
     Rcpp::Named("smoothed_probabilities") = smooth_probs
   );
 }