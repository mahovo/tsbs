#' Hidden Markov Model (HMM) Bootstrap for Multivariate Time Series
#'
#' Fits a Gaussian Hidden Markov Model (HMM) to a multivariate time series
#'   and generates bootstrap replicates by resampling regime-specific blocks.
#'
#' @param x Numeric vector representing the time series.
#' @param n_boot Length of bootstrap series.
#' @param num_states Integer number of hidden states for the HMM.
#' @param num_blocks Integer number of blocks to sample for each bootstrap
#'   replicate.
#' @param num_boots Integer number of bootstrap replicates to generate.
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores.
#' @param verbose Logical. If TRUE, prints HMM fitting information and warnings.
#'   Defaults to FALSE.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits a Gaussian HMM to `x` using `depmixS4::depmix()` and 
#'     `depmixS4::fit()`.
#'   \item Uses Viterbi decoding (`posterior(fit, type = "viterbi")$state`)
#'     to assign each observation to a state.
#'   \item Samples contiguous blocks of observations belonging to each state.
#' }
#' 
#' If `n_boot` is set, the last block will be truncated when necessary to match 
#' the length (`n_boot`) of the bootstrap series. This is the only way to ensure 
#' equal length of all bootstrap series, as the length of each block is random. 
#' If `n_boot` is not set, `num_blocks` must be set, and the length of each 
#' bootstrap series will be determined by the number of blocks and the random 
#' lengths of the individual blocks for that particular series. This almost 
#' certainly results in bootstrap series of different lengths.  
#'   
#' For multivariate series (matrices or data frames), the function fits a single
#' HMM where all variables are assumed to depend on the same underlying hidden
#' state sequence. The returned bootstrap samples are matrices with the same
#' number of columns as the input `x`.  
#'   
#' Hidden Markov Model definition:  
#' 
#' - \eqn{T}: sequence length  
#' - \eqn{K}: number of hidden states  
#' - \eqn{\mathbf{X} = (X_1, \dots, X_T)}: observed sequence  
#' - \eqn{\mathbf{S} = (S_1, \dots, S_T)}: hidden (latent) state sequence  
#' - \eqn{\pi_i = \mathbb{P}(S_1 = i)}: initial state distribution  
#' - \eqn{A = [a_{ij}], \text{ where } a_{ij} = \mathbb{P}(S_{t+1} = j \mid S_t = i)}: transition matrix  
#' - \eqn{b_j(x_t) = \mathbb{P}(X_t = x_t \mid S_t = j)}: output probability
#' 
#' Joint probability of the observations and the hidden states:  
#' 
#' \eqn{\mathbb{P}(\mathbf{X}, \mathbf{S}) = \pi_{S_1} b_{S_1}(X_1) \prod_{t=2}^{T} a_{S_{t-1} S_t} b_{S_t}(X_t)}
#'   
#' Marginal probability of the observed data is obtained by summing over all 
#' possible hidden state sequences:  
#' 
#' \eqn{\mathbb{P}(\mathbf{X}) = \sum_{\mathbf{S}} \mathbb{P}(\mathbf{X}, \mathbf{S})}
#' 
#' (Beware of the "double use of data" problem: The bootstrap procedure relies 
#'  on regime classification, but the regimes themselves are estimated from the 
#'  same data and depend on the parameters being resampled.)
#'   
#' @return A list of numeric vectors, each one a bootstrap replicate.
#' 
#' @references
#' Holst, U., Lindgren, G., Holst, J. and Thuvesholmen, M. (1994), Recursive 
#'   Estimation In Switching Autoregressions With A Markov Regime. Journal of 
#'   Time Series Analysis, 15: 489-506. 
#'   [https://doi.org/10.1111/j.1467-9892.1994.tb00206.x](https://doi.org/10.1111/j.1467-9892.1994.tb00206.x)
#' 
#' @examples
#' \donttest{
#' # Requires depmixS4 package
#' if (requireNamespace("depmixS4", quietly = TRUE)) {
#'   
#'   ## Example 1: Univariate time series with regime switching
#'   set.seed(123)
#'   # Generate data with two regimes
#'   n1 <- 100
#'   n2 <- 100
#'   regime1 <- rnorm(n1, mean = 0, sd = 1)
#'   regime2 <- rnorm(n2, mean = 5, sd = 2)
#'   x_univar <- c(regime1, regime2)
#'   
#'   # Generate bootstrap samples
#'   boot_samples <- hmm_bootstrap(
#'     x = x_univar,
#'     n_boot = 150,
#'     num_states = 2,
#'     num_boots = 100
#'   )
#'   
#'   # Inspect results
#'   length(boot_samples)  # 100 bootstrap replicates
#'   dim(boot_samples[[1]])  # Each is a 150 x 1 matrix
#'   
#'   # Compare original and bootstrap distributions
#'   original_mean <- mean(x_univar)
#'   boot_means <- sapply(boot_samples, mean)
#'   
#'   hist(boot_means, main = "Bootstrap Distribution of Mean",
#'        xlab = "Mean", col = "lightblue")
#'   abline(v = original_mean, col = "red", lwd = 2, lty = 2)
#'   
#'   
#'   ## Example 2: Multivariate time series
#'   set.seed(456)
#'   n <- 200
#'   # Two correlated series with regime switching
#'   states_true <- rep(1:2, each = n/2)
#'   x1 <- ifelse(states_true == 1, 
#'                rnorm(n, 0, 1), 
#'                rnorm(n, 4, 2))
#'   x2 <- 0.7 * x1 + rnorm(n, 0, 0.5)
#'   x_multivar <- cbind(x1, x2)
#'   
#'   boot_samples_mv <- hmm_bootstrap(
#'     x = x_multivar,
#'     n_boot = 180,
#'     num_states = 2,
#'     num_boots = 50
#'   )
#'   
#'   dim(boot_samples_mv[[1]])  # 180 x 2 matrix
#'   
#'   # Compute bootstrap correlation estimates
#'   boot_cors <- sapply(boot_samples_mv, function(b) cor(b[,1], b[,2]))
#'   original_cor <- cor(x_multivar[,1], x_multivar[,2])
#'   
#'   hist(boot_cors, main = "Bootstrap Distribution of Correlation",
#'        xlab = "Correlation", col = "lightgreen")
#'   abline(v = original_cor, col = "red", lwd = 2, lty = 2)
#'   
#'   
#'   ## Example 3: Variable-length bootstrap samples
#'   set.seed(789)
#'   x_ar <- arima.sim(n = 150, list(ar = 0.8))
#'   
#'   # Don't specify n_boot to get variable-length samples
#'   boot_samples_var <- hmm_bootstrap(
#'     x = x_ar,
#'     n_boot = NULL,  # Variable length
#'     num_states = 2,
#'     num_blocks = 15,
#'     num_boots = 20
#'   )
#'   
#'   # Check lengths vary
#'   sample_lengths <- sapply(boot_samples_var, nrow)
#'   summary(sample_lengths)
#'   
#'   
#'   ## Example 4: Using verbose mode for diagnostics
#'   set.seed(321)
#'   x_diag <- rnorm(100)
#'   
#'   boot_samples_verbose <- hmm_bootstrap(
#'     x = x_diag,
#'     n_boot = 80,
#'     num_states = 3,
#'     num_boots = 5,
#'     verbose = TRUE  # Print diagnostic information
#'   )
#'   
#'   
#'   ## Example 5: Bootstrap confidence intervals
#'   set.seed(654)
#'   # Data with heteroskedasticity
#'   n <- 200
#'   x_hetero <- numeric(n)
#'   for (i in 1:n) {
#'     sigma <- ifelse(i <= n/2, 1, 3)
#'     x_hetero[i] <- rnorm(1, mean = 2, sd = sigma)
#'   }
#'   
#'   boot_samples_ci <- hmm_bootstrap(
#'     x = x_hetero,
#'     n_boot = 200,
#'     num_states = 2,
#'     num_boots = 500
#'   )
#'   
#'   # Compute bootstrap confidence interval for the mean
#'   boot_means_ci <- sapply(boot_samples_ci, mean)
#'   ci_95 <- quantile(boot_means_ci, c(0.025, 0.975))
#'   
#'   cat("95% Bootstrap CI for mean:", ci_95[1], "to", ci_95[2], "\n")
#'   cat("Original mean:", mean(x_hetero), "\n")
#' }
#' }
#'   
#' @export
hmm_bootstrap <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_blocks = NULL,
    num_boots = 100,
    parallel = FALSE,
    num_cores = 1L,
    verbose = FALSE
) {
  
  ## ---- Dependency Check ----
  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    stop("depmixS4 package required for HMM bootstrap. Install with: install.packages('depmixS4')",
         call. = FALSE)
  }
  
  ## ---- Data Preparation ----
  x <- as.matrix(x)
  n <- nrow(x)
  k <- ncol(x)
  
  ## Basic sanity checks (detailed validation happens in tsbs())
  if (n < num_states * 3) {
    warning("Series length (", n, ") is very short relative to num_states (", num_states, 
            "). HMM fitting may be unreliable.", call. = FALSE)
  }
  
  ## ---- Create Formula for depmixS4 ----
  df <- as.data.frame(x)
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  
  # Build formula based on dimensionality
  if (k == 1) {
    formula <- as.formula("V1 ~ 1")
    family_spec <- list(gaussian())
  } else {
    # For multivariate: each variable gets its own response
    formula_list <- lapply(colnames(df), function(v) as.formula(paste(v, "~ 1")))
    # Replicate gaussian() family for each variable
    family_spec <- replicate(k, gaussian(), simplify = FALSE)
  }
  
  ## ---- Fit HMM ----
  if (verbose) {
    message("Fitting ", num_states, "-state Gaussian HMM to ", k, 
            "-dimensional series of length ", n, "...")
  }
  
  model <- tryCatch({
    if (k == 1) {
      depmixS4::depmix(
        formula, 
        data = df, 
        nstates = num_states, 
        family = family_spec[[1]]
      )
    } else {
      depmixS4::depmix(
        formula_list, 
        data = df, 
        nstates = num_states, 
        family = family_spec
      )
    }
  }, error = function(e) {
    stop("Failed to create HMM model specification. ",
         "This may be due to:\n",
         "  - Incompatible data structure for depmixS4\n",
         "  - Issues with the formula specification\n",
         "  - Data contains invalid values (NA, Inf, NaN)\n",
         "Original error: ", e$message, call. = FALSE)
  })
  
  fit <- tryCatch({
    depmixS4::fit(model, verbose = verbose)
  }, error = function(e) {
    stop("HMM fitting failed. This could be due to:\n",
         "  - Convergence issues (try different starting values or simpler model)\n",
         "  - Data not suitable for ", num_states, "-state HMM (try reducing num_states)\n",
         "  - Insufficient data (need more observations relative to num_states)\n",
         "  - Data quality issues (constant values, extreme outliers, perfect multicollinearity)\n",
         "Original error: ", e$message, call. = FALSE)
  }, warning = function(w) {
    if (verbose) {
      warning("HMM fitting produced a warning: ", w$message, call. = FALSE)
    }
    ## Continue with the fit despite warnings
    suppressWarnings(depmixS4::fit(model, verbose = verbose))
  })
  
  ## Check convergence
  if (!is.null(fit@message) && verbose) {
    message("HMM fitting message: ", fit@message)
  }
  
  ## ---- Extract State Sequence ----
  states <- tryCatch({
    depmixS4::posterior(fit, type = "viterbi")$state
  }, error = function(e) {
    stop("Failed to extract state sequence from fitted HMM. ",
         "Original error: ", e$message, call. = FALSE)
  })
  
  ## Validate state sequence
  if (length(states) != n) {
    stop("State sequence length (", length(states), 
         ") does not match data length (", n, "). ",
         "This is an internal error.", call. = FALSE)
  }
  
  if (length(unique(states)) < num_states) {
    warning("Only ", length(unique(states)), " of ", num_states, 
            " states were observed in the Viterbi sequence. ",
            "Some states may be too rare or the model may not fit well.",
            call. = FALSE)
  }
  
  if (verbose) {
    state_counts <- table(states)
    message("State distribution: ", 
            paste(names(state_counts), "=", state_counts, collapse = ", "))
  }
  
  ## ---- Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating ", num_boots, " bootstrap samples...")
  }
  
  bootstrap_samples <- tryCatch({
    .sample_blocks(
      x = x, 
      n_boot = n_boot, 
      num_blocks = num_blocks, 
      states = states, 
      num_boots = num_boots,
      parallel = parallel,
      num_cores = num_cores
    )
  }, error = function(e) {
    stop("Bootstrap sampling failed. Original error: ", e$message, call. = FALSE)
  })
  
  if (verbose) {
    message("Bootstrap complete. Generated ", length(bootstrap_samples), " samples.")
    if (!is.null(n_boot)) {
      message("Each sample has length ", n_boot)
    } else {
      sample_lengths <- sapply(bootstrap_samples, nrow)
      message("Sample lengths: min=", min(sample_lengths), 
              ", max=", max(sample_lengths), 
              ", mean=", round(mean(sample_lengths), 1))
    }
  }
  
  bootstrap_samples
}



#' Stationary Bootstrap for a 2-State MS-VAR(1) Model
#'
#' This function first fits a 2-state Markov-Switching Vector Autoregressive
#' (MS-VAR) model of order 1 to the provided multivariate time series data.
#' It then uses the estimated state sequence to perform a stationary bootstrap,
#' generating resampled time series that preserve the state-dependent properties
#' of the original data.
#' 
#' For a stationary bootstrap based on a more general \eqn{n}-state MS-VECTOR
#' ARIMA(\eqn{p, d, q})-GARCH model see [ms_varma_garch_bs()].
#'
#' @param x A numeric matrix or data frame where rows are observations and
#'   columns are the time series variables.
#' @param n_boot An integer specifying the length of each bootstrapped series.
#'   If NULL (the default), the length of the original series is used.
#' @param num_blocks An integer specifying the number of blocks to sample for
#'   the bootstrap. Defaults to 100.
#' @param num_boots An integer specifying the total number of bootstrap samples
#'   to generate. Defaults to 100.
#' @param parallel A logical value indicating whether to use parallel processing
#'   for generating bootstrap samples. Defaults to FALSE.
#' @param num_cores An integer specifying the number of cores to use for
#'   parallel processing. Only used if `parallel` is TRUE. Defaults to 1.
#'   
#' @details
#' - \eqn{y_t \in \mathbb{R}^K} be a **$K$-dimensional multivariate response vector** at time \eqn{t}
#' - \eqn{S_t \in {1, \dots, M}} be a **latent Markov chain** with \eqn{M} discrete regimes
#' - \eqn{p} be the **lag order** of the VAR model
#' 
#' \eqn{
#'   y_t = \mu^{(S_t)} + \sum_{i=1}^{p} A_i^{(S_t)} y_{t-i} + \varepsilon_t, \quad \varepsilon_t \sim \mathcal{N}(0, \Sigma^{(S_t)})
#' }
#' 
#' Where:
#' - \eqn{\mu^{(S_t)} \in \mathbb{R}^K} is the regime-specific intercept vector
#' - \eqn{A_i^{(S_t)} \in \mathbb{R}^{K \times K}} are the **regime-specific autoregressive coefficient matrices**
#' - \eqn{\Sigma^{(S_t)} \in \mathbb{R}^{K \times K}} is the regime-specific error covariance matrix
#' 
#' Note: The bootstrap helper function requires the 'foreach' and 'doParallel'
#' packages if `parallel = TRUE`.
#'
#' @return A list of bootstrapped time series matrices.
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' T_obs <- 250
#' y1 <- arima.sim(model = list(ar = 0.7), n = T_obs)
#' y2 <- 0.5 * y1 + arima.sim(model = list(ar = 0.3), n = T_obs)
#' sample_data <- cbind(y1, y2)
#'
#' # Run the bootstrap function (assuming fit_msvar is loaded)
#' # bootstrap_results <- msvar_bootstrap(sample_data, num_boots = 5)
#'
#' # View results
#' # str(bootstrap_results)
#'
#' @export
msvar_bootstrap <- function(
    x,
    n_boot = NULL,
    num_blocks = NULL,
    num_boots = 100,
    parallel = FALSE,
    num_cores = 1L
) {
  
  ## ---- 1. Fit the MS-VAR(1) Model ----
  ## The fit_msvar function handles data validation (e.g., converting to matrix)
  ## and all internal model fitting steps.
  ## It is hardcoded for a 2-state, VAR(1) model, simplifying this wrapper.
  message("Fitting the MS-VAR(1) model...")
  ms_model <- fit_msvar(x)
  
  ## ---- 2. Determine the Most Likely State Sequence ----
  ## The output from fit_msvar contains the smoothed probabilities for each state.
  ## We find the most likely state for each time point (row).
  smoothed_probs <- ms_model$smoothed_probabilities
  state_seq <- apply(smoothed_probs, 1, which.max)
  
  ## The state sequence is one observation shorter than the original series
  ## because the VAR model requires an initial lag. We prepend the first state
  ## to align the sequence with the original data's length.
  state_seq_aligned <- c(state_seq[1], state_seq)
  
  ## ---- 3. Generate Bootstrap Samples ----
  ## This part calls a helper function that performs the actual stationary
  ## block bootstrap based on the state sequence.
  message("Generating bootstrap samples...")
  
  ## Ensure x is a matrix for the helper function
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  bootstrap_samples <- .sample_blocks(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = state_seq_aligned, ## Use 'states' to match the helper function's 
                                ## parameter
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores
  )
  
  return(bootstrap_samples)
}


#' Stationary Bootstrap for a General MS-VARMA-GARCH Model
#'
#' Fits a flexible \eqn{n}-state Markov-Switching Vector ARIMA\eqn{(p, d, q)}-
#' GARCH model and then uses the estimated state sequence to perform a stationary 
#' block bootstrap. This generates resampled time series that preserve the 
#' state-dependent properties of the original data.
#'
#' @param x A numeric matrix where rows are observations and columns are the time
#'   series variables.
#' @param n_boot An integer specifying the length of each bootstrapped series.
#'   If NULL (the default), the length of the original series is used.
#' @param num_blocks An integer specifying the number of blocks to sample for
#'   the bootstrap. Defaults to 100.
#' @param num_boots An integer specifying the total number of bootstrap samples
#'   to generate. Defaults to 100.
#' @param M An integer specifying the number of states in the Markov chain.
#' @param d An integer specifying the order of differencing for the ARIMA model.
#' @param spec A list of model specifications, one for each of the M states.
#' @param model_type A character string, either "univariate" or "multivariate".
#' @param control A list of control parameters for the EM algorithm.
#' @param return_fit If `TRUE`, `tsbs()` will return model fit when 
#' `bs_type = "ms_varma_garch"`. Default is `return_fit = FALSE`. If 
#'  `bs_type = "ms_varma_garch"`, `return_fit = TRUE` and 
#'  `collect_diagnostics = TRUE` diagnostics can be extracted from
#'  `result$fit$diagnostics`. See \code{vignette("Diagnostics", package = "tsbs")}.
#' @param parallel A logical value indicating whether to use parallel processing.
#' @param num_cores An integer specifying the number of cores for parallel processing.
#' @param collect_diagnostics Logical. Collect diagnostics or not.
#' @param verbose Logical. If TRUE, print detailed diagnostic information during 
#'   estimation. Default is FALSE.
#' @param verbose_file Character string specifying path to file for verbose output.
#'   If NULL (default), verbose output goes to console. If specified, all verbose
#'   output is written to this file instead. Only used if verbose = TRUE.
#' 
#' @details
#' The fitted model is defined as:  
#' 
#' Let \eqn{y_t} be the \eqn{k \times 1} vector of observations at time \eqn{t}.
#' The model assumes that the data-generating process is governed by a latent
#' (unobserved) state variable, \eqn{S_t}, which follows a first-order Markov
#' chain with \eqn{M} states.
#'
#' \enumerate{
#'   \item \strong{State Process}: The evolution of the state is described by the
#'   \eqn{M \times M} transition probability matrix \eqn{P}, where the element
#'   \eqn{p_{ij}} is the probability of transitioning from state \eqn{i} to state \eqn{j}:
#'   \deqn{p_{ij} = P(S_t = j | S_{t-1} = i)}
#'   The matrix \eqn{P} is structured such that \eqn{P_{ij} = p_{ij}}, and its
#'   rows sum to one: \eqn{\sum_{j=1}^{M} p_{ij} = 1} for all \eqn{i=1, \dots, M}.
#'
#'   \item \strong{Observation Process}: Conditional on the system being in state
#'   \eqn{S_t = j}, each of the \eqn{k} time series, \eqn{y_{i,t}} for
#'   \eqn{i=1, \dots, k}, is assumed to follow an independent
#'   ARIMA(\eqn{p_j, d_j, q_j})-GARCH(\eqn{q'_j, p'_j}) process. The parameters
#'   for both the mean and variance equations are specific to the state \eqn{j}.
#'   \itemize{
#'     \item \strong{Mean Equation (ARIMA)}:
#'       \deqn{\phi_j(L)(1-L)^{d_j} (y_{i,t} - \mu_j) = \theta_j(L) \varepsilon_{i,t}}
#'       where \eqn{\phi_j(L)} and \eqn{\theta_j(L)} are the AR and MA lag
#'       polynomials, \eqn{\mu_j} is the mean, and \eqn{\varepsilon_{i,t}} is the
#'       innovation term, all specific to state \eqn{j}.
#'     \item \strong{Variance Equation (GARCH)}: The innovations have a conditional
#'       variance \eqn{\sigma_{i,t}^2} that evolves according to:
#'       \deqn{\varepsilon_{i,t} = \sigma_{i,t} z_{i,t}, \quad z_{i,t} \sim \mathcal{N}(0,1)}
#'       \deqn{\sigma_{i,t}^2 = \omega_j + \sum_{l=1}^{q'_j} \alpha_{j,l} \varepsilon_{i,t-l}^2 + \sum_{l=1}^{p'_j} \beta_{j,l} \sigma_{i,t-l}^2}
#'       where \eqn{\omega_j, \alpha_{j,l}, \beta_{j,l}} are the GARCH parameters for state \eqn{j}.
#'   }
#' }
#' Let \eqn{\Psi_j = \{\mu_j, \phi_j, \theta_j, \omega_j, \alpha_j, \beta_j\}} be
#' the complete set of ARIMA-GARCH parameters for state \eqn{j}, and let
#' \eqn{\Psi = \{\Psi_1, \dots, \Psi_M, P\}} be the full parameter set for the
#' entire model. The EM algorithm using a Hamilton Filter & Kim Smoother for the 
#' E-step is used to find the Maximum Likelihood Estimate (MLE) of \eqn{\Psi}.
#' 
#' The tsbs package uses different optimization strategies for DCC models
#' depending on the model order:
#' 
#' \strong{DCC(1,1) - Reparameterized Optimization}
#' 
#' For the common DCC(1,1) case, we use a reparameterization that transforms

#' the constrained problem into an unconstrained one:
#' \itemize{
#'   \item Original: \eqn{\alpha \in (0,1), \beta \in (0,1), \alpha + \beta < 1}
#'   \item Reparameterized: \eqn{persistence \in (0,1), ratio \in (0,1)}
#' }
#' 
#' where \eqn{persistence = \alpha + \beta} and \eqn{ratio = \alpha/(\alpha + \beta)}.
#' 
#' This eliminates the stationarity constraint since \eqn{\alpha + \beta = persistence < 1}
#' is automatically satisfied by the box constraint on persistence.
#' 
#' Benefits:
#' \itemize{
#'   \item No penalty function discontinuities
#'   \item More stable optimization near high-persistence regions
#'   \item Eliminates "dcc_penalty" warnings from optimizer exploration
#' }
#' 
#' \strong{DCC(p,q) with max(p,q) > 1 - Penalty Method}
#' 
#' For higher-order DCC models, reparameterization becomes significantly more
#' complex (requiring softmax distributions over parameter vectors). We therefore
#' use the standard penalty method:
#' \itemize{
#'   \item Box constraints: \eqn{\alpha_j, \beta_j \in (\epsilon, 1-\epsilon)}
#'   \item Stationarity enforced via penalty when \eqn{\sum \alpha + \sum \beta \geq 1}
#' }
#' 
#' This may result in optimizer instability warnings for models with high
#' persistence.
#' 
#' @return A list of bootstrapped time series matrices.
#' 
#' @references
#' Natatou Moutari, D. et al. (2021). Dependence Modeling and Risk Assessment 
#'   of a Financial Portfolio with ARMA-APARCH-EVT models based on HACs. 
#'   [arXiv:2105.09473](http://arxiv.org/abs/2105.09473)
#'
#' @export
ms_varma_garch_bs <- function(
    x,
    n_boot = NULL,
    num_blocks = 100,
    num_boots = 100,
    M,
    d = 0,
    spec,
    model_type = c("univariate", "multivariate"),
    control = list(),
    parallel = FALSE,
    num_cores = 1L,
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    verbose_file = NULL
) {
  
  ## ---- 1. Input Validation ----
  model_type <- match.arg(model_type)
  if (is.null(num_boots) || !is.numeric(num_boots) || num_boots < 1) {
    stop("num_boots must be a positive integer.", call. = FALSE)
  }
  if (is.null(n_boot) && is.null(num_blocks)) {
    stop("Must provide a valid value for either n_boot or num_blocks", call. = FALSE)
  }
  
  ## ---- 2. Fit the MS-VARMA-GARCH Model ----
  ## This function handles all data validation and model fitting.
  ms_model <- fit_ms_varma_garch(
    y = x,
    M = M,
    d = d,
    spec = spec,
    model_type = model_type,
    control = control,
    # parallel = parallel,
    # num_cores = num_cores,
    collect_diagnostics = collect_diagnostics,
    verbose = verbose,
    verbose_file = verbose_file
  )
  
  ## ---- 3. Determine the Most Likely State Sequence ----
  ## The smoothed_probabilities from the fit object are already aligned
  ## with the original data 'x', but may have leading NAs if d > 0.
  smoothed_probs <- ms_model$smoothed_probabilities
  state_seq <- apply(smoothed_probs, 1, which.max)
  
  ## Handle potential leading NAs from differencing by forward-filling the first valid state.
  if (anyNA(state_seq)) {
    first_valid_state_idx <- min(which(!is.na(state_seq)))
    first_valid_state <- state_seq[first_valid_state_idx]
    state_seq[1:(first_valid_state_idx - 1)] <- first_valid_state
  }
  
  ## ---- 4. Generate Bootstrap Samples ----
  ## This part calls the existing helper function that performs the actual
  ## stationary block bootstrap based on the state sequence.
  message("Generating bootstrap samples...")
  
  ## Ensure x is a matrix for the helper function
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  bootstrap_samples <- .sample_blocks(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = state_seq,
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores
  )
  
  if (return_fit) {
    list(
      bootstrap_series = bootstrap_samples,
      fit = ms_model
    )
  } else {
    bootstrap_samples
  }
}



#' Wild Bootstrap for Time Series Residuals
#'
#' Generates wild bootstrap replicates of a vector or matrix of residuals
#' by multiplying each observation by a random Rademacher weight (+1 or -1).
#'
#' @param x Numeric vector or matrix of residuals.
#' @param num_boots Integer number of bootstrap replicates.
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores.
#'
#' @return A list of numeric matrices, each one a wild bootstrap replicate.
#'
#' @details
#' The wild bootstrap is often used to resample regression or model residuals
#' when heteroskedasticity or other non-i.i.d. errors are present.
#' Each replicate is constructed by multiplying every observation by +1 or -1,
#' where the signs are drawn randomly with equal probability.
#' 
#' @references
#' A. Colin Cameron & Jonah B. Gelbach & Douglas L. Miller, 2008. 
#'   "Bootstrap-Based Improvements for Inference with Clustered Errors", 
#'   The Review of Economics and Statistics, MIT Press, vol. 90(3), pages 
#'   414-427, August.
#'
#' @examples
#' set.seed(123)
#' resids <- rnorm(100)
#' boot_reps <- wild_bootstrap(resids, num_boots = 5)
#' length(boot_reps)           # 5 replicates
#' dim(boot_reps[[1]])         # 100 x 1 matrix
#'
#' @export
wild_bootstrap <- function(
    x, 
    num_boots = 100,
    parallel = FALSE,
    num_cores = 2
    #wild_type = "Rademacher" ## No other options currently implemented
  ) {
  if (!is.numeric(x)) stop("`x` must be numeric.")
  if (!is.numeric(num_boots) || num_boots < 1) stop("`num_boots` must be a positive integer.")
  
  ## Coerce vector to a column matrix
  if (is.vector(x)) {
    x <- matrix(x, ncol = 1)
  }
  
  n <- nrow(x)
  

  ## ---- Parallel Backend Setup ----
  
  ## The `%dopar%` operator from foreach is special and needs to be imported
  ## or defined. Define it locally if the package is found.
  `%dopar%` <- if (parallel && num_cores > 1) {
      foreach::`%dopar%`
    } else {
      foreach::`%do%`
    }
  
  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) || 
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for parallel execution.", 
      call. = FALSE)
    }
    if (is.null(num_cores) || num_cores < 1) {
      stop(
        paste0("To run in parallel, you must specify 'num_cores'.
       The number of cores on your machine is ", 
               as.character(parallel::detectCores()), "."),
        call. = FALSE
      )
    }
    if (num_cores > 1) {
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      `%dopar%` <- foreach::`%dopar%`
    } else {
      ## Prevent a warning message from being issued if the ⁠%dopar%⁠ function is
      ## alled and no parallel backend has been registered.
      foreach::registerDoSEQ()
    }
  } else {
    foreach::registerDoSEQ()
  }
  
  
  ## ---- Bootstrap ----
  
  # bootstrap_samples <- vector("list", num_boots)
  # 
  # for (b in seq_len(num_boots)) {
  #   # Rademacher weights (+1 or -1)
  #   v <- sample(c(-1, 1), size = n, replace = TRUE)
  #   # Multiply each column by v
  #   bootstrap_samples[[b]] <- x * v
  # }
  
  bootstrap_samples <- foreach::foreach(b = seq_len(num_boots)) %dopar% {
    # Rademacher weights (+1 or -1)
    v <- sample(c(-1, 1), size = n, replace = TRUE)
    # Multiply each observation by the random weight v
    x * v
  }
  
  bootstrap_samples
}



#' Helper Function for Multivariate Stationary Block Bootstrap
#' @param x The original multivariate time series as a matrix.
#' @param n_boot The target length for each bootstrap sample.
#' @param num_blocks The number of blocks to sample.
#' @param states A vector of integer states corresponding to each row of x.
#' @param num_boots The number of bootstrap series to generate.
#' @param parallel Logical, whether to use parallel processing.
#' @param num_cores The number of cores for parallel processing.
#' @return A list of bootstrapped series.
#' @keywords internal
.sample_blocks <- function(
    x,
    n_boot,
    num_blocks,
    states,
    num_boots,
    parallel = FALSE,
    num_cores = 1L
) {
  
  ## ---- Identify blocks ----
  r <- rle(states)
  ends <- cumsum(r$lengths)
  ## The state blocks are identified from the original series length,
  ## so we need to account for the rows lost to lagging.
  n_orig <- nrow(x)
  n_states <- length(states)
  offset <- n_orig - n_states
  
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(
    function(s, e) list(start = s + offset, end = e + offset),
    starts, ends
  )
  
  get_block <- function(b) {
    x[b$start:b$end, , drop = FALSE]
  }
  
  ## ---- Parallel Backend Setup ----
  `%dopar%` <- foreach::`%do%`
  if (parallel) {
    if (!requireNamespace("foreach", quietly = TRUE) || !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' are required for parallel execution.", call. = FALSE)
    }
    if (is.null(num_cores) || num_cores < 1) {
      stop("To run in parallel, you must specify a positive 'num_cores'.", call. = FALSE)
    }
    if (num_cores > 1) {
      cl <- parallel::makeCluster(num_cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      `%dopar%` <- foreach::`%dopar%`
    }
  }
  if (is.null(n_boot) && is.null(num_blocks)) {
    stop("Must provide a valid value for either n_boot or num_blocks")
  }
  if (!is.null(n_boot) && !is.null(num_blocks)) {
    warning("`num_blocks` is ignored when `n_boot` is set.")
  }
  
  sampled_series <- foreach::foreach(i = seq_len(num_boots)) %dopar% {
    if (!is.null(n_boot)) {
      ## Initialize an empty matrix and use rbind
      bootstrap_series <- x[0, , drop = FALSE]
      while (nrow(bootstrap_series) < n_boot) {
        sampled_block <- sample(blocks, 1, replace = TRUE)[[1]]
        bootstrap_series <- rbind(bootstrap_series, get_block(sampled_block))
      }
      bootstrap_series[seq_len(n_boot), , drop = FALSE]
    } else {
      ## Sample blocks and combine with do.call(rbind, ...)
      bootstrap_blocks <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
      do.call(rbind, bootstrap_blocks)
    }
  }
  
  sampled_series
}