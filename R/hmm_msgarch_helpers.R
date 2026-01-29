#' @title MSGARCH Helper Functions for HMM Bootstrap
#' @description Internal functions for regime-switching bootstrap using MSGARCH.
#'
#' This module provides helper functions that wrap the MSGARCH package to enable
#' semi-parametric bootstrap for time series with regime-switching behavior and
#' non-Gaussian conditional distributions (skew Student-t, Student-t, GED, etc.).
#'
#' @details
#' ## Literature Background
#'
#' The implementation combines several established techniques:
#'
#' - **Markov-switching GARCH**: Haas, Mittnik & Paolella (2004), implemented
#'   via the MSGARCH package (Ardia et al., 2019)
#' - **Skew Student-t distribution**: Fernández & Steel (1998) transformation
#' - **Semi-parametric bootstrap**: State sequences are simulated from the fitted
#'   Markov chain (parametric), while innovations are resampled from empirical
#'   state-specific pools (nonparametric)
#'
#' @references
#' Ardia, D., Bluteau, K., Boudt, K., Catania, L., & Trottier, D.-A. (2019).
#' Markov-Switching GARCH Models in R: The MSGARCH Package.
#' Journal of Statistical Software, 91(4), 1-38. doi:10.18637/jss.v091.i04
#'
#' Fernández, C., & Steel, M. F. (1998). On Bayesian modeling of fat tails
#' and skewness. Journal of the American Statistical Association, 93(441), 359-371.
#'
#' Haas, M., Mittnik, S., & Paolella, M. S. (2004). A New Approach to Markov-
#' Switching GARCH Models. Journal of Financial Econometrics, 2, 493-530.
#'
#' Hamilton, J. D. (1989). A New Approach to the Economic Analysis of
#' Nonstationary Time Series and the Business Cycle. Econometrica, 57(2), 357-384.
#'
#' @name hmm_msgarch_helpers
#' @keywords internal
NULL


#' Fit MSGARCH Model
#'
#' Thin wrapper around MSGARCH::FitML with sensible defaults for bootstrap use.
#' Provides consistent interface and error handling for the tsbs package.
#'
#' @param y Numeric vector of returns (univariate).
#' @param n_states Integer, number of regimes (default 2).
#' @param variance_model Character, GARCH specification: one of
#'   \code{"sGARCH"}, \code{"eGARCH"}, \code{"gjrGARCH"}, \code{"tGARCH"}.
#'   Default is \code{"sGARCH"}.
#' @param distribution Character, conditional distribution: one of
#'   \code{"norm"}, \code{"std"}, \code{"ged"}, \code{"snorm"}, \code{"sstd"},
#'   \code{"sged"}. Default is \code{"sstd"} (skew Student-t).
#' @param ... Additional arguments passed to \code{MSGARCH::FitML}.
#'
#' @return An MSGARCH fit object (class \code{MSGARCH_ML_FIT}).
#'
#' @details
#' This function creates a Markov-switching GARCH specification using the
#' Haas et al. (2004) approach, which avoids the path-dependency problem by
#' defining K separate GARCH processes that evolve in parallel.
#'
#' The skewed distributions use the Fernández & Steel (1998) transformation,
#' which introduces a skewness parameter \eqn{\xi > 0} where \eqn{\xi = 1}
#' corresponds to a symmetric distribution.
#'
#' @seealso \code{\link[MSGARCH]{CreateSpec}}, \code{\link[MSGARCH]{FitML}}
#'
#' @keywords internal
fit_msgarch_model <- function(
    y,
    n_states = 2L,
    variance_model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH"),
    distribution = c("sstd", "std", "norm", "snorm", "sged", "ged"),
    ...
) {
  
  ## ---- Input Validation ----
  
  # Check MSGARCH availability
  if (!requireNamespace("MSGARCH", quietly = TRUE)) {
    stop("Package 'MSGARCH' is required for non-Gaussian HMM bootstrap. ",
         "Install with: install.packages('MSGARCH')", call. = FALSE)
  }
  
  # Validate y
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector.", call. = FALSE)
  }
  if (length(y) < 50) {
    warning("Series length < 50 may lead to unreliable MSGARCH estimation.",
            call. = FALSE)
  }
  if (any(!is.finite(y))) {
    stop("'y' contains non-finite values (NA, NaN, Inf).", call. = FALSE)
  }
  
  # Validate n_states
  n_states <- as.integer(n_states)
  if (n_states < 2L) {
    stop("'n_states' must be at least 2.", call. = FALSE)
  }
  
  if (n_states > 5L) {
    warning("More than 5 states may lead to estimation difficulties.",
            call. = FALSE)
  }
  
  # Match arguments
  variance_model <- match.arg(variance_model)
  distribution <- match.arg(distribution)
  
  ## ---- Create MSGARCH Specification ----
  
  # Map variance model names to MSGARCH format
  msgarch_model <- switch(variance_model,
                          "sGARCH" = "sGARCH",
                          "eGARCH" = "eGARCH",
                          "gjrGARCH" = "gjrGARCH",
                          "tGARCH" = "tGARCH",
                          stop("Unknown variance model: ", variance_model, call. = FALSE)
  )
  
  # Create specification
  spec <- MSGARCH::CreateSpec(
    variance.spec = list(model = rep(msgarch_model, n_states)),
    distribution.spec = list(distribution = rep(distribution, n_states)),
    switch.spec = list(do.mix = FALSE)
  )
  
  ## ---- Fit Model ----
  
  fit <- tryCatch({
    MSGARCH::FitML(spec = spec, data = y, ...)
  }, error = function(e) {
    stop("MSGARCH model fitting failed: ", e$message,
         "\nConsider trying different starting values or fewer states.",
         call. = FALSE)
  })
  
  # Check convergence
  if (!is.null(fit$convergence) && fit$convergence != 0) {
    warning("MSGARCH optimization may not have converged (code: ",
            fit$convergence, ").", call. = FALSE)
  }
  
  return(fit)
}


#' Extract State Information from MSGARCH Fit
#'
#' Extracts decoded state sequence, transition matrix, and state probabilities
#' from a fitted MSGARCH model.
#'
#' @param fit An MSGARCH fit object (output from \code{fit_msgarch_model} or
#'   \code{MSGARCH::FitML}).
#' @param y Numeric vector, the original data used for fitting. Required for
#'   computing residuals.
#'
#' @return A list containing:
#'   \describe{
#'     \item{states}{Integer vector of Viterbi-decoded state assignments.}
#'     \item{transition_matrix}{K x K matrix of transition probabilities.}
#'     \item{smoothed_probs}{T x K matrix of smoothed state probabilities.}
#'     \item{n_states}{Integer, number of states K.}
#'     \item{coefficients}{Named vector of estimated model parameters.}
#'     \item{state_durations}{Numeric vector of expected state durations
#'       (1 / (1 - p_ii) for each state i).}
#'   }
#'
#' @keywords internal
extract_msgarch_states <- function(fit, y) {
  
  ## ---- Input Validation ----
  
  if (!inherits(fit, "MSGARCH_ML_FIT")) {
    stop("'fit' must be an MSGARCH_ML_FIT object.", call. = FALSE)
  }
  
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector.", call. = FALSE)
  }
  
  ## ---- Extract State Information ----
  
  # Get state probabilities and Viterbi path
  state_probs <- MSGARCH::State(fit)
  
  # Viterbi decoding gives most likely state sequence
  # Viterbi is returned as a matrix/array, convert to vector
  viterbi <- as.vector(state_probs$Viterbi)
  T_len <- length(viterbi)
  
  # Smoothed probabilities
  # State() returns SmoothProb as a 3D array with dimensions [T+1, n_draws, K]
  # where:
  #   - First dimension is time (T+1 because it includes initial state prob)
  #   - Second dimension is number of draws (1 for ML estimation)
  #   - Third dimension is the K states
  # We need a [T, K] matrix, so we:
  #   1. Drop the first time point (initial state prob before any observations)
  #   2. Take the first draw (index 1 in second dimension)
  #   3. Keep all states (third dimension)
  smoothed_raw <- state_probs$SmoothProb
  n_states <- dim(smoothed_raw)[3]
  
  # Extract [T, K] matrix: drop first row, take first draw, all states
  smoothed <- matrix(smoothed_raw[-1, 1, ], nrow = T_len, ncol = n_states)
  
  ## ---- Get Transition Matrix ----
  
  trans_mat <- MSGARCH::TransMat(fit)
  
  # Ensure proper matrix format
  if (!is.matrix(trans_mat)) {
    trans_mat <- as.matrix(trans_mat)
  }
  
  ## ---- Compute State Durations ----
  
  n_states <- nrow(trans_mat)
  state_durations <- numeric(n_states)
  for (k in 1:n_states) {
    p_stay <- trans_mat[k, k]
    # Expected duration = 1 / (1 - p_stay) if p_stay < 1
    if (p_stay < 1) {
      state_durations[k] <- 1 / (1 - p_stay)
    } else {
      state_durations[k] <- Inf
    }
  }
  
  ## ---- Package Results ----
  
  list(
    states = viterbi,
    transition_matrix = trans_mat,
    smoothed_probs = smoothed,
    n_states = n_states,
    coefficients = fit$par,
    state_durations = state_durations
  )
}


#' Extract State-Specific Innovation Pools
#'
#' Builds pools of time indices for each regime based on Viterbi decoding.
#' These pools are used for synchronized sampling in the bootstrap procedure.
#'
#' @param fit An MSGARCH fit object.
#' @param y Numeric vector, the original data.
#' @param state_info Output from \code{extract_msgarch_states}.
#'
#' @return A list of length \code{n_states}, where each element contains:
#'   \describe{
#'     \item{indices}{Integer vector of time indices assigned to this state.}
#'     \item{n}{Integer, count of observations in this state.}
#'     \item{proportion}{Numeric, proportion of total observations.}
#'   }
#'
#' @details
#' For multivariate bootstrap with synchronized sampling, we sample entire
#' time indices (not individual residuals) from each state pool. This preserves
#' the cross-sectional dependence structure that existed at each time point.
#'
#' @keywords internal
extract_state_pools <- function(fit, y, state_info) {
  
  ## ---- Input Validation ----
  
  if (!is.list(state_info) || is.null(state_info$states)) {
    stop("'state_info' must be output from extract_msgarch_states().",
         call. = FALSE)
  }
  
  n_states <- state_info$n_states
  states <- state_info$states
  T_len <- length(y)
  
  # Check alignment
  if (length(states) != T_len) {
    stop("Length mismatch: states (", length(states),
         ") vs y (", T_len, ").", call. = FALSE)
  }
  
  ## ---- Build Pools ----
  
  pools <- vector("list", n_states)
  
  for (k in 1:n_states) {
    idx <- which(states == k)
    pools[[k]] <- list(
      indices = idx,
      n = length(idx),
      proportion = length(idx) / T_len
    )
  }
  
  # Check for empty pools (model found fewer effective states than requested)
  pool_sizes <- sapply(pools, function(p) p$n)
  empty_states <- which(pool_sizes == 0)
  
  if (length(empty_states) > 0) {
    stop("Viterbi decoding assigned no observations to state(s): ",
         paste(empty_states, collapse = ", "), ". ",
         "The model effectively has fewer states than specified. ",
         "Consider reducing 'num_states' or using different data.",
         call. = FALSE)
  }
  
  # Warn if any pool is very small
  min_pool_size <- min(pool_sizes)
  if (min_pool_size < 10) {
    warning("State ", which.min(pool_sizes),
            " has only ", min_pool_size, " observations. ",
            "Bootstrap samples may have limited diversity.", call. = FALSE)
  }
  
  return(pools)
}


#' Compute Stationary Distribution of Markov Chain
#'
#' Computes the stationary (ergodic) distribution of a finite-state Markov chain
#' by solving \eqn{\pi P = \pi} subject to \eqn{\sum_k \pi_k = 1}.
#'
#' @param trans_mat Square matrix of transition probabilities. Row i, column j
#'   gives P(S_{t+1} = j | S_t = i). Rows must sum to 1.
#'
#' @return Numeric vector of length K giving the stationary probabilities.
#'
#' @details
#' Solves the system \eqn{(P' - I) \pi = 0} with the constraint \eqn{\sum \pi = 1}
#' by replacing one equation with the normalization constraint.
#'
#' @keywords internal
compute_stationary_dist <- function(trans_mat) {
  
  ## ---- Input Validation ----
  
  if (!is.matrix(trans_mat)) {
    stop("'trans_mat' must be a matrix.", call. = FALSE)
  }
  
  K <- nrow(trans_mat)
  if (K != ncol(trans_mat)) {
    stop("'trans_mat' must be square.", call. = FALSE)
  }
  
  # Check rows sum to 1 (within tolerance)
  row_sums <- rowSums(trans_mat)
  if (any(abs(row_sums - 1) > 1e-6)) {
    warning("Transition matrix rows do not sum to 1. Normalizing.",
            call. = FALSE)
    trans_mat <- trans_mat / row_sums
  }
  
  ## ---- Solve for Stationary Distribution ----
  
  # System: (P' - I) * pi = 0, with sum(pi) = 1
  # Replace last row with normalization constraint
  A <- t(trans_mat) - diag(K)
  A[K, ] <- 1
  b <- c(rep(0, K - 1), 1)
  
  pi_stat <- tryCatch({
    solve(A, b)
  }, error = function(e) {
    # Fallback: use eigenvalue decomposition
    # Stationary distribution is the left eigenvector with eigenvalue 1
    ev <- eigen(t(trans_mat))
    idx <- which.min(abs(ev$values - 1))
    pi_raw <- Re(ev$vectors[, idx])
    pi_raw / sum(pi_raw)
  })
  
  # Ensure non-negative (numerical issues can cause tiny negatives)
  pi_stat <- pmax(pi_stat, 0)
  pi_stat <- pi_stat / sum(pi_stat)
  
  return(pi_stat)
}


#' Simulate Markov Chain
#'
#' Generates a realization of a discrete-state Markov chain given transition
#' probabilities and initial distribution.
#'
#' @param T_len Integer, length of sequence to generate.
#' @param trans_mat Square matrix of transition probabilities.
#' @param init_dist Numeric vector of initial state probabilities. If NULL,
#'   uses the stationary distribution.
#'
#' @return Integer vector of length T_len containing state assignments (1 to K).
#'
#' @keywords internal
simulate_markov_chain <- function(T_len, trans_mat, init_dist = NULL) {
  
  ## ---- Input Validation ----
  
  T_len <- as.integer(T_len)
  if (T_len < 1L) {
    stop("'T_len' must be at least 1.", call. = FALSE)
  }
  
  K <- nrow(trans_mat)
  
  if (is.null(init_dist)) {
    init_dist <- compute_stationary_dist(trans_mat)
  }
  
  if (length(init_dist) != K) {
    stop("'init_dist' length must match number of states.", call. = FALSE)
  }
  
  ## ---- Simulate Chain ----
  
  states <- integer(T_len)
  
  # Initial state
  states[1L] <- sample.int(K, size = 1L, prob = init_dist)
  
  # Subsequent states
  if (T_len > 1L) {
    for (t in 2L:T_len) {
      states[t] <- sample.int(K, size = 1L, prob = trans_mat[states[t - 1L], ])
    }
  }
  
  return(states)
}


#' Sample from State-Specific Pools with Synchronized Sampling
#'
#' For each position in a simulated state sequence, samples a time index
#' from the original data where that state occurred. Supports both individual
#' (iid) sampling and micro-block sampling within states.
#'
#' @param sim_states Integer vector of simulated state sequence (length T).
#' @param pools State-specific pools from \code{extract_state_pools}.
#' @param micro_block_length Integer, length of micro-blocks for within-state
#'   sampling. Use 1 for iid sampling (default), >1 to preserve local dependence.
#'
#' @return Integer vector of length T containing sampled time indices from
#'   the original data.
#'
#' @details
#' When \code{micro_block_length = 1}, each time point is sampled independently
#' from the appropriate state pool.
#'
#' When \code{micro_block_length > 1}, the function identifies "runs" of
#' consecutive time points in the same state and fills them with micro-blocks
#' sampled from that state's pool. This preserves some local autocorrelation
#' in innovations while still allowing the Markov chain to drive regime dynamics.
#'
#' For multivariate data, the same time indices are used across all variables,
#' which preserves cross-sectional dependence.
#'
#' @keywords internal
sample_from_pools_sync <- function(sim_states, pools, micro_block_length = 1L) {
  
  ## ---- Input Validation ----
  
  T_len <- length(sim_states)
  if (T_len < 1L) {
    stop("'sim_states' must have length >= 1.", call. = FALSE)
  }
  
  micro_block_length <- as.integer(micro_block_length)
  if (micro_block_length < 1L) {
    stop("'micro_block_length' must be at least 1.", call. = FALSE)
  }
  
  n_states <- length(pools)
  if (any(sim_states < 1L) || any(sim_states > n_states)) {
    stop("'sim_states' contains invalid state values.", call. = FALSE)
  }
  
  # Check pools have indices
  for (k in 1:n_states) {
    if (length(pools[[k]]$indices) == 0) {
      stop("Pool for state ", k, " is empty.", call. = FALSE)
    }
  }
  
  ## ---- Sample Indices ----
  
  sampled_idx <- integer(T_len)
  
  if (micro_block_length == 1L) {
    ## Simple iid sampling within each state
    for (t in 1:T_len) {
      k <- sim_states[t]
      pool_idx <- pools[[k]]$indices
      sampled_idx[t] <- pool_idx[sample.int(length(pool_idx), size = 1L)]
    }
    
  } else {
    ## Micro-block sampling
    t <- 1L
    while (t <= T_len) {
      k <- sim_states[t]
      pool_idx <- pools[[k]]$indices
      n_pool <- length(pool_idx)
      
      # Find length of current run in same state
      run_end <- t
      while (run_end < T_len && sim_states[run_end + 1L] == k) {
        run_end <- run_end + 1L
      }
      run_length <- run_end - t + 1L
      
      # Fill run with micro-blocks
      filled <- 0L
      while (filled < run_length) {
        remaining <- run_length - filled
        block_len <- min(micro_block_length, remaining)
        
        # Find valid starting points for a block of this length
        # Need consecutive indices in the pool
        if (block_len > 1L && n_pool >= block_len) {
          # Find consecutive sequences in pool_idx
          valid_starts <- .find_consecutive_starts(pool_idx, block_len)
          if (length(valid_starts) == 0L) {
            # Fall back to allowing any start with wrapping
            valid_starts <- pool_idx
            block_len <- 1L
          }
        } else {
          valid_starts <- pool_idx
          block_len <- 1L
        }
        
        # Sample starting point
        start_idx <- valid_starts[sample.int(length(valid_starts), size = 1L)]
        
        if (block_len > 1L) {
          # Get consecutive indices
          block_indices <- start_idx:(start_idx + block_len - 1L)
        } else {
          block_indices <- start_idx
        }
        
        # Assign to output
        assign_range <- (t + filled):(t + filled + length(block_indices) - 1L)
        sampled_idx[assign_range] <- block_indices
        filled <- filled + length(block_indices)
      }
      
      t <- run_end + 1L
    }
  }
  
  return(sampled_idx)
}


#' Find Starting Points for Consecutive Sequences
#'
#' Helper function to find indices in a pool that can serve as starting points
#' for blocks of a specified length with consecutive values.
#'
#' @param pool_idx Integer vector of available indices.
#' @param block_len Integer, required block length.
#'
#' @return Integer vector of valid starting points.
#'
#' @keywords internal
.find_consecutive_starts <- function(pool_idx, block_len) {
  
  if (length(pool_idx) < block_len) {
    return(integer(0))
  }
  
  # Sort pool indices
  sorted_idx <- sort(pool_idx)
  n <- length(sorted_idx)
  
  # Find runs of consecutive integers
  diffs <- diff(sorted_idx)
  run_starts <- c(1L, which(diffs != 1L) + 1L)
  run_ends <- c(which(diffs != 1L), n)
  run_lengths <- run_ends - run_starts + 1L
  
  # Collect valid starting points
  valid_starts <- integer(0)
  for (i in seq_along(run_starts)) {
    if (run_lengths[i] >= block_len) {
      # Can start at any position in this run that leaves room for block_len
      n_valid <- run_lengths[i] - block_len + 1L
      start_positions <- sorted_idx[run_starts[i]:(run_starts[i] + n_valid - 1L)]
      valid_starts <- c(valid_starts, start_positions)
    }
  }
  
  return(valid_starts)
}


#' Get Series for Regime Identification (Multivariate Case)
#'
#' For multivariate data, determines which series (or aggregate) to use
#' for fitting the regime model.
#'
#' @param y_mat Numeric matrix (T x N) of returns.
#' @param regime_basis Character or integer specifying how to identify regimes:
#'   \describe{
#'     \item{"market"}{Use equal-weighted average (rowMeans).}
#'     \item{"first_pc"}{Use first principal component.}
#'     \item{integer}{Use the specified column index.}
#'   }
#'
#' @return Numeric vector of length T for regime identification.
#'
#' @keywords internal
get_regime_series <- function(y_mat, regime_basis = "market") {
  
  ## ---- Input Validation ----
  
  if (!is.matrix(y_mat)) {
    y_mat <- as.matrix(y_mat)
  }
  
  T_len <- nrow(y_mat)
  N <- ncol(y_mat)
  
  ## ---- Extract Series ----
  
  if (is.numeric(regime_basis) && length(regime_basis) == 1L) {
    # Use specified column
    col_idx <- as.integer(regime_basis)
    if (col_idx < 1L || col_idx > N) {
      stop("'regime_basis' column index ", col_idx,
           " is out of range [1, ", N, "].", call. = FALSE)
    }
    return(y_mat[, col_idx])
    
  } else if (is.character(regime_basis)) {
    regime_basis <- tolower(regime_basis)
    
    if (regime_basis == "market") {
      # Equal-weighted aggregate
      return(rowMeans(y_mat))
      
    } else if (regime_basis %in% c("first_pc", "firstpc", "pc1")) {
      # First principal component
      if (N < 2L) {
        warning("Cannot compute PC with single column. Using column directly.",
                call. = FALSE)
        return(y_mat[, 1L])
      }
      pc <- prcomp(y_mat, center = TRUE, scale. = TRUE)
      return(pc$x[, 1L])
      
    } else {
      stop("Unknown regime_basis: '", regime_basis, "'. ",
           "Use 'market', 'first_pc', or a column index.", call. = FALSE)
    }
    
  } else {
    stop("'regime_basis' must be 'market', 'first_pc', or an integer column index.",
         call. = FALSE)
  }
}


#' Generate Bootstrap Samples Using MSGARCH
#'
#' Core bootstrap generation function. Simulates new state sequences from the
#' fitted Markov chain and samples observations from state-specific empirical
#' pools (semi-parametric bootstrap).
#'
#' @param fit MSGARCH fit object from \code{fit_msgarch_model}.
#' @param y Original data. Can be a vector (univariate) or matrix (multivariate).
#'   For multivariate data with synchronized sampling, the regime model is fit
#'   to an aggregate, and full cross-sections are sampled together.
#' @param state_info Output from \code{extract_msgarch_states}.
#' @param pools Output from \code{extract_state_pools}.
#' @param n_boot Integer, number of bootstrap replicates (default 1000).
#' @param micro_block_length Integer, block length for within-state sampling.
#'   Use 1 for iid sampling (default), >1 to preserve local dependence.
#' @param sync_sampling Logical. If TRUE and y is multivariate, sample the same
#'   time indices across all assets to preserve cross-sectional dependence.
#'   Default is TRUE.
#' @param seed Integer, random seed for reproducibility. Default is NULL.
#'
#' @return A list containing:
#'   \describe{
#'     \item{samples}{3D array of bootstrap samples with dimensions
#'       (T x N x n_boot) where T is series length and N is number of variables.}
#'     \item{states}{Matrix of simulated state sequences (T x n_boot).}
#'     \item{sampled_indices}{Matrix of sampled time indices (T x n_boot).}
#'     \item{state_info}{The input state_info for reference.}
#'     \item{pools}{The input pools for reference.}
#'     \item{params}{List of bootstrap parameters used.}
#'   }
#'
#' @details
#' The bootstrap procedure is semi-parametric:
#' \itemize{
#'   \item \strong{Parametric component}: State sequences are simulated from
#'     the fitted Markov chain with transition matrix P.
#'   \item \strong{Nonparametric component}: Observations are resampled from
#'     empirical state-specific pools, preserving the actual distributional
#'     characteristics without parametric assumptions.
#' }
#'
#' For multivariate data, the \code{sync_sampling = TRUE} option ensures that
#' the same time index is sampled for all variables, preserving whatever
#' cross-sectional dependence existed at that moment in the original data.
#' This is simpler than copula-based approaches and automatically captures
#' the empirical dependence structure.
#'
#' @seealso \code{\link{fit_msgarch_model}}, \code{\link{extract_msgarch_states}},
#'   \code{\link{extract_state_pools}}
#'
#' @keywords internal
generate_msgarch_bootstrap <- function(
    fit,
    y,
    state_info,
    pools,
    n_boot = 1000L,
    micro_block_length = 1L,
    sync_sampling = TRUE,
    seed = NULL
) {
  
  ## ---- Input Validation ----
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Coerce y to matrix
  y_mat <- as.matrix(y)
  T_len <- nrow(y_mat)
  N <- ncol(y_mat)
  
  n_boot <- as.integer(n_boot)
  if (n_boot < 1L) {
    stop("'n_boot' must be at least 1.", call. = FALSE)
  }
  
  micro_block_length <- as.integer(micro_block_length)
  if (micro_block_length < 1L) {
    stop("'micro_block_length' must be at least 1.", call. = FALSE)
  }
  
  n_states <- state_info$n_states
  trans_mat <- state_info$transition_matrix
  
  ## ---- Compute Stationary Distribution ----
  
  stat_dist <- compute_stationary_dist(trans_mat)
  
  ## ---- Storage Allocation ----
  
  boot_samples <- array(NA_real_, dim = c(T_len, N, n_boot))
  boot_states <- matrix(NA_integer_, nrow = T_len, ncol = n_boot)
  boot_indices <- matrix(NA_integer_, nrow = T_len, ncol = n_boot)
  
  ## ---- Generate Bootstrap Replicates ----
  
  for (b in 1:n_boot) {
    
    # Step 1: Simulate new state sequence from Markov chain
    sim_states <- simulate_markov_chain(T_len, trans_mat, stat_dist)
    boot_states[, b] <- sim_states
    
    # Step 2: Sample time indices from state-specific pools
    if (sync_sampling || N == 1L) {
      # Synchronized sampling: same time index for all variables
      sampled_idx <- sample_from_pools_sync(sim_states, pools, micro_block_length)
      boot_indices[, b] <- sampled_idx
      boot_samples[, , b] <- y_mat[sampled_idx, , drop = FALSE]
    } else {
      # Independent sampling per variable (advanced option, not typically used)
      # This would require per-variable state info, which we don't have
      # Fall back to synchronized
      warning("Independent per-variable sampling not implemented. ",
              "Using synchronized sampling.", call. = FALSE)
      sampled_idx <- sample_from_pools_sync(sim_states, pools, micro_block_length)
      boot_indices[, b] <- sampled_idx
      boot_samples[, , b] <- y_mat[sampled_idx, , drop = FALSE]
    }
  }
  
  ## ---- Return Results ----
  
  list(
    samples = boot_samples,
    states = boot_states,
    sampled_indices = boot_indices,
    state_info = state_info,
    pools = pools,
    params = list(
      n_boot = n_boot,
      micro_block_length = micro_block_length,
      sync_sampling = sync_sampling,
      seed = seed
    )
  )
}


# =============================================================================
# RAW RETURNS HMM (NO GARCH) - Using Standalone EM Algorithm
# =============================================================================

#' Fit HMM with Skew Student-t Emissions (No GARCH)
#'
#' Fits a Hidden Markov Model with skew Student-t distributed emissions
#' directly to returns, without the GARCH volatility layer. Uses an EM
#' algorithm with numerical optimization for the M-step.
#'
#' @param y Numeric vector of returns.
#' @param n_states Integer, number of hidden states (default 2).
#' @param distribution Character, emission distribution. Currently supports
#'   \code{"sstd"} (skew Student-t, default), \code{"std"} (Student-t),
#'   and \code{"norm"} (Gaussian).
#' @param max_iter Integer, maximum EM iterations (default 200).
#' @param tol Numeric, convergence tolerance for log-likelihood (default 1e-6).
#' @param verbose Logical, print progress (default FALSE).
#'
#' @return A list with class \code{"hmm_sstd_fit"} containing:
#'   \describe{
#'     \item{transition_matrix}{K x K transition probabilities.}
#'     \item{initial_probs}{Length K initial state distribution.}
#'     \item{state_params}{List of distribution parameters per state:
#'       mu (mean), sigma (sd), nu (df for t), xi (skewness for sstd).}
#'     \item{states}{Viterbi-decoded state sequence.}
#'     \item{smoothed_probs}{T x K smoothed state probabilities.}
#'     \item{loglik}{Final log-likelihood.}
#'     \item{n_iter}{Number of EM iterations.}
#'     \item{converged}{Logical, whether EM converged.}
#'   }
#'
#' @details
#' This provides an alternative to MSGARCH when you want regime-switching
#' based purely on distributional characteristics (mean, variance, skewness,
#' kurtosis) without modeling GARCH dynamics.
#'
#' The EM algorithm alternates between:
#' \itemize{
#'   \item \strong{E-step}: Forward-backward algorithm to compute state
#'     probabilities given current parameters.
#'   \item \strong{M-step}: Update distribution parameters by maximizing
#'     the expected complete-data log-likelihood.
#' }
#'
#' The skew Student-t density uses the Fernández & Steel (1998) parameterization
#' as implemented in \code{fGarch::dsstd}.
#'
#' @seealso \code{\link{fit_msgarch_model}} for GARCH-based regime switching.
#'
#' @keywords internal
fit_hmm_sstd <- function(
    y,
    n_states = 2L,
    distribution = c("sstd", "std", "norm"),
    max_iter = 200L,
    tol = 1e-6,
    verbose = FALSE
) {
  
  ## ---- Input Validation ----
  
  if (!is.numeric(y) || !is.vector(y)) {
    stop("'y' must be a numeric vector.", call. = FALSE)
  }
  
  T_len <- length(y)
  if (T_len < 30) {
    warning("Series length < 30 may lead to unreliable HMM estimation.",
            call. = FALSE)
  }
  
  n_states <- as.integer(n_states)
  if (n_states < 2L) {
    stop("'n_states' must be at least 2.", call. = FALSE)
  }
  
  distribution <- match.arg(distribution)
  
  # Check fGarch for sstd/std
  if (distribution %in% c("sstd", "std")) {
    if (!requireNamespace("fGarch", quietly = TRUE)) {
      stop("Package 'fGarch' is required for '", distribution, "' distribution. ",
           "Install with: install.packages('fGarch')", call. = FALSE)
    }
  }
  
  ## ---- Initialize Parameters ----
  
  # Use k-means for initial state assignment
  
  km <- kmeans(y, centers = n_states, nstart = 10)
  init_states <- km$cluster
  
  # Initialize transition matrix with persistence
  trans_mat <- matrix(0.1 / (n_states - 1), nrow = n_states, ncol = n_states)
  diag(trans_mat) <- 0.9
  
  # Initialize emission parameters from clustered data
  state_params <- vector("list", n_states)
  for (k in 1:n_states) {
    y_k <- y[init_states == k]
    if (length(y_k) < 5) {
      # Fall back if cluster is too small
      y_k <- y[sample.int(T_len, min(20, T_len))]
    }
    
    state_params[[k]] <- list(
      mu = mean(y_k),
      sigma = sd(y_k),
      nu = 10,      # Moderate df
      xi = 1        # Symmetric initially
    )
  }
  
  # Initial distribution (stationary)
  init_probs <- compute_stationary_dist(trans_mat)
  
  ## ---- EM Algorithm ----
  
  loglik_prev <- -Inf
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    
    ## ---- E-Step: Forward-Backward ----
    
    fb_result <- .forward_backward_hmm(y, trans_mat, init_probs,
                                       state_params, distribution)
    
    gamma <- fb_result$gamma      # T x K smoothed probs
    xi_mat <- fb_result$xi        # T-1 x K x K transition probs
    loglik <- fb_result$loglik
    
    # Check convergence
    if (abs(loglik - loglik_prev) < tol) {
      converged <- TRUE
      if (verbose) message("EM converged at iteration ", iter)
      break
    }
    loglik_prev <- loglik
    
    if (verbose && iter %% 10 == 0) {
      message("Iteration ", iter, ": loglik = ", round(loglik, 4))
    }
    
    ## ---- M-Step: Update Parameters ----
    
    # Update transition matrix
    for (i in 1:n_states) {
      for (j in 1:n_states) {
        trans_mat[i, j] <- sum(xi_mat[, i, j]) / sum(gamma[-T_len, i])
      }
    }
    # Normalize rows
    trans_mat <- trans_mat / rowSums(trans_mat)
    
    # Update initial distribution
    init_probs <- gamma[1, ]
    init_probs <- init_probs / sum(init_probs)
    
    # Update emission parameters
    for (k in 1:n_states) {
      state_params[[k]] <- .update_emission_params(
        y, gamma[, k], state_params[[k]], distribution
      )
    }
  }
  
  if (!converged && verbose) {
    warning("EM did not converge within ", max_iter, " iterations.",
            call. = FALSE)
  }
  
  ## ---- Viterbi Decoding ----
  
  states <- .viterbi_hmm(y, trans_mat, init_probs, state_params, distribution)
  
  ## ---- Return Results ----
  
  result <- list(
    transition_matrix = trans_mat,
    initial_probs = init_probs,
    state_params = state_params,
    states = states,
    smoothed_probs = gamma,
    loglik = loglik,
    n_iter = iter,
    converged = converged,
    distribution = distribution,
    n_states = n_states
  )
  
  class(result) <- "hmm_sstd_fit"
  return(result)
}


#' Forward-Backward Algorithm for HMM
#'
#' Computes smoothed state probabilities using the forward-backward algorithm.
#'
#' @param y Numeric vector of observations.
#' @param trans_mat Transition matrix.
#' @param init_probs Initial state probabilities.
#' @param state_params List of emission parameters per state.
#' @param distribution Emission distribution type.
#'
#' @return List with gamma (smoothed probs), xi (transition probs), loglik.
#'
#' @keywords internal
.forward_backward_hmm <- function(
    y,
    trans_mat,
    init_probs,
    state_params,
    distribution
) {
  
  T_len <- length(y)
  K <- length(init_probs)
  
  # Compute emission probabilities
  B <- matrix(0, nrow = T_len, ncol = K)
  for (k in 1:K) {
    B[, k] <- .emission_density(y, state_params[[k]], distribution)
  }
  
  # Forward pass (scaled)
  alpha <- matrix(0, nrow = T_len, ncol = K)
  scale <- numeric(T_len)
  
  alpha[1, ] <- init_probs * B[1, ]
  scale[1] <- sum(alpha[1, ])
  alpha[1, ] <- alpha[1, ] / scale[1]
  
  for (t in 2:T_len) {
    for (j in 1:K) {
      alpha[t, j] <- sum(alpha[t-1, ] * trans_mat[, j]) * B[t, j]
    }
    scale[t] <- sum(alpha[t, ])
    if (scale[t] > 0) {
      alpha[t, ] <- alpha[t, ] / scale[t]
    } else {
      # Underflow protection
      alpha[t, ] <- 1 / K
      scale[t] <- .Machine$double.xmin
    }
  }
  
  # Backward pass
  beta <- matrix(0, nrow = T_len, ncol = K)
  beta[T_len, ] <- 1
  
  for (t in (T_len - 1):1) {
    for (i in 1:K) {
      beta[t, i] <- sum(trans_mat[i, ] * B[t+1, ] * beta[t+1, ])
    }
    if (scale[t+1] > 0) {
      beta[t, ] <- beta[t, ] / scale[t+1]
    }
  }
  
  # Smoothed state probabilities (gamma)
  gamma <- alpha * beta
  gamma <- gamma / rowSums(gamma)
  
  # Transition probabilities (xi)
  xi <- array(0, dim = c(T_len - 1, K, K))
  for (t in 1:(T_len - 1)) {
    denom <- 0
    for (i in 1:K) {
      for (j in 1:K) {
        xi[t, i, j] <- alpha[t, i] * trans_mat[i, j] * B[t+1, j] * beta[t+1, j]
        denom <- denom + xi[t, i, j]
      }
    }
    if (denom > 0) {
      xi[t, , ] <- xi[t, , ] / denom
    }
  }
  
  # Log-likelihood
  loglik <- sum(log(scale))
  
  list(gamma = gamma, xi = xi, loglik = loglik)
}


#' Emission Density Evaluation
#'
#' Evaluates the emission density for each observation given state parameters.
#'
#' @param y Numeric vector of observations.
#' @param params List with mu, sigma, nu (for t), xi (for sstd).
#' @param distribution Character, distribution type.
#'
#' @return Numeric vector of density values.
#'
#' @keywords internal
.emission_density <- function(y, params, distribution) {
  
  mu <- params$mu
  sigma <- params$sigma
  
  switch(distribution,
         "norm" = dnorm(y, mean = mu, sd = sigma),
         "std" = {
           nu <- params$nu
           # Standardize and use fGarch
           z <- (y - mu) / sigma
           fGarch::dstd(z, mean = 0, sd = 1, nu = nu) / sigma
         },
         "sstd" = {
           nu <- params$nu
           xi <- params$xi
           # fGarch::dsstd handles location/scale internally
           fGarch::dsstd(y, mean = mu, sd = sigma, nu = nu, xi = xi)
         },
         stop("Unknown distribution: ", distribution)
  )
}


#' Update Emission Parameters (M-step)
#'
#' Updates emission distribution parameters using weighted MLE.
#'
#' @param y Numeric vector of observations.
#' @param weights Numeric vector of state responsibilities (gamma_k).
#' @param params Current parameter values.
#' @param distribution Distribution type.
#'
#' @return Updated parameter list.
#'
#' @keywords internal
.update_emission_params <- function(y, weights, params, distribution) {
  
  # Normalize weights
  w <- weights / sum(weights)
  
  # Weighted mean and sd
  mu_new <- sum(w * y)
  sigma_new <- sqrt(sum(w * (y - mu_new)^2))
  sigma_new <- max(sigma_new, 1e-6)  # Ensure positive
  
  new_params <- list(mu = mu_new, sigma = sigma_new)
  
  if (distribution == "norm") {
    return(new_params)
  }
  
  # For t and sstd, optimize nu (and xi for sstd) via numerical optimization
  if (distribution == "std") {
    # Optimize nu
    opt_result <- tryCatch({
      optim(
        par = params$nu,
        fn = function(nu) {
          if (nu <= 2) return(1e10)
          z <- (y - mu_new) / sigma_new
          -sum(weights * log(fGarch::dstd(z, mean = 0, sd = 1, nu = nu) / sigma_new + 1e-300))
        },
        method = "Brent",
        lower = 2.1,
        upper = 100
      )
    }, error = function(e) list(par = params$nu))
    
    new_params$nu <- opt_result$par
    
  } else if (distribution == "sstd") {
    # Optimize nu and xi jointly
    opt_result <- tryCatch({
      optim(
        par = c(params$nu, params$xi),
        fn = function(p) {
          nu <- p[1]
          xi <- p[2]
          if (nu <= 2 || xi <= 0) return(1e10)
          dens <- fGarch::dsstd(y, mean = mu_new, sd = sigma_new, nu = nu, xi = xi)
          -sum(weights * log(dens + 1e-300))
        },
        method = "L-BFGS-B",
        lower = c(2.1, 0.1),
        upper = c(100, 10)
      )
    }, error = function(e) list(par = c(params$nu, params$xi)))
    
    new_params$nu <- opt_result$par[1]
    new_params$xi <- opt_result$par[2]
  }
  
  return(new_params)
}


#' Viterbi Algorithm for HMM
#'
#' Finds the most likely state sequence using the Viterbi algorithm.
#'
#' @param y Numeric vector of observations.
#' @param trans_mat Transition matrix.
#' @param init_probs Initial state probabilities.
#' @param state_params Emission parameters per state.
#' @param distribution Emission distribution type.
#'
#' @return Integer vector of most likely states.
#'
#' @keywords internal
.viterbi_hmm <- function(
    y,
    trans_mat,
    init_probs,
    state_params,
    distribution
) {
  
  T_len <- length(y)
  K <- length(init_probs)
  
  # Compute log emission probabilities
  log_B <- matrix(0, nrow = T_len, ncol = K)
  for (k in 1:K) {
    dens <- .emission_density(y, state_params[[k]], distribution)
    log_B[, k] <- log(dens + 1e-300)
  }
  
  log_trans <- log(trans_mat + 1e-300)
  log_init <- log(init_probs + 1e-300)
  
  # Forward pass
  delta <- matrix(-Inf, nrow = T_len, ncol = K)
  psi <- matrix(0L, nrow = T_len, ncol = K)
  
  delta[1, ] <- log_init + log_B[1, ]
  
  for (t in 2:T_len) {
    for (j in 1:K) {
      candidates <- delta[t-1, ] + log_trans[, j]
      psi[t, j] <- which.max(candidates)
      delta[t, j] <- max(candidates) + log_B[t, j]
    }
  }
  
  # Backtrack
  states <- integer(T_len)
  states[T_len] <- which.max(delta[T_len, ])
  
  for (t in (T_len - 1):1) {
    states[t] <- psi[t + 1, states[t + 1]]
  }
  
  return(states)
}


#' Extract State Info from HMM-SSTD Fit
#'
#' Wrapper to extract state information from a \code{hmm_sstd_fit} object
#' in the same format as \code{extract_msgarch_states}.
#'
#' @param fit Object of class \code{hmm_sstd_fit}.
#' @param y Original data (for compatibility, not used).
#'
#' @return List compatible with \code{extract_msgarch_states} output.
#'
#' @keywords internal
extract_hmm_sstd_states <- function(fit, y = NULL) {
  
  if (!inherits(fit, "hmm_sstd_fit")) {
    stop("'fit' must be an hmm_sstd_fit object.", call. = FALSE)
  }
  
  trans_mat <- fit$transition_matrix
  n_states <- fit$n_states
  
  # Compute expected state durations
  state_durations <- numeric(n_states)
  for (k in 1:n_states) {
    p_stay <- trans_mat[k, k]
    if (p_stay < 1) {
      state_durations[k] <- 1 / (1 - p_stay)
    } else {
      state_durations[k] <- Inf
    }
  }
  
  list(
    states = fit$states,
    transition_matrix = trans_mat,
    smoothed_probs = fit$smoothed_probs,
    n_states = n_states,
    coefficients = unlist(fit$state_params),
    state_durations = state_durations
  )
}