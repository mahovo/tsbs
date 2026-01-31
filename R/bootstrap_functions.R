#' Hidden Markov Model (HMM) Bootstrap for Multivariate Time Series
#' 
#' Fits a Gaussian Hidden Markov Model (HMM) to a multivariate time series
#'   and generates bootstrap replicates by resampling regime-specific blocks.
#' 
#' @param x Numeric vector or matrix representing the time series.
#' @param n_boot Length of bootstrap series.
#' @param num_states Integer number of hidden states for the HMM.
#' @param num_blocks Integer number of blocks to sample for each bootstrap
#'   replicate.
#' @param num_boots Integer number of bootstrap replicates to generate.
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores.
#' @param return_fit Logical. If TRUE, returns the fitted HMM model along with
#'   bootstrap samples. Default is FALSE.
#' @param collect_diagnostics Logical. If TRUE, collects detailed diagnostic
#'   information including regime composition of each bootstrap replicate.
#'   Default is FALSE.
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
#' If `n_boot` is set, the last block will be trimmed when necessary.
#' If `n_boot` is not set, and `num_blocks` is set, the length of each
#' bootstrap series will be determined by the number of blocks and the random
#' lengths of the individual blocks for that particular series.
#' If neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
#' number of rows in `x` and the last block will be trimmed when necessary.
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
#' When `collect_diagnostics = TRUE`, the function records:
#' \itemize{
#'   \item Original state sequence from Viterbi decoding
#'   \item State sequence for each bootstrap replicate
#'   \item Block information (lengths, source positions)
#' }
#' This information can be used with `plot_regime_composition()` to visualize
#' how the bootstrap samples are composed from different regimes.
#' 
#' @return
#' If `return_fit = FALSE` and `collect_diagnostics = FALSE`:
#'   A list of bootstrap replicate matrices.
#' 
#' If `return_fit = TRUE` or `collect_diagnostics = TRUE`:
#'   A list containing:
#'   \describe{
#'     \item{bootstrap_series}{List of bootstrap replicate matrices}
#'     \item{fit}{(if return_fit = TRUE) The fitted depmixS4 model object}
#'     \item{states}{The Viterbi state sequence for the original data}
#'     \item{diagnostics}{(if collect_diagnostics = TRUE) A tsbs_diagnostics object}
#'   }
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
#' 
#' 
#'   ## Example 6: Basic usage
#'   set.seed(123)
#'   x <- matrix(rnorm(500), ncol = 2)
#'   boot_samples <- hmm_bootstrap(x, n_boot = 400, num_states = 2, num_boots = 50)
#' 
#'   ## With diagnostics for regime visualization
#'   result <- hmm_bootstrap(
#'     x, n_boot = 400, num_states = 2, num_boots = 10,
#'     collect_diagnostics = TRUE, return_fit = TRUE
#'   )
#' 
#'   ## Plot regime composition
#'   plot_regime_composition(result, x)
#' }
#' }
#' 
#' @export

hmm_bootstrap_old <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_blocks = NULL,
    num_boots = 100,
    parallel = FALSE,
    num_cores = 1L,
    return_fit = FALSE,
    collect_diagnostics = FALSE,
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

  ## Basic sanity checks
  if (n < num_states * 3) {
    warning("Series length (", n, ") is very short relative to num_states (", num_states,
            "). HMM fitting may be unreliable.", call. = FALSE)
  }

  ## ---- Create Formula for depmixS4 ----
  df <- as.data.frame(x)
  colnames(df) <- paste0("V", seq_len(ncol(df)))

  if (k == 1) {
    formula <- as.formula("V1 ~ 1")
    family_spec <- list(gaussian())
  } else {
    formula_list <- lapply(colnames(df), function(v) as.formula(paste(v, "~ 1")))
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
         "Original error: ", e$message, call. = FALSE)
  })

  fit <- tryCatch({
    depmixS4::fit(model, verbose = verbose)
  }, error = function(e) {
    stop("HMM fitting failed. Original error: ", e$message, call. = FALSE)
  }, warning = function(w) {
    if (verbose) {
      warning("HMM fitting produced a warning: ", w$message, call. = FALSE)
    }
    suppressWarnings(depmixS4::fit(model, verbose = verbose))
  })

  ## ---- Extract State Sequence ----
  states <- tryCatch({
    depmixS4::posterior(fit, type = "viterbi")$state
  }, error = function(e) {
    stop("Failed to extract state sequence from fitted HMM. ",
         "Original error: ", e$message, call. = FALSE)
  })

  if (length(unique(states)) < num_states && verbose) {
    warning("Only ", length(unique(states)), " of ", num_states,
            " states were observed in the Viterbi sequence.", call. = FALSE)
  }

  if (verbose) {
    state_counts <- table(states)
    message("State distribution: ",
            paste(names(state_counts), "=", state_counts, collapse = ", "))
  }

  ## ---- Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "hmm",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )

    ## Record original series stats
    diagnostics <- record_original_stats(diagnostics, x)

    ## Record regime information
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = states,
      num_states = num_states,
      state_labels = state_labels
    )

    ## Store HMM-specific config
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      num_blocks = num_blocks
    )
  }

  ## ---- Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating ", num_boots, " bootstrap samples...")
  }

  ## Use the enhanced sampling function that tracks regime composition
  sample_result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = states,
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores,
    collect_diagnostics = collect_diagnostics
  )

  bootstrap_samples <- sample_result$samples

  ## ---- Record Per-Replicate Diagnostics ----
  if (collect_diagnostics && !is.null(sample_result$replicate_info)) {
    for (i in seq_len(num_boots)) {
      rep_info <- sample_result$replicate_info[[i]]

      ## Record block information
      if (!is.null(rep_info$block_lengths)) {
        diagnostics <- record_blocks(
          diagnostics,
          replicate_idx = i,
          block_lengths = rep_info$block_lengths,
          start_positions = rep_info$block_starts,
          block_type = "regime"
        )
      }

      ## Record regime composition (including source indices for probability lookup)
      if (!is.null(rep_info$states)) {
        diagnostics <- record_replicate_regimes(
          diagnostics,
          replicate_idx = i,
          replicate_states = rep_info$states,
          block_states = rep_info$block_states,
          source_indices = rep_info$source_indices
        )
      }

      ## Record series statistics
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }

  if (verbose) {
    message("Bootstrap complete. Generated ", length(bootstrap_samples), " samples.")
  }

  ## ---- Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }

  result <- list(bootstrap_series = bootstrap_samples)

  if (return_fit) {
    result$fit <- fit
    result$states <- states

    ## Extract and store smoothed probabilities for plotting
    ## type = "viterbi" - Returns the most likely state sequence (used for $state)
    ## type = "smoothing" - Returns smoothed state probabilities P(S_t | all data) 
    ## (used for $S1, $S2, etc.)
    result$smoothed_probabilities <- tryCatch({
      post <- depmixS4::posterior(fit, type = "smoothing")
      prob_cols <- grep("^S", names(post), value = TRUE)
      as.matrix(post[, prob_cols])
    }, error = function(e) NULL)
  }

  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }

  return(result)
}

#' Hidden Markov Model (HMM) Bootstrap for Multivariate Time Series
#'
#' Fits a Hidden Markov Model (HMM) to a multivariate time series and generates
#' bootstrap replicates by resampling regime-specific blocks. Supports both
#' Gaussian emissions (via depmixS4) and non-Gaussian emissions including
#' skew Student-t distributions (via MSGARCH or standalone EM).
#'
#' @param x Numeric vector or matrix representing the time series.
#' @param n_boot Integer, length of each bootstrap series. If NULL, defaults to
#'   the length of the original series.
#' @param num_states Integer, number of hidden states for the HMM. Default is 2.
#' @param num_blocks Integer, number of blocks to sample for each bootstrap
#'   replicate. Only used when \code{distribution = "gaussian"}.
#' @param num_boots Integer, number of bootstrap replicates to generate.
#'   Default is 100.
#' @param distribution Character, emission distribution for the HMM. One of:
#'   \describe{
#'     \item{"gaussian"}{(Default) Gaussian emissions via depmixS4.
#'       This is the original behavior.}
#'     \item{"sstd"}{Skew Student-t via MSGARCH (with GARCH dynamics).}
#'     \item{"std"}{Student-t via MSGARCH (with GARCH dynamics).}
#'     \item{"snorm"}{Skew normal via MSGARCH (with GARCH dynamics).}
#'     \item{"sged"}{Skew GED via MSGARCH (with GARCH dynamics).}
#'     \item{"norm"}{Normal via MSGARCH (with GARCH dynamics).}
#'     \item{"ged"}{GED via MSGARCH (with GARCH dynamics).}
#'     \item{"sstd_raw"}{Skew Student-t without GARCH (standalone HMM-EM).}
#'     \item{"std_raw"}{Student-t without GARCH (standalone HMM-EM).}
#'     \item{"norm_raw"}{Normal without GARCH (standalone HMM-EM, equivalent
#'       to depmixS4 but using our EM implementation).}
#'   }
#' @param variance_model Character, GARCH specification when using MSGARCH-based
#'   distributions. One of \code{"sGARCH"} (default), \code{"eGARCH"},
#'   \code{"gjrGARCH"}, or \code{"tGARCH"}. Ignored for Gaussian and raw
#'   distributions.
#' @param micro_block_length Integer, block length for within-state sampling
#'   when using MSGARCH-based bootstrap. Use 1 (default) for iid sampling
#'   within states, or >1 to preserve some local dependence. Ignored for
#'   Gaussian distribution which uses the traditional block sampling.
#' @param regime_basis For multivariate data with MSGARCH: how to identify
#'   regimes. One of \code{"market"} (equal-weighted average, default),
#'   \code{"first_pc"} (first principal component), or an integer column index.
#'   Ignored for Gaussian distribution.
#' @param parallel Logical, parallelize computation? Default is FALSE.
#' @param num_cores Integer, number of cores for parallel processing.
#' @param return_fit Logical. If TRUE, returns the fitted model along with
#'   bootstrap samples. Default is FALSE.
#' @param collect_diagnostics Logical. If TRUE, collects detailed diagnostic
#'   information including regime composition. Default is FALSE.
#' @param verbose Logical. If TRUE, prints fitting information. Default is FALSE.
#' @param seed Integer, random seed for reproducibility. Default is NULL.
#' @param ... Additional arguments passed to the underlying fitting functions.
#'
#' @details
#' ## Method Selection
#'
#' When \code{distribution = "gaussian"}, uses depmixS4 for HMM fitting with
#' Gaussian emissions. This is the original behavior and uses block resampling
#' within identified regimes.
#'
#' When \code{distribution} is one of the MSGARCH distributions ("sstd", "std",
#' etc.), uses the MSGARCH package for Markov-switching GARCH estimation. The
#' bootstrap is then semi-parametric: state sequences are simulated from the
#' fitted Markov chain, while observations are resampled from state-specific
#' empirical pools.
#'
#' When \code{distribution} ends with "_raw" (e.g., "sstd_raw"), fits an HMM
#' directly to the returns without the GARCH volatility layer, using a
#' standalone EM algorithm with the specified emission distribution.
#'
#' ## Literature Background
#'
#' The MSGARCH-based implementation combines several established techniques:
#'
#' \itemize{
#'   \item \strong{Markov-switching GARCH}: Haas, Mittnik & Paolella (2004),
#'     implemented via the MSGARCH package (Ardia et al., 2019)
#'   \item \strong{Skew Student-t distribution}: Fernández & Steel (1998)
#'     transformation
#'   \item \strong{Semi-parametric bootstrap}: State sequences are simulated
#'     from the fitted Markov chain (parametric), while innovations are
#'     resampled from empirical state-specific pools (nonparametric)
#' }
#'
#' ## Multivariate Handling
#'
#' For multivariate data:
#' \itemize{
#'   \item With Gaussian: fits a joint HMM where all variables depend on the
#'     same hidden state sequence.
#'   \item With MSGARCH: fits the regime model to an aggregate series (see
#'     \code{regime_basis}), then samples full cross-sections synchronously
#'     to preserve empirical dependence structure.
#' }
#'
#' @return
#' If \code{return_fit = FALSE} and \code{collect_diagnostics = FALSE}:
#'   A list of bootstrap replicate matrices.
#'
#' If \code{return_fit = TRUE} or \code{collect_diagnostics = TRUE}:
#'   A list containing:
#'   \describe{
#'     \item{bootstrap_series}{List of bootstrap replicate matrices}
#'     \item{fit}{(if return_fit = TRUE) The fitted model object}
#'     \item{states}{The Viterbi state sequence for the original data}
#'     \item{smoothed_probabilities}{State probabilities matrix}
#'     \item{diagnostics}{(if collect_diagnostics = TRUE) A tsbs_diagnostics
#'       object}
#'     \item{method}{Character indicating which method was used}
#'   }
#'
#' @references
#' Ardia, D., Bluteau, K., Boudt, K., Catania, L., & Trottier, D.-A. (2019).
#' Markov-Switching GARCH Models in R: The MSGARCH Package.
#' Journal of Statistical Software, 91(4), 1-38. doi:10.18637/jss.v091.i04
#'
#' Fernández, C., & Steel, M. F. (1998). On Bayesian modeling of fat tails
#' and skewness. Journal of the American Statistical Association, 93(441),
#' 359-371.
#'
#' Haas, M., Mittnik, S., & Paolella, M. S. (2004). A New Approach to Markov-
#' Switching GARCH Models. Journal of Financial Econometrics, 2, 493-530.
#'
#' Hamilton, J. D. (1989). A New Approach to the Economic Analysis of
#' Nonstationary Time Series and the Business Cycle. Econometrica, 57(2),
#' 357-384.
#'
#' Holst, U., Lindgren, G., Holst, J. and Thuvesholmen, M. (1994), Recursive
#' Estimation In Switching Autoregressions With A Markov Regime. Journal of
#' Time Series Analysis, 15: 489-506.
#'
#' @examples
#' \donttest{
#' ## Example 1: Traditional Gaussian HMM bootstrap (original behavior)
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol = 2)
#' boot_gaussian <- hmm_bootstrap(x, n_boot = 400, num_states = 2, num_boots = 50)
#'
#' ## Example 2: Skew Student-t with MSGARCH (requires MSGARCH package)
#' if (requireNamespace("MSGARCH", quietly = TRUE)) {
#'   # Univariate
#'   y <- rnorm(300)
#'   boot_sstd <- hmm_bootstrap(y, num_states = 2, num_boots = 50,
#'                               distribution = "sstd")
#'
#'   # With diagnostics
#'   result <- hmm_bootstrap(y, num_states = 2, num_boots = 20,
#'                            distribution = "sstd",
#'                            return_fit = TRUE,
#'                            collect_diagnostics = TRUE)
#' }
#'
#' ## Example 3: Raw returns HMM without GARCH (requires fGarch)
#' if (requireNamespace("fGarch", quietly = TRUE)) {
#'   boot_raw <- hmm_bootstrap(y, num_states = 2, num_boots = 50,
#'                              distribution = "sstd_raw")
#' }
#' }
#'
#' @seealso \code{\link{msvar_bootstrap}}, \code{\link{ms_varma_garch_bs}}
#'
#' @export
hmm_bootstrap <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_blocks = NULL,
    num_boots = 100,
    distribution = c("gaussian", "sstd", "std", "snorm", "sged", "norm", "ged",
                     "sstd_raw", "std_raw", "norm_raw"),
    variance_model = c("sGARCH", "eGARCH", "gjrGARCH", "tGARCH"),
    micro_block_length = 1L,
    regime_basis = "market",
    parallel = FALSE,
    num_cores = 1L,
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    seed = NULL,
    ...
) {
  
  ## ---- Match Arguments ----
  distribution <- match.arg(distribution)
  variance_model <- match.arg(variance_model)
  
  ## ---- Dispatch Based on Distribution ----
  if (distribution == "gaussian") {
    ## Original Gaussian HMM via depmixS4
    .hmm_bootstrap_gaussian(
      x = x,
      n_boot = n_boot,
      num_states = num_states,
      num_blocks = num_blocks,
      num_boots = num_boots,
      parallel = parallel,
      num_cores = num_cores,
      return_fit = return_fit,
      collect_diagnostics = collect_diagnostics,
      verbose = verbose,
      ...
    )
    
  } else if (grepl("_raw$", distribution)) {
    ## Raw returns HMM (no GARCH) via standalone EM
    raw_dist <- sub("_raw$", "", distribution)
    .hmm_bootstrap_raw(
      x = x,
      n_boot = n_boot,
      num_states = num_states,
      num_boots = num_boots,
      distribution = raw_dist,
      micro_block_length = micro_block_length,
      regime_basis = regime_basis,
      return_fit = return_fit,
      collect_diagnostics = collect_diagnostics,
      verbose = verbose,
      seed = seed,
      ...
    )
    
  } else {
    ## MSGARCH-based bootstrap (with GARCH dynamics)
    .hmm_bootstrap_msgarch(
      x = x,
      n_boot = n_boot,
      num_states = num_states,
      num_boots = num_boots,
      distribution = distribution,
      variance_model = variance_model,
      micro_block_length = micro_block_length,
      regime_basis = regime_basis,
      return_fit = return_fit,
      collect_diagnostics = collect_diagnostics,
      verbose = verbose,
      seed = seed,
      ...
    )
  }
}


#' Gaussian HMM Bootstran
#'
#' Internal function implementing the original Gaussian HMM bootstrap via
#' depmixS4. This preserves the original behavior for backward compatibility.
#'
#' @inheritParams hmm_bootstrap
#' @keywords internal
.hmm_bootstrap_gaussian <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_blocks = NULL,
    num_boots = 100,
    parallel = FALSE,
    num_cores = 1L,
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    ...
) {
  
  ## ---- Dependency Check ----
  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    stop("Package 'depmixS4' is required for Gaussian HMM bootstrap. ",
         "Install with: install.packages('depmixS4')", call. = FALSE)
  }
  
  ## ---- Data Preparation ----
  x <- as.matrix(x)
  n <- nrow(x)
  k <- ncol(x)
  
  ## Basic sanity checks
  if (n < num_states * 3) {
    warning("Series length (", n, ") is very short relative to num_states (",
            num_states, "). HMM fitting may be unreliable.", call. = FALSE)
  }
  
  ## ---- Create Formula for depmixS4 ----
  df <- as.data.frame(x)
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  
  if (k == 1) {
    formula <- as.formula("V1 ~ 1")
    family_spec <- list(gaussian())
  } else {
    formula_list <- lapply(colnames(df), function(v) as.formula(paste(v, "~ 1")))
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
         "Original error: ", e$message, call. = FALSE)
  })
  
  fit <- tryCatch({
    depmixS4::fit(model, verbose = verbose)
  }, error = function(e) {
    stop("HMM fitting failed. Original error: ", e$message, call. = FALSE)
  }, warning = function(w) {
    if (verbose) {
      warning("HMM fitting produced a warning: ", w$message, call. = FALSE)
    }
    suppressWarnings(depmixS4::fit(model, verbose = verbose))
  })
  
  ## ---- Extract State Sequence ----
  states <- tryCatch({
    depmixS4::posterior(fit, type = "viterbi")$state
  }, error = function(e) {
    stop("Failed to extract state sequence from fitted HMM. ",
         "Original error: ", e$message, call. = FALSE)
  })
  
  if (length(unique(states)) < num_states && verbose) {
    warning("Only ", length(unique(states)), " of ", num_states,
            " states were observed in the Viterbi sequence.", call. = FALSE)
  }
  
  if (verbose) {
    state_counts <- table(states)
    message("State distribution: ",
            paste(names(state_counts), "=", state_counts, collapse = ", "))
  }
  
  ## ---- Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "hmm",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    diagnostics <- record_original_stats(diagnostics, x)
    
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = states,
      num_states = num_states,
      state_labels = state_labels
    )
    
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      num_blocks = num_blocks,
      distribution = "gaussian"
    )
  }
  
  ## ---- Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating ", num_boots, " bootstrap samples...")
  }
  
  sample_result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = states,
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores,
    collect_diagnostics = collect_diagnostics
  )
  
  bootstrap_samples <- sample_result$samples
  
  ## ---- Record Per-Replicate Diagnostics ----
  if (collect_diagnostics && !is.null(sample_result$replicate_info)) {
    for (i in seq_len(num_boots)) {
      rep_info <- sample_result$replicate_info[[i]]
      
      if (!is.null(rep_info$block_lengths)) {
        diagnostics <- record_blocks(
          diagnostics,
          replicate_idx = i,
          block_lengths = rep_info$block_lengths,
          start_positions = rep_info$block_starts,
          block_type = "regime"
        )
      }
      
      if (!is.null(rep_info$states)) {
        diagnostics <- record_replicate_regimes(
          diagnostics,
          replicate_idx = i,
          replicate_states = rep_info$states,
          block_states = rep_info$block_states,
          source_indices = rep_info$source_indices
        )
      }
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }
  
  if (verbose) {
    message("Bootstrap complete. Generated ", length(bootstrap_samples), " samples.")
  }
  
  ## ---- Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  result <- list(bootstrap_series = bootstrap_samples)
  result$method <- "hmm_gaussian"
  
  if (return_fit) {
    result$fit <- fit
    result$states <- states
    
    result$smoothed_probabilities <- tryCatch({
      post <- depmixS4::posterior(fit, type = "smoothing")
      prob_cols <- grep("^S", names(post), value = TRUE)
      as.matrix(post[, prob_cols])
    }, error = function(e) NULL)
  }
  
  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }
  
  return(result)
}


#' MSGARCH-based HMM Bootstrap
#'
#' Internal function implementing HMM bootstrap using MSGARCH for non-Gaussian
#' conditional distributions with GARCH dynamics.
#'
#' @inheritParams hmm_bootstrap
#' @keywords internal
.hmm_bootstrap_msgarch <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_boots = 100,
    distribution = "sstd",
    variance_model = "sGARCH",
    micro_block_length = 1L,
    regime_basis = "market",
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    seed = NULL,
    ...
) {
  
  ## ---- Dependency Check ----
  if (!requireNamespace("MSGARCH", quietly = TRUE)) {
    stop("Package 'MSGARCH' is required for '", distribution, "' distribution. ",
         "Install with: install.packages('MSGARCH')", call. = FALSE)
  }
  
  ## ---- Set Seed ----
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## ---- Data Preparation ----
  x_mat <- as.matrix(x)
  n <- nrow(x_mat)
  k <- ncol(x_mat)
  
  if (is.null(n_boot)) {
    n_boot <- n
  }
  
  ## Basic sanity checks
  if (n < 50) {
    warning("Series length (", n, ") is short. ",
            "MSGARCH estimation may be unreliable.", call. = FALSE)
  }
  
  ## ---- Determine Regime Series ----
  if (k > 1) {
    regime_series <- get_regime_series(x_mat, regime_basis)
    if (verbose) {
      message("Multivariate data: using '", regime_basis,
              "' for regime identification.")
    }
  } else {
    regime_series <- x_mat[, 1]
  }
  
  ## ---- Fit MSGARCH Model ----
  if (verbose) {
    message("Fitting ", num_states, "-state MSGARCH(", variance_model,
            ") with ", distribution, " distribution...")
  }
  
  fit <- fit_msgarch_model(
    y = regime_series,
    n_states = num_states,
    variance_model = variance_model,
    distribution = distribution,
    ...
  )
  
  ## ---- Extract State Information ----
  state_info <- extract_msgarch_states(fit, regime_series)
  
  if (verbose) {
    message("States identified. Expected durations: ",
            paste(round(state_info$state_durations, 1), collapse = ", "))
  }
  
  ## ---- Build State Pools ----
  pools <- extract_state_pools(fit, regime_series, state_info)
  
  if (verbose) {
    pool_sizes <- sapply(pools, function(p) p$n)
    message("State pool sizes: ", paste(pool_sizes, collapse = ", "))
  }
  
  ## ---- Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "hmm_msgarch",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    diagnostics <- record_original_stats(diagnostics, x_mat)
    
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = state_info$states,
      num_states = num_states,
      state_labels = state_labels
    )
    
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      distribution = distribution,
      variance_model = variance_model,
      micro_block_length = micro_block_length,
      regime_basis = regime_basis
    )
  }
  
  ## ---- Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating ", num_boots, " bootstrap samples...")
  }
  
  boot_result <- generate_msgarch_bootstrap(
    fit = fit,
    y = x_mat,
    state_info = state_info,
    pools = pools,
    n_boot = num_boots,
    micro_block_length = micro_block_length,
    sync_sampling = TRUE,
    seed = NULL  # Already set above if provided
  )
  
  ## Convert 3D array to list of matrices for consistency
  bootstrap_samples <- lapply(1:num_boots, function(b) {
    boot_result$samples[, , b, drop = FALSE]
  })
  ## Drop the third dimension
  bootstrap_samples <- lapply(bootstrap_samples, function(m) {
    matrix(m, nrow = dim(m)[1], ncol = dim(m)[2])
  })
  
  ## Trim to n_boot if different from original length
  if (n_boot != n) {
    bootstrap_samples <- lapply(bootstrap_samples, function(m) {
      m[seq_len(min(n_boot, nrow(m))), , drop = FALSE]
    })
  }
  
  ## ---- Record Per-Replicate Diagnostics ----
  if (collect_diagnostics) {
    for (i in seq_len(num_boots)) {
      diagnostics <- record_replicate_regimes(
        diagnostics,
        replicate_idx = i,
        replicate_states = boot_result$states[, i],
        block_states = NULL,
        source_indices = boot_result$sampled_indices[, i]
      )
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }
  
  if (verbose) {
    message("Bootstrap complete. Generated ", length(bootstrap_samples), " samples.")
  }
  
  ## ---- Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  result <- list(bootstrap_series = bootstrap_samples)
  result$method <- "hmm_msgarch"
  
  if (return_fit) {
    result$fit <- fit
    result$states <- state_info$states
    result$smoothed_probabilities <- state_info$smoothed_probs
    result$transition_matrix <- state_info$transition_matrix
    result$state_info <- state_info
    result$pools <- pools
  }
  
  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }
  
  return(result)
}


#' Raw Returns HMM Bootstrap (No GARCH)
#'
#' Internal function implementing HMM bootstrap using standalone EM algorithm
#' with non-Gaussian emissions but without GARCH dynamics.
#'
#' @inheritParams hmm_bootstrap
#' @keywords internal
.hmm_bootstrap_raw <- function(
    x,
    n_boot = NULL,
    num_states = 2,
    num_boots = 100,
    distribution = "sstd",
    micro_block_length = 1L,
    regime_basis = "market",
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = FALSE,
    seed = NULL,
    ...
) {
  
  ## ---- Set Seed ----
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ## ---- Data Preparation ----
  x_mat <- as.matrix(x)
  n <- nrow(x_mat)
  k <- ncol(x_mat)
  
  if (is.null(n_boot)) {
    n_boot <- n
  }
  
  ## ---- Determine Regime Series ----
  if (k > 1) {
    regime_series <- get_regime_series(x_mat, regime_basis)
    if (verbose) {
      message("Multivariate data: using '", regime_basis,
              "' for regime identification.")
    }
  } else {
    regime_series <- x_mat[, 1]
  }
  
  ## ---- Fit HMM with sstd Emissions ----
  if (verbose) {
    message("Fitting ", num_states, "-state HMM with ", distribution,
            " emissions (no GARCH)...")
  }
  
  fit <- fit_hmm_sstd(
    y = regime_series,
    n_states = num_states,
    distribution = distribution,
    verbose = verbose,
    ...
  )
  
  if (verbose && fit$converged) {
    message("EM converged after ", fit$n_iter, " iterations.")
  }
  
  ## ---- Extract State Information ----
  state_info <- extract_hmm_sstd_states(fit)
  
  if (verbose) {
    message("States identified. Expected durations: ",
            paste(round(state_info$state_durations, 1), collapse = ", "))
  }
  
  ## ---- Build State Pools ----
  # Create pools structure compatible with extract_state_pools output
  pools <- vector("list", num_states)
  for (kk in 1:num_states) {
    idx <- which(state_info$states == kk)
    pools[[kk]] <- list(
      indices = idx,
      n = length(idx),
      proportion = length(idx) / n
    )
  }
  
  if (verbose) {
    pool_sizes <- sapply(pools, function(p) p$n)
    message("State pool sizes: ", paste(pool_sizes, collapse = ", "))
  }
  
  ## ---- Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "hmm_sstd_raw",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    diagnostics <- record_original_stats(diagnostics, x_mat)
    
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = state_info$states,
      num_states = num_states,
      state_labels = state_labels
    )
    
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      distribution = paste0(distribution, "_raw"),
      micro_block_length = micro_block_length,
      regime_basis = regime_basis
    )
  }
  
  ## ---- Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating ", num_boots, " bootstrap samples...")
  }
  
  ## Compute stationary distribution
  stat_dist <- compute_stationary_dist(state_info$transition_matrix)
  
  ## Storage
  bootstrap_samples <- vector("list", num_boots)
  boot_states <- matrix(NA_integer_, nrow = n, ncol = num_boots)
  boot_indices <- matrix(NA_integer_, nrow = n, ncol = num_boots)
  
  for (b in 1:num_boots) {
    # Simulate state sequence
    sim_states <- simulate_markov_chain(n, state_info$transition_matrix, stat_dist)
    boot_states[, b] <- sim_states
    
    # Sample from pools
    sampled_idx <- sample_from_pools_sync(sim_states, pools, micro_block_length)
    boot_indices[, b] <- sampled_idx
    
    # Extract bootstrap sample
    boot_mat <- x_mat[sampled_idx, , drop = FALSE]
    
    # Trim if needed
    if (n_boot != n) {
      boot_mat <- boot_mat[seq_len(min(n_boot, nrow(boot_mat))), , drop = FALSE]
    }
    
    bootstrap_samples[[b]] <- boot_mat
  }
  
  ## ---- Record Per-Replicate Diagnostics ----
  if (collect_diagnostics) {
    for (i in seq_len(num_boots)) {
      diagnostics <- record_replicate_regimes(
        diagnostics,
        replicate_idx = i,
        replicate_states = boot_states[seq_len(nrow(bootstrap_samples[[i]])), i],
        block_states = NULL,
        source_indices = boot_indices[seq_len(nrow(bootstrap_samples[[i]])), i]
      )
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }
  
  if (verbose) {
    message("Bootstrap complete. Generated ", length(bootstrap_samples), " samples.")
  }
  
  ## ---- Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  result <- list(bootstrap_series = bootstrap_samples)
  result$method <- "hmm_sstd_raw"
  
  if (return_fit) {
    result$fit <- fit
    result$states <- state_info$states
    result$smoothed_probabilities <- state_info$smoothed_probs
    result$transition_matrix <- state_info$transition_matrix
    result$state_params <- fit$state_params
  }
  
  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }
  
  return(result)
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
#' @param return_fit Logical. If TRUE, returns the fitted MS-VAR model along with
#'   bootstrap samples. Default is FALSE.
#' @param collect_diagnostics Logical. If TRUE, collects detailed diagnostic
#'   information including regime composition of each bootstrap replicate.
#'   Default is FALSE.
#' @param verbose Logical. If TRUE, prints fitting information. Default is TRUE.
#'   
#' @details
#' This function:
#' \itemize{
#'   \item Fits a 2-state MS-VAR(1) model using \code{fit_msvar()}
#'   \item Uses smoothed probabilities to determine the most likely state sequence
#'   \item Samples contiguous blocks of observations within each regime
#' }
#' 
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
#' If `n_boot` is set, the last block will be trimmed when necessary.
#' If `n_boot` is not set, and `num_blocks` is set, the length of each 
#' bootstrap series will be determined by the number of blocks and the random 
#' lengths of the individual blocks for that particular series.
#' If neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
#' number of rows in `x` and the last block will be trimmed when necessary.
#' 
#' When \code{collect_diagnostics = TRUE}, the function records:
#' \itemize{
#'   \item Original state sequence from smoothed probabilities
#'   \item State sequence for each bootstrap replicate
#'   \item Block information (lengths, source positions)
#'   \item Source indices for probability reconstruction
#' }
#'
#' @return 
#' If \code{return_fit = FALSE} and \code{collect_diagnostics = FALSE}: 
#'   A list of bootstrap replicate matrices.
#'   
#' Otherwise, a list containing:
#' \describe{
#'   \item{bootstrap_series}{List of bootstrap replicate matrices}
#'   \item{fit}{(if return_fit = TRUE) The fitted MS-VAR model object}
#'   \item{states}{The state sequence for the original data}
#'   \item{smoothed_probabilities}{(if return_fit = TRUE) Matrix of smoothed state probabilities}
#'   \item{diagnostics}{(if collect_diagnostics = TRUE) A tsbs_diagnostics object}
#' }
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' T_obs <- 250
#' y1 <- arima.sim(model = list(ar = 0.7), n = T_obs)
#' y2 <- 0.5 * y1 + arima.sim(model = list(ar = 0.3), n = T_obs)
#' sample_data <- cbind(y1, y2)
#'
#' # Basic bootstrap
#' boot_samples <- msvar_bootstrap(sample_data, num_boots = 50)
#'
#' # With diagnostics for visualization
#' result <- msvar_bootstrap(
#'   sample_data, 
#'   num_boots = 10,
#'   return_fit = TRUE,
#'   collect_diagnostics = TRUE
#' )
#' plot_regime_composition(result, sample_data)
#' }
#'
#' @export
msvar_bootstrap <- function(
    x,
    n_boot = NULL,
    num_blocks = NULL,
    num_boots = 100,
    parallel = FALSE,
    num_cores = 1L,
    return_fit = FALSE,
    collect_diagnostics = FALSE,
    verbose = TRUE
) {
  
  ## ---- Data Preparation ----
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  n <- nrow(x)
  k <- ncol(x)
  
  ## ---- 1. Fit the MS-VAR(1) Model ----
  if (verbose) {
    message("Fitting the MS-VAR(1) model...")
  }
  ms_model <- fit_msvar(x)
  
  ## ---- 2. Determine the Most Likely State Sequence ----
  smoothed_probs <- ms_model$smoothed_probabilities
  state_seq <- apply(smoothed_probs, 1, which.max)
  
  ## The state sequence is one observation shorter than the original series
  ## because the VAR model requires an initial lag. Prepend the first state.
  state_seq_aligned <- c(state_seq[1], state_seq)
  
  ## Also align smoothed probabilities (prepend first row)
  smoothed_probs_aligned <- rbind(smoothed_probs[1, , drop = FALSE], smoothed_probs)
  
  num_states <- ncol(smoothed_probs)
  
  ## ---- 3. Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "msvar",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    ## Record original series stats
    diagnostics <- record_original_stats(diagnostics, x)
    
    ## Record regime information
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = state_seq_aligned,
      num_states = num_states,
      state_labels = state_labels
    )
    
    ## Store config
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      num_blocks = num_blocks
    )
  }
  
  ## ---- 4. Generate Bootstrap Samples ----
  if (verbose) {
    message("Generating bootstrap samples...")
  }
  
  ## Use the enhanced sampling function that tracks regime composition
  sample_result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = state_seq_aligned,
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores,
    collect_diagnostics = collect_diagnostics
  )
  
  bootstrap_samples <- sample_result$samples
  
  ## ---- 5. Record Per-Replicate Diagnostics ----
  if (collect_diagnostics && !is.null(sample_result$replicate_info)) {
    for (i in seq_len(num_boots)) {
      rep_info <- sample_result$replicate_info[[i]]
      
      if (!is.null(rep_info$block_lengths)) {
        diagnostics <- record_blocks(
          diagnostics,
          replicate_idx = i,
          block_lengths = rep_info$block_lengths,
          start_positions = rep_info$block_starts,
          block_type = "regime"
        )
      }
      
      if (!is.null(rep_info$states)) {
        diagnostics <- record_replicate_regimes(
          diagnostics,
          replicate_idx = i,
          replicate_states = rep_info$states,
          block_states = rep_info$block_states,
          source_indices = rep_info$source_indices
        )
      }
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }
  
  ## ---- 6. Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  result <- list(bootstrap_series = bootstrap_samples)
  
  if (return_fit) {
    result$fit <- ms_model
    result$states <- state_seq_aligned
    result$smoothed_probabilities <- smoothed_probs_aligned
  }
  
  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }
  
  return(result)
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
#'   Default is `NULL`.
#' @param num_blocks An integer specifying the number of blocks to sample for
#'   the bootstrap. Defaults to 100.
#' @param num_boots An integer specifying the total number of bootstrapped series
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
#' 
#' If `n_boot` is set, the last block will be trimmed when necessary.
#' If `n_boot` is not set, and `num_blocks` is set, the length of each 
#' bootstrap series will be determined by the number of blocks and the random 
#' lengths of the individual blocks for that particular series.
#' If neither `n_boot` nor `num_blocks` is set, `n_boot` will default to the
#' number of rows in `x` and the last block will be trimmed when necessary.
#' 
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
  
  ## Redundant precaution: Ensure x is a matrix.
  ## At this point x is supposed to have been coerced into a matrix by tsbs()
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  
  n <- nrow(x)
  k <- ncol(x)
  
  ## Ensure number of bootstrapped series is valid
  if (is.null(num_boots) || !is.numeric(num_boots) || num_boots < 1) {
    stop("num_boots must be a positive integer.", call. = FALSE)
  }
  
  ## If n_boot is not set, num_blocks must be set, and the length of each 
  ## bootstrap series will be determined by the number of blocks and the random 
  ## lengths of the individual blocks for that particular series.
  ## If neither of n_boot and num_blocks are set, n_boot will be set to the 
  ## number of rows in x
  if (is.null(n_boot) && is.null(num_blocks)) {
    #stop("Must provide a valid value for either n_boot or num_blocks", call. = FALSE)
    n_boot <- n
  }
  if (!is.null(n_boot) && !is.null(num_blocks)) {
    warning("`num_blocks` is ignored when `n_boot` is set.")
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
  
  ## Handle potential leading NAs from differencing by forward-filling the first 
  ## valid state.
  if (anyNA(state_seq)) {
    first_valid_state_idx <- min(which(!is.na(state_seq)))
    first_valid_state <- state_seq[first_valid_state_idx]
    state_seq[1:(first_valid_state_idx - 1)] <- first_valid_state
    
    ## Also fill NAs in smoothed_probs for consistency
    for (i in 1:(first_valid_state_idx - 1)) {
      smoothed_probs[i, ] <- smoothed_probs[first_valid_state_idx, ]
    }
  }
  
  num_states <- M
  
  ## ---- 4. Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = "ms_varma_garch",
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    ## Record original series stats
    diagnostics <- record_original_stats(diagnostics, x)
    
    ## Record regime information
    state_labels <- paste("State", seq_len(num_states))
    diagnostics <- record_regime_info(
      diagnostics,
      original_states = state_seq,
      num_states = num_states,
      state_labels = state_labels
    )
    
    ## Store config
    diagnostics$config <- list(
      num_states = num_states,
      n_boot = n_boot,
      num_blocks = num_blocks,
      model_type = model_type,
      d = d
    )
    
    ## If ms_model has its own diagnostics, merge them
    if (!is.null(ms_model$diagnostics)) {
      diagnostics$method_specific$ms_model_diagnostics <- ms_model$diagnostics
    }
  }
  
  ## ---- 5. Generate Bootstrap Samples ----
  ## This part calls the existing helper function that performs the actual
  ## stationary block bootstrap based on the state sequence.
  message("Generating bootstrap samples...")

  sample_result <- .sample_blocks_with_diagnostics(
    x = x,
    n_boot = n_boot,
    num_blocks = num_blocks,
    states = state_seq,
    num_boots = num_boots,
    parallel = parallel,
    num_cores = num_cores,
    collect_diagnostics = collect_diagnostics
  )
  
  bootstrap_samples <- sample_result$samples
  
  ## ---- 6. Record Per-Replicate Diagnostics ----
  if (collect_diagnostics && !is.null(sample_result$replicate_info)) {
    for (i in seq_len(num_boots)) {
      rep_info <- sample_result$replicate_info[[i]]
      
      if (!is.null(rep_info$block_lengths)) {
        diagnostics <- record_blocks(
          diagnostics,
          replicate_idx = i,
          block_lengths = rep_info$block_lengths,
          start_positions = rep_info$block_starts,
          block_type = "regime"
        )
      }
      
      if (!is.null(rep_info$states)) {
        diagnostics <- record_replicate_regimes(
          diagnostics,
          replicate_idx = i,
          replicate_states = rep_info$states,
          block_states = rep_info$block_states,
          source_indices = rep_info$source_indices
        )
      }
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = i,
        bootstrap_matrix = bootstrap_samples[[i]]
      )
    }
  }
  
  ## ---- 7. Return Results ----
  if (!return_fit && !collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  result <- list(bootstrap_series = bootstrap_samples)
  
  if (return_fit) {
    result$fit <- ms_model
    result$states <- state_seq
    result$smoothed_probabilities <- smoothed_probs
  }
  
  if (collect_diagnostics) {
    result$diagnostics <- diagnostics
  }
  
  return(result)
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
  
  ## If n_boot is not set, num_blocks must be set, and the length of each 
  ## bootstrap series will be determined by the number of blocks and the random 
  ## lengths of the individual blocks for that particular series.
  ## If neither of n_boot and num_blocks are set, n_boot will be set to the 
  ## number of rows in x
  if (is.null(n_boot) && is.null(num_blocks)) {
    #stop("Must provide a valid value for either n_boot or num_blocks")
    n_boot <- nrow(x) 
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


#' Enhanced Block Sampling with Diagnostic Tracking
#'
#' Internal function that performs regime-based block sampling while optionally
#' tracking diagnostic information about each bootstrap replicate.
#'
#' @param x Original data matrix.
#' @param n_boot Target bootstrap length.
#' @param num_blocks Number of blocks to sample.
#' @param states State sequence for original data.
#' @param num_boots Number of bootstrap replicates.
#' @param parallel Use parallel processing.
#' @param num_cores Number of cores.
#' @param collect_diagnostics Track diagnostic info.
#'
#' @return List with samples and optional replicate_info.
#' @keywords internal
.sample_blocks_with_diagnostics <- function(
    x,
    n_boot,
    num_blocks,
    states,
    num_boots,
    parallel = FALSE,
    num_cores = 1L,
    collect_diagnostics = FALSE
) {
  
  ## If n_boot is not set, num_blocks must be set, and the length of each 
  ## bootstrap series will be determined by the number of blocks and the random 
  ## lengths of the individual blocks for that particular series.
  ## If neither of n_boot and num_blocks are set, n_boot will be set to the 
  ## number of rows in x
  if (is.null(n_boot) && is.null(num_blocks)) {
    #stop("Must provide a valid value for either n_boot or num_blocks")
    n_boot <- nrow(x) 
  }
  if (!is.null(n_boot) && !is.null(num_blocks)) {
    warning("`num_blocks` is ignored when `n_boot` is set.")
  }
  
  ## ---- Identify blocks ----
  r <- rle(states)
  ends <- cumsum(r$lengths)
  n_orig <- nrow(x)
  n_states <- length(states)
  offset <- n_orig - n_states
  
  starts <- c(1, head(ends, -1) + 1)
  block_states <- r$values
  
  blocks <- Map(
    function(s, e, st) list(start = s + offset, end = e + offset, state = st),
    starts, ends, block_states
  )
  
  get_block <- function(b) {
    x[b$start:b$end, , drop = FALSE]
  }
  
  ## ---- Sampling Function ----
  sample_one_replicate <- function(idx) {
    if (!is.null(n_boot)) {
      ## Fixed-length sampling
      bootstrap_series <- x[0, , drop = FALSE]
      sampled_block_indices <- integer(0)
      source_row_indices <- integer(0)  # Track original row indices
      
      while (nrow(bootstrap_series) < n_boot) {
        block_idx <- sample(length(blocks), 1)
        sampled_block_indices <- c(sampled_block_indices, block_idx)
        b <- blocks[[block_idx]]
        source_row_indices <- c(source_row_indices, b$start:b$end)
        bootstrap_series <- rbind(bootstrap_series, get_block(b))
      }
      
      final_length <- n_boot
      bootstrap_series <- bootstrap_series[seq_len(final_length), , drop = FALSE]
      source_row_indices <- source_row_indices[seq_len(final_length)]
      
    } else {
      ## Variable-length sampling (num_blocks determines length)
      sampled_block_indices <- sample(length(blocks), num_blocks, replace = TRUE)
      bootstrap_blocks <- lapply(sampled_block_indices, function(i) get_block(blocks[[i]]))
      bootstrap_series <- do.call(rbind, bootstrap_blocks)
      
      ## Track source indices
      source_row_indices <- unlist(lapply(sampled_block_indices, function(i) {
        b <- blocks[[i]]
        b$start:b$end
      }))
      
      final_length <- nrow(bootstrap_series)
    }
    
    ## Compute replicate state sequence
    if (collect_diagnostics) {
      rep_states <- integer(0)
      block_lengths <- integer(0)
      block_starts <- integer(0)
      block_state_vec <- integer(0)
      
      for (bi in sampled_block_indices) {
        b <- blocks[[bi]]
        block_len <- b$end - b$start + 1
        rep_states <- c(rep_states, rep(b$state, block_len))
        block_lengths <- c(block_lengths, block_len)
        block_starts <- c(block_starts, b$start)
        block_state_vec <- c(block_state_vec, b$state)
      }
      
      ## Trim to match series length
      rep_states <- rep_states[seq_len(final_length)]
      
      list(
        series = bootstrap_series,
        info = list(
          states = rep_states,
          block_lengths = block_lengths,
          block_starts = block_starts,
          block_states = block_state_vec,
          source_indices = source_row_indices  # For probability lookup
        )
      )
    } else {
      list(series = bootstrap_series, info = NULL)
    }
  }
  
  ## ---- Run Sampling ----
  if (parallel && num_cores > 1) {
    if (!requireNamespace("foreach", quietly = TRUE) || 
        !requireNamespace("doParallel", quietly = TRUE)) {
      stop("Packages 'foreach' and 'doParallel' required for parallel execution.")
    }
    
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    `%dopar%` <- foreach::`%dopar%`
    results <- foreach::foreach(i = seq_len(num_boots)) %dopar% {
      sample_one_replicate(i)
    }
  } else {
    results <- lapply(seq_len(num_boots), sample_one_replicate)
  }
  
  ## ---- Extract Results ----
  samples <- lapply(results, `[[`, "series")
  replicate_info <- if (collect_diagnostics) {
    lapply(results, `[[`, "info")
  } else {
    NULL
  }
  
  list(samples = samples, replicate_info = replicate_info)
}