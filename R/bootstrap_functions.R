#' Hidden Markov Model (HMM) Bootstrap for Time Series
#'
#' Fits a Gaussian Hidden Markov Model (HMM) to a univariate time series
#' and generates bootstrap replicates by resampling regime-specific blocks.
#'
#' @param x Numeric vector representing the time series.
#' @param n_boot Length of bootstrap series.
#' @param num_states Integer number of hidden states for the HMM.
#' @param num_blocks Integer number of blocks to sample for each bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates to generate.
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores.
#'
#' @return A list of numeric vectors, each one a bootstrap replicate.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits a Gaussian HMM to `x` using `depmixS4::depmix()` and `depmixS4::fit()`.
#'   \item Uses Viterbi decoding (`posterior(fit, type = "viterbi")$state`)
#'     to assign each observation to a state.
#'   \item Samples contiguous blocks of observations belonging to each state.
#' }
#' If `n_boot` is set, the last block will be truncated when necessary to match 
#'   the length (`n_boot`) of the bootstrap series. This is the only way to 
#'   ensure equal length of all bootstrap series, as the length of each block is 
#'   random. If `n_boot` is not set, `num_blocks` must be set, and the length of 
#'   each bootstrap series will be determined by the number of blocks and the 
#'   random lengths of the individual blocks for that particular series. Note 
#'   that this almost certainly results in bootstrap series of different lengths.
#'
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 200, list(ar = 0.5))
#' hmm_samples <- hmm_bootstrap(x, num_states = 2, num_blocks = 10, num_boots = 5)
#' length(hmm_samples)     # 5 bootstrap replicates
#' length(hmm_samples[[1]]) # length of one bootstrap series
#'
#' @importFrom depmixS4 depmix fit posterior
#' @importFrom stats gaussian
#' @export
hmm_bootstrap <- function(
    x, 
    n_boot = NULL, 
    num_states = 2, 
    num_blocks = 100, 
    num_boots = 100,
    parallel = FALSE, 
    num_cores = 1L
  ) {
  df <- data.frame(return = x)
  
  ## Fit HMM
  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    stop("depmixS4 package required for HMM bootstrap.")
  }
  
  model <- depmixS4::depmix(return ~ 1, data = df, nstates = num_states, family = gaussian())
  fit <- depmixS4::fit(model, verbose = FALSE)
  states <- depmixS4::posterior(fit, type = "viterbi")$state
  
  ## Sample bootstraps
  .sample_blocks(
    x, 
    n_boot, 
    num_blocks, 
    states, 
    num_boots,
    parallel = parallel, 
    num_cores = num_cores
  )
}


#' Markov-Switching Autoregressive (MSAR) Bootstrap for Time Series
#'
#' Fits a Markov-switching autoregressive model (MSAR) to a univariate time series
#' and generates bootstrap replicates by resampling regime-specific blocks.
#'
#' @param x Numeric vector representing the time series.
#' @param ar_order Integer order of the autoregressive model.
#' @param num_states Integer number of regimes (hidden states) in the MSAR model.
#' @param n_boot Length of bootstrap series.
#' @param num_blocks Integer number of blocks to sample per bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates.
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores.
#'
#' @return A list of numeric vectors, each a bootstrap replicate.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Fits a Markov-switching autoregressive model using `MSwM::msmFit()` on an `lm()` fit of `x`.
#'   \item Uses Viterbi decoding (`ms_model@Fit@smoProb`) to classify each observation into states.
#'   \item Groups contiguous observations belonging to the same state into blocks.
#'   \item Samples these regime-specific blocks with replacement to generate synthetic series.
#' }
#' If `n_boot` is set, the last block will be truncated when necessary to match 
#'   the length (`n_boot`) of the bootstrap series. This is the only way to 
#'   ensure equal length of all bootstrap series, as the length of each block is 
#'   random. If `n_boot` is not set, `num_blocks` must be set, and the length of 
#'   each bootstrap series will be determined by the number of blocks and the 
#'   random lengths of the individual blocks for that particular series. Note 
#'   that this almost certainly results in bootstrap series of different lengths.
#'   Note that if `n_boot` and `num_blocks` are both set, `num_blocks` will be 
#'   ignored. 
#'
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 200, list(ar = 0.7))
#' msar_samples <- msar_bootstrap(x, ar_order = 1, num_states = 2, num_blocks = 10, num_boots = 5)
#' length(msar_samples)     # 5 replicates
#' length(msar_samples[[1]]) # length of one bootstrap series
#'
#' @importFrom MSwM msmFit
#' @importFrom stats lm
#' @export
msar_bootstrap <- function(
    x, 
    ar_order = 1, 
    num_states = 2, 
    n_boot = NULL, 
    num_blocks = 100, 
    num_boots = 100,
    parallel = FALSE, 
    num_cores = 1L
  ) {

  ## Prepare lagged design matrix
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package required for MSAR bootstrap.")
  }
  if (!requireNamespace("MSwM", quietly = TRUE)) {
    stop("MSwM package required for MSAR bootstrap.")
  }
  
  df <- data.frame(y = x)
  
  for (i in seq_len(ar_order)) df[[paste0("lag", i)]] <- dplyr::lag(x, i)
  df <- na.omit(df)
  
  ## Fit base AR model
  base_model <- lm(y ~ ., data = df)
  
  ## Fit MSAR model
  ms_model <- MSwM::msmFit(base_model, k = num_states, sw = rep(TRUE, length(coef(base_model)) + 1))
  
  ## State sequence
  states <- ms_model@Fit@smoProb
  state_seq <- apply(states, 1, which.max)
  
  ## Sample bootstraps
  .sample_blocks(
    x, 
    n_boot, 
    num_blocks, 
    state_seq, 
    num_boots,
    parallel = parallel, 
    num_cores = num_cores
  )
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
#' @examples
#' set.seed(123)
#' resids <- rnorm(100)
#' boot_reps <- wild_bootstrap(resids, num_boots = 5)
#' length(boot_reps)           # 5 replicates
#' dim(boot_reps[[1]])         # 100 x 1 matrix
#'
#' @export
wild_bootstrap <- function(x, num_boots = 100) {
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


#' Sample bootstrap blocks
#' 
#' Internal helper function.
#'
#' @param x x
#' @param n_boot Length of bootstrap series
#' @param num_blocks Number of blocks
#' @param states States (see tsbs() help)
#' @param num_boots Number of bootstraps
#' @param parallel Boolean
#' @param num_cores Number of cores
#'
#' @returns
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#'
#' @examples
.sample_blocks <- function(
    x, 
    n_boot, 
    num_blocks, 
    states, 
    num_boots, 
    parallel = FALSE,
    num_cores = 1L
  ) {
  ## Get contiguous blocks per state
  r <- rle(states)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  blocks <- Map(
    function(s, e) list(start = s, end = e),
    starts, ends
  )
  
  ## Sample `num_blocks` blocks with replacement
  get_block <- function(b) x[b$start:b$end]
  
  
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
  
  ## Sample num_blocks blocks per bootstrap
  # sampled_series <- replicate(num_boots, {
  #   sampled_blocks <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
  #   unlist(sampled_blocks)
  # }, simplify = FALSE)
  
  # if(!is.null(n_boot)) {
  #   if(!is.null(num_blocks)) {
  #     warning("`num_blocks` is ingored when `n_boot` is set.")
  #   }
  #   sampled_series <- replicate(num_boots, {
  #     bootstrap_series <- numeric(0)
  #     while(length(bootstrap_series) < n_boot) {
  #       sampled_block <- sample(blocks, 1, replace = TRUE)[[1]]
  #       bootstrap_series <- c(bootstrap_series, get_block(sampled_block))
  #     }
  #     bootstrap_series <- bootstrap_series[seq_len(n_boot)]
  #   }, simplify = FALSE)
  # } else if(!is.null(num_blocks)) {
  #   sampled_series <- replicate(num_boots, {
  #     bootstrap_series <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
  #     unlist(bootstrap_series)
  #   }, simplify = FALSE)
  # } else {
  #   stop("Must provide valid value for either n_boot or num_blocks")
  # }
  
  ## ---- Parallelized sampling ----
  
  if (is.null(n_boot) && is.null(num_blocks)) {
    stop("Must provide a valid value for either n_boot or num_blocks")
  }
  if (!is.null(n_boot) && !is.null(num_blocks)) {
    warning("`num_blocks` is ignored when `n_boot` is set.")
  }
  
  sampled_series <- foreach::foreach(i = seq_len(num_boots)) %dopar% {
    if (!is.null(n_boot)) {
      ## Generate series to a fixed length n_boot
      bootstrap_series <- numeric(0)
      while (length(bootstrap_series) < n_boot) {
        sampled_block <- sample(blocks, 1, replace = TRUE)[[1]]
        bootstrap_series <- c(bootstrap_series, get_block(sampled_block))
      }
      bootstrap_series[seq_len(n_boot)]
    } else {
      # Generate series from a fixed number of blocks
      bootstrap_series <- lapply(sample(blocks, num_blocks, replace = TRUE), get_block)
      unlist(bootstrap_series)
    }
  }

  sampled_series
}