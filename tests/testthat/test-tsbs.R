#### ______________________________________________________________________ ####
#### PART 1: tsbs() End-to-End Pipeline Tests                              ####

context("tsbs End-to-End Pipeline")

## -----------------------------------------------------------------------------
## Setup: devFixtures for pipeline tests
## -----------------------------------------------------------------------------

## Small multivariate data for fast tests
set.seed(999)
y_mv_small <- simulate_dcc_garch(
  n = 150, k = 2,
  omega = c(0.05, 0.08),
  alpha_garch = c(0.08, 0.10),
  beta_garch = c(0.85, 0.82),
  alpha_dcc = 0.05,
  beta_dcc = 0.90,
  seed = 999
)
colnames(y_mv_small) <- c("Asset1", "Asset2")

## 3-asset data
set.seed(888)
y_mv_3asset <- simulate_dcc_garch(
  n = 150, k = 3,
  omega = c(0.05, 0.06, 0.07),
  alpha_garch = c(0.08, 0.09, 0.10),
  beta_garch = c(0.85, 0.84, 0.83),
  alpha_dcc = 0.04,
  beta_dcc = 0.92,
  seed = 888
)
colnames(y_mv_3asset) <- c("A", "B", "C")

## Spec generators
make_2state_dcc_spec <- function(k) {
  spec_uni <- list(model = "garch", garch_order = c(1, 1), distribution = "norm")
  
  list(
    ## State 1
    list(
      var_order = 1,
      garch_spec_fun = "dcc_modelspec",
      distribution = "mvn",
      garch_spec_args = list(
        dcc_order = c(1, 1),
        dynamics = "dcc",
        garch_model = list(
          univariate = replicate(k, spec_uni, simplify = FALSE)
        )
      ),
      start_pars = list(
        var_pars = rep(0, k * (1 + k)),
        garch_pars = replicate(k, list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85), simplify = FALSE),
        dcc_pars = list(alpha_1 = 0.03, beta_1 = 0.92)
      )
    ),
    ## State 2
    list(
      var_order = 1,
      garch_spec_fun = "dcc_modelspec",
      distribution = "mvn",
      garch_spec_args = list(
        dcc_order = c(1, 1),
        dynamics = "dcc",
        garch_model = list(
          univariate = replicate(k, spec_uni, simplify = FALSE)
        )
      ),
      start_pars = list(
        var_pars = rep(0, k * (1 + k)),
        garch_pars = replicate(k, list(omega = 0.10, alpha1 = 0.12, beta1 = 0.78), simplify = FALSE),
        dcc_pars = list(alpha_1 = 0.06, beta_1 = 0.88)
      )
    )
  )
}

spec_2asset <- make_2state_dcc_spec(2)
spec_3asset <- make_2state_dcc_spec(3)


## -----------------------------------------------------------------------------
## 1a: Basic tsbs() with ms_varma_garch returns correct structure
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch returns correct structure", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_type(result, "list")
  expect_named(result, c("bootstrap_series", "func_outs", "func_out_means"))
  expect_length(result$bootstrap_series, 3)
  expect_true(all(sapply(result$bootstrap_series, is.matrix)))
  expect_equal(ncol(result$bootstrap_series[[1]]), 2)
})


## -----------------------------------------------------------------------------
## 1b: tsbs() with return_fit = TRUE
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch with return_fit = TRUE includes fit object", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = TRUE,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_named(result, c("bootstrap_series", "func_outs", "func_out_means", "fit"))
  expect_type(result$fit, "list")
  expect_true("model_fits" %in% names(result$fit))
  expect_true("P" %in% names(result$fit))
  expect_true("log_likelihood" %in% names(result$fit))
  expect_true("smoothed_probabilities" %in% names(result$fit))
})


## -----------------------------------------------------------------------------
## 1c: tsbs() with return_fit = TRUE and collect_diagnostics = TRUE
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch with diagnostics includes diagnostics in fit", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = TRUE,
    collect_diagnostics = TRUE,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_true("diagnostics" %in% names(result$fit))
  expect_s3_class(result$fit$diagnostics, "ms_diagnostics")
})


## -----------------------------------------------------------------------------
## 1d: tsbs() with custom func (apply_func_to = "all")
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch applies custom func to each replicate", {
  skip_on_cran()
  
  ## Simple function that returns column means
  my_func <- function(df) {
    colMeans(df)
  }
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 5,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    func = my_func,
    apply_func_to = "all",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_length(result$func_outs, 5)
  expect_true(all(sapply(result$func_outs, length) == 2))
  expect_true(is.numeric(result$func_out_means))
  expect_length(result$func_out_means, 2)
})


## -----------------------------------------------------------------------------
## 1e: tsbs() with portfolio optimization func
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch works with portfolio optimization func", {
  skip_on_cran()
  
  ## Simple min variance (no quadprog dependency)
  simple_min_var <- function(returns) {
    Sigma <- cov(returns)
    n <- ncol(returns)
    vols <- sqrt(diag(Sigma))
    weights <- (1 / vols) / sum(1 / vols)
    names(weights) <- colnames(returns)
    weights
  }
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 5,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    func = simple_min_var,
    apply_func_to = "all",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  ## Each func_out should be weights that sum to 1
  for (i in seq_along(result$func_outs)) {
    expect_equal(sum(result$func_outs[[i]]), 1, tolerance = 1e-10)
    expect_true(all(result$func_outs[[i]] >= 0))
  }
  
  ## func_out_means should also sum to 1
  expect_equal(sum(result$func_out_means), 1, tolerance = 1e-10)
})


## -----------------------------------------------------------------------------
## 1f: tsbs() with 3 assets
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch works with 3 assets", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_3asset,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_3asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_length(result$bootstrap_series, 3)
  expect_equal(ncol(result$bootstrap_series[[1]]), 3)
})


## -----------------------------------------------------------------------------
## 1g: tsbs() bootstrap series have correct length
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch bootstrap series have correct length", {
  skip_on_cran()
  
  num_blocks <- 10
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = num_blocks,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  ## Each bootstrap series should have approximately num_blocks * avg_block_length rows
  ## The exact length depends on state sequence, so just check it's reasonable
  for (i in seq_along(result$bootstrap_series)) {
    expect_gt(nrow(result$bootstrap_series[[i]]), 0)
    expect_lte(nrow(result$bootstrap_series[[i]]), nrow(y_mv_small) * 2)
  }
})


## -----------------------------------------------------------------------------
## 1h: tsbs() with n_boot truncates correctly
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch respects n_boot truncation", {
  skip_on_cran()
  
  n_boot <- 50
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    n_boot = n_boot,
    num_boots = 3,
    num_blocks = 20,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  for (i in seq_along(result$bootstrap_series)) {
    expect_equal(nrow(result$bootstrap_series[[i]]), n_boot)
  }
})


## -----------------------------------------------------------------------------
## 1i: tsbs() return_fit = FALSE doesn't include fit
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch with return_fit = FALSE excludes fit", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = FALSE,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_false("fit" %in% names(result))
})


## -----------------------------------------------------------------------------
## 1j: tsbs() with func that returns scalar
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch works with scalar-returning func", {
  skip_on_cran()
  
  ## Function that returns portfolio volatility
  port_vol <- function(returns) {
    w <- rep(1/ncol(returns), ncol(returns))
    port_ret <- returns %*% w
    sd(port_ret)
  }
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 5,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    func = port_vol,
    apply_func_to = "all",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_true(all(sapply(result$func_outs, is.numeric)))
  expect_true(all(sapply(result$func_outs, length) == 1))
  expect_length(result$func_out_means, 1)
})



#### ______________________________________________________________________ ####
#### PART 2: ms_varma_garch_bs() Direct Tests                              ####

context("ms_varma_garch_bs Direct Tests")

## -----------------------------------------------------------------------------
## 2a: ms_varma_garch_bs returns correct structure (return_fit = FALSE)
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs returns list of matrices when return_fit = FALSE", {
  skip_on_cran()
  
  result <- ms_varma_garch_bs(
    x = y_mv_small,
    num_blocks = 5,
    num_boots = 3,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = FALSE,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_type(result, "list")
  expect_length(result, 3)
  expect_true(all(sapply(result, is.matrix)))
})


## -----------------------------------------------------------------------------
## 2b: ms_varma_garch_bs returns correct structure (return_fit = TRUE)
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs returns list with fit when return_fit = TRUE", {
  skip_on_cran()
  
  result <- ms_varma_garch_bs(
    x = y_mv_small,
    num_blocks = 5,
    num_boots = 3,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = TRUE,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_type(result, "list")
  expect_named(result, c("bootstrap_series", "fit"))
  expect_length(result$bootstrap_series, 3)
  expect_type(result$fit, "list")
})


## -----------------------------------------------------------------------------
## 2c: ms_varma_garch_bs generates valid state sequences
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs generates bootstrap from fitted states", {
  skip_on_cran()
  
  result <- ms_varma_garch_bs(
    x = y_mv_small,
    num_blocks = 10,
    num_boots = 5,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = TRUE,
    control = list(max_iter = 10, tol = 0.1)
  )
  
  ## Check that smoothed probabilities are valid
  probs <- result$fit$smoothed_probabilities
  expect_equal(ncol(probs), 2)
  expect_true(all(probs >= 0 & probs <= 1))
  expect_true(all(abs(rowSums(probs) - 1) < 1e-6))
})


## -----------------------------------------------------------------------------
## 2d: ms_varma_garch_bs with parallel = TRUE
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs works with parallel = TRUE", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  result <- ms_varma_garch_bs(
    x = y_mv_small,
    num_blocks = 5,
    num_boots = 3,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    return_fit = FALSE,
    parallel = TRUE,
    num_cores = 2,
    control = list(max_iter = 5, tol = 0.5)
  )
  
  expect_length(result, 3)
  expect_true(all(sapply(result, is.matrix)))
})


## -----------------------------------------------------------------------------
## 2e: ms_varma_garch_bs verbose output works
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs verbose output doesn't error", {
  skip_on_cran()
  
  temp_file <- tempfile(fileext = ".txt")
  
  expect_no_error({
    result <- ms_varma_garch_bs(
      x = y_mv_small,
      num_blocks = 5,
      num_boots = 2,
      M = 2,
      spec = spec_2asset,
      model_type = "multivariate",
      return_fit = FALSE,
      verbose = TRUE,
      verbose_file = temp_file,
      control = list(max_iter = 3, tol = 0.5)
    )
  })
  
  ## Check verbose file was created and has content
  expect_true(file.exists(temp_file))
  expect_gt(file.size(temp_file), 0)
  
  unlink(temp_file)
})


#### ______________________________________________________________________ ####
#### PART 3: .sample_blocks() Tests                                        ####

context("Bootstrap Block Sampling")

## -----------------------------------------------------------------------------
## 3a: .sample_blocks produces valid output
## -----------------------------------------------------------------------------
test_that(".sample_blocks produces list of matrices", {
  skip_on_cran()
  
  ## First need to fit to get state sequence
  fit <- fit_ms_varma_garch(
    y = y_mv_small,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  ## Get most likely state sequence
  states <- apply(fit$smoothed_probabilities, 1, which.max)
  
  ## Sample blocks
  result <- .sample_blocks(
    x = y_mv_small,
    n_boot = NULL,
    num_blocks = 10,
    states = states,
    num_boots = 5,
    parallel = FALSE,
    num_cores = 1
  )
  
  expect_type(result, "list")
  expect_length(result, 5)
  expect_true(all(sapply(result, is.matrix)))
  expect_true(all(sapply(result, function(m) ncol(m) == 2)))
})


## -----------------------------------------------------------------------------
## 3b: .sample_blocks respects n_boot
## -----------------------------------------------------------------------------
test_that(".sample_blocks respects n_boot parameter", {
  skip_on_cran()
  
  fit <- fit_ms_varma_garch(
    y = y_mv_small,
    M = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  states <- apply(fit$smoothed_probabilities, 1, which.max)
  
  n_boot <- 75
  
  result <- .sample_blocks(
    x = y_mv_small,
    n_boot = n_boot,
    num_blocks = 20,
    states = states,
    num_boots = 3,
    parallel = FALSE,
    num_cores = 1
  )
  
  for (i in seq_along(result)) {
    expect_equal(nrow(result[[i]]), n_boot)
  }
})


#### ______________________________________________________________________ ####
#### PART 4: Error Handling and Edge Cases                                 ####

context("Error Handling and Edge Cases")

## -----------------------------------------------------------------------------
## 4a: tsbs() with invalid bs_type for return_fit
## -----------------------------------------------------------------------------
test_that("tsbs return_fit is ignored for non-ms_varma_garch types", {
  skip_on_cran()
  
  ## For moving block bootstrap, return_fit should be ignored
  result <- tsbs(
    x = y_mv_small,
    bs_type = "moving",
    block_length = 10,
    num_blocks = 5,
    num_boots = 3,
    return_fit = TRUE  # Should be ignored
  )
  
  expect_false("fit" %in% names(result))
})


## -----------------------------------------------------------------------------
## 4b: tsbs() handles func that errors gracefully
## -----------------------------------------------------------------------------
test_that("tsbs handles func errors appropriately", {
  skip_on_cran()
  
  ## Function that sometimes errors
  bad_func <- function(df) {
    if (nrow(df) < 50) stop("Not enough data")
    colMeans(df)
  }
  
  ## With short bootstrap samples, this should error
  expect_error(
    tsbs(
      x = y_mv_small,
      bs_type = "ms_varma_garch",
      n_boot = 30,  # Short samples
      num_boots = 3,
      num_blocks = 5,
      num_states = 2,
      spec = spec_2asset,
      model_type = "multivariate",
      func = bad_func,
      apply_func_to = "all",
      control = list(max_iter = 5, tol = 0.5)
    )
  )
})


## -----------------------------------------------------------------------------
## 4c: ms_varma_garch_bs with very short data
## -----------------------------------------------------------------------------
test_that("ms_varma_garch_bs handles short data gracefully", {
  skip_on_cran()
  
  y_short <- y_mv_small[1:50, ]
  
  ## Should work but possibly with warnings
  result <- suppressWarnings(
    ms_varma_garch_bs(
      x = y_short,
      num_blocks = 5,
      num_boots = 2,
      M = 2,
      spec = spec_2asset,
      model_type = "multivariate",
      return_fit = FALSE,
      control = list(max_iter = 3, tol = 1)
    )
  )
  
  expect_type(result, "list")
  expect_length(result, 2)
})


## -----------------------------------------------------------------------------
## 4d: tsbs() with single state (M = 1)
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch handles single state correctly", {
  skip_on_cran()
  
  ## Single state spec
  spec_1state <- list(spec_2asset[[1]])
  
  expect_error(
    result <- tsbs(
      x = y_mv_small,
      bs_type = "ms_varma_garch",
      num_boots = 3,
      num_blocks = 5,
      num_states = 1,
      spec = spec_1state,
      model_type = "multivariate",
      control = list(max_iter = 5, tol = 0.5)
    )
  )
})



#### ______________________________________________________________________ ####
#### PART 5: Consistency and Reproducibility                               ####

context("Consistency and Reproducibility")

## -----------------------------------------------------------------------------
## 5a: tsbs() results are reproducible with seed
## -----------------------------------------------------------------------------
test_that("tsbs ms_varma_garch is reproducible with set.seed", {
  skip_on_cran()
  
  set.seed(123)
  result1 <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  set.seed(123)
  result2 <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 3,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  ## Bootstrap series should be identical
  for (i in 1:3) {
    expect_equal(result1$bootstrap_series[[i]], result2$bootstrap_series[[i]])
  }
})


## -----------------------------------------------------------------------------
## 5b: func_out_means is correctly computed
## -----------------------------------------------------------------------------
test_that("tsbs func_out_means equals mean of func_outs", {
  skip_on_cran()
  
  result <- tsbs(
    x = y_mv_small,
    bs_type = "ms_varma_garch",
    num_boots = 5,
    num_blocks = 5,
    num_states = 2,
    spec = spec_2asset,
    model_type = "multivariate",
    func = colMeans,
    apply_func_to = "all",
    control = list(max_iter = 5, tol = 0.5)
  )
  
  ## Manually compute mean of func_outs
  manual_mean <- Reduce(`+`, result$func_outs) / length(result$func_outs)
  
  expect_equal(result$func_out_means, manual_mean)
})