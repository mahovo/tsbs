#' Advanced Bootstrap for Time Series
#'
#' Generates bootstrap replicates of multivariate or (multiple) univariate time
#' series. Supports moving and stationary block, HMM, MSVAR, MS VARMA GARCH and 
#' wild bootstrap types.
#'
#' @param x Numeric vector, matrix, data frame or time series observations 
#'   (rows = time points, cols = variables).
#' @param n_boot Integer, optional desired length of each bootstrap replicate.
#'   See details below.
#' @param block_length Integer, length of each block. If `NULL`, an automatic 
#'   heuristic is used: The method first calculates the average absolute 
#'   first-order autocorrelation (\eqn{\rho_1}) across all time series columns. 
#'   A candidate block length is then calculated based on this persistence 
#'   measure using the formula \eqn{⌊10/(1−\rho_1)⌋}. The final block length is 
#'   constrained to be between 5 and the square root of the series length. 
#'   For stationary bootstrap, `block_length` is the expected block length, when 
#'   `p_method="1/n"`.
#' @param bs_type Bootstrap type. Character string: One of `"moving"`, 
#'   `"stationary"`, `"hmm"`, `"msvar"`, `"ms_varma_garch"`, or `"wild"`. 
#'   `"msvar"` is a lightweight special case of `"ms_varma_garch"`. See details 
#'   below.
#' @param block_type Block type. Character string: One of `"non-overlapping"`, 
#'   `"overlapping"`, or `"tapered"`. Only affects `bs_type="moving"` and 
#'   `bs_type="stationary"`. `block_type="tapered"` will smooth out transitions
#'   between blocks. It follows that `block_type="tapered"` can not be used when
#'   `block_length=1`.
#' @param taper_type Tapering window function. Character. One of `"cosine"`, 
#'   `"bartlett"`, or `"tukey"`. Only affects `block_type="tapered"`. See
#'   details below.
#' @param tukey_alpha numeric, alpha parameter for `taper_type = "tukey"`. 
#'   \eqn{\alpha \in [0, 1]} is the fraction of the window length tapered at 
#'   each end, using indices \eqn{0 \dots n−1}, where \eqn{n} is the block length.
#' @param num_blocks Integer number of blocks per bootstrap replicate.
#' @param num_boots Integer number of bootstrap replicates.
#' @param func A function to apply to each bootstrap replicate. The function's
#'   expected input depends on the `apply_func_to` parameter:
#'   \itemize{
#'     \item If `apply_func_to = "cols"`: `func` receives a single numeric 
#'       vector (one column at a time) and should return a scalar or numeric 
#'       vector. Examples: `mean`, `sd`, `quantile`.
#'     \item If `apply_func_to = "all"`: `func` receives the entire bootstrap 
#'       replicate as a numeric matrix (rows = time points, columns = variables)
#'       and can return:
#'       \itemize{
#'         \item A numeric scalar or vector (e.g., portfolio weights)
#'         \item A named list of numeric values (e.g., multiple summary 
#'           statistics)
#'       }
#'       Examples: portfolio optimization functions, functions returning 
#'       multiple statistics.
#'   }
#'   The function outputs are collected in `func_outs` (one entry per 
#'   replicate), and `func_out_means` contains the element-wise average across
#'   replicates. For list-returning functions, each named element is averaged
#'   separately. Default is `mean`.
#' @param apply_func_to Character string: `"cols"` to apply columnwise or `"all"` 
#'   to apply on the full data frame or matrix.
#' @param p_method Character string to choose method for stationary bootstrap 
#'   parameter: `"1/n"`, `"plugin"`, or `"cv"`.
#' @param p numeric \eqn{p \in (0, 1)}. Probability parameter for the geometric 
#'   block length (used in Stationary Bootstrap). 
#' @param overlap Logical indicating if overlapping blocks are allowed.
#' @param num_states Integer number of states for HMM, MSVAR, or 
#'        MS-VARMA-GARCH models.
#' @param d Integer differencing order for the MS-VARMA-GARCH model.
#' @param spec A list defining the state-specific models for the MS-VARMA-GARCH
#'   bootstrap (bs_type = "ms_varma_garch"). This argument is required for
#'   this bootstrap type and must be a list of length `num_states`. Each element
#'   of the list is itself a list specifying the model for that state. The
#'   structure depends on whether the model is univariate or multivariate.
#'   Fitting relies on the tsgarch or tsmarch packages, respectively. See 
#'   details below.
#' @param model_type Character string for MS-VARMA-GARCH: `"univariate"` or 
#'   `"multivariate"`.
#' @param control A list of control parameters for the MS-VARMA-GARCH EM 
#'   algorithm.
#' @param model_func Model-fitting function for cross-validation. 
#'   See [k_fold_cv_ts].
#' @param score_func Scoring function for cross-validation.
#' @param stationary_max_percentile Stationary max percentile.
#' @param stationary_max_fraction_of_n Stationary max fraction of n.
#' @param return_diagnostics If `TRUE`, returns bootstrap diagnostics including
#'   block composition, length statistics, and quality metrics. Diagnostics are
#'   returned in `result$diagnostics` as a `tsbs_diagnostics` object. Use 
#'   `summary()` and `plot()` methods on the diagnostics for detailed analysis.
#'   Default is `FALSE`.
#' @param return_fit If `TRUE`, `tsbs()` will return model fit. Default is
#'   `return_fit = FALSE`. If `return_fit = TRUE` and 
#'   `bs_type = "ms_varma_garch"`, diagnostics can be extracted from 
#'   `result$fit$diagnostics`. See `?fit_ms_varma_garch`.
#' @param fail_mode Character string. One of `"predictably"` (development/
#'   debugging - fail fast on any unexpected behavior) or `"gracefully"` 
#'   (production/robust pipelines - continue despite validation errors). 
#' @param parallel Parallelize computation? `TRUE` or `FALSE`.
#' @param num_cores Number of cores when `parallel=TRUE`.
#' 
#' @details 
#' ## Bootstrap types
#' 
#' See documentation for the individual bootstrap functions: 
#' * For moving or stationary block bootstraps, see [blockBootstrap()]
#' * [hmm_bootstrap()]  
#' * [msvar_bootstrap()]
#' * [ms_varma_garch_bs()]
#' * [wild_bootstrap()]
#' 
#' ### `bs_type="moving"`
#' If `n_boot` is set, the last block will be truncated when necessary to
#'   match the length (`n_boot`) of the bootstrap series. If `n_boot` is not
#'   set, `block_length` and `num_blocks` must be set, and `n_boot` will
#'   automatically be set to `block_length * num_blocks`.
#' 
#' ### `bs_type="stationary"`, `bs_type="hmm"`, `bs_type="msvar"` If `n_boot` is
#' set, the last block will be truncated when necessary to match the length
#' (`n_boot`) of the bootstrap series. This is the only way to ensure equal
#' length of all bootstrap series, as the length of each block is random. If
#' `n_boot` is not set, `num_blocks` must be set, and the length of each
#' bootstrap series will be determined by the number of blocks and the lengths
#' of the individual blocks for that particular series. Truncation ensures equal
#' length of bootstrapped series for bootstrap types with random block length.
#' For stationary bootstrap, `block_length` is the expected block length, when
#' `p_method="1/n"`.
#' 
#' ### `bs_type="wild"`
#' `n_boot`, `block_length` and `num_blocks` are ignored. The length of the
#'   bootstrap series is always the same as the original series.
#'   
#' ### `bs_type = "ms_varma_garch"`
#' 
#' `spec` list structure:
#'
#' \itemize{
#' \item \strong{For univariate models (one-column x):}
#' Each element spec\[\[j\]\] must be a list with the following components. For a
#' complete list of all possible arguments (e.g., for different GARCH model
#' flavors), please refer to the documentation for the tsgarch package
#' (`?tsgarch` and
#' `vignette("garch_models", package = "tsgarch")`).
#'   \itemize{
#'     \item \code{arma_order}: A numeric vector c(p,q) for the ARMA(p,q) order.
#'     \item \code{garch_model}: A character string for the GARCH model type
#'     (e.g., "garch", "egarch").
#'     \item \code{garch_order}: A numeric vector c(g,h) for the GARCH(g,h) order.
#'     \item \code{distribution}: A character string for the conditional
#'     distribution.
#'     \itemize{
#'        \item `"norm"` = Normal
#'        \item `"snorm"` = Skew normal
#'        \item `"std"` = Student t
#'        \item `"sstd"` = Skew Student
#'        \item `"ged"`  = Generalized error
#'        \item `"ged"`  = Skew generalized error
#'        \item `"ghyp"`  = Generalized hyperbolic
#'        \item `"ghst"`  = Generalized hyperbolic skew Student
#'        \item `"jsu"`  = Johnson reparameterized SU
#'     }
#'     See \url{https://www.nopredict.com/packages/tsdistributions}.
#'     \item \code{start_pars}: A list containing the starting parameters for
#'     the optimization, with two named elements:
#'       \itemize{
#'         \item \code{arma_pars}: A named numeric vector for the ARMA parameters
#'         (e.g., c(ar1 = 0.5)).
#'         \item \code{garch_pars}: A named list for the GARCH parameters (e.g.,
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)).
#'         \item \code{dist_pars}: A named list for the starting values for
#'         univariate dist params (e.g. list(skew = 1.0, shape = 6.0)).
#'       }
#'   }
#'
#' \item \strong{For multivariate models (multi-column x):}  
#' Each element spec\[\[j\]\] must be a list with the following components. For a
#' complete list of all possible arguments, please refer to the documentation
#' for the relevant tsmarch specification function. For a DCC model see
#' `?dcc_modelspec.tsgarch.multi_estimate`, for a Copula GARCH model see
#' `?cgarch_modelspec.tsgarch.multi_estimate`. (See also `?tsmarch`,
#' `vignette("tsmarch_demo", package = "tsmarch")` and
#' `vignette("feasible_multivariate_garch", package = "tsmarch")`.)
#'   \itemize{
#'     \item \code{var_order}: An integer p for the VAR(p) order.
#'     \item \code{garch_spec_fun}: A character string with the name of the
#'     tsmarch specification function to use (e.g., "dcc_modelspec",
#'     "gogarch_modelspec").
#'     \item \code{distribution}: The multivariate distribution. Valid
#'          choices are "mvn" (for multivariate normal) and "mvt" (for
#'          multivariate Student's t).
#'     \item \code{garch_spec_args}: A list of arguments to pass to the function
#'     specified in garch_spec_fun. Key arguments include:
#'       \itemize{
#'          \item For models like DCC or CGARCH, this list must also contain a
#'          \code{garch_model} definition for the underlying univariate GARCH
#'          fits.
#'          Note: With the exception of the Copula model, the marginal
#'          distributions of the univariate GARCH models should always be
#'          Normal, irrespective of whether a multivariate Normal or Student is
#'          chosen as the DCC model distribution. There are no checks performed
#'          for this and it is up to the user to ensure that this is the case.
#'       }
#'     \item \code{start_pars}: A list containing the starting parameters, with
#'     two named elements:
#'       \itemize{
#'         \item \code{var_pars}: A numeric vector for the VAR parameters,
#'         stacked column by column (intercept, lags for eq1, intercept, lags
#'         for eq2, ...).
#'         \item \code{garch_pars}: A list of named lists for the multivariate 
#'         GARCH parameters, e.g. 
#'         \code{replicate(k, list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85), simplify = FALSE)}.
#'         \item \code{dcc pars}: A named list, e.g. 
#'         \code{list(alpha_1 = 0.05, beta_1 = 0.90)}.
#'         \item \code{dist pars}: A named list, e.g.
#'         \code{list(shape = 8.0)}.
#'       }
#'   }
#' } 
#' 
#' ### DCC (Dynamic Conditional Correlation) Models
#'
#' DCC models capture time-varying correlations between multiple series by 
#' modeling the correlation dynamics directly on standardized residuals. This
#' is the standard approach for multivariate GARCH modeling in the tsbs package.
#'
#' To use DCC, set `garch_spec_fun = "dcc_modelspec"` in your spec.
#' See [dcc_modelspec()] for the tsmarch specification function.
#'
#' **DCC-specific `garch_spec_args`:**
#' \itemize{
#'   \item `dcc_order`: Integer vector c(p, q) for DCC(p, q) order. Default is
#'     c(1, 1). Note: DCC(1,1) uses analytical gradients; higher orders use
#'     numerical differentiation.
#'   \item `dynamics`: Character. Correlation dynamics type:
#'     \itemize{
#'       \item `"constant"`: Time-invariant correlation (reduces to CCC model)
#'       \item `"dcc"`: Standard DCC dynamics
#'       \item `"adcc"`: Asymmetric DCC (negative shocks have larger impact)
#'     }
#'   \item `distribution`: Character. Multivariate distribution:
#'     \itemize{
#'       \item `"mvn"`: Multivariate Normal
#'       \item `"mvt"`: Multivariate Student-t (adds shape parameter)
#'     }
#' }
#'
#' **Important**: For DCC models, the univariate GARCH distributions should
#' always be Normal (`distribution = "norm"`), regardless of whether the
#' multivariate distribution is MVN or MVT. The multivariate distribution
#' applies to the joint model, not the marginals.
#'
#' **Example DCC specification:**
#' ```
#' spec_dcc <- list(
#'   list(
#'     var_order = 1,
#'     garch_spec_fun = "dcc_modelspec",
#'     distribution = "mvn",
#'     garch_spec_args = list(
#'       dcc_order = c(1, 1),
#'       dynamics = "dcc",
#'       distribution = "mvn",
#'       garch_model = list(univariate = list(
#'         list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
#'         list(model = "garch", garch_order = c(1, 1), distribution = "norm")
#'       ))
#'     ),
#'     start_pars = list(
#'       var_pars = rep(0, 6),
#'       garch_pars = list(
#'         list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85),
#'         list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85)
#'       ),
#'       dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
#'       dist_pars = NULL
#'     )
#'   )
#' )
#' ```
#'
#' For ADCC (asymmetric correlation dynamics), use `dynamics = "adcc"` and
#' include gamma parameters in `dcc_pars`:
#' ```
#' dcc_pars = list(alpha_1 = 0.05, gamma_1 = 0.03, beta_1 = 0.90)
#' ```
#'
#' @seealso
#' \itemize{
#'   \item [dcc_modelspec()] for creating DCC specifications
#'   \item [estimate_garch_weighted_dcc()] for the weighted estimation method
#'   \item `vignette("cgarch_vs_dcc", package = "tsbs")` for model comparison
#' }
#'
#'
#' ### Copula GARCH (CGARCH) Models
#'
#' Copula GARCH models provide an alternative to DCC for modeling dynamic
#' correlations in multivariate time series. While DCC directly models
#' correlation dynamics, CGARCH separates marginal distributions from the
#' dependence structure using copulas.
#'
#' To use CGARCH instead of DCC, set `garch_spec_fun = "cgarch_modelspec"` in
#' your spec. See [cgarch_modelspec()] for the tsmarch specification function.
#'
#' **Key differences from DCC:**
#' \itemize{
#'   \item Uses Probability Integral Transform (PIT) to transform standardized
#'     residuals to uniform margins before modeling dependence
#'   \item Supports multiple transformation methods: "parametric", "empirical",
#'     or "spd" (semi-parametric distribution)
#'   \item Copula distribution can be "mvn" (Gaussian) or "mvt" (Student-t)
#'   \item More flexible tail dependence modeling with Student-t copula
#' }
#'
#' **CGARCH-specific `garch_spec_args`:**
#' \itemize{
#'   \item `transformation`: Character. Method for PIT transformation:
#'     \itemize{
#'       \item `"parametric"`: Uses fitted univariate distribution CDFs
#'       \item `"empirical"`: Uses empirical CDFs (non-parametric)
#'       \item `"spd"`: Semi-parametric distribution combining kernel density
#'         in the center with parametric tails
#'     }
#'   \item `copula`: Character. Copula distribution: `"mvn"` or `"mvt"`
#'   \item `dynamics`: Character. Correlation dynamics type:
#'     \itemize{
#'       \item `"constant"`: Time-invariant correlation
#'       \item `"dcc"`: Dynamic Conditional Correlation
#'       \item `"adcc"`: Asymmetric DCC (captures leverage effects where
#'         negative returns have larger impact on correlation)
#'     }
#'   \item `dcc_order`: Integer vector c(p, q) for DCC/ADCC order
#' }
#'
#' **Example CGARCH specification:**
#' ```
#' spec_cgarch <- list(
#'   list(
#'     var_order = 1,
#'     garch_spec_fun = "cgarch_modelspec",
#'     distribution = "mvn",
#'     garch_spec_args = list(
#'       dcc_order = c(1, 1),
#'       dynamics = "dcc",
#'       transformation = "parametric",
#'       copula = "mvn",
#'       garch_model = list(univariate = list(
#'         list(model = "garch", garch_order = c(1, 1), distribution = "norm"),
#'         list(model = "garch", garch_order = c(1, 1), distribution = "norm")
#'       ))
#'     ),
#'     start_pars = list(
#'       var_pars = rep(0, 6),
#'       garch_pars = list(
#'         list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85),
#'         list(omega = 0.05, alpha1 = 0.08, beta1 = 0.85)
#'       ),
#'       dcc_pars = list(alpha_1 = 0.05, beta_1 = 0.90),
#'       dist_pars = NULL
#'     )
#'   )
#' )
#' ```
#'
#' For ADCC (asymmetric correlation dynamics), use `dynamics = "adcc"` and
#' include gamma parameters:
#' ```
#' garch_spec_args = list(
#'   dynamics = "adcc",
#'   ...
#' ),
#' start_pars = list(
#'   ...
#'   dcc_pars = list(alpha_1 = 0.05, gamma_1 = 0.02, beta_1 = 0.90)
#' )
#' ```
#'
#' @seealso
#' \itemize{
#'   \item [dcc_modelspec()] for creating DCC specifications
#'   \item [cgarch_modelspec()] for creating CGARCH specifications
#'   \item [estimate_garch_weighted_dcc()] for DCC weighted estimation
#'   \item [estimate_garch_weighted_cgarch()] for CGARCH weighted estimation
#'   \item [compute_pit_transform()] for PIT transformation details
#'   \item [adcc_recursion()] for asymmetric DCC dynamics
#'   \item `vignette("cgarch_vs_dcc", package = "tsbs")` for model comparison
#'   \item `vignette("Diagnostics", package = "tsbs")` for diagnostic system
#' }
#' 
#' **GOGARCH-specific `garch_spec_args`:**
#' \itemize{
#'   \item `model`: Character. GARCH model type for components. Default `"garch"`.
#'   \item `order`: Integer vector c(p, q) for GARCH order. Default `c(1, 1)`.
#'   \item `ica`: Character. ICA algorithm:
#'     \itemize{
#'       \item `"radical"`: RADICAL algorithm (recommended, robust to outliers)
#'       \item `"fastica"`: FastICA algorithm (currently not supported by tsbs, 
#'       waiting for tsmarch v1.0.1 to appear on CRAN...)
#'     }
#'   \item `components`: Integer. Number of ICA components to extract.
#'     Default equals number of series. Set lower for dimension reduction.
#'   \item `lambda_range`: Numeric vector c(min, max) for GH lambda parameter.
#'   \item `shape_range`: Numeric vector c(min, max) for GH shape parameter.
#' }
#'
#' **GOGARCH Distributions:**
#' \itemize{
#'   \item `"norm"`: Normal distribution (simplest, fastest)
#'   \item `"nig"`: Normal Inverse Gaussian (captures skewness and kurtosis)
#'   \item `"gh"`: Generalized Hyperbolic (most flexible, includes NIG as 
#'   special case)
#' }
#'
#' **Example GOGARCH specification (Normal):**
#' ```
#' spec_gogarch <- list(
#'   list(
#'     var_order = 0,
#'     garch_spec_fun = "gogarch_modelspec",
#'     distribution = "norm",
#'     garch_spec_args = list(
#'       model = "garch",
#'       order = c(1, 1),
#'       ica = "radical",
#'       components = 2
#'     ),
#'     start_pars = list(
#'       var_pars = NULL,
#'       garch_pars = list(
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
#'       ),
#'       dist_pars = NULL
#'     )
#'   )
#' )
#' ```
#'
#' **Example GOGARCH specification (NIG distribution):**
#' ```
#' spec_gogarch_nig <- list(
#'   list(
#'     var_order = 0,
#'     garch_spec_fun = "gogarch_modelspec",
#'     distribution = "nig",
#'     garch_spec_args = list(
#'       model = "garch",
#'       order = c(1, 1),
#'       ica = "radical",
#'       components = 2
#'     ),
#'     start_pars = list(
#'       var_pars = NULL,
#'       garch_pars = list(
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
#'       ),
#'       dist_pars = list(skew = 0, shape = 1)
#'     )
#'   )
#' )
#' ```
#'
#' **Example GOGARCH specification (GH distribution):**
#' ```
#' spec_gogarch_gh <- list(
#'   list(
#'     var_order = 0,
#'     garch_spec_fun = "gogarch_modelspec",
#'     distribution = "gh",
#'     garch_spec_args = list(
#'       model = "garch",
#'       order = c(1, 1),
#'       ica = "radical",
#'       components = 2,
#'       lambda_range = c(-5, 5),
#'       shape_range = c(0.1, 25)
#'     ),
#'     start_pars = list(
#'       var_pars = NULL,
#'       garch_pars = list(
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8),
#'         list(omega = 0.1, alpha1 = 0.1, beta1 = 0.8)
#'       ),
#'       dist_pars = list(skew = 0, shape = 1, lambda = -0.5)
#'     )
#'   )
#' )
#' ```
#'
#' @seealso
#' \itemize{
#'   \item [gogarch_modelspec()] for creating GOGARCH specifications
#'   \item [estimate_garch_weighted_gogarch()] for GOGARCH weighted estimation
#'   \item [radical()] for the RADICAL ICA algorithm
#'   \item [compute_gogarch_loglik_ms()] for GOGARCH log-likelihood computation
#'   \item `vignette("cgarch_vs_dcc", package = "tsbs")` for model comparison
#' }
#' 
#' ### Diagnostics for MS VARMA GARCH 
#' Diagnostic System: A comprehensive diagnostic system monitors
#' convergence and numerical stability during estimation.
#' Additional parameters when \code{bs_type
#' = "ms_varma_garch"}
#' \itemize{
#'  \item \code{collect_diagnostics} Logical. Enable diagnostics with
#'    \code{collect_diagnostics = TRUE}. Diagnostics can only be extracted if
#'    \code{return_fit = TRUE}.
#'  \item \code{verbose} Logical. If TRUE, print detailed diagnostic information
#'    during estimation. Default is FALSE.
#'  \item \code{verbose_file} Character string specifying path to file for
#'    verbose output. If NULL (default), verbose output goes to console. If
#'    specified, all verbose output is written to this file instead. Only used if
#'    verbose = TRUE.
#' }
#' See details in \code{vignette("Diagnostics", package = "tsbs")}.
#'   
#' ## `taper_type` when `block_type="tapered"`
#' 
#' For block length \eqn{n}, and index \eqn{i} within the block, and \eqn{0 \leq i \leq n},
#' - `"cosine"`: Hann window 
#'     \deqn{w(i) = \dfrac{1}{2} \left(1 - \cos \left( \dfrac{2 \pi i}{n - i} \right ) \right)}
#' - `"bartlett"`: Triangular window  
#'     \deqn{w(i) = 1 - \left | \dfrac{i - (n - 1) / 2}{(n - 1) / 2} \right |,\ \ \ \alpha \in [0, 1]}
#' - `"tukey"`: Cosine tapered window. At \eqn{\alpha = 0} it becomes 
#'     rectangular, and at \eqn{\alpha = 1} it becomes a Hann window. 
#'     
#'     We implement the Tukey (tapered cosine) window using the convention 
#'     \eqn{i = 0, \dots, n-1} for a window of length \eqn{n} and taper 
#'     parameter \eqn{\alpha \in [0,1]}, consistent with common numerical 
#'     libraries such as MATLAB and SciPy (see MATLAB tukeywin, 
#'     SciPy scipy.signal.windows.tukey). After index shift and reflection this 
#'     is algebraically equivalent to the classical definition in Harris (1978), 
#'     which instead centers the index symmetrically around zero.  
#'     
#'     Let \eqn{N = n - 1} be the last index of the block, and let the running 
#'     index be \eqn{i = 0, 1, \dots, N}.
#'
#'     The window weights \eqn{w[i]} are defined piecewise as follows:
#'
#'     1. Left taper region (\eqn{0 \le i < \alpha N / 2}): 
#'     \eqn{
#'       w[i] = \frac{1}{2} \left[ 1 + \cos\left( \pi \left( \frac{2 i}{\alpha N} - 1 \right) \right) \right]
#'     }
#'
#'     2. Central (non-tapered) region (\eqn{\alpha N / 2 \le i \le N (1 - \alpha/2)}):
#'     \eqn{
#'       w[i] = 1
#'     }
#'
#'     3. Right taper region (\eqn{N (1 - \alpha/2) < i \le N}):
#'     \eqn{
#'       w[i] = \frac{1}{2} \left[ 1 + \cos\left( \pi \left( \frac{2 i}{\alpha N} - \frac{2}{\alpha} + 1 \right) \right) \right]
#'     }  
#'     
#' (Harris, 1978)
#'
#' @return A list containing:
#' \describe{
#'   \item{bootstrap_series}{List of bootstrap replicate matrices.}
#'   \item{func_outs}{List of computed function outputs for each replicate.}
#'   \item{func_out_means}{Mean of the computed outputs across replicates.}
#'   \item{fit}{If `return_fit = TRUE`, a list containing model fit.}
#'   \item{diagnostics}{If `return_diagnostics = TRUE`, a `tsbs_diagnostics` 
#'     object containing block composition, statistics, and quality metrics. 
#'     See \code{\link{compute_bootstrap_diagnostics}}.}
#' }
#' 
#' @references
#' Harris, F. J. (1978). "On the Use of Windows for 
#'   Harmonic Analysis with the Discrete Fourier Transform." (Proceedings of the 
#'   IEEE, 66(1), pp. 66f)
#'   
#' @examples
#' set.seed(123)
#' x <- arima.sim(n = 100, list(ar = 0.8))
#' result <- tsbs(
#'   x = as.matrix(x),
#'   block_length = 10,
#'   bs_type = "stationary",
#'   num_blocks = 5,
#'   num_boots = 10,
#'   func = mean,
#'   apply_func_to = "cols"
#' )
#' print(result$func_out_means)
#' 
#' @useDynLib tsbs, .registration = TRUE
#' @export
#' @md
tsbs <- function(
    x,
    n_boot = NULL,
    block_length = NULL,
    bs_type = c("moving", "stationary", "hmm", "msvar", "ms_varma_garch", "wild"),
    block_type = c("overlapping", "non-overlapping", "tapered"),
    taper_type = c("cosine", "bartlett", "tukey"),
    tukey_alpha = 0.5, ## Only applies to block_type="tapered" with 
                       ## taper_type="tukey"
    num_blocks = NULL,
    num_boots = 100L,
    func = mean,
    apply_func_to = c("cols", "all"),
    p_method = c("1/n", "plugin", "cv"),
    p = NULL,
    overlap = TRUE,
    num_states = 2L,
    d = 0,
    spec = NULL,
    model_type = c("univariate", "multivariate"),
    control = list(),
    parallel = FALSE,
    num_cores = 1L,
    model_func = default_model_func,
    score_func = mse,
    stationary_max_percentile = 0.99,
    stationary_max_fraction_of_n = 0.10,
    return_fit = FALSE,
    return_diagnostics = FALSE,
    fail_mode = c("predictably", "gracefully"),
    ...
    #collect_diagnostics = FALSE,
    #verbose = FALSE,
    #verbose_file = NULL
) {
  
  ## ---- Validation ----
  ## match.arg() picks the first element in the vector as default
  bs_type <- match.arg(bs_type)
  block_type <- match.arg(block_type)
  taper_type <- match.arg(taper_type)
  apply_func_to <- match.arg(apply_func_to)
  p_method <- match.arg(p_method)
  model_type <- match.arg(model_type)
  fail_mode <- match.arg(fail_mode)
  
  # if(.is_invalid_data(x, fail_mode = fail_mode))
  #   stop("No valid x value provided.")
  # if (is.vector(x) || is.data.frame(x) || is.ts(x)) x <- as.matrix(x)
  
  x <- .coerce_and_validate_data(x)
  
  n <- nrow(x)
  
  if (!is.function(func))
    stop("`func` must be a valid function.")
  
  ## Note: If NULL, value is calculated automatically below
  if (.is_invalid_count(n_boot, fail_mode = fail_mode))
    stop("`n_boot` must be a positive integer or NULL.")
  if (.is_invalid_count(block_length, fail_mode = fail_mode))
    stop("`block_length` must be a positive integer or NULL.")
  if (.is_invalid_count(num_blocks, fail_mode = fail_mode))
    stop("`num_blocks` must be a positive integer or NULL.")
  
  ## Need to provide either n_boot or num_blocks.
  ## If missing, block_length will be calculated automatically.
  ## n_boot = block_length * num_blocks.
  if(.is_invalid_count(n_boot, fail_mode = fail_mode) && 
     .is_invalid_count(num_blocks, fail_mode = fail_mode)) {
    stop("Must provide either n_boot or num_blocks.")
  } 
  
  ## Commented out section:
  ## This is done in blockBootstrap.cpp and does not apply to bootstrap methods
  ## which rely on .sample_blocks() or .sample_blocks_with diagnostics().
  # if(.is_invalid_count(n_boot, fail_mode = fail_mode)) {
  #   if (is.null(block_length)) {block_length <- compute_default_block_length(x)}
  #   n_boot <- num_blocks * block_length
  # }
  # if(.is_invalid_count(num_blocks, fail_mode = fail_mode)) {
  #   if (is.null(block_length)) {block_length <- compute_default_block_length(x)}
  #   num_blocks <- n_boot / block_length
  # }
  
  
  ## Fails if NULL. Value is not calculated automatically.
  if(.is_invalid_count(num_states, allow_null = FALSE, fail_mode = fail_mode))
    stop("`num_states` must be a positive integer.")
  if(.is_invalid_count(num_boots, allow_null = FALSE, fail_mode = fail_mode))
    stop("`num_boots` must be a positive integer.")
  
  ## Validate p
  if(.not_in_unit_interval(
      p, 
      allow_null = TRUE, 
      interval_type = "open", ## open unit interval, (0, 1)
      fail_mode = fail_mode
    ) 
  ) {
    stop("`p` must be a single number in (0,1) or NULL.")
  }
  
  ## Validate tukey_alpha
  if(block_type == "tapered" && taper_type == "tukey") {
    if(
      .not_in_unit_interval(
        tukey_alpha, 
        allow_null = FALSE, 
        interval_type = "closed", ## closed unit interval, [0, 1]
        fail_mode = fail_mode
      )
    ) {
      stop("`tukey_alpha` must be a single number in [0,1] or NULL.")
    }
  }

  ## ---- Bootstrap ----
  
  .blockBootstrap <- function() {
    blockBootstrap(
      x, 
      n_boot, 
      block_length, 
      bs_type,
      block_type,
      taper_type,
      tukey_alpha,
      num_blocks, 
      num_boots, 
      p, 
      stationary_max_percentile,
      stationary_max_fraction_of_n
    )
  }
  
  bootstrap_series <- switch(
    bs_type,
    moving = {
      .blockBootstrap()
    },
    stationary = {
      if(is.null(p)) {p <- .estimate_p(x, p_method, block_length, model_func, score_func)}
      .blockBootstrap()
    },
    hmm = hmm_bootstrap(
      x, 
      n_boot = n_boot,
      num_states = num_states, 
      num_blocks = num_blocks, 
      num_boots = num_boots, 
      parallel = parallel, 
      num_cores = num_cores
    ),
    msvar = msvar_bootstrap(
      x, 
      n_boot = n_boot,
      num_blocks = num_blocks, 
      num_boots = num_boots, 
      parallel = parallel, 
      num_cores = num_cores
    ),
    ms_varma_garch = ms_varma_garch_bs(
      x = x,
      n_boot = n_boot,
      num_blocks = num_blocks,
      num_boots = num_boots,
      M = num_states,
      d = d,
      spec = spec,
      model_type = model_type,
      return_fit = return_fit,
      control = control,
      parallel = parallel,
      num_cores = num_cores,
      ...
      # collect_diagnostics = collect_diagnostics,
      # verbose = verbose,
      # verbose_file = verbose_file
    ),
    wild = wild_bootstrap(
      x, 
      num_boots, 
      parallel = parallel, 
      num_cores = num_cores
    ),
    stop("Unsupported bootstrap type.")
  )
  
  ## Unpack if ms_varma_garch with return_fit
  fit_object <- NULL
  if (bs_type == "ms_varma_garch" && return_fit) {
    fit_object <- bootstrap_series$fit
    bootstrap_series <- bootstrap_series$bootstrap_series
  }
  
  func_outs <- lapply(bootstrap_series, function(sampled) {
    if (apply_func_to == "all") {
      #func(as.data.frame(sampled))
      func(sampled)
    } else {
      apply(sampled, 2, func)
    }
  })
  
  #func_out_means <- Reduce(`+`, func_outs) / length(func_outs)
  if (is.list(func_outs[[1]]) && !is.data.frame(func_outs[[1]])) {
    ## Average each element of the list across bootstrap replicates
    func_out_means <- lapply(names(func_outs[[1]]), function(nm) {
      vals <- lapply(func_outs, `[[`, nm)
      tryCatch(Reduce(`+`, vals) / length(vals), error = function(e) NULL)
    })
    names(func_out_means) <- names(func_outs[[1]])
  } else {
    func_out_means <- Reduce(`+`, func_outs) / length(func_outs)
  }
  
  ## ---- Diagnostics Collection ----
  diagnostics_object <- NULL
  if (return_diagnostics) {
    
    # Build configuration list
    config <- list(
      bs_type = bs_type,
      block_type = block_type,
      taper_type = if (block_type == "tapered") taper_type else NA,
      tukey_alpha = if (block_type == "tapered" && taper_type == "tukey") tukey_alpha else NA,
      block_length = block_length,
      n_boot = n_boot,
      num_blocks = num_blocks,
      num_boots = num_boots,
      p = p,
      p_method = if (bs_type == "stationary") p_method else NA,
      num_states = if (bs_type %in% c("hmm", "msvar", "ms_varma_garch")) num_states else NA
    )
    
    # Compute diagnostics from bootstrap series
    diagnostics_object <- compute_bootstrap_diagnostics(
      bootstrap_series = bootstrap_series,
      original_data = x,
      bs_type = bs_type,
      config = config
    )
    
    # Add method-specific diagnostics
    if (bs_type == "stationary") {
      diagnostics_object <- record_method_specific(
        diagnostics_object,
        p_estimated = p,
        expected_block_length = if (!is.null(p)) 1/p else NA
      )
    }
    
    if (bs_type %in% c("hmm", "msvar", "ms_varma_garch") && !is.null(fit_object)) {
      # Extract state sequence if available
      if (!is.null(fit_object$smoothed_probabilities)) {
        state_seq <- apply(fit_object$smoothed_probabilities, 1, which.max)
        state_counts <- table(state_seq)
        diagnostics_object <- record_method_specific(
          diagnostics_object,
          state_counts = as.list(state_counts),
          transition_matrix = if (!is.null(fit_object$P)) fit_object$P else NA
        )
      }
    }
  }

  ## Build result
  result <- list(
    bootstrap_series = bootstrap_series,
    func_outs = func_outs,
    func_out_means = func_out_means
  )
  
  if (return_fit && !is.null(fit_object)) {
    result$fit <- fit_object
  }
  
  if (return_diagnostics && !is.null(diagnostics_object)) {
    result$diagnostics <- diagnostics_object
  }
  
  result
}


#' Check if expression or variable exists and is evaluable
#'
#' DEPRECATED. Only called by .is_invalid_data(), which is DEPRECATED.
#'
#' Helper function that distinguishes between: 1. Expressions that evaluate
#' successfully (e.g., numeric(0), 1+1) 2. Undefined variables (e.g.,
#' nonexistent_var) 3. Expressions that fail to evaluate for other reasons
#'
#' @param x The expression/variable to check (as substitute() result)
#' @param parent_env The parent environment to check for variable existence
#'
#' @returns TRUE if the expression/variable is invalid/non-existent, FALSE
#'   otherwise
.check_expression_validity <- function(expr_sub, parent_env = parent.frame()) {
  
  ## First, try to evaluate the expression in the parent environment
  tryCatch({
    ## Try to evaluate the substituted expression
    eval(expr_sub, envir = parent_env)
    return(FALSE)  ## Expression evaluated successfully
  }, error = function(e) {
    ## Evaluation failed - check if it's a simple undefined variable
    var_name <- deparse(expr_sub)
    
    ## Check if this looks like a simple variable name (not a complex expression)
    is_simple_name <- grepl("^[a-zA-Z][a-zA-Z0-9_.]*$", var_name)
    
    if (is_simple_name && !exists(var_name, where = parent_env)) {
      return(TRUE)  ## Simple variable name that doesn't exist
    } else {
      ## Either it's a complex expression that failed, or a variable that exists 
      ## but evaluation failed for other reasons. Re-throw the original error.
      stop(e)
    }
  })
}


#' @title Coerce, Validate, and Prepare Time Series Data
#' @description A robust helper function that takes a user-provided object,
#'   validates it, and coerces it into a clean, numeric matrix suitable for
#'   downstream analysis. It provides informative error messages or fails
#'   gracefully based on the specified mode.
#'
#' @param x The user-provided data (e.g., vector, matrix, data.frame, ts, xts,
#'   zoo).
#' @param fail_mode How to handle validation errors: "predictably" (the default)
#'   will stop with an informative error. "gracefully" will return NULL.
#' @return A numeric matrix if the input is valid, otherwise stops or returns
#'   NULL.
#' @keywords internal
.coerce_and_validate_data <- function(x, fail_mode = c("predictably", "gracefully")) {
  fail_mode <- match.arg(fail_mode)
  x_name <- deparse(substitute(x))
  
  ## --- Internal handler for failing predictably or gracefully ---
  handle_failure <- function(message) {
    if (fail_mode == "predictably") {
      stop(message, call. = FALSE)
    } else {
      return(NULL)
    }
  }
  
  ## --- 1. Initial NULL check ---
  if (is.null(x)) {
    return(handle_failure(paste0("Input '", x_name, "' is NULL.")))
  }
  
  ## --- 2. Coercion based on type ---
  data_matrix <- NULL
  tryCatch({
    if (inherits(x, c("xts", "zoo"))) {
      data_matrix <- zoo::coredata(x)
    } else if (is.matrix(x)) {
      data_matrix <- x
    } else if (is.data.frame(x)) {
      data_matrix <- as.matrix(x)
    } else if (is.ts(x)) {
      data_matrix <- as.matrix(x)
    } else if (is.vector(x) && !is.list(x)) { ## Explicitly exclude lists.
                                              ## A vector in R is either an 
                                              ## atomic vector i.e., one of the 
                                              ## atomic types, or of type 
                                              ## or mode list or expression.
      data_matrix <- as.matrix(x)
    } else {
      ## This block now correctly catches lists and other unsupported types
      return(handle_failure(paste0("Unsupported data type for '", x_name, "': '", class(x)[1], "'.")))
    }
  }, error = function(e) {
    ## This catch is a final safety net for unexpected coercion errors
    data_matrix <<- handle_failure(paste0("Failed to coerce '", x_name, "' to a matrix: ", e$message))
  })
  
  ## If coercion failed gracefully, data_matrix will be NULL.
  if (is.null(data_matrix)) return(NULL)
  
  
  ## --- 3. Universal Validation on the resulting matrix ---
  if (NROW(data_matrix) < 1 || NCOL(data_matrix) < 1) {
    return(handle_failure(paste0("Input '", x_name, "' must not be empty.")))
  }
  if (!is.numeric(data_matrix)) {
    return(handle_failure(paste0("Input '", x_name, "' must be numeric.")))
  }
  if (any(!is.finite(data_matrix))) {
    return(handle_failure(paste0("Input '", x_name, "' contains non-finite values (NA, NaN, Inf).")))
  }

  
  ## --- 4. Success ---
  return(data_matrix)
}


#' Invalid input data
#'
#' DEPRECATED
#'
#' Returns TRUE if input data is not valid, FALSE otherwise.
#' 
#' @param x x, vector, matrix, data.frame or time series
#' @param allow_null NULL value may be allowed when functions calculates value
#'   of x automatically.
#' @param fail_mode How to handle validation errors: "predictably" (fail fast) 
#'   or "gracefully" (return FALSE on error)
#'
#' @returns boolean
#' @export
.is_invalid_data <- function(x, allow_null = TRUE, fail_mode = c("predictably", "gracefully")) {
  fail_mode <- match.arg(fail_mode)
  
  ## Check if the expression/variable itself is valid and evaluable
  if (.check_expression_validity(substitute(x), parent.frame())) {
    return(TRUE)
  }
  
  ## Handle NULL values based on the allow_null parameter
  if (is.null(x)) {
    return(!allow_null)
  }
  
  ## --- Unified Validation Logic with fail_mode ---
  
  ## Define the core coercion logic in a helper function
  perform_coercion <- function(input) {
    if (inherits(input, c("xts", "zoo"))) {
      return(zoo::coredata(input))
    } else if (inherits(input, c("ts", "data.frame", "matrix"))) {
      return(as.matrix(input))
    } else if (is.vector(input)) {
      return(as.matrix(input))
    } else {
      ## Return NULL for unsupported types, which signals an invalid state
      return(NULL)
    }
  }
  
  data_to_check <- NULL
  if (fail_mode == "gracefully") {
    ## Graceful mode: Catch any errors during coercion and treat as invalid
    data_to_check <- tryCatch(perform_coercion(x), error = function(e) NULL)
  } else {
    ## Predictable mode: Let any coercion errors bubble up
    data_to_check <- perform_coercion(x)
  }
  
  ## If coercion failed or the type was unsupported, the data is invalid.
  if (is.null(data_to_check)) {
    return(TRUE)
  }
  
  ## Now, perform a single, universal check on the resulting matrix.
  check_expression <- {
    !is.numeric(data_to_check) ||
      any(is.na(data_to_check)) ||
      any(!is.finite(data_to_check))
  }
  
  is_invalid <- FALSE
  if (fail_mode == "gracefully") {
    ## In graceful mode, if the check itself throws an error,
    ## we treat the data as invalid.
    is_invalid <- tryCatch(check_expression, error = function(e) TRUE)
  } else {
    ## In predictable mode, let any error from the check bubble up.
    is_invalid <- check_expression
  }
  
  if (is_invalid) {
    return(TRUE)
  }
  
  ## If we reach here, all checks passed
  return(FALSE)
}

# .is_invalid_data <- function(x, allow_null = TRUE, fail_mode = c("predictably", "gracefully")) {
#   
#   fail_mode <- match.arg(fail_mode)
#   
#   ## Check if expression/variable is valid and evaluable
#   if (.check_expression_validity(substitute(x), parent.frame())) {
#     return(TRUE)
#   }
#   
#   ## Handle NULL values based on allow_null parameter
#   if (is.null(x)) {
#     return(!allow_null)
#   }
#   
#   ## Helper function to handle errors based on fail_mode
#   safe_check <- function(expr) {
#     if (fail_mode == "predictably") {
#       ## Let errors bubble up
#       expr
#     } else {
#       ## Handle errors gracefully
#       tryCatch(expr, error = function(e) FALSE)
#     }
#   }
#   
#   ## Check if x is numeric (for vectors, matrices, or data frames)
#   if (is.vector(x) || is.matrix(x)) {
#     if (!is.numeric(x)) return(TRUE)
#   } else if (is.data.frame(x)) {
#     ## Check if all columns are numeric and data frame is not empty
#     if (!all(sapply(x, is.numeric)) || nrow(x) < 1 || ncol(x) < 1) {
#       return(TRUE)
#     }
#     
#     ## Check for NA values in data frame
#     if (safe_check(any(is.na(x)))) return(TRUE)
#     
#     ## Check for non-finite values in data frame
#     if (safe_check(any(!sapply(x, function(col) all(is.finite(col)))))) return(TRUE)
#   } else {
#     ## For other object types, consider them invalid
#     return(TRUE)
#   }
#   
#   ## Handle remaining checks based on object type
#   if (is.data.frame(x)) {
#     ## Data frame specific checks (already handled NA and finite above)
#     if (safe_check(any(sapply(x, is.null)))) return(TRUE)
#   } else {
#     ## Vector and matrix checks - ultra-consolidated
#     if (safe_check({
#       any(is.null(x)) || any(is.na(x)) || (is.numeric(x) && any(!is.finite(x)))
#     })) return(TRUE)
#   }
#   
#   ## If we reach here, all checks passed
#   return(FALSE)
# }



#' Invalid input count
#'
#' Returns TRUE if input count is not valid, FALSE otherwise.
#' 
#' @param n Count value to validate
#' @param allow_null NULL value may be allowed when functions calculates value
#'   of n automatically.
#' @param fail_mode How to handle validation errors: "predictably" (fail fast) 
#'   or "gracefully" (return FALSE on error)
#'
#' @returns boolean
#' @export
.is_invalid_count <- function(n, allow_null = TRUE, fail_mode = c("predictably", "gracefully")) {
  
  fail_mode <- match.arg(fail_mode)
  
  ## Check if expression/variable is valid and evaluable
  if (.check_expression_validity(substitute(n), parent.frame())) {
    return(TRUE)
  }
  
  ## Handle NULL values based on allow_null parameter
  if (is.null(n)) {
    return(!allow_null)
  }
  
  ## Helper function to handle errors based on fail_mode
  safe_check <- function(expr) {
    if (fail_mode == "predictably") {
      expr
    } else {
      tryCatch(expr, error = function(e) FALSE)
    }
  }
  
  ## Check if n is numeric and has length 1
  if (!is.numeric(n) || length(n) != 1) {
    return(TRUE)
  }
  
  ## NA, non-finite, non-integer, non-positive
  if (safe_check({
    is.na(n) || !is.finite(n) || n <= 0 || 
      ## Use trunc(n) rather than as.integer(n) which may return NA due to
      ## integer overflow.
      ## trunc(n) never fails - it just truncates to the integer part.
      ## .Machine$integer.max is the largest representable integer 
      ## (usually 2,147,483,647)
      n != trunc(n) || n > .Machine$integer.max
  })) {
    return(TRUE)
  }
  
  ## If we reach here, all checks passed
  return(FALSE)
}


#' Value not in unit interval
#'
#' Returns TRUE if `p` is not in the unit interval, FALSE otherwise.
#' 
#' If 
#' 
#' @param p Value to validate, e.g. a percentage (decimal fraction) .
#' @param allow_null boolean, allow NULL?
#' @param interval_type One of `"open"` or `"closed"`. Indicates if the valid
#'   interval is an open or closed unit interval, i.e. \eqn{p \in (0, 1)} or
#'   \eqn{p \in [0, 1]} respectively.
#' @param fail_mode How to handle validation errors: "predictably" (fail fast) 
#'   or "gracefully" (return FALSE on error)
#'
#' @returns boolean
#' @export
.not_in_unit_interval <- function(
    p, 
    allow_null = TRUE, 
    interval_type = c("open", "closed"),
    fail_mode = c("predictably", "gracefully")
  ) {
  
  interval_type <- match.arg(interval_type)
  fail_mode <- match.arg(fail_mode)
  
  ## Check if expression/variable is valid and evaluable
  if (.check_expression_validity(substitute(p), parent.frame())) {
    return(TRUE)
  }
  
  ## Handle NULL values based on allow_null parameter
  if (is.null(p)) {
    return(!allow_null)
  }
  
  ## Helper function to handle errors based on fail_mode
  safe_check <- function(expr) {
    if (fail_mode == "predictably") {
      expr
    } else {
      tryCatch(expr, error = function(e) FALSE)
    }
  }
  
  ## Check if p is numeric and has length 1
  if (!is.numeric(p) || length(p) != 1) {
    return(TRUE)
  }
  
  ## NA, non-finite, out of valid range (0,1)
  invalid_expr <- switch(
    interval_type,
    open   = is.na(p) || !is.finite(p) || p <= 0 || p >= 1,
    closed = is.na(p) || !is.finite(p) || p < 0  || p > 1
  )
  
  if (safe_check(invalid_expr)) {
    return(TRUE)
  }
  
  ## If we reach here, all checks passed
  return(FALSE)
}


#' Estimate geometric distribution parameter p for stationary bootstrap
#'
#' @param x x
#' @param p_method p estimation method
#' @param block_length Expected block length
#' @param model_func Model function for k fold cross validation
#' @param score_func Score function for k fold cross validation
#'
#' @returns p
.estimate_p <- function(x, p_method, block_length, model_func, score_func) {
  if (p_method == "1/n") {
    p <- 1 / if (is.null(block_length)) compute_default_block_length(x) else block_length
  } else if (p_method == "plugin") {
    ## acf[2,,1] is the autocorrelation for the 1st order (index 2) the first 
    ## variable (col 1).
    #ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,1]
    
    ## acf[2,,] is the vector of autocorrelations for the 1st order (index 2) 
    ## of each variable (if multivariate).
    ac1 <- acf(x, lag.max = 1, plot = FALSE)$acf[2,,]
    
    ## The larger the average 1st order autocorrelation, the smaller the p
    p <- 1 - abs(mean(ac1, na.rm = TRUE))
  } else if (p_method == "cv") {
    if (is.null(model_func) || is.null(score_func))
      stop("For p_method = 'cross validation', provide `model_func` and `score_func`.")
    p <- k_fold_cv_ts(x, x[,1], model_func, score_func)
  }
}