## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
## CGARCH (Copula GARCH) Documentation for tsbs Package
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
##
## This file contains roxygen2 documentation to be incorporated into the tsbs
## package. The documentation covers:
##   1. Updates to tsbs() main help page (additions to existing documentation)
##   2. Internal function documentation (linked via @seealso)
##   3. A summary vignette section for cgarch vs dcc
##
## = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =


## ---------------------------------------------------------------------------
## SECTION 1: Additions to tsbs() Documentation
## ---------------------------------------------------------------------------
##
## The following should be ADDED to the existing tsbs.R documentation in the
## appropriate sections.

## In the @details section, after the DCC multivariate documentation, ADD:
##
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
#'   \item [cgarch_modelspec()] for creating CGARCH specifications
#'   \item [estimate_garch_weighted_cgarch()] for the weighted estimation method
#'   \item [compute_pit_transform()] for PIT transformation details
#'   \item [adcc_recursion()] for asymmetric DCC dynamics
#'   \item `vignette("cgarch_vs_dcc", package = "tsbs")` for model comparison
#' }


## ---------------------------------------------------------------------------
## SECTION 2: Internal Function Documentation
## ---------------------------------------------------------------------------

#' @title Estimate Copula GARCH Parameters with Weighted Likelihood
#' @name estimate_garch_weighted_cgarch
#'
#' @description
#' Two-stage weighted maximum likelihood estimation for Copula GARCH models
#' within the MS-VARMA-GARCH framework. This is the main workhorse function
#' called during the EM algorithm's M-step for multivariate models using
#' copula-based dependence structures.
#'
#' @details
#' ## Two-Stage Estimation
#'
#' Stage 1 estimates univariate GARCH parameters for each series using
#' weighted quasi-maximum likelihood. Stage 2 estimates the copula dependence
#' parameters (DCC/ADCC dynamics and distribution shape) given the Stage 1
#' volatilities.
#'
#' ## Estimation Flow
#'
#' 1. **Stage 1 (Marginals)**: For each series i = 1, ..., k:
#'    - Fit univariate GARCH(p,q) with weighted likelihood

#'    - Extract conditional volatility σ_{i,t}
#'    - Compute standardized residuals z_{i,t} = ε_{i,t} / σ_{i,t}
#'
#' 2. **PIT Transformation**: Transform z to uniform margins U:
#'    - Parametric: U_{i,t} = F_i(z_{i,t}; θ_i)
#'    - Empirical: U_{i,t} = F̂_n(z_{i,t})
#'    - SPD: Combines kernel density center with parametric tails
#'
#' 3. **Copula Residuals**: Transform U to copula space:
#'    - MVN copula: Z_{i,t} = Φ^{-1}(U_{i,t})
#'    - MVT copula: Z_{i,t} = t_ν^{-1}(U_{i,t}) / √(ν/(ν-2))
#'
#' 4. **Stage 2 (Dependence)**: Optimize copula parameters:
#'    - DCC(1,1): Uses analytical gradient with reparameterization
#'    - DCC(p,q): Uses numerical gradient
#'    - ADCC: Uses numerical gradient with asymmetry parameter
#'
#' ## Reparameterization for DCC(1,1)
#'
#' For numerical stability, DCC(1,1) parameters are transformed:
#' - ψ = logit(α + β) (persistence in logit space)
#' - φ = log(α / β) (ratio in log space)
#'
#' This maps the constrained region {α > 0, β > 0, α + β < 1} to
#' the unconstrained R² space.
#'
#' @param residuals Matrix of residuals from conditional mean model (T × k)
#' @param weights Vector of observation weights from E-step (length T)
#' @param spec Model specification list containing:
#'   \describe{
#'     \item{garch_spec_args}{List with transformation, copula, dynamics, etc.}
#'     \item{start_pars}{Starting parameter values}
#'   }
#' @param state Integer state index (for diagnostics)
#' @param diagnostics Optional diagnostics collector object
#' @param iteration Current EM iteration (for diagnostics)
#' @param verbose Logical; print progress information
#'
#' @return A list containing:
#' \describe{
#'   \item{garch_pars}{List of k lists, each with univariate GARCH parameters}
#'   \item{dcc_pars}{Named list with DCC/ADCC parameters (alpha_1, beta_1,
#'     optionally gamma_1 for ADCC)}
#'   \item{dist_pars}{Named list with distribution parameters (shape for MVT)}
#'   \item{weighted_ll}{Scalar weighted log-likelihood}
#'   \item{warnings}{List of any warnings generated during estimation}
#' }
#'
#' @seealso
#' [estimate_copula_parameters_weighted()] for Stage 2 estimation details
#' [compute_pit_transform()] for PIT transformation
#' [copula_nll()] for copula likelihood computation
#'
#' @keywords internal
NULL


#' @title Estimate Copula Parameters with Weighted Likelihood
#' @name estimate_copula_parameters_weighted
#'
#' @description
#' Stage 2 of copula GARCH estimation: estimates the copula dependence
#' parameters given Stage 1 univariate GARCH fits.
#'
#' @details
#' This function handles:
#' - PIT transformation of standardized residuals
#' - Computation of copula residuals
#' - Optimization of DCC/ADCC parameters
#' - Distribution shape estimation for MVT copula
#'
#' The optimization strategy depends on the model:
#' - **DCC(1,1)**: Analytical gradient with (ψ, φ) reparameterization
#' - **DCC(p,q)**: Numerical gradient with box constraints
#' - **ADCC**: Numerical gradient with (α, γ, β) parameterization
#'
#' @param residuals Matrix of raw residuals (T × k)
#' @param weights Observation weights (length T)
#' @param garch_pars Stage 1 univariate GARCH parameter estimates
#' @param dcc_start_pars Starting values for DCC parameters
#' @param dist_start_pars Starting values for distribution parameters
#' @param spec Full model specification
#' @param transformation PIT transformation type
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param dynamics Correlation dynamics ("constant", "dcc", or "adcc")
#' @param verbose Print progress
#'
#' @return List with dcc_pars, dist_pars, weighted_ll, warnings
#'
#' @seealso [copula_gradient()] for analytical gradients
#' @keywords internal
NULL


#' @title Probability Integral Transform for Copula Models
#' @name compute_pit_transform
#'
#' @description
#' Transforms standardized residuals to uniform margins using various methods.
#' This is a critical step in copula models that separates the marginal
#' distributions from the dependence structure.
#'
#' @details
#' ## Transformation Methods
#'
#' **Parametric (`transformation = "parametric"`)**
#'
#' Uses the fitted univariate distribution CDF:
#' U_{i,t} = F_i(z_{i,t}; θ_i)
#'
#' where F_i is the CDF implied by the univariate GARCH model's distribution
#' (normal, Student-t, GED, etc.) and θ_i are the estimated distribution
#' parameters.
#'
#' **Empirical (`transformation = "empirical"`)**
#'
#' Uses the empirical CDF:
#' U_{i,t} = (1/T) Σ_{s=1}^T I(z_{i,s} ≤ z_{i,t})
#'
#' Scaled by T/(T+1) to avoid exact 0 and 1 values.
#'
#' **Semi-Parametric Distribution (`transformation = "spd"`)**
#'
#' Combines kernel density estimation in the interior with parametric tails:
#' - Interior (10th-90th percentile): Kernel density with Sheather-Jones
#'   bandwidth
#' - Tails: Empirical quantile-based scaling
#'
#' This provides robustness to misspecification while maintaining smooth
#' tail behavior.
#'
#' @param std_residuals Matrix of standardized residuals (T × k)
#' @param uni_fit_list List of univariate GARCH fit objects
#' @param transformation Character: "parametric", "empirical", or "spd"
#' @param copula_dist Copula distribution (for scaling if needed)
#' @param dist_pars Distribution parameters (shape for MVT)
#'
#' @return Matrix of uniform values (T × k), bounded in (0, 1)
#'
#' @seealso [fit_spd_transform()] for SPD details
#' @keywords internal
NULL


#' @title Compute Copula Residuals from Uniform Margins
#' @name compute_copula_residuals
#'
#' @description
#' Transforms uniform margins to copula-space residuals using the inverse
#' copula marginal CDF.
#'
#' @details
#' For Gaussian copula (MVN):
#' Z_{i,t} = Φ^{-1}(U_{i,t})
#'
#' For Student-t copula (MVT):
#' Z_{i,t} = t_ν^{-1}(U_{i,t}) / √(ν/(ν-2))
#'
#' The scaling for MVT ensures unit variance.
#'
#' @param u_matrix Matrix of uniform values (T × k)
#' @param copula_dist Copula distribution ("mvn" or "mvt")
#' @param dist_pars Distribution parameters (shape for MVT)
#'
#' @return Matrix of copula residuals (T × k)
#'
#' @keywords internal
NULL


#' @title Copula Negative Log-Likelihood for DCC(1,1)
#' @name copula_nll
#'
#' @description
#' Computes the weighted copula negative log-likelihood for DCC(1,1) models.
#'
#' @details
#' ## Gaussian Copula (MVN)
#'
#' The copula log-likelihood contribution at time t is:
#' ℓ_t = -1/2 (log|R_t| + z_t' R_t^{-1} z_t - z_t' z_t)
#'
#' where R_t is the dynamic correlation matrix from DCC recursion.
#'
#' ## Student-t Copula (MVT)
#'
#' ℓ_t = log Γ((k+ν)/2) - log Γ(ν/2) - (k/2) log(π(ν-2))
#'       - (1/2) log|R_t|
#'       - ((ν+k)/2) log(1 + q_t/(ν-2))
#'       - Σ_i log f_t(z_{i,t}; ν)
#'
#' where q_t = z_t' R_t^{-1} z_t is the Mahalanobis distance and f_t is the
#' univariate Student-t density.
#'
#' ## Reparameterization
#'
#' When `use_reparam = TRUE`, the function expects parameters in the
#' unconstrained (ψ, φ) space rather than (α, β).
#'
#' @param params Parameter vector
#' @param z_matrix Copula residuals (T × k)
#' @param weights Observation weights
#' @param Qbar Unconditional covariance matrix
#' @param copula_dist "mvn" or "mvt"
#' @param use_reparam Use reparameterized space
#'
#' @return Scalar negative log-likelihood
#'
#' @seealso [copula_gradient()] for the gradient
#' @keywords internal
NULL


#' @title Analytical Copula Gradient for DCC(1,1)
#' @name copula_gradient
#'
#' @description
#' Computes the analytical gradient of the copula negative log-likelihood
#' with respect to DCC parameters.
#'
#' @details
#' The gradient is computed using the chain rule through:
#' 1. Forward pass: DCC recursion with gradient storage
#' 2. Backward pass: Backpropagate through R → Q normalization
#' 3. Chain rule: Transform gradient to (ψ, φ) space if reparameterized
#'
#' For MVT copula, the shape gradient is computed numerically since the
#' analytical form involves digamma functions.
#'
#' @param params Parameter vector
#' @param z_matrix Copula residuals
#' @param weights Observation weights
#' @param Qbar Unconditional covariance
#' @param copula_dist Copula distribution
#' @param use_reparam Use reparameterized space
#'
#' @return Gradient vector (same length as params)
#'
#' @keywords internal
NULL


#' @title ADCC (Asymmetric DCC) Recursion
#' @name adcc_recursion
#'
#' @description
#' Computes Q and R matrices for ADCC(1,1) model, which captures asymmetric
#' response to positive vs. negative shocks.
#'
#' @details
#' The ADCC recursion is:
#' Q_t = Ω + α(z_{t-1} z_{t-1}') + γ(n_{t-1} n_{t-1}') + β Q_{t-1}
#'
#' where:
#' - n_t = z_t ⊙ I(z_t < 0) (element-wise negative indicator)
#' - Ω = (1 - α - β) Q̄ - γ N̄
#' - N̄ is the unconditional covariance of negative shocks
#'
#' The correlation matrix is:
#' R_t = D_t^{-1} Q_t D_t^{-1}
#'
#' where D_t = diag(√diag(Q_t)).
#'
#' ## Stationarity Constraint
#'
#' For ADCC stationarity: α + β + δγ < 1
#' where δ ≈ 0.5 (data-dependent scaling factor).
#'
#' @param std_resid Standardized residuals (T × k)
#' @param Qbar Unconditional covariance matrix
#' @param alpha DCC alpha (shock response)
#' @param gamma ADCC gamma (asymmetry parameter)
#' @param beta DCC beta (persistence)
#' @param Nbar Unconditional covariance of negative shocks (computed if NULL)
#'
#' @return List with success, Q, R arrays, and Nbar
#'
#' @seealso [adcc_copula_nll()] for ADCC likelihood
#' @keywords internal
NULL


#' @title Semi-Parametric Distribution Transformation
#' @name fit_spd_transform
#'
#' @description
#' Fits a Semi-Parametric Distribution (SPD) for PIT transformation.
#'
#' @details
#' SPD combines:
#' - **Interior**: Kernel density estimation with Sheather-Jones bandwidth
#' - **Tails**: Empirical quantile-based scaling (or GPD if tsdistributions
#'   package is available)
#'
#' This approach is more robust to marginal distribution misspecification
#' than purely parametric transforms, while maintaining smooth behavior
#' compared to purely empirical transforms.
#'
#' The default thresholds are 10th and 90th percentiles.
#'
#' @param z Standardized residuals for one series
#' @param uni_fit Univariate GARCH fit (for parametric tails)
#' @param dist_pars Distribution parameters
#' @param lower_threshold Lower quantile for tail cutoff (default 0.1)
#' @param upper_threshold Upper quantile for tail cutoff (default 0.9)
#'
#' @return List with u (uniform values) and spd_model (if available)
#'
#' @keywords internal
NULL


## ---------------------------------------------------------------------------
## SECTION 3: Comparison Summary (for vignette reference)
## ---------------------------------------------------------------------------

#' @title CGARCH vs DCC Comparison
#' @name cgarch_vs_dcc
#'
#' @description
#' Summary of differences between Copula GARCH (CGARCH) and DCC models.
#'
#' @details
#' ## Model Structure
#'
#' | Aspect | DCC | CGARCH |
#' |--------|-----|--------|
#' | Marginals | Normal | Flexible (via PIT) |
#' | Dependence | Direct correlation | Copula-based |
#' | Tail dependence | Limited | Student-t copula |
#' | Computation | Faster | More flexible |
#'
#' ## When to Use CGARCH
#'
#' Consider CGARCH when:
#' - Non-normal marginal distributions are important
#' - Tail dependence is a concern (use MVT copula)
#' - Flexibility in marginal specification is needed
#' - Semi-parametric transforms are preferred
#'
#' ## When to Use DCC
#'
#' DCC may be preferable when:
#' - Computational speed is critical
#' - Normal margins are reasonable
#' - Simpler model interpretation is desired
#' - Analytical gradients for all parameters are needed
#'
#' ## Parameter Correspondence
#'
#' Both models use α (shock response) and β (persistence) for dynamics.
#' ADCC adds γ (asymmetry) in both frameworks.
#'
#' The shape parameter ν in MVT copula controls tail dependence:
#' - Lower ν → stronger tail dependence
#' - ν → ∞ → Gaussian copula (no tail dependence)
#'
#' @seealso `vignette("cgarch_vs_dcc", package = "tsbs")`
#' @keywords internal
NULL
