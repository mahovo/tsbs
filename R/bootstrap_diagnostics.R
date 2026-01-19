#' @title Bootstrap Diagnostics System for tsbs Package
#' @description Comprehensive diagnostic collection and reporting for all 
#'   bootstrap types supported by tsbs(). Provides insights into how bootstrap
#'   series are composed, block statistics, and quality metrics.
#' @name bootstrap_diagnostics
NULL


# =============================================================================
# Diagnostic Collector Creation
# =============================================================================

#' Create Bootstrap Diagnostic Collector
#'
#' Initializes a diagnostic collector object for storing bootstrap metadata
#' and statistics during the bootstrap procedure.
#'
#' @param bs_type Character string specifying bootstrap type.
#' @param n_original Integer, length of original series.
#' @param n_vars Integer, number of variables/columns.
#' @param num_boots Integer, number of bootstrap replicates.
#'
#' @return An object of class \code{tsbs_diagnostics}.
#' @keywords internal
#' @export
create_bootstrap_diagnostics <- function(
    bs_type,
    n_original,
    n_vars,
    num_boots
) {
  
  diagnostics <- list(
    # Metadata
    meta = list(
      bs_type = bs_type,
      n_original = n_original,
      n_vars = n_vars,
      num_boots = num_boots,
      timestamp = Sys.time()
    ),
    
    # Block composition (for block-based methods)
    blocks = list(
      # Per-replicate block info: list of data frames
      replicate_blocks = vector("list", num_boots),
      # Aggregated block length distribution
      all_block_lengths = integer(0),
      # Block starting positions
      all_start_positions = integer(0)
    ),
    
    # Bootstrap series statistics
    series_stats = list(
      # Per-replicate: means, sds, autocorrelations
      replicate_means = vector("list", num_boots),
      replicate_sds = vector("list", num_boots),
      replicate_ac1 = vector("list", num_boots),  # lag-1 autocorrelation
      replicate_lengths = integer(num_boots)
    ),
    
    # Original series reference statistics
    original_stats = list(
      means = NULL,
      sds = NULL,
      ac1 = NULL
    ),
    
    # Method-specific diagnostics
    method_specific = list(),
    
    # Configuration used
    config = list()
  )
  
  class(diagnostics) <- "tsbs_diagnostics"
  return(diagnostics)
}


# =============================================================================
# Recording Functions
# =============================================================================

#' Record Block Information for a Bootstrap Replicate
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param replicate_idx Integer, which bootstrap replicate (1-indexed).
#' @param block_lengths Integer vector of block lengths used.
#' @param start_positions Integer vector of starting positions in original series.
#' @param block_type Character, type of block (e.g., "overlapping").
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_blocks <- function(
    diagnostics,
    replicate_idx,
    block_lengths,
    start_positions,
    block_type = NA_character_
) {
  
  block_info <- data.frame(
    block_num = seq_along(block_lengths),
    length = block_lengths,
    start_pos = start_positions,
    block_type = block_type,
    stringsAsFactors = FALSE
  )
  
  diagnostics$blocks$replicate_blocks[[replicate_idx]] <- block_info
  diagnostics$blocks$all_block_lengths <- c(
    diagnostics$blocks$all_block_lengths, 
    block_lengths
  )
  diagnostics$blocks$all_start_positions <- c(
    diagnostics$blocks$all_start_positions, 
    start_positions
  )
  
  return(diagnostics)
}


#' Record Statistics for a Bootstrap Replicate
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param replicate_idx Integer, which bootstrap replicate.
#' @param bootstrap_matrix The bootstrap series matrix.
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_replicate_stats <- function(diagnostics, replicate_idx, bootstrap_matrix) {
  
  bootstrap_matrix <- as.matrix(bootstrap_matrix)
  n <- nrow(bootstrap_matrix)
  k <- ncol(bootstrap_matrix)
  
  # Column means
  diagnostics$series_stats$replicate_means[[replicate_idx]] <- colMeans(bootstrap_matrix)
  

  # Column standard deviations
  diagnostics$series_stats$replicate_sds[[replicate_idx]] <- apply(bootstrap_matrix, 2, sd)
  
  # Lag-1 autocorrelation for each column
  ac1_vec <- numeric(k)
  for (j in seq_len(k)) {
    if (n > 1) {
      ac1_vec[j] <- tryCatch(
        acf(bootstrap_matrix[, j], lag.max = 1, plot = FALSE)$acf[2, 1, 1],
        error = function(e) NA_real_
      )
    } else {
      ac1_vec[j] <- NA_real_
    }
  }
  diagnostics$series_stats$replicate_ac1[[replicate_idx]] <- ac1_vec
  
  # Series length
  diagnostics$series_stats$replicate_lengths[replicate_idx] <- n
  
  return(diagnostics)
}


#' Record Original Series Statistics
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param x Original data matrix.
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_original_stats <- function(diagnostics, x) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  k <- ncol(x)
  
  diagnostics$original_stats$means <- colMeans(x)
  diagnostics$original_stats$sds <- apply(x, 2, sd)
  
  ac1_vec <- numeric(k)
  for (j in seq_len(k)) {
    if (n > 1) {
      ac1_vec[j] <- tryCatch(
        acf(x[, j], lag.max = 1, plot = FALSE)$acf[2, 1, 1],
        error = function(e) NA_real_
      )
    } else {
      ac1_vec[j] <- NA_real_
    }
  }
  diagnostics$original_stats$ac1 <- ac1_vec
  
  return(diagnostics)
}


#' Record Method-Specific Diagnostics
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param ... Named method-specific items.
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_method_specific <- function(diagnostics, ...) {
  items <- list(...)
  diagnostics$method_specific <- modifyList(diagnostics$method_specific, items)
  return(diagnostics)
}


# =============================================================================
# Post-Processing: Compute from Bootstrap Series
# =============================================================================

#' Compute Diagnostics from Bootstrap Series
#'
#' After bootstrap series have been generated, this function computes diagnostic
#' statistics from the actual series. Use this when block-level tracking is not
#' available (e.g., when using C++ backend).
#'
#' @param bootstrap_series List of bootstrap replicate matrices.
#' @param original_data Original data matrix.
#' @param bs_type Character string specifying bootstrap type.
#' @param config Named list of configuration parameters used.
#'
#' @return A \code{tsbs_diagnostics} object.
#' @export
compute_bootstrap_diagnostics <- function(
    bootstrap_series,
    original_data,
    bs_type,
    config = list()
) {
  
  original_data <- as.matrix(original_data)
  n_original <- nrow(original_data)
  n_vars <- ncol(original_data)
  num_boots <- length(bootstrap_series)
  

  # Create collector
  diagnostics <- create_bootstrap_diagnostics(
    bs_type = bs_type,
    n_original = n_original,
    n_vars = n_vars,
    num_boots = num_boots
  )
  
  # Record original stats
  diagnostics <- record_original_stats(diagnostics, original_data)
  
  # Record config
  diagnostics$config <- config
  
  # Record replicate stats
  for (i in seq_len(num_boots)) {
    diagnostics <- record_replicate_stats(diagnostics, i, bootstrap_series[[i]])
  }
  
  # Compute block length estimates for block-based methods
  if (bs_type %in% c("moving", "stationary")) {
    diagnostics <- estimate_block_lengths_from_series(
      diagnostics, 
      bootstrap_series, 
      original_data
    )
  }
  
  return(diagnostics)
}


#' Estimate Block Lengths from Bootstrap Series
#'
#' For block bootstrap methods, attempts to identify block boundaries by 
#' detecting discontinuities in the bootstrap series relative to the original.
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param bootstrap_series List of bootstrap replicate matrices.
#' @param original_data Original data matrix.
#'
#' @return Updated diagnostics object with estimated block information.
#' @keywords internal
estimate_block_lengths_from_series <- function(
    diagnostics, 
    bootstrap_series, 
    original_data
) {
  
  # This is a heuristic approach since we don't have direct block tracking

  # from the C++ implementation
  
  original_data <- as.matrix(original_data)
  n_orig <- nrow(original_data)
  
  all_estimated_lengths <- integer(0)
  
  for (b in seq_along(bootstrap_series)) {
    boot_mat <- as.matrix(bootstrap_series[[b]])
    n_boot <- nrow(boot_mat)
    
    if (n_boot < 2) next
    
    # Find potential block boundaries by detecting jumps
    # A new block starts when the bootstrap series "jumps" to a different

    # position in the original series
    
    # For each row, find which original row it matches best
    matched_rows <- integer(n_boot)
    for (i in seq_len(n_boot)) {
      # Euclidean distance to all original rows
      dists <- rowSums((sweep(original_data, 2, boot_mat[i, ]))^2)
      matched_rows[i] <- which.min(dists)
    }
    
    # Detect block boundaries: when consecutive rows don't follow sequentially
    block_starts <- 1L
    for (i in 2:n_boot) {
      expected_next <- (matched_rows[i - 1] %% n_orig) + 1
      if (matched_rows[i] != expected_next) {
        block_starts <- c(block_starts, i)
      }
    }
    
    # Compute block lengths
    block_ends <- c(block_starts[-1] - 1, n_boot)
    block_lengths <- block_ends - block_starts + 1
    
    all_estimated_lengths <- c(all_estimated_lengths, block_lengths)
    
    # Record for this replicate
    diagnostics$blocks$replicate_blocks[[b]] <- data.frame(
      block_num = seq_along(block_lengths),
      length = block_lengths,
      start_pos = matched_rows[block_starts],
      block_type = "estimated",
      stringsAsFactors = FALSE
    )
  }
  
  diagnostics$blocks$all_block_lengths <- all_estimated_lengths
  
  return(diagnostics)
}


# =============================================================================
# Summary Method
# =============================================================================

#' Summary Method for Bootstrap Diagnostics
#'
#' @param object A \code{tsbs_diagnostics} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the object. Called for side effects (printing).
#' @export
summary.tsbs_diagnostics <- function(object, ...) {
  
  cat("========================================\n")
  cat("  tsbs Bootstrap Diagnostics Summary\n")
  cat("========================================\n\n")
  
  # --- Metadata ---
  cat("BOOTSTRAP CONFIGURATION:\n")
  cat("  Bootstrap type:", object$meta$bs_type, "\n")
  cat("  Original series length:", object$meta$n_original, "\n")
  cat("  Number of variables:", object$meta$n_vars, "\n")
  cat("  Number of replicates:", object$meta$num_boots, "\n")
  cat("  Generated:", format(object$meta$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
  
  # --- Configuration ---
  if (length(object$config) > 0) {
    cat("\nCONFIGURATION PARAMETERS:\n")
    for (nm in names(object$config)) {
      val <- object$config[[nm]]
      if (is.null(val)) {
        cat("  ", nm, ": NULL (auto-computed)\n", sep = "")
      } else {
        cat("  ", nm, ": ", as.character(val), "\n", sep = "")
      }
    }
  }
  
  # --- Bootstrap Series Lengths ---
  lengths <- object$series_stats$replicate_lengths
  if (length(lengths) > 0 && !all(is.na(lengths))) {
    cat("\nBOOTSTRAP SERIES LENGTHS:\n")
    cat("  Min:", min(lengths, na.rm = TRUE), "\n")
    cat("  Max:", max(lengths, na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(lengths, na.rm = TRUE), 2), "\n")
    if (min(lengths) != max(lengths)) {
      cat("  SD:", round(sd(lengths, na.rm = TRUE), 2), "\n")
      cat("  Note: Variable-length bootstrap series detected\n")
    }
  }
  
  # --- Block Statistics ---
  if (object$meta$bs_type %in% c("moving", "stationary", "hmm", "msvar")) {
    block_lengths <- object$blocks$all_block_lengths
    if (length(block_lengths) > 0) {
      cat("\nBLOCK LENGTH STATISTICS:\n")
      cat("  Total blocks sampled:", length(block_lengths), "\n")
      cat("  Mean block length:", round(mean(block_lengths), 2), "\n")
      cat("  SD block length:", round(sd(block_lengths), 2), "\n")
      cat("  Min block length:", min(block_lengths), "\n")
      cat("  Max block length:", max(block_lengths), "\n")
      cat("  Median block length:", median(block_lengths), "\n")
      
      # Quartiles
      qs <- quantile(block_lengths, probs = c(0.25, 0.75))
      cat("  25th percentile:", qs[1], "\n")
      cat("  75th percentile:", qs[2], "\n")
      
      # Blocks per replicate
      blocks_per_rep <- sapply(object$blocks$replicate_blocks, function(x) {
        if (is.null(x)) NA_integer_ else nrow(x)
      })
      blocks_per_rep <- blocks_per_rep[!is.na(blocks_per_rep)]
      if (length(blocks_per_rep) > 0) {
        cat("  Mean blocks per replicate:", round(mean(blocks_per_rep), 2), "\n")
      }
    }
  }
  
  # --- Original vs Bootstrap Comparison ---
  if (!is.null(object$original_stats$means)) {
    cat("\nORIGINAL vs BOOTSTRAP STATISTICS:\n")
    
    # Means comparison
    orig_means <- object$original_stats$means
    boot_means <- do.call(rbind, object$series_stats$replicate_means)
    if (!is.null(boot_means) && nrow(boot_means) > 0) {
      boot_mean_avg <- colMeans(boot_means, na.rm = TRUE)
      boot_mean_sd <- apply(boot_means, 2, sd, na.rm = TRUE)
      
      cat("\n  MEANS:\n")
      for (j in seq_along(orig_means)) {
        cat("    Variable", j, ":\n")
        cat("      Original:", round(orig_means[j], 4), "\n")
        cat("      Bootstrap avg:", round(boot_mean_avg[j], 4), "\n")
        cat("      Bootstrap SD:", round(boot_mean_sd[j], 4), "\n")
        bias <- boot_mean_avg[j] - orig_means[j]
        cat("      Bias:", round(bias, 4), "\n")
      }
    }
    
    # Autocorrelation comparison
    orig_ac1 <- object$original_stats$ac1
    boot_ac1 <- do.call(rbind, object$series_stats$replicate_ac1)
    if (!is.null(boot_ac1) && nrow(boot_ac1) > 0 && !all(is.na(boot_ac1))) {
      boot_ac1_avg <- colMeans(boot_ac1, na.rm = TRUE)
      boot_ac1_sd <- apply(boot_ac1, 2, sd, na.rm = TRUE)
      
      cat("\n  LAG-1 AUTOCORRELATION:\n")
      for (j in seq_along(orig_ac1)) {
        cat("    Variable", j, ":\n")
        cat("      Original:", round(orig_ac1[j], 4), "\n")
        cat("      Bootstrap avg:", round(boot_ac1_avg[j], 4), "\n")
        cat("      Bootstrap SD:", round(boot_ac1_sd[j], 4), "\n")
      }
    }
  }
  
  # --- Method-Specific Information ---
  if (length(object$method_specific) > 0) {
    cat("\nMETHOD-SPECIFIC DIAGNOSTICS:\n")
    for (nm in names(object$method_specific)) {
      val <- object$method_specific[[nm]]
      if (is.atomic(val) && length(val) == 1) {
        cat("  ", nm, ": ", as.character(val), "\n", sep = "")
      } else if (is.atomic(val) && length(val) <= 5) {
        cat("  ", nm, ": ", paste(round(val, 4), collapse = ", "), "\n", sep = "")
      } else {
        cat("  ", nm, ": [", class(val)[1], " of length ", length(val), "]\n", sep = "")
      }
    }
  }
  
  cat("\n========================================\n")
  
  invisible(object)
}


# =============================================================================
# Print Method
# =============================================================================

#' Print Method for Bootstrap Diagnostics
#'
#' @param x A \code{tsbs_diagnostics} object.
#' @param ... Additional arguments passed to summary.
#'
#' @return Invisibly returns the object.
#' @export
print.tsbs_diagnostics <- function(x, ...) {
  cat("tsbs Bootstrap Diagnostics\n")
  cat("  Type:", x$meta$bs_type, "\n")
  cat("  Replicates:", x$meta$num_boots, "\n")
  cat("  Original series: ", x$meta$n_original, " x ", x$meta$n_vars, "\n", sep = "")
  cat("\nUse summary() for detailed diagnostics.\n")
  cat("Use plot() for visualizations.\n")
  invisible(x)
}


# =============================================================================
# Plot Method
# =============================================================================

#' Plot Method for Bootstrap Diagnostics
#'
#' Creates diagnostic visualizations for bootstrap analysis.
#'
#' @param x A \code{tsbs_diagnostics} object.
#' @param type Character string specifying plot type. One of:
#'   \describe{
#'     \item{\code{"all"}}{Produce all applicable plots (default).}
#'     \item{\code{"block_lengths"}}{Histogram of block lengths.}
#'     \item{\code{"start_positions"}}{Histogram of block starting positions.}
#'     \item{\code{"means_comparison"}}{Compare original vs bootstrap means.}
#'     \item{\code{"acf_comparison"}}{Compare original vs bootstrap autocorrelation.}
#'     \item{\code{"length_distribution"}}{Distribution of bootstrap series lengths.}
#'   }
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns list of ggplot objects if ggplot2 available,
#'   otherwise uses base R graphics.
#' @export
plot.tsbs_diagnostics <- function(
    x, 
    type = c("all", "block_lengths", "start_positions", "means_comparison", 
             "acf_comparison", "length_distribution"),
    ...
) {
  
  type <- match.arg(type)
  use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  
  plots <- list()
  
  # Block length distribution
  if (type %in% c("all", "block_lengths")) {
    if (length(x$blocks$all_block_lengths) > 0) {
      plots$block_lengths <- plot_block_lengths(x, use_ggplot)
    }
  }
  
  # Start position distribution
  if (type %in% c("all", "start_positions")) {
    if (length(x$blocks$all_start_positions) > 0) {
      plots$start_positions <- plot_start_positions(x, use_ggplot)
    }
  }
  
  # Means comparison
  if (type %in% c("all", "means_comparison")) {
    if (!is.null(x$original_stats$means)) {
      plots$means_comparison <- plot_means_comparison(x, use_ggplot)
    }
  }
  
  # ACF comparison
  if (type %in% c("all", "acf_comparison")) {
    if (!is.null(x$original_stats$ac1)) {
      plots$acf_comparison <- plot_acf_comparison(x, use_ggplot)
    }
  }
  
  # Length distribution
  if (type %in% c("all", "length_distribution")) {
    lengths <- x$series_stats$replicate_lengths
    lengths <- lengths[!is.na(lengths)]
    if (length(lengths) > 0 && length(unique(lengths)) > 1) {
      plots$length_distribution <- plot_length_distribution(x, use_ggplot)
    } else if (type == "length_distribution") {
      if (length(lengths) == 0) {
        message("No series length data available.")
      } else {
        message("All ", length(lengths), " bootstrap series have the same length (", 
                unique(lengths), ") - no distribution to plot.")
      }
    }
  }
  
  if (length(plots) == 0 && type == "all") {
    message("No diagnostic plots available for this bootstrap configuration.")
  }
  
  invisible(NULL)
}


#' @keywords internal
plot_block_lengths <- function(x, use_ggplot = TRUE) {
  
  block_lengths <- x$blocks$all_block_lengths
  
  if (use_ggplot) {
    df <- data.frame(length = block_lengths)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = length)) +
      ggplot2::geom_histogram(bins = min(30, length(unique(block_lengths))),
                              fill = "steelblue", color = "white", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = mean(block_lengths), 
                          linetype = "dashed", color = "red", linewidth = 1) +
      ggplot2::labs(
        title = "Distribution of Block Lengths",
        subtitle = paste0("Mean = ", round(mean(block_lengths), 2), 
                          ", SD = ", round(sd(block_lengths), 2)),
        x = "Block Length",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
    return(p)
  } else {
    hist(block_lengths, 
         main = "Distribution of Block Lengths",
         xlab = "Block Length",
         col = "steelblue",
         border = "white")
    abline(v = mean(block_lengths), col = "red", lwd = 2, lty = 2)
    legend("topright", 
           legend = paste("Mean =", round(mean(block_lengths), 2)),
           col = "red", lty = 2, lwd = 2)
    return(invisible(NULL))
  }
}


#' @keywords internal
plot_start_positions <- function(x, use_ggplot = TRUE) {
  
  start_pos <- x$blocks$all_start_positions
  n_orig <- x$meta$n_original
  
  if (use_ggplot) {
    df <- data.frame(position = start_pos)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = position)) +
      ggplot2::geom_histogram(bins = min(50, n_orig / 2),
                              fill = "darkgreen", color = "white", alpha = 0.7) +
      ggplot2::labs(
        title = "Distribution of Block Starting Positions",
        subtitle = paste0("Original series length = ", n_orig),
        x = "Starting Position in Original Series",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
    return(p)
  } else {
    hist(start_pos,
         main = "Distribution of Block Starting Positions",
         xlab = "Starting Position",
         col = "darkgreen",
         border = "white")
    return(invisible(NULL))
  }
}


#' @keywords internal
plot_means_comparison <- function(x, use_ggplot = TRUE) {
  
  orig_means <- x$original_stats$means
  boot_means <- do.call(rbind, x$series_stats$replicate_means)
  
  if (is.null(boot_means) || nrow(boot_means) == 0) {
    message("No bootstrap means available for plotting")
    return(invisible(NULL))
  }
  
  n_vars <- length(orig_means)
  
  if (use_ggplot) {
    # Reshape for ggplot
    df_list <- lapply(seq_len(n_vars), function(j) {
      data.frame(
        variable = paste0("V", j),
        value = boot_means[, j],
        type = "Bootstrap"
      )
    })
    df_boot <- do.call(rbind, df_list)
    
    df_orig <- data.frame(
      variable = paste0("V", seq_len(n_vars)),
      value = orig_means,
      type = "Original"
    )
    
    p <- ggplot2::ggplot(df_boot, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", 
                              color = "white", alpha = 0.7) +
      ggplot2::geom_vline(data = df_orig, 
                          ggplot2::aes(xintercept = value),
                          color = "red", linewidth = 1, linetype = "dashed") +
      ggplot2::facet_wrap(~ variable, scales = "free") +
      ggplot2::labs(
        title = "Bootstrap Distribution of Means",
        subtitle = "Red dashed line = original mean",
        x = "Mean",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
    return(p)
  } else {
    par(mfrow = c(ceiling(sqrt(n_vars)), ceiling(n_vars / ceiling(sqrt(n_vars)))))
    for (j in seq_len(n_vars)) {
      hist(boot_means[, j],
           main = paste("Variable", j, "- Mean"),
           xlab = "Bootstrap Mean",
           col = "steelblue",
           border = "white")
      abline(v = orig_means[j], col = "red", lwd = 2, lty = 2)
    }
    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }
}


#' @keywords internal
plot_acf_comparison <- function(x, use_ggplot = TRUE) {
  
  orig_ac1 <- x$original_stats$ac1
  boot_ac1 <- do.call(rbind, x$series_stats$replicate_ac1)
  
  if (is.null(boot_ac1) || nrow(boot_ac1) == 0 || all(is.na(boot_ac1))) {
    message("No bootstrap autocorrelations available for plotting")
    return(invisible(NULL))
  }
  
  n_vars <- length(orig_ac1)
  
  if (use_ggplot) {
    df_list <- lapply(seq_len(n_vars), function(j) {
      data.frame(
        variable = paste0("V", j),
        value = boot_ac1[, j]
      )
    })
    df_boot <- do.call(rbind, df_list)
    
    df_orig <- data.frame(
      variable = paste0("V", seq_len(n_vars)),
      value = orig_ac1
    )
    
    p <- ggplot2::ggplot(df_boot, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = 30, fill = "purple", 
                              color = "white", alpha = 0.7) +
      ggplot2::geom_vline(data = df_orig,
                          ggplot2::aes(xintercept = value),
                          color = "red", linewidth = 1, linetype = "dashed") +
      ggplot2::facet_wrap(~ variable, scales = "free") +
      ggplot2::labs(
        title = "Bootstrap Distribution of Lag-1 Autocorrelation",
        subtitle = "Red dashed line = original AC(1)",
        x = "Lag-1 Autocorrelation",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
    return(p)
  } else {
    par(mfrow = c(ceiling(sqrt(n_vars)), ceiling(n_vars / ceiling(sqrt(n_vars)))))
    for (j in seq_len(n_vars)) {
      ac_vals <- boot_ac1[, j]
      ac_vals <- ac_vals[!is.na(ac_vals)]
      if (length(ac_vals) > 0) {
        hist(ac_vals,
             main = paste("Variable", j, "- AC(1)"),
             xlab = "Lag-1 Autocorrelation",
             col = "purple",
             border = "white")
        if (!is.na(orig_ac1[j])) {
          abline(v = orig_ac1[j], col = "red", lwd = 2, lty = 2)
        }
      }
    }
    par(mfrow = c(1, 1))
    return(invisible(NULL))
  }
}


#' @keywords internal
plot_length_distribution <- function(x, use_ggplot = TRUE) {
  
  lengths <- x$series_stats$replicate_lengths
  lengths <- lengths[!is.na(lengths)]
  
  if (length(lengths) == 0) {
    message("No length data available")
    return(invisible(NULL))
  }
  
  if (use_ggplot) {
    df <- data.frame(length = lengths)
    p <- ggplot2::ggplot(df, ggplot2::aes(x = length)) +
      ggplot2::geom_histogram(bins = min(30, length(unique(lengths))),
                              fill = "orange", color = "white", alpha = 0.7) +
      ggplot2::labs(
        title = "Distribution of Bootstrap Series Lengths",
        subtitle = paste0("Mean = ", round(mean(lengths), 2),
                          ", Range = [", min(lengths), ", ", max(lengths), "]"),
        x = "Series Length",
        y = "Count"
      ) +
      ggplot2::theme_minimal()
    print(p)
    return(p)
  } else {
    hist(lengths,
         main = "Distribution of Bootstrap Series Lengths",
         xlab = "Series Length",
         col = "orange",
         border = "white")
    return(invisible(NULL))
  }
}


# =============================================================================
# Extraction Functions
# =============================================================================

#' Extract Block Information from Diagnostics
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param replicate Optional integer specifying which replicate. If NULL, 
#'   returns combined data for all replicates.
#'
#' @return A data frame with block information.
#' @export
extract_blocks <- function(diagnostics, replicate = NULL) {
  
  if (!is.null(replicate)) {
    return(diagnostics$blocks$replicate_blocks[[replicate]])
  }
  
  # Combine all replicates
  all_blocks <- lapply(seq_along(diagnostics$blocks$replicate_blocks), function(i) {
    block_df <- diagnostics$blocks$replicate_blocks[[i]]
    if (!is.null(block_df)) {
      block_df$replicate <- i
      block_df
    } else {
      NULL
    }
  })
  
  do.call(rbind, all_blocks)
}


#' Extract Summary Statistics from Diagnostics
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#'
#' @return A list with summary statistics.
#' @export
extract_summary_stats <- function(diagnostics) {
  
  boot_means <- do.call(rbind, diagnostics$series_stats$replicate_means)
  boot_sds <- do.call(rbind, diagnostics$series_stats$replicate_sds)
  boot_ac1 <- do.call(rbind, diagnostics$series_stats$replicate_ac1)
  
  list(
    original = diagnostics$original_stats,
    bootstrap = list(
      means = if (!is.null(boot_means)) {
        list(
          mean = colMeans(boot_means, na.rm = TRUE),
          sd = apply(boot_means, 2, sd, na.rm = TRUE),
          quantiles = apply(boot_means, 2, quantile, 
                            probs = c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
        )
      },
      sds = if (!is.null(boot_sds)) {
        list(
          mean = colMeans(boot_sds, na.rm = TRUE),
          sd = apply(boot_sds, 2, sd, na.rm = TRUE)
        )
      },
      ac1 = if (!is.null(boot_ac1)) {
        list(
          mean = colMeans(boot_ac1, na.rm = TRUE),
          sd = apply(boot_ac1, 2, sd, na.rm = TRUE)
        )
      }
    ),
    block_lengths = if (length(diagnostics$blocks$all_block_lengths) > 0) {
      list(
        mean = mean(diagnostics$blocks$all_block_lengths),
        sd = sd(diagnostics$blocks$all_block_lengths),
        min = min(diagnostics$blocks$all_block_lengths),
        max = max(diagnostics$blocks$all_block_lengths),
        median = median(diagnostics$blocks$all_block_lengths)
      )
    }
  )
}


#' Convert Diagnostics to Data Frame
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param what Character specifying what to convert: "blocks", "stats", or "all".
#'
#' @return A data frame.
#' @export
as.data.frame.tsbs_diagnostics <- function(x, ..., what = c("stats", "blocks", "all")) {
  
  what <- match.arg(what)
  
  if (what == "blocks") {
    return(extract_blocks(x))
  }
  
  if (what == "stats") {
    # Replicate-level statistics
    n <- x$meta$num_boots
    
    df <- data.frame(
      replicate = seq_len(n),
      length = x$series_stats$replicate_lengths
    )
    
    # Add means for each variable
    boot_means <- do.call(rbind, x$series_stats$replicate_means)
    if (!is.null(boot_means)) {
      colnames(boot_means) <- paste0("mean_V", seq_len(ncol(boot_means)))
      df <- cbind(df, boot_means)
    }
    
    # Add sds for each variable
    boot_sds <- do.call(rbind, x$series_stats$replicate_sds)
    if (!is.null(boot_sds)) {
      colnames(boot_sds) <- paste0("sd_V", seq_len(ncol(boot_sds)))
      df <- cbind(df, boot_sds)
    }
    
    # Add ac1 for each variable
    boot_ac1 <- do.call(rbind, x$series_stats$replicate_ac1)
    if (!is.null(boot_ac1)) {
      colnames(boot_ac1) <- paste0("ac1_V", seq_len(ncol(boot_ac1)))
      df <- cbind(df, boot_ac1)
    }
    
    return(df)
  }
  
  # what == "all" - combine blocks and stats
  stats_df <- as.data.frame.tsbs_diagnostics(x, what = "stats")
  blocks_df <- extract_blocks(x)
  
  if (!is.null(blocks_df)) {
    # Summarize blocks per replicate
    block_summary <- aggregate(
      length ~ replicate,
      data = blocks_df,
      FUN = function(x) c(n_blocks = length(x), 
                          mean_length = mean(x), 
                          sd_length = sd(x))
    )
    block_summary <- do.call(data.frame, block_summary)
    names(block_summary) <- c("replicate", "n_blocks", "mean_block_length", "sd_block_length")
    
    stats_df <- merge(stats_df, block_summary, by = "replicate", all.x = TRUE)
  }
  
  stats_df
}
