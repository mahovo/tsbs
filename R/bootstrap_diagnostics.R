#' @title Bootstrap Diagnostics System for tsbs Package
#' @description Comprehensive diagnostic collection and reporting for all 
#'   bootstrap types supported by tsbs(). Provides insights into how bootstrap
#'   series are composed, block statistics, and quality metrics.
#' @name bootstrap_diagnostics
NULL


# Diagnostic Collector Creation ================================================

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


# Recording Functions ==========================================================

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


# Post-Processing: Compute from Bootstrap Series ===============================

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


# Summary Method ===============================================================

#' Summary Method for Bootstrap Diagnostics
#'
#' @param object A \code{tsbs_diagnostics} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns the object. Called for side effects (printing).
#' @method summary tsbs_diagnostics
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


# Print Method =================================================================

#' Print Method for Bootstrap Diagnostics
#'
#' @param x A \code{tsbs_diagnostics} object.
#' @param ... Additional arguments passed to summary.
#'
#' @return Invisibly returns the object.
#' @method print tsbs_diagnostics
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


# Plot Method ==================================================================

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
#' @method plot tsbs_diagnostics
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


# Extraction Functions =========================================================

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


#' Summarize Bootstrap Functional Outputs
#'
#' Computes summary statistics (mean, SD, confidence intervals) for bootstrap
#' functional outputs such as portfolio weights or other derived quantities.
#'
#' @param func_outs List of functional outputs from \code{tsbs()}, or a matrix
#'   where each row is a bootstrap replicate and columns are output dimensions.
#' @param names Optional character vector of names for output dimensions.
#' @param probs Numeric vector of probabilities for quantile calculation.
#'   Defaults to \code{c(0.025, 0.975)} for 95% confidence intervals.
#'
#' @return A data frame with columns:
#'   \item{Name}{Dimension name}
#'   \item{Mean}{Bootstrap mean}
#'   \item{SD}{Bootstrap standard deviation}
#'   \item{CI_Lower}{Lower confidence bound}
#'   \item{CI_Upper}{Upper confidence bound}
#'   \item{CI_Width}{Width of confidence interval}
#'
#' @examples
#' \dontrun{
#' # After running tsbs() with a portfolio function
#' result <- tsbs(x, bs_type = "ms_varma_garch", func = risk_parity_portfolio, ...)
#' 
#' # Summarize the bootstrap weights
#' summary_df <- summarize_func_outs(result$func_outs, names = c("SPY", "EFA", "BND"))
#' print(summary_df)
#' }
#'
#' @export
summarize_func_outs <- function(func_outs, names = NULL, probs = c(0.025, 0.975)) {
  
  
  ## Convert list of vectors/matrices to matrix form
  if (is.list(func_outs) && !is.data.frame(func_outs)) {
    ## Handle various output formats from tsbs
    boot_mat <- tryCatch({
      do.call(rbind, lapply(func_outs, function(w) {
        if (is.matrix(w)) as.vector(t(w)) else as.vector(w)
      }))
    }, error = function(e) NULL)
    
    if (is.null(boot_mat)) {
      warning("Could not convert func_outs to matrix format")
      return(NULL)
    }
  } else if (is.matrix(func_outs) || is.data.frame(func_outs)) {
    boot_mat <- as.matrix(func_outs)
  } else {
    warning("func_outs must be a list, matrix, or data frame")
    return(NULL)
  }
  
  ## Validate
  
  if (nrow(boot_mat) < 2) {
    warning("Need at least 2 bootstrap replicates for summary statistics")
    return(NULL)
  }
  
  n_dims <- ncol(boot_mat)
  
  ## Set names if not provided
  if (is.null(names)) {
    if (!is.null(colnames(boot_mat))) {
      names <- colnames(boot_mat)
    } else {
      names <- paste0("V", seq_len(n_dims))
    }
  }
  
  if (length(names) != n_dims) {
    warning("Length of 'names' does not match number of dimensions")
    names <- paste0("V", seq_len(n_dims))
  }
  
  ## Compute statistics
  boot_mean <- colMeans(boot_mat, na.rm = TRUE)
  boot_sd <- apply(boot_mat, 2, sd, na.rm = TRUE)
  boot_quantiles <- apply(boot_mat, 2, quantile, probs = probs, na.rm = TRUE)
  
  data.frame(
    Name = names,
    Mean = round(boot_mean, 4),
    SD = round(boot_sd, 4),
    CI_Lower = round(boot_quantiles[1, ], 4),
    CI_Upper = round(boot_quantiles[2, ], 4),
    CI_Width = round(boot_quantiles[2, ] - boot_quantiles[1, ], 4),
    stringsAsFactors = FALSE
  )
}


#' Compute Coefficient of Variation for Bootstrap Outputs
#'
#' Calculates the coefficient of variation (CV) for each dimension of bootstrap
#' functional outputs, providing a measure of estimation stability.
#'
#' @param func_outs List of functional outputs from \code{tsbs()}, or a matrix.
#' @param names Optional character vector of names for output dimensions.
#' @param cv_thresholds Named numeric vector with thresholds for stability
#'   classification. Defaults to \code{c(Stable = 0.3, Moderate = 0.6)}.
#'
#' @return A data frame with columns:
#'   \item{Name}{Dimension name}
#'   \item{Mean}{Bootstrap mean}
#'   \item{SD}{Bootstrap standard deviation}
#'   \item{CV}{Coefficient of variation (SD/Mean)}
#'   \item{Stability}{Stability classification based on CV thresholds}
#'
#' @details
#' The coefficient of variation (CV) is defined as SD/Mean. Lower CV values
#' indicate more stable estimates. Default thresholds classify estimates as:
#' \itemize{
#'   \item \strong{Stable}: CV < 0.3
#'   \item \strong{Moderate}: 0.3 <= CV < 0.6
#'   \item \strong{Unstable}: CV >= 0.6
#' }
#'
#' @examples
#' \dontrun{
#' # After running tsbs() with a portfolio function
#' result <- tsbs(x, bs_type = "ms_varma_garch", func = risk_parity_portfolio, ...)
#' 
#' # Assess stability of bootstrap weights
#' stability_df <- compute_func_out_cv(result$func_outs, names = c("SPY", "EFA", "BND"))
#' print(stability_df)
#' }
#' @export
compute_func_out_cv <- function(func_outs, names = NULL, 
                                cv_thresholds = c(Stable = 0.3, Moderate = 0.6)) {
  
  ## Convert to matrix (reuse logic from summarize_func_outs)
  if (is.list(func_outs) && !is.data.frame(func_outs)) {
    boot_mat <- tryCatch({
      do.call(rbind, lapply(func_outs, function(w) {
        if (is.matrix(w)) as.vector(t(w)) else as.vector(w)
      }))
    }, error = function(e) NULL)
    
    if (is.null(boot_mat)) {
      warning("Could not convert func_outs to matrix format")
      return(NULL)
    }
  } else if (is.matrix(func_outs) || is.data.frame(func_outs)) {
    boot_mat <- as.matrix(func_outs)
  } else {
    warning("func_outs must be a list, matrix, or data frame")
    return(NULL)
  }
  
  if (nrow(boot_mat) < 2) {
    warning("Need at least 2 bootstrap replicates for CV calculation")
    return(NULL)
  }
  
  n_dims <- ncol(boot_mat)
  
  ## Set names
  if (is.null(names)) {
    if (!is.null(colnames(boot_mat))) {
      names <- colnames(boot_mat)
    } else {
      names <- paste0("V", seq_len(n_dims))
    }
  }
  
  ## Compute CV
  boot_mean <- colMeans(boot_mat, na.rm = TRUE)
  boot_sd <- apply(boot_mat, 2, sd, na.rm = TRUE)
  cv <- boot_sd / (abs(boot_mean) + 1e-8)  # Add small constant to avoid division by zero
  
  ## Classify stability
  stability <- sapply(cv, function(x) {
    if (x < cv_thresholds["Stable"]) {
      "Stable"
    } else if (x < cv_thresholds["Moderate"]) {
      "Moderate"
    } else {
      "Unstable"
    }
  })
  
  data.frame(
    Name = names,
    Mean = round(boot_mean, 4),
    SD = round(boot_sd, 4),
    CV = round(cv, 3),
    Stability = stability,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


#' Extract Bootstrap Output Matrix
#'
#' Converts the list of functional outputs from \code{tsbs()} into a matrix
#' format suitable for further analysis.
#'
#' @param func_outs List of functional outputs from \code{tsbs()}.
#'
#' @return A matrix where each row is a bootstrap replicate and each column
#'   is an output dimension. Returns NULL if conversion fails.
#'
#' @export
extract_func_out_matrix <- function(func_outs) {
  
  if (!is.list(func_outs)) {
    if (is.matrix(func_outs) || is.data.frame(func_outs)) {
      return(as.matrix(func_outs))
    }
    warning("func_outs must be a list, matrix, or data frame")
    return(NULL)
  }
  
  tryCatch({
    do.call(rbind, lapply(func_outs, function(w) {
      if (is.matrix(w)) as.vector(t(w)) else as.vector(w)
    }))
  }, error = function(e) {
    warning("Could not convert func_outs to matrix: ", e$message)
    NULL
  })
}


#' Compute Robust Estimates from Bootstrap Outputs
#'
#' Computes various robust estimators (mean, median, winsorized mean,
#' conservative quantile) from bootstrap functional outputs.
#'
#' @param func_outs List of functional outputs from \code{tsbs()}, or a matrix.
#' @param names Optional character vector of names for output dimensions.
#' @param point_est Optional numeric vector of point estimates to include
#'   in the comparison.
#' @param trim Trim proportion for winsorized mean. Defaults to 0.1.
#' @param conservative_quantile Quantile for conservative estimate. Defaults
#'   to 0.25.
#'
#' @return A data frame with columns for each robust estimator.
#'
#' @details
#' This function computes:
#' \itemize{
#'   \item \strong{Point}: Original point estimate (if provided)
#'   \item \strong{Boot_Mean}: Bootstrap mean
#'   \item \strong{Boot_Median}: Bootstrap median
#'   \item \strong{Winsorized}: Winsorized mean (trimmed mean)
#'   \item \strong{Conservative}: Lower quantile estimate
#' }
#'
#' For portfolio weights, the conservative estimate is renormalized to sum to 1.
#'
#' @export
compute_robust_estimates <- function(
    func_outs, 
    names = NULL, 
    point_est = NULL,
    trim = 0.1, 
    conservative_quantile = 0.25
  ) {
  
  ## Convert to matrix
  boot_mat <- extract_func_out_matrix(func_outs)
  
  if (is.null(boot_mat) || nrow(boot_mat) < 2) {
    warning("Insufficient bootstrap data for robust estimates")
    return(NULL)
  }
  
  n_dims <- ncol(boot_mat)
  
  ## Set names
  if (is.null(names)) {
    if (!is.null(colnames(boot_mat))) {
      names <- colnames(boot_mat)
    } else {
      names <- paste0("V", seq_len(n_dims))
    }
  }
  
  ## Compute estimators
  boot_mean <- colMeans(boot_mat, na.rm = TRUE)
  boot_median <- apply(boot_mat, 2, median, na.rm = TRUE)
  boot_winsor <- apply(boot_mat, 2, function(x) mean(x, trim = trim, na.rm = TRUE))
  boot_conservative <- apply(boot_mat, 2, quantile, probs = conservative_quantile, na.rm = TRUE)
  
  ## Renormalize conservative estimate if it looks like weights (all positive, sum ~ 1)
  ## Handle NA values gracefully
  if (!any(is.na(boot_conservative)) && 
      all(boot_conservative >= 0, na.rm = TRUE) && 
      abs(sum(boot_mean, na.rm = TRUE) - 1) < 0.1) {
    boot_conservative <- boot_conservative / sum(boot_conservative)
  }
  
  ## Build result
  result <- data.frame(
    Name = names,
    Boot_Mean = round(boot_mean, 3),
    Boot_Median = round(boot_median, 3),
    Winsorized = round(boot_winsor, 3),
    Conservative = round(boot_conservative, 3),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  ## Add point estimate if provided
  if (!is.null(point_est)) {
    if (length(point_est) == n_dims) {
      result <- cbind(
        data.frame(Name = names, Point = round(point_est, 3), stringsAsFactors = FALSE),
        result[, -1]  # Remove duplicate Name column
      )
    } else {
      warning("point_est length does not match number of dimensions")
    }
  }
  
  result
}


# Bootstrap Benchmarking Utilities =============================================

#' Benchmark Bootstrap Methods
#'
#' Compares the runtime of different \code{tsbs()} configurations by varying
#' a single parameter (e.g., series length, number of assets, or number of
#' bootstrap replicates).
#'
#' @param x Numeric matrix of input data for bootstrapping.
#' @param setups A named list of \code{tsbs()} argument lists. Each element
#'   should be a list of arguments to pass to \code{tsbs()}, excluding \code{x}
#'   and the varying parameter.
#' @param vary Character string specifying which parameter to vary. One of:
#'   \describe{
#'     \item{\code{"n_boot"}}{Vary the length of bootstrap series}
#'     \item{\code{"num_boots"}}{Vary the number of bootstrap replicates}
#'     \item{\code{"n_assets"}}{Vary the number of assets (columns)}
#'   }
#' @param values Numeric vector of values for the varying parameter.
#' @param times Integer number of times to repeat each benchmark for averaging.
#'   Default is 3.
#' @param verbose Logical. If TRUE, print progress messages. Default is TRUE.
#'
#' @return An object of class \code{tsbs_benchmark} containing:
#'   \item{results}{Data frame with columns: Setup, Parameter, Value, Time, Rep}
#'   \item{summary}{Data frame with mean and sd of times by Setup and Value}
#'   \item{vary}{The parameter that was varied}
#'   \item{setups}{The setup configurations used}
#'
#' @details
#' This function is useful for:
#' \itemize{
#'   \item Comparing computational cost of different bootstrap methods
#'   \item Understanding how runtime scales with data size
#'   \item Choosing appropriate settings for production use
#' }
#'
#' @examples
#' \dontrun{
#' # Create test data
#' set.seed(123)
#' x <- matrix(rnorm(500 * 3), ncol = 3)
#'
#' # Define setups to compare
#' setups <- list(
#'   "Plain" = list(bs_type = "moving", block_length = 1),
#'   "Moving" = list(bs_type = "moving", block_length = 5),
#'   "Stationary" = list(bs_type = "stationary")
#' )
#'
#' # Benchmark varying number of replicates
#' bench <- benchmark_tsbs(
#'   x = x,
#'   setups = setups,
#'   vary = "num_boots",
#'   values = c(10, 25, 50, 100),
#'   times = 3
#' )
#'
#' # View results
#' print(bench)
#' plot(bench)
#' }
#'
#' @seealso \code{\link{plot.tsbs_benchmark}} for visualization.
#' @export
benchmark_tsbs <- function(
    x,
    setups,
    vary = c("num_boots", "n_boot", "n_assets"),
    values,
    times = 3,
    verbose = TRUE
) {
  
  vary <- match.arg(vary)
  
  ## Validate inputs
  if (!is.list(setups) || is.null(names(setups))) {
    stop("'setups' must be a named list of tsbs() argument lists")
  }
  
  if (!is.numeric(values) || length(values) < 2) {
    stop("'values' must be a numeric vector with at least 2 elements")
  }
  
  x <- as.matrix(x)
  n_obs <- nrow(x)
  n_vars <- ncol(x)
  
  ## Storage for results
  results <- data.frame(
    Setup = character(),
    Parameter = character(),
    Value = numeric(),
    Time = numeric(),
    Rep = integer(),
    stringsAsFactors = FALSE
  )
  
  ## Total iterations for progress
  total_iter <- length(setups) * length(values) * times
  current_iter <- 0
  
  ## Run benchmarks
  for (setup_name in names(setups)) {
    setup_args <- setups[[setup_name]]
    
    for (val in values) {
      for (rep in seq_len(times)) {
        current_iter <- current_iter + 1
        
        if (verbose) {
          cat(sprintf("\r[%d/%d] %s: %s = %g (rep %d/%d)    ",
                      current_iter, total_iter, setup_name, vary, val, rep, times))
        }
        
        ## Prepare data and arguments based on vary parameter
        if (vary == "n_assets") {
          ## Subset columns
          if (val > n_vars) {
            warning(sprintf("Skipping n_assets = %d (only %d available)", val, n_vars))
            next
          }
          x_use <- x[, 1:val, drop = FALSE]
          args <- c(list(x = x_use), setup_args)
        } else if (vary == "n_boot") {
          x_use <- x
          args <- c(list(x = x_use, n_boot = val), setup_args)
        } else if (vary == "num_boots") {
          x_use <- x
          args <- c(list(x = x_use, num_boots = val), setup_args)
        }
        
        ## Time the call
        start_time <- Sys.time()
        run_success <- tryCatch({
          result <- do.call(tsbs, args)
          TRUE
        }, error = function(e) {
          if (verbose) {
            message("\n  Error in ", setup_name, ": ", conditionMessage(e))
          }
          FALSE
        })
        end_time <- Sys.time()
        
        elapsed <- if (isTRUE(run_success)) {
          as.numeric(difftime(end_time, start_time, units = "secs"))
        } else {
          NA_real_
        }
        
        ## Store result
        results <- rbind(results, data.frame(
          Setup = setup_name,
          Parameter = vary,
          Value = val,
          Time = elapsed,
          Rep = rep,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  if (verbose) cat("\n")
  
  ## Compute summary statistics
  summary_df <- aggregate(
    Time ~ Setup + Value,
    data = results,
    FUN = function(x) c(Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE))
  )
  summary_df <- do.call(data.frame, summary_df)
  names(summary_df) <- c("Setup", "Value", "Mean_Time", "SD_Time")
  
  ## Create result object
  result <- list(
    results = results,
    summary = summary_df,
    vary = vary,
    values = values,
    setups = setups,
    times = times,
    data_dim = c(n_obs = n_obs, n_vars = n_vars)
  )
  
  class(result) <- "tsbs_benchmark"
  return(result)
}


#' Print Method for tsbs_benchmark
#'
#' @param x A \code{tsbs_benchmark} object.
#' @param ... Additional arguments (unused).
#'
#' @export
print.tsbs_benchmark <- function(x, ...) {
  
  cat("tsbs Benchmark Results\n")
  cat("======================\n\n")
  
  cat("Varying parameter:", x$vary, "\n")
  cat("Values tested:", paste(x$values, collapse = ", "), "\n")
  cat("Repetitions per test:", x$times, "\n")
  cat("Input data dimensions:", x$data_dim["n_obs"], "x", x$data_dim["n_vars"], "\n\n")
  
  cat("Setups compared:\n")
  for (name in names(x$setups)) {
    setup <- x$setups[[name]]
    cat("  -", name, ": bs_type =", setup$bs_type)
    if (!is.null(setup$block_length)) cat(", block_length =", setup$block_length)
    cat("\n")
  }
  
  cat("\nMean runtime (seconds) by setup and", x$vary, ":\n\n")
  
  ## Reshape summary for display
  summary_wide <- reshape(
    x$summary[, c("Setup", "Value", "Mean_Time")],
    idvar = "Value",
    timevar = "Setup",
    direction = "wide"
  )
  names(summary_wide) <- gsub("Mean_Time\\.", "", names(summary_wide))
  
  print(summary_wide, row.names = FALSE, digits = 3)
  
  invisible(x)
}


#' Plot Method for tsbs_benchmark
#'
#' Creates a visualization of benchmark results showing how runtime scales
#' with the varying parameter across different bootstrap setups.
#'
#' @param x A \code{tsbs_benchmark} object.
#' @param log_scale Logical. If TRUE, use log scale for y-axis. Default is FALSE.
#' @param show_points Logical. If TRUE, show individual timing points. Default is TRUE.
#' @param show_ribbon Logical. If TRUE, show +/- 1 SD ribbon. Default is TRUE.
#' @param ... Additional arguments passed to ggplot.
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot.tsbs_benchmark <- function(
    x,
    log_scale = FALSE,
    show_points = TRUE,
    show_ribbon = TRUE,
    ...
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting benchmarks")
  }
  
  ## Merge individual results with summary for plotting
  plot_data <- merge(x$results, x$summary, by = c("Setup", "Value"))
  
  ## Determine axis labels
  x_label <- switch(
    x$vary,
    "num_boots" = "Number of Bootstrap Replicates",
    "n_boot" = "Bootstrap Series Length",
    "n_assets" = "Number of Assets"
  )
  
  ## Build plot
  p <- ggplot2::ggplot(x$summary, ggplot2::aes(x = Value, y = Mean_Time, color = Setup)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3)
  
  ## Add ribbon for SD
  if (show_ribbon && any(!is.na(x$summary$SD_Time))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = Mean_Time - SD_Time, ymax = Mean_Time + SD_Time, fill = Setup),
      alpha = 0.2,
      color = NA
    )
  }
  
  ## Add individual points
  if (show_points) {
    p <- p + ggplot2::geom_point(
      data = x$results,
      ggplot2::aes(x = Value, y = Time),
      alpha = 0.3,
      size = 1.5
    )
  }
  
  ## Styling
  p <- p +
    ggplot2::labs(
      title = "Bootstrap Method Runtime Comparison",
      subtitle = paste("Varying:", x_label),
      x = x_label,
      y = "Runtime (seconds)",
      color = "Method",
      fill = "Method"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  ## Log scale if requested
  if (log_scale) {
    p <- p + ggplot2::scale_y_log10()
  }
  
  print(p)
  invisible(p)
}


#' Summary Method for tsbs_benchmark
#'
#' @param object A \code{tsbs_benchmark} object.
#' @param ... Additional arguments (unused).
#'
#' @return The summary data frame (invisibly).
#' @export
summary.tsbs_benchmark <- function(object, ...) {
  
  cat("Benchmark Summary: Runtime by", object$vary, "\n\n")
  
  ## Add relative speed column (relative to fastest)
  summary_df <- object$summary
  
  for (val in unique(summary_df$Value)) {
    idx <- summary_df$Value == val
    min_time <- min(summary_df$Mean_Time[idx], na.rm = TRUE)
    summary_df$Relative[idx] <- summary_df$Mean_Time[idx] / min_time
  }
  
  summary_df$Mean_Time <- round(summary_df$Mean_Time, 4)
  summary_df$SD_Time <- round(summary_df$SD_Time, 4)
  summary_df$Relative <- round(summary_df$Relative, 2)
  
  print(summary_df, row.names = FALSE)
  
  ## Identify fastest method overall
  avg_by_setup <- aggregate(Mean_Time ~ Setup, data = object$summary, FUN = mean, na.rm = TRUE)
  fastest <- avg_by_setup$Setup[which.min(avg_by_setup$Mean_Time)]
  cat("\nFastest method (on average):", fastest, "\n")
  
  invisible(summary_df)
}


# Regime diagnostics ===========================================================

#' Record Regime/State Information for Bootstrap Diagnostics
#'
#' Records the state sequence from the original data and the state composition
#' of each bootstrap replicate. This is used by HMM and MS-VARMA-GARCH bootstrap
#' methods to track how regimes are sampled.
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param original_states Integer vector of states for the original series.
#' @param num_states Integer, total number of possible states.
#' @param state_labels Optional character vector of state labels.
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_regime_info <- function(
    diagnostics,
    original_states,
    num_states,
    state_labels = NULL
) {
  
  if (is.null(state_labels)) {
    state_labels <- paste("State", seq_len(num_states))
  }
  
  diagnostics$method_specific$regime_info <- list(
    original_states = original_states,
    num_states = num_states,
    state_labels = state_labels,
    original_state_counts = table(factor(original_states, levels = seq_len(num_states))),
    # Per-replicate state sequences (filled in by record_replicate_regimes)
    replicate_states = vector("list", diagnostics$meta$num_boots),
    # Block-to-state mapping for each replicate
    replicate_block_states = vector("list", diagnostics$meta$num_boots)
  )
  
  return(diagnostics)
}


#' Record Regime Composition for a Bootstrap Replicate
#'
#' Records the state sequence and block-state mapping for a single bootstrap
#' replicate.
#'
#' @param diagnostics A \code{tsbs_diagnostics} object.
#' @param replicate_idx Integer, which bootstrap replicate (1-indexed).
#' @param replicate_states Integer vector of states for this replicate.
#' @param block_states Optional integer vector of states for each block used.
#' @param block_source_indices Optional integer vector of which original block
#'   each bootstrap block came from.
#' @param source_indices Optional integer vector mapping each bootstrap time point
#'   to its source index in the original series. Used for probability lookup.
#'
#' @return Updated diagnostics object.
#' @keywords internal
#' @export
record_replicate_regimes <- function(
    diagnostics,
    replicate_idx,
    replicate_states,
    block_states = NULL,
    block_source_indices = NULL,
    source_indices = NULL
) {
  
  if (is.null(diagnostics$method_specific$regime_info)) {
    warning("record_regime_info() should be called before record_replicate_regimes()")
    return(diagnostics)
  }
  
  diagnostics$method_specific$regime_info$replicate_states[[replicate_idx]] <- replicate_states
  
  if (!is.null(block_states) || !is.null(source_indices)) {
    diagnostics$method_specific$regime_info$replicate_block_states[[replicate_idx]] <- list(
      states = block_states,
      block_source_indices = block_source_indices,
      source_indices = source_indices  # Per-observation source mapping
    )
  }
  
  return(diagnostics)
}


#' Plot Regime Composition of Bootstrap Series
#'
#' Creates a visualization comparing the original series with one or more
#' bootstrap replicates, showing the regime structure. Can display either
#' discrete regime bands or smoothed state probabilities.
#'
#' @param tsbs_result A list returned by \code{tsbs()} with 
#'   \code{collect_diagnostics = TRUE}, or a \code{tsbs_diagnostics} object.
#' @param original_data The original data matrix used for bootstrapping.
#' @param replicate_idx Integer or vector of integers specifying which bootstrap
#'   replicate(s) to plot. Default is 1 (first replicate).
#' @param series_idx Integer specifying which column/series to plot if 
#'   multivariate. Default is 1.
#' @param show_prices Logical. If TRUE and data appears to be returns, convert
#'   to cumulative prices for visualization. Default is TRUE.
#' @param initial_price Numeric. Starting price for cumulative calculation.
#'   Default is 100.
#' @param show_probabilities Logical. If TRUE, overlay smoothed state 
#'   probabilities on the plot instead of (or in addition to) discrete bands.
#'   Requires that smoothed probabilities are available in the diagnostics.
#'   Default is FALSE.
#' @param probability_style Character. How to display probabilities:
#'   \itemize{
#'     \item \code{"ribbon"}: Stacked ribbons showing probability of each state
#'     \item \code{"line"}: Line plot of probability for each state
#'     \item \code{"bands_alpha"}: Regime bands with alpha proportional to probability
#'   }
#'   Default is "ribbon".
#' @param show_bands Logical. If TRUE, show discrete regime bands. If 
#'   \code{show_probabilities = TRUE} and \code{show_bands = TRUE}, both are shown.
#'   Default is TRUE when \code{show_probabilities = FALSE}.
#' @param state_colors Optional named vector of colors for each state.
#' @param title Optional plot title.
#'
#' @return A ggplot object (invisibly).
#'
#' @details
#' This function creates a faceted plot with:
#' \itemize{
#'   \item Top panel: Original series with regime information
#'   \item Bottom panel(s): Bootstrap replicate(s) with their regime information
#' }
#'
#' When \code{show_probabilities = TRUE}, the smoothed state probabilities from
#' the fitted model (HMM or MS-VARMA-GARCH) are displayed. This shows the 
#' model's uncertainty about regime membership at each time point, which is
#' more informative than the hard Viterbi assignments.
#'
#' For bootstrap replicates, the "probabilities" shown are actually the 
#' probabilities from the original series at the source time points of each
#' resampled block. This reveals how the bootstrap combines observations from
#' different regime-certainty periods.
#'
#' @examples
#' \dontrun{
#' # Run HMM bootstrap with diagnostics
#' result <- tsbs(
#'   x = returns_data,
#'   bs_type = "hmm",
#'   num_states = 2,
#'   num_boots = 10,
#'   collect_diagnostics = TRUE,
#'   return_fit = TRUE
#' )
#'
#' # Plot with discrete regime bands
#' plot_regime_composition(result, returns_data)
#'
#' # Plot with smoothed probabilities as ribbons
#' plot_regime_composition(result, returns_data, show_probabilities = TRUE)
#'
#' # Plot with probability lines
#' plot_regime_composition(result, returns_data, 
#'                         show_probabilities = TRUE, 
#'                         probability_style = "line")
#' }
#'
#' @export
plot_regime_composition <- function(
    tsbs_result,
    original_data,
    replicate_idx = 1,
    series_idx = 1,
    show_prices = TRUE,
    initial_price = 100,
    show_probabilities = FALSE,
    probability_style = c("ribbon", "line", "bands_alpha"),
    show_bands = NULL,
    state_colors = NULL,
    title = NULL
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_regime_composition()")
  }
  
  probability_style <- match.arg(probability_style)
  
  ## Default show_bands: TRUE if not showing probabilities
  if (is.null(show_bands)) {
    show_bands <- !show_probabilities
  }
  
  ## Extract diagnostics object and other components
  if (inherits(tsbs_result, "tsbs_diagnostics")) {
    diagnostics <- tsbs_result
    bootstrap_series <- NULL
    smoothed_probs <- NULL
    fit_object <- NULL
  } else if (is.list(tsbs_result)) {
    diagnostics <- tsbs_result$diagnostics
    bootstrap_series <- tsbs_result$bootstrap_series
    smoothed_probs <- tsbs_result$smoothed_probabilities
    fit_object <- tsbs_result$fit
  } else {
    stop("tsbs_result must be a tsbs_diagnostics object or a tsbs() result with diagnostics")
  }
  
  ## Check for regime info
  if (is.null(diagnostics) || is.null(diagnostics$method_specific$regime_info)) {
    stop("No regime information found in diagnostics. ",
         "This plot is only available for HMM and MS-VARMA-GARCH bootstrap types.")
  }
  
  regime_info <- diagnostics$method_specific$regime_info
  
  ## Try to extract smoothed probabilities if requested but not directly available
  if (show_probabilities && is.null(smoothed_probs)) {
    smoothed_probs <- extract_smoothed_probabilities(fit_object, regime_info)
    if (is.null(smoothed_probs)) {
      warning("Smoothed probabilities not available. Falling back to discrete bands.")
      show_probabilities <- FALSE
      show_bands <- TRUE
    }
  }
  
  ## Validate replicate_idx
  replicate_idx <- as.integer(replicate_idx)
  if (any(replicate_idx < 1) || any(replicate_idx > diagnostics$meta$num_boots)) {
    stop("replicate_idx must be between 1 and ", diagnostics$meta$num_boots)
  }
  
  ## Prepare original data
  original_data <- as.matrix(original_data)
  if (series_idx > ncol(original_data)) {
    stop("series_idx (", series_idx, ") exceeds number of columns (", ncol(original_data), ")")
  }
  
  orig_values <- original_data[, series_idx]
  orig_states <- regime_info$original_states
  n_orig <- length(orig_values)
  
  ## Convert to prices if requested
  if (show_prices) {
    scale_factor <- if (max(abs(orig_values), na.rm = TRUE) > 1) 100 else 1
    orig_prices <- initial_price * exp(cumsum(orig_values / scale_factor))
  } else {
    orig_prices <- orig_values
  }
  
  ## Set up state colors
  num_states <- regime_info$num_states
  state_labels <- regime_info$state_labels
  
  if (is.null(state_colors)) {
    default_colors <- c(
      "lightblue", "salmon", "lightgreen", "plum", 
      "lightyellow", "lightcyan", "peachpuff", "lavender"
    )
    state_colors <- setNames(default_colors[seq_len(num_states)], state_labels)
  }
  
  ## Build plot data for series
  plot_data <- data.frame(
    Index = seq_len(n_orig),
    Value = orig_prices,
    Series = "Original",
    stringsAsFactors = FALSE
  )
  
  ## Build regime bands for original (if needed)
  regime_bands <- NULL
  if (show_bands) {
    regime_bands <- create_regime_bands_df(orig_states, "Original")
  }
  
  ## Build probability data for original (if needed)
  prob_data <- NULL
  if (show_probabilities && !is.null(smoothed_probs)) {
    prob_data <- create_probability_df(smoothed_probs, "Original", state_labels)
  }
  
  ## Add bootstrap replicates
  for (rep_idx in replicate_idx) {
    rep_states <- regime_info$replicate_states[[rep_idx]]
    
    if (is.null(rep_states)) {
      warning("No regime data for replicate ", rep_idx)
      next
    }
    
    ## Get bootstrap series values
    if (!is.null(bootstrap_series) && length(bootstrap_series) >= rep_idx) {
      boot_mat <- as.matrix(bootstrap_series[[rep_idx]])
      boot_values <- boot_mat[, min(series_idx, ncol(boot_mat))]
    } else {
      boot_values <- rep(NA, length(rep_states))
    }
    
    ## Convert to prices
    if (show_prices && !all(is.na(boot_values))) {
      scale_factor <- if (max(abs(boot_values), na.rm = TRUE) > 1) 100 else 1
      boot_prices <- initial_price * exp(cumsum(boot_values / scale_factor))
    } else {
      boot_prices <- boot_values
    }
    
    series_label <- paste("Bootstrap", rep_idx)
    
    plot_data <- rbind(plot_data, data.frame(
      Index = seq_along(boot_prices),
      Value = boot_prices,
      Series = series_label,
      stringsAsFactors = FALSE
    ))
    
    if (show_bands) {
      regime_bands <- rbind(regime_bands, create_regime_bands_df(rep_states, series_label))
    }
    
    ## For bootstrap, reconstruct probabilities from source blocks
    if (show_probabilities && !is.null(smoothed_probs)) {
      boot_probs <- reconstruct_bootstrap_probabilities(
        regime_info, rep_idx, smoothed_probs
      )
      if (!is.null(boot_probs)) {
        prob_data <- rbind(prob_data, 
                           create_probability_df(boot_probs, series_label, state_labels))
      }
    }
  }
  
  ## Map regime values to labels (for bands)
  if (!is.null(regime_bands)) {
    regime_bands$State_Label <- state_labels[regime_bands$State]
  }
  
  ## Set factor order for Series
  series_order <- c("Original", paste("Bootstrap", replicate_idx))
  plot_data$Series <- factor(plot_data$Series, levels = series_order)
  if (!is.null(regime_bands)) {
    regime_bands$Series <- factor(regime_bands$Series, levels = series_order)
  }
  if (!is.null(prob_data)) {
    prob_data$Series <- factor(prob_data$Series, levels = series_order)
  }
  
  ## Build plot
  if (is.null(title)) {
    prob_text <- if (show_probabilities) " with Smoothed Probabilities" else ""
    title <- paste0("Regime Composition: ", diagnostics$meta$bs_type, " Bootstrap", prob_text)
  }
  
  p <- ggplot2::ggplot()
  
  ## Add regime bands if requested
  if (show_bands && !is.null(regime_bands)) {
    p <- p + ggplot2::geom_rect(
      data = regime_bands,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = State_Label),
      alpha = 0.3
    )
  }
  
  ## Add probability visualization if requested
  if (show_probabilities && !is.null(prob_data)) {
    p <- add_probability_layer(p, prob_data, probability_style, state_colors, plot_data)
  }
  
  ## Add series lines
  p <- p + ggplot2::geom_line(
    data = plot_data,
    ggplot2::aes(x = Index, y = Value),
    color = "black",
    linewidth = 0.5
  )
  
  ## Facet by series
  p <- p + ggplot2::facet_wrap(~ Series, ncol = 1, scales = "free_y")
  
  ## Colors for bands
  if (show_bands) {
    p <- p + ggplot2::scale_fill_manual(values = state_colors, name = "Regime")
  }
  
  ## Labels and theme
  subtitle_text <- if (show_probabilities) {
    "Shading intensity shows regime probability; darker = higher confidence"
  } else {
    "Colored bands show regime membership; bootstrap resamples blocks within regimes"
  }
  
  p <- p + ggplot2::labs(
    title = title,
    subtitle = subtitle_text,
    x = "Time Index",
    y = if (show_prices) "Cumulative Value" else "Value"
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  print(p)
  invisible(p)
}


#' Extract Smoothed Probabilities from Fitted Model
#'
#' Attempts to extract smoothed state probabilities from a fitted HMM,
#' MS-VAR, or MS-VARMA-GARCH model object.
#'
#' @param fit_object Fitted model object (depmixS4, ms_var, or ms_varma_garch fit).
#' @param regime_info Regime info from diagnostics.
#'
#' @return Matrix of smoothed probabilities (n x num_states) or NULL.
#' @keywords internal
extract_smoothed_probabilities <- function(fit_object, regime_info) {
  
  if (is.null(fit_object)) {
    return(NULL)
  }
  
  ## Try depmixS4 (HMM)
  if (inherits(fit_object, "depmix.fitted")) {
    return(tryCatch({
      post <- depmixS4::posterior(fit_object)
      ## posterior() returns data frame with 'state' and probability columns
      prob_cols <- grep("^S", names(post), value = TRUE)
      as.matrix(post[, prob_cols])
    }, error = function(e) NULL))
  }
  
  ## Try MS-VARMA-GARCH (has smoothed_probabilities directly)
  if (is.list(fit_object) && "smoothed_probabilities" %in% names(fit_object)) {
    probs <- fit_object$smoothed_probabilities
    if (is.matrix(probs) || is.data.frame(probs)) {
      return(as.matrix(probs))
    }
  }
  
  ## Try MSGARCH package models
  if (inherits(fit_object, "MSGARCH_ML_FIT") || inherits(fit_object, "MSGARCH_MCMC_FIT")) {
    return(tryCatch({
      ## MSGARCH uses State() function for smoothed probabilities
      if (requireNamespace("MSGARCH", quietly = TRUE)) {
        state_probs <- MSGARCH::State(fit_object)
        if (is.list(state_probs) && "SmoothProb" %in% names(state_probs)) {
          return(state_probs$SmoothProb)
        }
      }
      NULL
    }, error = function(e) NULL))
  }
  
  ## Try MSwM package models
  if (inherits(fit_object, "MSM.lm") || inherits(fit_object, "MSM.glm")) {
    return(tryCatch({
      ## MSwM stores smoothed probabilities in @Fit@smoProb
      if (methods::hasMethod("slot", class(fit_object))) {
        fit_slot <- methods::slot(fit_object, "Fit")
        if ("smoProb" %in% slotNames(fit_slot)) {
          return(methods::slot(fit_slot, "smoProb"))
        }
      }
      NULL
    }, error = function(e) NULL))
  }
  
  NULL
}


#' Create Probability Data Frame for Plotting
#'
#' Converts a probability matrix to long format suitable for ggplot.
#'
#' @param probs Matrix of probabilities (n x num_states).
#' @param series_label Character label for this series.
#' @param state_labels Character vector of state names.
#'
#' @return Data frame with columns: Index, State, Probability, Series.
#' @keywords internal
create_probability_df <- function(probs, series_label, state_labels) {
  
  n <- nrow(probs)
  num_states <- ncol(probs)
  
  ## Long format
  df_list <- lapply(seq_len(num_states), function(s) {
    data.frame(
      Index = seq_len(n),
      State = state_labels[s],
      Probability = probs[, s],
      Series = series_label,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, df_list)
}


#' Reconstruct Bootstrap Probabilities from Source Blocks
#'
#' For a bootstrap replicate, reconstruct the "probabilities" by looking up
#' the original probabilities at the source time points of each resampled block.
#'
#' @param regime_info Regime info from diagnostics.
#' @param replicate_idx Which replicate.
#' @param original_probs Original smoothed probability matrix.
#'
#' @return Matrix of reconstructed probabilities or NULL.
#' @keywords internal
reconstruct_bootstrap_probabilities <- function(regime_info, replicate_idx, original_probs) {
  
  block_info <- regime_info$replicate_block_states[[replicate_idx]]
  
  ## Best case: we have per-observation source indices
  if (!is.null(block_info) && !is.null(block_info$source_indices)) {
    source_idx <- block_info$source_indices
    
    ## Validate indices
    n_orig <- nrow(original_probs)
    valid_idx <- source_idx >= 1 & source_idx <= n_orig
    
    if (all(valid_idx)) {
      return(original_probs[source_idx, , drop = FALSE])
    } else {
      ## Some indices out of range - handle gracefully
      n_boot <- length(source_idx)
      num_states <- ncol(original_probs)
      probs <- matrix(NA_real_, nrow = n_boot, ncol = num_states)
      probs[valid_idx, ] <- original_probs[source_idx[valid_idx], , drop = FALSE]
      return(probs)
    }
  }
  
  ## Fallback: use state assignments to create pseudo-probabilities (one-hot)
  rep_states <- regime_info$replicate_states[[replicate_idx]]
  if (is.null(rep_states)) return(NULL)
  
  n <- length(rep_states)
  num_states <- ncol(original_probs)
  
  ## Create one-hot encoding from states
  probs <- matrix(0, nrow = n, ncol = num_states)
  for (i in seq_len(n)) {
    if (!is.na(rep_states[i]) && rep_states[i] >= 1 && rep_states[i] <= num_states) {
      probs[i, rep_states[i]] <- 1
    }
  }
  
  probs
}


#' Add Probability Visualization Layer to Plot
#'
#' Adds the appropriate ggplot layer for probability visualization based on style.
#'
#' @param p Existing ggplot object.
#' @param prob_data Probability data frame.
#' @param style One of "ribbon", "line", "bands_alpha".
#' @param state_colors Named color vector.
#' @param plot_data Series plot data (for y-axis scaling).
#'
#' @return Updated ggplot object.
#' @keywords internal
add_probability_layer <- function(p, prob_data, style, state_colors, plot_data) {
  
  if (style == "line") {
    ## Line plot of probabilities (on secondary axis conceptually)
    ## Scale probabilities to fit in plot range
    p <- p + ggplot2::geom_line(
      data = prob_data,
      ggplot2::aes(x = Index, y = Probability, color = State),
      linewidth = 0.8,
      alpha = 0.7
    ) +
      ggplot2::scale_color_manual(values = state_colors, name = "State Probability")
    
  } else if (style == "ribbon") {
    ## Stacked ribbon showing cumulative probabilities
    ## Need to compute cumulative probabilities for stacking
    prob_data <- compute_stacked_probabilities(prob_data)
    
    p <- p + ggplot2::geom_ribbon(
      data = prob_data,
      ggplot2::aes(x = Index, ymin = ymin, ymax = ymax, fill = State),
      alpha = 0.5
    ) +
      ggplot2::scale_fill_manual(values = state_colors, name = "Regime")
    
  } else if (style == "bands_alpha") {
    ## Regime bands with alpha proportional to probability
    ## Convert to band format with probability-based alpha
    bands_alpha <- create_probability_bands(prob_data)
    
    p <- p + ggplot2::geom_rect(
      data = bands_alpha,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, 
                   fill = State, alpha = Probability)
    ) +
      ggplot2::scale_fill_manual(values = state_colors, name = "Regime") +
      ggplot2::scale_alpha_continuous(range = c(0.1, 0.5), guide = "none")
  }
  
  p
}


#' Compute Stacked Probabilities for Ribbon Plot
#'
#' Computes cumulative probabilities for stacked ribbon visualization.
#'
#' @param prob_data Probability data frame in long format.
#'
#' @return Data frame with ymin and ymax columns added.
#' @keywords internal
compute_stacked_probabilities <- function(prob_data) {
  
  ## Split by Series and Index, compute cumulative
  prob_data <- prob_data[order(prob_data$Series, prob_data$Index, prob_data$State), ]
  
  result_list <- lapply(split(prob_data, list(prob_data$Series, prob_data$Index)), function(df) {
    df <- df[order(df$State), ]
    cum_prob <- cumsum(df$Probability)
    df$ymax <- cum_prob
    df$ymin <- c(0, head(cum_prob, -1))
    df
  })
  
  do.call(rbind, result_list)
}


#' Create Probability-Weighted Bands
#'
#' Creates band data with probability for alpha mapping.
#'
#' @param prob_data Probability data frame.
#'
#' @return Data frame suitable for geom_rect with alpha mapping.
#' @keywords internal
create_probability_bands <- function(prob_data) {
  
  ## For each state, create bands where probability > threshold
  ## Simplified: just use point-wise probability
  
  ## Aggregate to get dominant state per time point with its probability
  prob_wide <- reshape(
    prob_data,
    direction = "wide",
    idvar = c("Index", "Series"),
    timevar = "State",
    v.names = "Probability"
  )
  
  ## For now, return simple per-index bands
  ## This is a simplified version - could be enhanced
  bands <- data.frame(
    xmin = prob_data$Index - 0.5,
    xmax = prob_data$Index + 0.5,
    State = prob_data$State,
    Probability = prob_data$Probability,
    Series = prob_data$Series,
    stringsAsFactors = FALSE
  )
  
  bands
}


#' Create Regime Bands Data Frame
#'
#' Helper function to convert a state sequence into a data frame of regime
#' bands suitable for geom_rect().
#'
#' @param states Integer vector of state assignments.
#' @param series_label Character string identifying this series.
#'
#' @return Data frame with columns: xmin, xmax, State, Series.
#' @keywords internal
create_regime_bands_df <- function(states, series_label) {
  
  rle_states <- rle(states)
  ends <- cumsum(rle_states$lengths)
  starts <- c(1, head(ends, -1) + 1)
  
  data.frame(
    xmin = starts,
    xmax = ends,
    State = rle_states$values,
    Series = series_label,
    stringsAsFactors = FALSE
  )
}


#' Summary of Regime Bootstrap Diagnostics
#'
#' Prints a summary of regime-related diagnostics including state distributions
#' in the original data and across bootstrap replicates.
#'
#' @param diagnostics A \code{tsbs_diagnostics} object with regime information.
#'
#' @return Invisibly returns a list with summary statistics.
#' @keywords internal
#' @export
summarize_regime_diagnostics <- function(diagnostics) {
  
  if (is.null(diagnostics$method_specific$regime_info)) {
    cat("No regime information available.\n")
    return(invisible(NULL))
  }
  
  regime_info <- diagnostics$method_specific$regime_info
  
  cat("=== Regime Bootstrap Diagnostics ===\n\n")
  
  cat("Number of states:", regime_info$num_states, "\n")
  cat("State labels:", paste(regime_info$state_labels, collapse = ", "), "\n\n")
  
  cat("Original series state distribution:\n")
  orig_counts <- regime_info$original_state_counts
  orig_props <- prop.table(orig_counts)
  for (i in seq_along(orig_counts)) {
    cat(sprintf("  %s: %d observations (%.1f%%)\n",
                regime_info$state_labels[i],
                orig_counts[i],
                orig_props[i] * 100))
  }
  
  ## Compute bootstrap state distribution summary
  replicate_state_props <- lapply(regime_info$replicate_states, function(states) {
    if (is.null(states)) return(NULL)
    table(factor(states, levels = seq_len(regime_info$num_states))) / length(states)
  })
  replicate_state_props <- replicate_state_props[!sapply(replicate_state_props, is.null)]
  
  if (length(replicate_state_props) > 0) {
    cat("\nBootstrap replicate state distribution (mean across replicates):\n")
    mean_props <- Reduce(`+`, replicate_state_props) / length(replicate_state_props)
    for (i in seq_along(mean_props)) {
      cat(sprintf("  %s: %.1f%% (original: %.1f%%)\n",
                  regime_info$state_labels[i],
                  mean_props[i] * 100,
                  orig_props[i] * 100))
    }
  }
  
  cat("\n")
  
  invisible(list(
    original_counts = orig_counts,
    original_props = orig_props,
    bootstrap_mean_props = if (exists("mean_props")) mean_props else NULL
  ))
}


# Block Bootstrap Diagnostics for Moving and Stationary Bootstrap ==============
# These functions provide diagnostic tracking for non-regime-based bootstrap
# methods (Moving Block Bootstrap and Stationary Bootstrap).
#
# Key diagnostics for block bootstrap:
# - Block lengths (fixed for moving, random for stationary)
# - Block starting positions (where each block was sampled from)
# - Coverage analysis (which parts of original series are represented)
# - Overlap patterns

#' Block Bootstrap with Diagnostics
#'
#' Performs Moving Block Bootstrap or Stationary Bootstrap while collecting
#' detailed diagnostic information about block composition.
#'
#' @param x Numeric matrix with rows as time points, columns as variables.
#' @param n_boot Integer. Desired length of bootstrap series. If NULL, uses
#'   length of original series.
#' @param block_length Integer. Block length for moving bootstrap, or expected
#'   block length for stationary bootstrap. If NULL, computed automatically.
#' @param bs_type Character. Either "moving" or "stationary".
#' @param block_type Character. One of "overlapping", "non-overlapping", or
#'   "tapered".
#' @param num_boots Integer. Number of bootstrap replicates.
#' @param p Numeric in (0,1). Probability parameter for geometric distribution
#'   in stationary bootstrap. Default 0.1.
#' @param collect_diagnostics Logical. If TRUE, collect block-level diagnostics.
#' @param taper_type Character. Taper window type if block_type = "tapered".
#' @param tukey_alpha Numeric. Alpha parameter for Tukey window.
#'
#' @return If collect_diagnostics = FALSE, a list of bootstrap matrices.
#'   If collect_diagnostics = TRUE, a list with:
#'   \describe{
#'     \item{bootstrap_series}{List of bootstrap matrices}
#'     \item{diagnostics}{A tsbs_diagnostics object with block information}
#'   }
#'
#' @details
#' This function provides an R implementation of block bootstrap that tracks
#' diagnostic information. For each bootstrap replicate, it records:
#' \itemize{
#'   \item Block lengths used
#'   \item Starting positions in original series
#'   \item Source indices mapping each bootstrap observation to original
#' }
#'
#' For the **Moving Block Bootstrap**, blocks have fixed length and are sampled
#' uniformly from all possible starting positions.
#'
#' For the **Stationary Bootstrap**, block lengths are drawn from a geometric
#' distribution with parameter p, giving expected block length 1/p.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(500), ncol = 2)
#'
#' # Moving block bootstrap with diagnostics
#' result <- blockBootstrap_with_diagnostics(
#'   x, block_length = 10, bs_type = "moving",
#'   num_boots = 50, collect_diagnostics = TRUE
#' )
#'
#' # View block length distribution
#' summary(result$diagnostics)
#'
#' # Plot coverage
#' plot_block_coverage(result$diagnostics)
#' }
#'
#' @export
blockBootstrap_with_diagnostics <- function(
    x,
    n_boot = NULL,
    block_length = NULL,
    bs_type = c("moving", "stationary"),
    block_type = c("overlapping", "non-overlapping"),
    num_boots = 100L,
    p = 0.1,
    collect_diagnostics = TRUE,
    taper_type = "cosine",
    tukey_alpha = 0.5
) {
  
  bs_type <- match.arg(bs_type)
  block_type <- match.arg(block_type)
  
  ## ---- Data Preparation ----
  x <- as.matrix(x)
  n <- nrow(x)
  k <- ncol(x)
  
  ## Set n_boot if not specified
  if (is.null(n_boot)) {
    n_boot <- n
  }
  
  ## Compute default block length if not specified
  if (is.null(block_length)) {
    block_length <- compute_default_block_length_r(x)
  }
  
  ## For stationary bootstrap, block_length is expected length (1/p)
  if (bs_type == "stationary" && is.null(p)) {
    p <- 1 / block_length
  }
  
  ## ---- Initialize Diagnostics ----
  diagnostics <- NULL
  if (collect_diagnostics) {
    diagnostics <- create_bootstrap_diagnostics(
      bs_type = bs_type,
      n_original = n,
      n_vars = k,
      num_boots = num_boots
    )
    
    ## Record original series stats
    diagnostics <- record_original_stats(diagnostics, x)
    
    ## Store config
    diagnostics$config <- list(
      block_length = block_length,
      block_type = block_type,
      n_boot = n_boot,
      p = if (bs_type == "stationary") p else NULL
    )
    
    ## Initialize block-specific tracking
    diagnostics$method_specific$block_info <- list(
      ## Per-replicate: list of block details
      replicate_blocks = vector("list", num_boots),
      ## Source indices for each replicate
      replicate_source_indices = vector("list", num_boots)
    )
  }
  
  ## ---- Generate Bootstrap Samples ----
  bootstrap_samples <- vector("list", num_boots)
  
  for (b in seq_len(num_boots)) {
    
    ## Storage for this replicate
    sample_matrix <- matrix(0, nrow = n_boot, ncol = k)
    pos <- 1
    
    ## Track blocks for diagnostics
    block_lengths_used <- integer(0)
    block_starts <- integer(0)
    source_indices <- integer(0)
    
    while (pos <= n_boot) {
      
      ## Determine block length
      if (bs_type == "moving") {
        current_block_len <- block_length
      } else {
        ## Stationary: geometric distribution
        current_block_len <- rgeom(1, p) + 1
        ## Cap at reasonable maximum
        max_len <- min(ceiling(qgeom(0.99, p) + 1), floor(0.1 * n))
        current_block_len <- min(current_block_len, max(1, max_len))
      }
      
      ## Determine starting position
      if (block_type == "overlapping") {
        max_start <- n - current_block_len + 1
        if (max_start < 1) max_start <- 1
        start_idx <- sample.int(max_start, 1)
      } else {
        ## Non-overlapping: sample from block boundaries
        num_complete_blocks <- n %/% current_block_len
        if (num_complete_blocks < 1) num_complete_blocks <- 1
        block_num <- sample.int(num_complete_blocks, 1)
        start_idx <- (block_num - 1) * current_block_len + 1
      }
      
      ## How many rows to copy (might be truncated at end)
      rows_to_copy <- min(current_block_len, n_boot - pos + 1)
      
      ## Copy block (with wraparound if needed)
      for (i in seq_len(rows_to_copy)) {
        idx <- ((start_idx - 1 + i - 1) %% n) + 1
        sample_matrix[pos + i - 1, ] <- x[idx, ]
        source_indices <- c(source_indices, idx)
      }
      
      ## Record block info
      block_lengths_used <- c(block_lengths_used, rows_to_copy)
      block_starts <- c(block_starts, start_idx)
      
      pos <- pos + rows_to_copy
    }
    
    bootstrap_samples[[b]] <- sample_matrix
    
    ## Record diagnostics for this replicate
    if (collect_diagnostics) {
      diagnostics <- record_blocks(
        diagnostics,
        replicate_idx = b,
        block_lengths = block_lengths_used,
        start_positions = block_starts,
        block_type = block_type
      )
      
      diagnostics$method_specific$block_info$replicate_source_indices[[b]] <- source_indices
      
      diagnostics <- record_replicate_stats(
        diagnostics,
        replicate_idx = b,
        bootstrap_matrix = sample_matrix
      )
    }
  }
  
  ## ---- Return Results ----
  if (!collect_diagnostics) {
    return(bootstrap_samples)
  }
  
  list(
    bootstrap_series = bootstrap_samples,
    diagnostics = diagnostics
  )
}


#' Compute Default Block Length (R version)
#'
#' Computes a default block length based on average lag-1 autocorrelation.
#'
#' @param x Numeric matrix.
#' @return Integer block length.
#' @keywords internal
compute_default_block_length_r <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  k <- ncol(x)
  
  ## Compute lag-1 autocorrelation for each column
  ac1 <- numeric(k)
  for (j in seq_len(k)) {
    ac1[j] <- tryCatch({
      abs(acf(x[, j], lag.max = 1, plot = FALSE)$acf[2, 1, 1])
    }, error = function(e) 0)
  }
  
  rho1 <- mean(ac1, na.rm = TRUE)
  if (is.na(rho1) || rho1 >= 1) rho1 <- 0.5
  
  candidate <- floor(10 / (1 - rho1))
  max(5L, min(candidate, as.integer(sqrt(n))))
}


#' Plot Block Coverage for Bootstrap Diagnostics
#'
#' Visualizes which parts of the original series are sampled across bootstrap
#' replicates, showing coverage patterns and potential gaps.
#'
#' @param diagnostics A tsbs_diagnostics object from block bootstrap.
#' @param type Character. Type of plot:
#'   \itemize{
#'     \item "heatmap": Heatmap showing coverage intensity at each time point
#'     \item "histogram": Histogram of block starting positions
#'     \item "blocks": Visual representation of blocks sampled per replicate
#'   }
#' @param max_replicates Integer. Maximum number of replicates to show in
#'   "blocks" plot. Default 20.
#'
#' @return A ggplot object (invisibly).
#'
#' @examples
#' \dontrun{
#' result <- blockBootstrap_with_diagnostics(x, bs_type = "moving",
#'                                            collect_diagnostics = TRUE)
#' plot_block_coverage(result$diagnostics, type = "heatmap")
#' plot_block_coverage(result$diagnostics, type = "histogram")
#' }
#'
#' @export
plot_block_coverage <- function(
    diagnostics,
    type = c("heatmap", "histogram", "blocks"),
    max_replicates = 20
) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_block_coverage()")
  }
  
  type <- match.arg(type)
  
  n_original <- diagnostics$meta$n_original
  num_boots <- diagnostics$meta$num_boots
  bs_type <- diagnostics$meta$bs_type
  
  if (type == "histogram") {
    ## Histogram of block starting positions
    starts <- diagnostics$blocks$all_start_positions
    
    df <- data.frame(start_pos = starts)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = start_pos)) +
      ggplot2::geom_histogram(bins = min(50, n_original / 5),
                              fill = "steelblue", color = "white", alpha = 0.7) +
      ggplot2::labs(
        title = paste("Block Starting Positions:", bs_type, "Bootstrap"),
        subtitle = paste(length(starts), "blocks across", num_boots, "replicates"),
        x = "Starting Position in Original Series",
        y = "Count"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::geom_vline(xintercept = c(1, n_original), 
                          linetype = "dashed", color = "red", alpha = 0.5)
    
  } else if (type == "heatmap") {
    ## Coverage heatmap: how often each time point is sampled
    source_indices <- diagnostics$method_specific$block_info$replicate_source_indices
    
    ## Count coverage at each position
    coverage <- tabulate(unlist(source_indices), nbins = n_original)
    
    df <- data.frame(
      position = seq_len(n_original),
      count = coverage
    )
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = position, y = 1, fill = count)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "Times\nSampled") +
      ggplot2::labs(
        title = paste("Coverage Heatmap:", bs_type, "Bootstrap"),
        subtitle = paste("How often each time point appears across", num_boots, "replicates"),
        x = "Position in Original Series",
        y = ""
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
    
  } else if (type == "blocks") {
    ## Visual representation of blocks per replicate
    block_data <- diagnostics$blocks$replicate_blocks
    
    n_show <- min(max_replicates, num_boots)
    
    plot_df <- do.call(rbind, lapply(seq_len(n_show), function(rep_idx) {
      blocks <- block_data[[rep_idx]]
      if (is.null(blocks) || nrow(blocks) == 0) return(NULL)
      
      data.frame(
        replicate = rep_idx,
        block_num = blocks$block_num,
        start = blocks$start_pos,
        end = blocks$start_pos + blocks$length - 1,
        length = blocks$length
      )
    }))
    
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      stop("No block data available for plotting")
    }
    
    p <- ggplot2::ggplot(plot_df) +
      ggplot2::geom_segment(
        ggplot2::aes(x = start, xend = end, 
                     y = replicate, yend = replicate,
                     color = factor(block_num %% 5)),
        linewidth = 2, alpha = 0.7
      ) +
      ggplot2::scale_y_reverse() +
      ggplot2::scale_color_brewer(palette = "Set2", guide = "none") +
      ggplot2::labs(
        title = paste("Block Composition:", bs_type, "Bootstrap"),
        subtitle = paste("First", n_show, "replicates; segments show sampled blocks"),
        x = "Position in Original Series",
        y = "Bootstrap Replicate"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::geom_vline(xintercept = c(1, n_original),
                          linetype = "dashed", color = "gray50", alpha = 0.5)
  }
  
  print(p)
  invisible(p)
}


#' Plot Block Length Distribution
#'
#' Visualizes the distribution of block lengths used in bootstrap replicates.
#' Particularly useful for Stationary Bootstrap where lengths are random.
#'
#' @param diagnostics A tsbs_diagnostics object from block bootstrap.
#' @param show_expected Logical. If TRUE and bs_type is "stationary", show
#'   the expected geometric distribution. Default TRUE.
#'
#' @return A ggplot object (invisibly).
#'
#' @export
plot_block_lengths <- function(diagnostics, show_expected = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_block_lengths()")
  }
  
  bs_type <- diagnostics$meta$bs_type
  block_lengths <- diagnostics$blocks$all_block_lengths
  
  df <- data.frame(length = block_lengths)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = length)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = max(10, length(unique(block_lengths))),
                            fill = "steelblue", color = "white", alpha = 0.7) +
    ggplot2::labs(
      title = paste("Block Length Distribution:", bs_type, "Bootstrap"),
      subtitle = paste("n =", length(block_lengths), "blocks"),
      x = "Block Length",
      y = "Density"
    ) +
    ggplot2::theme_minimal()
  
  ## Add expected distribution for stationary bootstrap
  if (show_expected && bs_type == "stationary" && !is.null(diagnostics$config$p)) {
    prob <- diagnostics$config$p
    max_len <- max(block_lengths)
    expected_df <- data.frame(
      length = seq_len(max_len),
      density = dgeom(seq_len(max_len) - 1, prob)
    )
    
    p <- p + ggplot2::geom_line(
      data = expected_df,
      ggplot2::aes(x = length, y = density),
      color = "red", linewidth = 1, linetype = "dashed"
    ) +
      ggplot2::annotate("text", x = max_len * 0.7, y = max(expected_df$density) * 0.8,
                        label = paste("Expected: Geom(p =", round(prob, 3), ")"),
                        color = "red", size = 3.5)
  }
  
  ## Add mean line
  mean_len <- mean(block_lengths)
  p <- p + ggplot2::geom_vline(xintercept = mean_len, 
                               color = "darkgreen", linewidth = 1) +
    ggplot2::annotate("text", x = mean_len, y = 0,
                      label = paste("Mean:", round(mean_len, 1)),
                      color = "darkgreen", vjust = -0.5, hjust = -0.1, size = 3.5)
  
  print(p)
  invisible(p)
}


#' Summary Method for Block Bootstrap Diagnostics
#'
#' Prints a summary of block bootstrap diagnostics.
#'
#' @param diagnostics A tsbs_diagnostics object.
#'
#' @return Invisibly returns a list with summary statistics.
#' @keywords internal
#' @export
summarize_block_diagnostics <- function(diagnostics) {
  
  cat("=== Block Bootstrap Diagnostics ===\n\n")
  
  cat("Bootstrap type:", diagnostics$meta$bs_type, "\n")
  cat("Original series length:", diagnostics$meta$n_original, "\n")
  cat("Number of variables:", diagnostics$meta$n_vars, "\n")
  cat("Number of replicates:", diagnostics$meta$num_boots, "\n\n")
  
  if (!is.null(diagnostics$config$block_length)) {
    cat("Configured block length:", diagnostics$config$block_length, "\n")
  }
  if (!is.null(diagnostics$config$p)) {
    cat("Geometric parameter p:", diagnostics$config$p, "\n")
    cat("Expected block length:", round(1/diagnostics$config$p, 1), "\n")
  }
  cat("Block type:", diagnostics$config$block_type, "\n\n")
  
  ## Block length statistics
  block_lengths <- diagnostics$blocks$all_block_lengths
  cat("Block length statistics:\n")
  cat("  Total blocks sampled:", length(block_lengths), "\n")
  cat("  Mean length:", round(mean(block_lengths), 2), "\n")
  cat("  SD length:", round(sd(block_lengths), 2), "\n")
  cat("  Min/Max:", min(block_lengths), "/", max(block_lengths), "\n\n")
  
  ## Coverage statistics
  if (!is.null(diagnostics$method_specific$block_info$replicate_source_indices)) {
    source_indices <- unlist(diagnostics$method_specific$block_info$replicate_source_indices)
    coverage <- tabulate(source_indices, nbins = diagnostics$meta$n_original)
    
    cat("Coverage statistics:\n")
    cat("  Time points never sampled:", sum(coverage == 0), "\n")
    cat("  Mean times sampled:", round(mean(coverage), 2), "\n")
    cat("  Max times sampled:", max(coverage), "\n")
    cat("  Coverage uniformity (CV):", round(sd(coverage) / mean(coverage), 3), "\n")
  }
  
  cat("\n")
  
  invisible(list(
    block_lengths = block_lengths,
    coverage = if (exists("coverage")) coverage else NULL
  ))
}
