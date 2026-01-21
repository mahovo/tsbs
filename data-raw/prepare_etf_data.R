## =============================================================================
## Data Preparation Script for tsbs Package
## =============================================================================
## This script downloads ETF data and saves it for inclusion in the tsbs package.
## Run this script to update the package data.
##
## The data will be saved as:
##   - data/etf_returns.rda  (R data file for lazy loading)
##   - data-raw/etf_returns.R (this script, for reproducibility)
## =============================================================================

library(quantmod)

## -----------------------------------------------------------------------------
## Download ETF Data
## -----------------------------------------------------------------------------

symbols <- c("SPY", "EFA", "BND", "GLD", "VNQ")
symbol_names <- c("US Equity", "Intl Equity", "US Bonds", "Gold", "REITs")

## Use a fixed date range for reproducibility
start_date <- "2018-01-01"
end_date <- "2024-12-31"

cat("Downloading ETF data for tsbs package...\n")
cat("Symbols:", paste(symbols, collapse = ", "), "\n")
cat("Date range:", start_date, "to", end_date, "\n\n")

prices_list <- lapply(symbols, function(sym) {
  cat("  Downloading", sym, "...")
  result <- tryCatch({
    data <- getSymbols(sym, src = "yahoo", from = start_date, to = end_date, 
                       auto.assign = FALSE)
    cat(" OK\n")
    data
  }, error = function(e) {
    cat(" FAILED:", e$message, "\n")
    NULL
  })
  result
})

names(prices_list) <- symbols

## Check for failures
valid_idx <- !sapply(prices_list, is.null)
if (sum(valid_idx) < length(symbols)) {
  warning("Some symbols failed to download: ", 
          paste(symbols[!valid_idx], collapse = ", "))
}

## -----------------------------------------------------------------------------
## Process Data
## -----------------------------------------------------------------------------

## Extract adjusted close prices
adj_close <- do.call(merge, lapply(prices_list[valid_idx], Ad))
colnames(adj_close) <- symbols[valid_idx]
adj_close <- na.omit(adj_close)

## Compute log returns (in percent)
returns_xts <- diff(log(adj_close)) * 100
returns_xts <- na.omit(returns_xts)

## Convert to matrix for the package
etf_returns <- as.matrix(coredata(returns_xts))
colnames(etf_returns) <- symbols[valid_idx]

## Store dates as an attribute
attr(etf_returns, "dates") <- index(returns_xts)
attr(etf_returns, "symbols") <- symbols[valid_idx]
attr(etf_returns, "symbol_names") <- symbol_names[valid_idx]
attr(etf_returns, "source") <- "Yahoo Finance via quantmod"
attr(etf_returns, "download_date") <- Sys.Date()
attr(etf_returns, "description") <- paste(
 "Daily log returns (in percent) for a diversified ETF portfolio.",
 "SPY = S&P 500, EFA = International Developed Markets,",
 "BND = US Aggregate Bonds, GLD = Gold, VNQ = US REITs."
)

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------

cat("\n")
cat("=== Data Summary ===\n")
cat("Observations:", nrow(etf_returns), "\n")
cat("Variables:", ncol(etf_returns), "\n")
cat("Date range:", as.character(attr(etf_returns, "dates")[1]), "to",
    as.character(attr(etf_returns, "dates")[nrow(etf_returns)]), "\n")
cat("\nAnnualized statistics:\n")
ann_stats <- rbind(
  "Return (%)" = colMeans(etf_returns) * 252,
  "Vol (%)" = apply(etf_returns, 2, sd) * sqrt(252),
  "Sharpe" = colMeans(etf_returns) / apply(etf_returns, 2, sd) * sqrt(252)
)
print(round(ann_stats, 2))

cat("\nCorrelation matrix:\n")
print(round(cor(etf_returns), 2))

## -----------------------------------------------------------------------------
## Save Data
## -----------------------------------------------------------------------------

## Save to data/ directory (for lazy loading in package)
save(etf_returns, file = "data/etf_returns.rda", compress = "xz")
cat("\nData saved to: data/etf_returns.rda\n")

## Also save the raw prices for the bootstrap illustration
etf_prices <- as.matrix(coredata(adj_close))
colnames(etf_prices) <- symbols[valid_idx]
attr(etf_prices, "dates") <- index(adj_close)
attr(etf_prices, "symbols") <- symbols[valid_idx]

save(etf_prices, file = "data/etf_prices.rda", compress = "xz")
cat("Data saved to: data/etf_prices.rda\n")

cat("\nDone!\n")
