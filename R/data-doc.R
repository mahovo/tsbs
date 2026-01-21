#' ETF Portfolio Returns Data
#'
#' Daily log returns (in percent) for a diversified ETF portfolio covering
#' major asset classes. This data is used in the portfolio optimization
#' vignette to demonstrate the tsbs bootstrap methods.
#'
#' @format A matrix with approximately 1,750 rows (trading days) and 5 columns:
#' \describe{
#'   \item{SPY}{S&P 500 ETF - US large-cap equity}
#'   \item{EFA}{iShares MSCI EAFE ETF - International developed markets equity}
#'   \item{BND}{Vanguard Total Bond Market ETF - US aggregate bonds}
#'   \item{GLD}{SPDR Gold Shares ETF - Gold}
#'   \item{VNQ}{Vanguard Real Estate ETF - US REITs}
#' }
#'
#' The matrix has the following attributes:
#' \describe{
#'   \item{dates}{Date vector of trading days}
#'   \item{symbols}{Character vector of ticker symbols}
#'   \item{symbol_names}{Character vector of asset class descriptions}
#'   \item{source}{Data source (Yahoo Finance)}
#'   \item{download_date}{Date when data was downloaded}
#'   \item{description}{Brief description of the dataset}
#' }
#'
#' @details
#' Returns are computed as \eqn{r_t = 100 \times \log(P_t / P_{t-1})} where
#' \eqn{P_t} is the adjusted closing price. The data spans from January 2018
#' to December 2024, covering various market conditions including the
#' COVID-19 crash and recovery, the 2022 bear market, and subsequent rally.
#'
#' @source Yahoo Finance via the \pkg{quantmod} package.
#'
#' @examples
#' # Load the data
#' data(etf_returns)
#'
#' # Check dimensions
#' dim(etf_returns)
#'
#' # View date range
#' dates <- attr(etf_returns, "dates")
#' cat("Date range:", as.character(dates[1]), "to", 
#'     as.character(dates[length(dates)]), "\n")
#'
#' # Compute annualized statistics
#' ann_return <- colMeans(etf_returns) * 252
#' ann_vol <- apply(etf_returns, 2, sd) * sqrt(252)
#' print(round(cbind(Return = ann_return, Vol = ann_vol), 2))
#'
#' # Correlation matrix
#' print(round(cor(etf_returns), 2))
#'
#' @seealso \code{\link{etf_prices}} for the corresponding price data.
#'
"etf_returns"


#' ETF Portfolio Price Data
#'
#' Daily adjusted closing prices for a diversified ETF portfolio. This is the
#' price data corresponding to \code{\link{etf_returns}}, useful for plotting
#' and visualization.
#'
#' @format A matrix with approximately 1,750 rows (trading days) and 5 columns:
#' \describe{
#'   \item{SPY}{S&P 500 ETF adjusted close price}
#'   \item{EFA}{iShares MSCI EAFE ETF adjusted close price}
#'   \item{BND}{Vanguard Total Bond Market ETF adjusted close price}
#'   \item{GLD}{SPDR Gold Shares ETF adjusted close price}
#'   \item{VNQ}{Vanguard Real Estate ETF adjusted close price}
#' }
#'
#' The matrix has the following attributes:
#' \describe{
#'   \item{dates}{Date vector of trading days}
#'   \item{symbols}{Character vector of ticker symbols}
#' }
#'
#' @source Yahoo Finance via the \pkg{quantmod} package.
#'
#' @examples
#' # Load the data
#' data(etf_prices)
#'
#' # Plot SPY prices
#' dates <- attr(etf_prices, "dates")
#' plot(dates, etf_prices[, "SPY"], type = "l",
#'      main = "SPY Price History", xlab = "Date", ylab = "Price ($)")
#'
#' @seealso \code{\link{etf_returns}} for the corresponding returns data.
#'
"etf_prices"
