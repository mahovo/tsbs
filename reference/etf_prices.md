# ETF Portfolio Price Data

Daily adjusted closing prices for a diversified ETF portfolio. This is
the price data corresponding to
[`etf_returns`](https://mahovo.github.io/tsbs/reference/etf_returns.md),
useful for plotting and visualization.

## Usage

``` r
etf_prices
```

## Format

A matrix with approximately 1,750 rows (trading days) and 5 columns:

- SPY:

  S&P 500 ETF adjusted close price

- EFA:

  iShares MSCI EAFE ETF adjusted close price

- BND:

  Vanguard Total Bond Market ETF adjusted close price

- GLD:

  SPDR Gold Shares ETF adjusted close price

- VNQ:

  Vanguard Real Estate ETF adjusted close price

The matrix has the following attributes:

- dates:

  Date vector of trading days

- symbols:

  Character vector of ticker symbols

## Source

Yahoo Finance via the quantmod package.

## See also

[`etf_returns`](https://mahovo.github.io/tsbs/reference/etf_returns.md)
for the corresponding returns data.

## Examples

``` r
# Load the data
data(etf_prices)

# Plot SPY prices
dates <- attr(etf_prices, "dates")
plot(dates, etf_prices[, "SPY"], type = "l",
     main = "SPY Price History", xlab = "Date", ylab = "Price ($)")

```
