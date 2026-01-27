# ETF Portfolio Returns Data

Daily log returns (in percent) for a diversified ETF portfolio covering
major asset classes. This data is used in the portfolio optimization
vignette to demonstrate the tsbs bootstrap methods.

## Usage

``` r
etf_returns
```

## Format

A matrix with approximately 1,750 rows (trading days) and 5 columns:

- SPY:

  S&P 500 ETF - US large-cap equity

- EFA:

  iShares MSCI EAFE ETF - International developed markets equity

- BND:

  Vanguard Total Bond Market ETF - US aggregate bonds

- GLD:

  SPDR Gold Shares ETF - Gold

- VNQ:

  Vanguard Real Estate ETF - US REITs

The matrix has the following attributes:

- dates:

  Date vector of trading days

- symbols:

  Character vector of ticker symbols

- symbol_names:

  Character vector of asset class descriptions

- source:

  Data source (Yahoo Finance)

- download_date:

  Date when data was downloaded

- description:

  Brief description of the dataset

## Source

Yahoo Finance via the quantmod package.

## Details

Returns are computed as \\r_t = 100 \times \log(P_t / P\_{t-1})\\ where
\\P_t\\ is the adjusted closing price. The data spans from January 2018
to December 2024, covering various market conditions including the
COVID-19 crash and recovery, the 2022 bear market, and subsequent rally.

## See also

[`etf_prices`](https://mahovo.github.io/tsbs/reference/etf_prices.md)
for the corresponding price data.

## Examples

``` r
# Load the data
data(etf_returns)

# Check dimensions
dim(etf_returns)
#> [1] 1759    5

# View date range
dates <- attr(etf_returns, "dates")
cat("Date range:", as.character(dates[1]), "to", 
    as.character(dates[length(dates)]), "\n")
#> Date range: 2018-01-03 to 2024-12-30 

# Compute annualized statistics
ann_return <- colMeans(etf_returns) * 252
ann_vol <- apply(etf_returns, 2, sd) * sqrt(252)
print(round(cbind(Return = ann_return, Vol = ann_vol), 2))
#>     Return   Vol
#> SPY  12.84 19.54
#> EFA   3.90 18.27
#> BND   1.06  6.16
#> GLD   9.37 14.33
#> VNQ   4.79 23.12

# Correlation matrix
print(round(cor(etf_returns), 2))
#>      SPY  EFA  BND  GLD  VNQ
#> SPY 1.00 0.87 0.16 0.11 0.76
#> EFA 0.87 1.00 0.20 0.21 0.71
#> BND 0.16 0.20 1.00 0.36 0.29
#> GLD 0.11 0.21 0.36 1.00 0.16
#> VNQ 0.76 0.71 0.29 0.16 1.00
```
