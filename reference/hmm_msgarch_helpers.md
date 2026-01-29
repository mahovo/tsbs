# MSGARCH Helper Functions for HMM Bootstrap

Internal functions for regime-switching bootstrap using MSGARCH.

This module provides helper functions that wrap the MSGARCH package to
enable semi-parametric bootstrap for time series with regime-switching
behavior and non-Gaussian conditional distributions (skew Student-t,
Student-t, GED, etc.).

## Details

### Literature Background

The implementation combines several established techniques:

- **Markov-switching GARCH**: Haas, Mittnik & Paolella (2004),
  implemented via the MSGARCH package (Ardia et al., 2019)

- **Skew Student-t distribution**: Fernández & Steel (1998)
  transformation

- **Semi-parametric bootstrap**: State sequences are simulated from the
  fitted Markov chain (parametric), while innovations are resampled from
  empirical state-specific pools (nonparametric)

## References

Ardia, D., Bluteau, K., Boudt, K., Catania, L., & Trottier, D.-A.
(2019). Markov-Switching GARCH Models in R: The MSGARCH Package. Journal
of Statistical Software, 91(4), 1-38. doi:10.18637/jss.v091.i04

Fernández, C., & Steel, M. F. (1998). On Bayesian modeling of fat tails
and skewness. Journal of the American Statistical Association, 93(441),
359-371.

Haas, M., Mittnik, S., & Paolella, M. S. (2004). A New Approach to
Markov- Switching GARCH Models. Journal of Financial Econometrics, 2,
493-530.

Hamilton, J. D. (1989). A New Approach to the Economic Analysis of
Nonstationary Time Series and the Business Cycle. Econometrica, 57(2),
357-384.
