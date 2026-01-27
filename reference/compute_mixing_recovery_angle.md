# Compute Mixing Matrix Recovery Angle

Compute angle between true and estimated mixing matrix subspaces. Uses
the Amari index or similar metric for ICA recovery assessment.

## Usage

``` r
compute_mixing_recovery_angle(A_true, A_est)
```

## Arguments

- A_true:

  True mixing matrix

- A_est:

  Estimated mixing matrix

## Value

Scalar angle metric (0 = perfect recovery, larger = worse)
