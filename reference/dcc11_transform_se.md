# Transform Standard Errors Between Parameterizations

Transform Standard Errors Between Parameterizations

## Usage

``` r
dcc11_transform_se(se_result, to_reparam = FALSE)
```

## Arguments

- se_result:

  Result from dcc11_standard_errors() or dcc11_hessian()

- to_reparam:

  Logical: transform TO (psi, phi)?

## Value

Updated se_result with transformed SEs
