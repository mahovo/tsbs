# List available demonstration Rmd files

These demonstrations are computationally intensive and not built as
vignettes. Use
[`get_demo_path()`](https://mahovo.github.io/tsbs/reference/get_demo_path.md)
to get the full path to a specific demo.

## Usage

``` r
list_package_demos()
```

## Value

Character vector of available demo Rmd file names

## Examples

``` r
list_package_demos()
#> [1] "hmm_bootstrap_demo.Rmd"               
#> [2] "ms-varma-garch_diagnostics_report.Rmd"
#> [3] "ms-varma-garch_inference_demo.Rmd"    
#> [4] "multivar_garch_comparison.Rmd"        
#> [5] "portfolio-optimization.Rmd"           
```
