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
#> [1] "ms-varma-garch_diagnostics_report.Rmd"
#> [2] "ms-varma-garch_inference_demo.Rmd"    
#> [3] "multivar_garch_comparison.Rmd"        
#> [4] "portfolio-optimization.Rmd"           
```
