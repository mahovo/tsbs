# Get path to a specific demo Rmd file

Get path to a specific demo Rmd file

## Usage

``` r
get_demo_path(name)
```

## Arguments

- name:

  Name of the demo file (with or without .Rmd extension)

## Value

Full path to the demo Rmd file

## Examples

``` r
# Get path to a demo
demo_path <- get_demo_path("example-demo")
#> Error in get_demo_path("example-demo"): Demo 'example-demo.Rmd' not found. Available demos:
#> - hmm_bootstrap_demo.Rmd
#> - ms-varma-garch_diagnostics_report.Rmd
#> - ms-varma-garch_inference_demo.Rmd
#> - multivar_garch_comparison.Rmd
#> - portfolio-optimization.Rmd
#> - portfolio-optimization_2.Rmd

# If you want to render it
# rmarkdown::render(demo_path, output_dir = tempdir())
```
