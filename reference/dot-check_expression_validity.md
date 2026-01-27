# Check if expression or variable exists and is evaluable

DEPRECATED. Only called by .is_invalid_data(), which is DEPRECATED.

## Usage

``` r
.check_expression_validity(expr_sub, parent_env = parent.frame())
```

## Arguments

- parent_env:

  The parent environment to check for variable existence

- x:

  The expression/variable to check (as substitute() result)

## Value

TRUE if the expression/variable is invalid/non-existent, FALSE otherwise

## Details

Helper function that distinguishes between: 1. Expressions that evaluate
successfully (e.g., numeric(0), 1+1) 2. Undefined variables (e.g.,
nonexistent_var) 3. Expressions that fail to evaluate for other reasons
