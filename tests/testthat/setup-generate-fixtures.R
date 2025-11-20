# Script to generate test fixtures
# Run this manually when setting up tests or when fixtures need updating
# 
# Usage:
#   source("tests/testthat/setup-generate-fixtures.R")
#   generate_all_fixtures()

#' Generate All Test Fixtures
#'
#' @description This function generates and caches all test fixtures.
#' It should be run manually when setting up the test suite or when
#' fixtures need to be regenerated.
#'
#' @param force If TRUE, regenerates all fixtures even if they exist
#' @return NULL (invisibly)
#' @examples
#' \dontrun{
#' generate_all_fixtures(force = FALSE)
#' }
generate_all_fixtures <- function(force = FALSE) {
  
  message("Generating test fixtures...")
  message("This may take several minutes as models are being estimated.")
  
  # Check required packages
  required_packages <- c("tsmarch", "tsgarch", "xts", "testthat")
  missing <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing) > 0) {
    stop("Missing required packages: ", paste(missing, collapse = ", "),
         "\nPlease install them first.")
  }
  
  # Create fixtures directory if it doesn't exist
  fixtures_dir <- file.path("tests", "testthat", "fixtures")
  if (!dir.exists(fixtures_dir)) {
    dir.create(fixtures_dir, recursive = TRUE)
    message("Created fixtures directory: ", fixtures_dir)
  }
  
  # List of fixtures to generate
  fixtures <- list(
    list(name = "DCC Dynamic (MVN)", 
         file = "dcc_dynamic_fit.rds",
         fun = create_test_dcc_dynamic),
    list(name = "DCC Constant (MVN)", 
         file = "dcc_constant_fit.rds",
         fun = create_test_dcc_constant),
    list(name = "DCC Dynamic (Student-t)", 
         file = "dcc_student_fit.rds",
         fun = create_test_dcc_student),
    list(name = "Copula-GARCH Dynamic", 
         file = "copula_dynamic_fit.rds",
         fun = create_test_copula_dynamic)
  )
  
  # Generate each fixture
  for (fixture in fixtures) {
    fixture_path <- file.path(fixtures_dir, fixture$file)
    
    if (!force && file.exists(fixture_path)) {
      message("  [SKIP] ", fixture$name, " (already exists)")
    } else {
      message("  [GENERATING] ", fixture$name, "...")
      tryCatch({
        fit <- fixture$fun(use_cache = FALSE)
        saveRDS(fit, fixture_path)
        
        # Verify the saved object
        test_load <- readRDS(fixture_path)
        ll <- as.numeric(logLik(test_load))
        
        message("    [OK] ", fixture$name, " (loglik = ", round(ll, 2), ")")
      }, error = function(e) {
        message("    [FAILED] ", fixture$name)
        message("    Error: ", e$message)
      })
    }
  }
  
  message("\nFixture generation complete!")
  message("Fixtures saved in: ", fixtures_dir)
  
  invisible(NULL)
}

#' Check Test Fixtures
#'
#' @description Checks whether all test fixtures exist and are loadable
#' @return A data.frame with fixture status
check_fixtures <- function() {
  
  fixtures_dir <- file.path("tests", "testthat", "fixtures")
  
  fixture_files <- c(
    "dcc_dynamic_fit.rds",
    "dcc_constant_fit.rds",
    "dcc_student_fit.rds",
    "copula_dynamic_fit.rds"
  )
  
  results <- data.frame(
    fixture = fixture_files,
    exists = FALSE,
    loadable = FALSE,
    class = NA_character_,
    loglik = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(fixture_files)) {
    fixture_path <- file.path(fixtures_dir, fixture_files[i])
    
    if (file.exists(fixture_path)) {
      results$exists[i] <- TRUE
      
      tryCatch({
        obj <- readRDS(fixture_path)
        results$loadable[i] <- TRUE
        results$class[i] <- paste(class(obj), collapse = ", ")
        results$loglik[i] <- as.numeric(logLik(obj))
      }, error = function(e) {
        results$loadable[i] <- FALSE
      })
    }
  }
  
  return(results)
}

# If running this script directly, generate fixtures
if (interactive() || identical(Sys.getenv("GENERATE_FIXTURES"), "TRUE")) {
  message("\n========================================")
  message("Test Fixture Generation")
  message("========================================\n")
  
  # Check current status
  message("Current fixture status:")
  status <- check_fixtures()
  print(status)
  message()
  
  # Ask user if they want to proceed
  if (interactive()) {
    response <- readline("Generate/update fixtures? (yes/no): ")
    if (tolower(response) %in% c("yes", "y")) {
      generate_all_fixtures(force = FALSE)
    } else {
      message("Fixture generation cancelled.")
    }
  } else {
    generate_all_fixtures(force = FALSE)
  }
  
  # Show updated status
  message("\nUpdated fixture status:")
  status <- check_fixtures()
  print(status)
}