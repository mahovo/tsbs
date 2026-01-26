#' List available demonstration Rmd files
#' 
#' These demonstrations are computationally intensive and not built as vignettes.
#' Use `get_demo_path()` to get the full path to a specific demo.
#' 
#' @return Character vector of available demo Rmd file names
#' @export
#' @examples
#' list_package_demos()
list_package_demos <- function() {
  demo_dir <- system.file("demos", package = "tsbs")
  if (demo_dir == "") {
    return(character(0))  # Return empty if no demos directory
  }
  files <- list.files(demo_dir, pattern = "\\.Rmd$", full.names = FALSE)
  return(files)
}

#' Get path to a specific demo Rmd file
#' 
#' @param name Name of the demo file (with or without .Rmd extension)
#' @return Full path to the demo Rmd file
#' @export
#' @examples
#' # Get path to a demo
#' demo_path <- get_demo_path("example-demo")
#' 
#' # If you want to render it
#' # rmarkdown::render(demo_path, output_dir = tempdir())
get_demo_path <- function(name) {
  # Add .Rmd extension if missing
  if (!grepl("\\.Rmd$", name)) {
    name <- paste0(name, ".Rmd")
  }
  
  path <- system.file("demos", name, package = "tsbs")
  
  if (path == "") {
    available <- list_package_demos()
    if (length(available) == 0) {
      stop("No demos found in the 'tsbs' package.")
    }
    stop("Demo '", name, "' not found. Available demos:\n",
         paste("-", available, collapse = "\n"))
  }
  
  return(path)
}

#' Open a demo Rmd file in RStudio or default editor
#' 
#' @param name Name of the demo file (with or without .Rmd extension)
#' @export
open_demo <- function(name) {
  if (!requireNamespace("rstudioapi", quietly = TRUE)) {
    stop("Please install the 'rstudioapi' package to use this function.")
  }
  
  path <- get_demo_path(name)
  
  if (rstudioapi::isAvailable()) {
    rstudioapi::navigateToFile(path)
  } else {
    message("Opening demo file: ", path)
    utils::file.show(path)
  }
  
  invisible(path)
}