.First.lib <- function(lib, pkg) {
  cat("Locfit for R 1.3.0.\n")
  cat("August 3, 2000.\n")
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}
