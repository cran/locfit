.First.lib <- function(lib, pkg) {
  cat("Locfit for R 1.2.0.\n")
  cat("2000-12-19.\n")
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}
