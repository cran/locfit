.First.lib <- function(lib, pkg) {
  cat("Locfit for R.\n")
  cat("August 3, 2000.  (Updated for R 1.7.0, March 21, 2003)\n")
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}
