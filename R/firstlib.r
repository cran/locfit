.First.lib <- function(lib, pkg) {
  cat("Locfit for R.\n")
  cat("August 3, 2000.\n")
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}
