.First.lib <- function(lib, pkg) {
  cat("Locfit for R 0.99.0.\n")
  cat("Feb 14, 2000.\n")
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}

print.gcvplot <- function(object, ...) plot.gcvplot(x=object,...)

provide(locfit)
