.First.lib <- function(lib, pkg) {
  library.dynam("locfit", pkg, lib)
  invisible(NULL)
}

provide(locfit)
