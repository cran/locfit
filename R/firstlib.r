## .First.lib <- function(lib, pkg) {
##  cat("Locfit for R-2.0.1.\n")
##  cat("December 16, 2004.\n")
##  library.dynam("locfit", pkg, lib)
##  invisible(NULL)
##}

.First.lib <- function(libname, pkgname) {
    ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                      fields=c("Version", "Date"))
    cat(paste(pkgname, ver[1], "\t", ver[2], "\n"))
    library.dynam("locfit", pkgname, libname)
}
