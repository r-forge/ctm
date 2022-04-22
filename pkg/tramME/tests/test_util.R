## Test utils: to be sourced from other test files

## Decide based on the package version number whether this is a CRAN version
if (length(strsplit(packageDescription("tramME")$Version, "\\.")[[1]]) > 3) {
  Sys.setenv("NOT_CRAN" = "true")
}

## Increase the default tolearnce levels for chkeq when run on CRAN
## to avoid "additional issues"
if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
  tol <- 1e4
} else {
  tol <- sqrt(.Machine$double.eps)
}

chkeq <- function(x, y, tolerance = tol, ...) {
  stopifnot(isTRUE(all.equal(x, y, tolerance = tolerance, ...)))
}

chkid <- function(x, y, ...) stopifnot(isTRUE(identical(x, y, ...)))

chkerr <- function(expr, em = NULL) {
  stopifnot(
    tryCatch(expr,
      error = function(e) {
        if (!is.null(em)) return(grepl(em, e))
        else return(TRUE)
      }
      )
  )
}

chkwarn <- function(expr, wm = NULL) {
  stopifnot(
    tryCatch(expr,
      warning = function(w) {
        if (!is.null(wm)) return(grepl(wm, w))
        else return(TRUE)
      }
      )
  )
}
