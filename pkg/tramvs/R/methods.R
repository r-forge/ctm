# S3 methods

#' Plot "tramvs" object
#' @param x object of class \code{"tramvs"}
#' @param which plotting either the regularization path (\code{"path"})
#'     or the information criterion against the support size (\code{"tune"}, default)
#' @param ... additional arguments to \code{plot()}
#' @method plot tramvs
#' @return Returns \code{invisible(NULL)}
#' @importFrom graphics par matplot text
#' @exportS3Method
plot.tramvs <- function(x, which = c("tune", "path"), ...) {
  which <- match.arg(which)
  opar <- par(no.readonly = TRUE)
  on.exit(opar)
  par(mar = c(5.1, 5.1, 4.1, 2.1), las = 1)
  switch(
    which,
    "tune" = plot(x$SIC, ... = ...),
    "path" = {
      coefs <- as.matrix(coef(x))
      supps <- as.numeric(colnames(coefs))
      matplot(supps, t(coefs), type = "l", ... = ..., xlab = "support",
              ylab = expression(hat(beta)[j]))
      text(max(supps) * 1.01, x$coefs[, ncol(coefs)],
           rownames(coefs), cex = 0.8)
    }
  )
  return(invisible(NULL))
}

#' Coef "tramvs"
#' @param object Object of class \code{"tramvs"}
#' @param best_only Wether to return the coefficients of the best model only
#'     (default: FALSE)
#' @param ... additional arguments to \code{coef()}
#' @return Vector (\code{best_only = TRUE}) or matrix (\code{best_only = FALSE})
#'     of coefficients
#' @method coef tramvs
#' @exportS3Method
coef.tramvs <- function(object, best_only = FALSE, ...) {
  if (best_only)
    return(coef(object$best_fit$mod, ... = ...))
  coefs <- as(do.call("cbind", lapply(object$all_fits, coef, ... = ...)),
              "sparseMatrix")
  colnames(coefs) <- seq_len(ncol(coefs))
  coefs
}

#' Predict "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... additional arguments to \code{predict.tram()}
#' @method predict tramvs
#' @return See \code{\link[tram]{predict.tram}}
#' @importFrom stats predict
#' @exportS3Method
predict.tramvs <- function(object, ...) {
  predict(object$best_fit$mod, ... = ...)
}

#' Simulate "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param nsim number of simulations
#' @param seed random seed for simulation
#' @param ... additional arguments to \code{simulate()}
#' @method simulate tramvs
#' @return See \code{\link[mlt]{simulate.mlt}}
#' @importFrom stats simulate
#' @exportS3Method
simulate.tramvs <- function(object, nsim = 1, seed = NULL, ...) {
  simulate(object$best_fit$mod, nsim = nsim, seed = seed, ... = ...)
}

#' Print "tramvs"
#' @param x object of class \code{"tramvs"}
#' @param ... ignored
#' @return \code{"tramvs"} object is returned invisibly
#' @method print tramvs
#' @exportS3Method
print.tramvs <- function(x, ...) {
  cat("\nL0-penalized tram:\n")
  print(x$best_fit$mod)
  cat("\nSIC:\n", min(x$SIC$SIC), "\n")
  cat("\nActive set:", x$best_fit$A, "\n\n")
  invisible(x)
}

#' Summary "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... ignored
#' @method summary tramvs
#' @return \code{"tramvs"} object is returned invisibly
#' @exportS3Method
summary.tramvs <- function(object, ...) {
  cat("\nL0-penalized tram:\n")
  print(object$best_fit$mod)
  cat("\nSIC:\n")
  print(object$SIC)
  cat("\n")
  cat("\nActive set:", object$best_fit$A, "\n\n")
  invisible(object)
}

#' logLik "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... additional arguments to \code{logLik()}
#' @method logLik tramvs
#' @return Numeric vector containing log-likelihood of best model,
#'      see \code{\link[tram]{logLik.tram}}
#' @exportS3Method
logLik.tramvs <- function(object, ...) {
  logLik(object$best_fit$mod, ... = ...)
}

#' AIC "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... additional arguments to \code{AIC()}
#' @method AIC tramvs
#' @return Numeric vector containing AIC of best model
#' @importFrom stats AIC
#' @exportS3Method
AIC.tramvs <- function(object, ...) {
  AIC(object$best_fit$mod)
}

#' SIC generic
#' @param object Model to compute SIC from
#' @param ... for methods compatibility only
#' @return Numeric vector (\code{best_only = TRUE}) or data.frame with SIC values
#' @export
SIC <- function(object, ...) {
  UseMethod("SIC", object)
}

#' SIC "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param best_only Wether to return the coefficients of the best model only
#'     (default: FALSE)
#' @param ... for methods compatibility only
#' @return Numeric vector (\code{best_only = TRUE}) or data.frame with SIC values
#' @method SIC tramvs
#' @exportS3Method
SIC.tramvs <- function(object, best_only = FALSE, ...) {
  if (best_only)
    return(min(object$SIC$SIC))
  object$SIC
}

#' Residuals "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... additional arguments to \code{residuals()}
#' @method residuals tramvs
#' @return Numeric vector containing residuals of best model,
#'      see \code{\link[tram]{residuals.tram}}
#' @exportS3Method
residuals.tramvs <- function(object, ...) {
  residuals(object$best_fit$mod, ... = ...)
}

#' Support "tramvs"
#' @param object object of class \code{"tramvs"}
#' @param ... ignored
#' @importFrom variables support
#' @method support tramvs
#' @return Character vector containing active set of best fit
#' @exportS3Method
support.tramvs <- function(object, ...) {
  object$best_fit$A
}

#' Coef "abess_tram"
#' @param object object of class \code{"tramvs"}
#' @param ... additional arguments to \code{coef()}
#' @method coef abess_tram
#' @return Named numeric vector containing coefficient estimates
#'      see \code{\link[tram]{coef.tram}}
#' @exportS3Method
coef.abess_tram <- function(object, ...) {
  coef(object$mod, ... = ...)
}
