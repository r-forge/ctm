#' S3 methods for class \code{"tramnet"}
#'
#' @rdname methods
#'
#' @param object Object of class \code{"tramnet"}.
#' @param parm Parameters to evaluate the log likelihood at.
#' @param w Optional vector of sample weights.
#' @param newdata Data to evaluate the log likelihood at.
#' @param add_penalty Whethr or not to return the penalized log-likelihood
#'     (default \code{add_penalty = FALSE}).
#' @param ... Additional arguments to \code{\link[mlt]{logLik.mlt}}
#'
#' @return Returns (potentially weighted \code{w}) log-likelihood based on
#' \code{object} evaluated at parameters \code{parm} and data \code{newdata}
#'
#' @exportS3Method logLik tramnet
#'
logLik.tramnet <- function(object,
                           parm = coef(object, tol = 0, with_baseline = TRUE),
                           w = NULL, newdata = NULL, add_penalty = FALSE, ...) {
  if (length(list(...)) > 0L)
    warning("additional arguments ignored")

  ctmobj <- .tramnet2ctm(object)
  newdata <- .get_tramnet_data(object, newdata)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  if (object$model$negative)
    parm <- .flip_sign(object, parm)
  ret <- logLik(mltobj, parm = parm, w = w, ...)
  pen <- if (add_penalty) .elnet(object) else 0
  attr(ret, "df") <- NA
  class(ret) <- "logLik"
  return(ret - pen)
}

#' @rdname methods
#'
#' @param with_baseline If \code{TRUE}, also prints coefficients for the
#'     baseline transformation.
#' @param tol Tolerance when an estimate should be considered 0 and not
#'     returned (default \code{tol = 1e-6}).
#' @param ... Additional arguments to \code{\link[mlt]{coef.mlt}}
#'
#' @return Numeric vector containing the model shift parameter estimates
#'
#' @exportS3Method coef tramnet
coef.tramnet <- function(object, with_baseline = FALSE, tol = 1e-6, ...) {
  if (length(list(...)) > 0L)
    warning("additional arguments ignored")
  beta <- c(object$beta[abs(object$beta) > tol])
  if (all(object$x == 0))
    beta <- NULL
  theta <- c(object$theta)
  names(theta) <- names(coef(as.mlt(object$model)))
  if (!with_baseline)
    return(beta)
  return(c(theta, beta))
}

#' @rdname methods
#'
#' @param as.lm See \code{\link[mlt]{coef.mlt}}
#'
#' @return Numeric vector containing the linear model shift parameter estimates
#'
#' @exportS3Method coef tramnet_Lm
coef.tramnet_Lm <- function(object, with_baseline = FALSE, tol = 1e-6,
                            as.lm = FALSE, ...) {
  class(object) <- class(object)[-1L]
  if (!as.lm)
    return(coef(object, with_baseline = with_baseline, tol = tol, ...))
  if (!is.null(object$stratacoef))
    stop("Cannot compute scaled coefficients with strata variables present")
  cf <- coef(object, with_baseline = TRUE, tol = 0, ...)
  cfx <- coef(object, with_baseline = FALSE, tol = 0, ...)
  cfs <- cf[!(names(cf) %in% names(cfx))]
  sd <- 1/cfs[names(cfs) != "(Intercept)"]
  ret <- c(-cf["(Intercept)"], cfx) * sd
  attr(ret, "scale") <- sd
  return(ret)
}


#' @rdname methods
#'
#' @return
#' Vector of predictions based on \code{object} evaluated at each row
#' of \code{newdata}
#'
#' @exportS3Method predict tramnet
#'
predict.tramnet <- function(object, newdata = NULL, ...) {
  newdata <- .get_tramnet_data(object, newdata)
  ctmobj <- .tramnet2ctm(object)
  predict(ctmobj, newdata = newdata, ...)
}

#' @rdname methods
#'
#' @param nsim Number of simulations, see \code{\link[mlt]{simulate.mlt}}.
#' @param seed Random seed, see \code{\link[mlt]{simulate.mlt}}.
#' @param bysim Return by simulation, see \code{\link[mlt]{simulate.mlt}}.
#'
#' @return Returns a \code{list} of \code{data.frames} containing parametric
#'     bootstrap samples of the response based on the data supplied in
#'     \code{newdata}
#'
#' @exportS3Method simulate tramnet
#'
simulate.tramnet <- function(object, nsim = 1, seed = NULL,
                             newdata = NULL,
                             bysim = TRUE, ...) {
  newdata <- .get_tramnet_data(object, newdata)
  ctmobj <- .tramnet2ctm(object)
  simulate(ctmobj, nsim = nsim, seed = seed,
           newdata = newdata, bysim = bysim, ...)
}

#' @rdname methods
#'
#' @param x Object of class \code{"tramnet"}.
#'
#' @return Matrix of score contributions w.r.t. model parameters evaluated at
#'     \code{parm}
#'
#' @exportS3Method estfun tramnet
#'
estfun.tramnet <- function(x, parm = coef(x, with_baseline = TRUE, tol = 0),
                           w = NULL, newdata = NULL, ...) {
  if (x$tuning_parm["lambda"] > 0)
    stop("Cannot compute the score for penalised parameters.")
  newdata <- .get_tramnet_data(x, newdata)
  ctmobj <- .tramnet2ctm(x)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  return(estfun(mltobj, parm = parm, w = w))
}

#' @rdname methods
#'
#' @return Returns a numeric vector of residuals for each row in \code{newdata}
#' @exportS3Method residuals tramnet
#'
residuals.tramnet <- function(object,
                              parm = coef(object, tol = 0,
                                          with_baseline = TRUE),
                              w = NULL, newdata = NULL, ...) {
  newdata <- .get_tramnet_data(object, newdata)
  ctmobj <- .tramnet2ctm(object)
  mltobj <- mlt(ctmobj, data = newdata, dofit = FALSE)
  return(residuals(object = mltobj, parm = parm,
                   w = w, newdata = newdata, ...))
}

#' Plot \code{"tramnet"} objects
#'
#' @param x Object of class \code{"tramnet"}.
#' @param newdata See \code{\link[mlt]{plot.ctm}}.
#' @param type See \code{\link[mlt]{plot.ctm}}.
#' @param q See \code{\link[mlt]{plot.ctm}}.
#' @param prob See \code{\link[mlt]{plot.ctm}}.
#' @param K See \code{\link[mlt]{plot.ctm}}.
#' @param col See \code{\link[mlt]{plot.ctm}}.
#' @param lty See \code{\link[mlt]{plot.ctm}}.
#' @param add See \code{\link[mlt]{plot.ctm}}.
#' @param ... Additional arguments passed to \code{\link[mlt]{plot.ctm}}.
#'
#' @exportS3Method plot tramnet
#'
plot.tramnet <- function(x, newdata = NULL,
                         type = c("distribution", "survivor", "density",
                                  "logdensity", "hazard", "loghazard",
                                  "cumhazard", "quantile", "trafo"),
                         q = NULL, prob = 1:(K - 1)/K, K = 50,
                         col = rgb(.1, .1, .1, .1), lty = 1, add = FALSE, ...) {
  newdata <- .get_tramnet_data(x, newdata)
  ctmobj <- .tramnet2ctm(x)
  plot(ctmobj, newdata = newdata, type = type, q = q, prob = prob, K = K,
       col = col, lty = lty, add = add, ...)
}

# print method for class "tramnet"

#' @rdname methods
#'
#' @return Object of class \code{"summary.tramnet"}.
#' @exportS3Method print tramnet
#'
print.tramnet <- function(x, ...) {
  print(summary(x, ...))
}

#' @rdname methods
#'
#' @return Object of class \code{"summary.tramnet"}.
#' @exportS3Method summary tramnet
#'
summary.tramnet <- function(object, ...) {
  tp <- object$model$tram
  if (object$tuning_parm[1] > 0)
    tp <- paste0("Regularized", tp)
  ret <- list(call = getCall(object), convergence = object$result$status,
              type = tp, logLik = logLik(object), coef = coef(object),
              sparsity = .tramnet_sparsity(object),
              tuning_parm = object$tuning_parm)
  class(ret) <- "summary.tramnet"
  ret
}

#' @rdname methods
#'
#' @param digits Number of digits to print.
#' @param ... Ignored.
#'
#' @return Invisible \code{x}.
#' @exportS3Method print summary.tramnet
#'
print.summary.tramnet <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nConvergence: ", x$convergence)
  cat("\nType: ", x$type)
  cat("\nLog-Likelihood: ", x$logLik)
  cat("\n")
  cat("\nCoefficients:\n")
  if (!is.null(x$coef)) {
    print(round(x$coef, digits = digits))
  } else {
    print(NULL)
  }
  cat("\nSparsity: ", x$sparsity, "\n")
  cat("\nTuning parameters:\n")
  print(round(x$tuning_parm, digits = digits))
  cat("\n\n")
  invisible(x)
}

# Helper functions

.tramnet2ctm <- function(object) {
  data <- .get_tramnet_data(object)
  cfx <- coef(object, with_baseline = TRUE, tol = 0)
  if (!is.null(object$model$model$model$binteracting)) {
    yBasis <- object$model$model$model$binteracting$iresponse
    iBasis <- object$model$model$model$binteracting$iinteracting
  } else {
    yBasis <- object$model$model$model$bresponse
    iBasis <- NULL
  }

  todistr <- switch(object$model$todistr$name,
                    "minimum extreme value" = "MinExtrVal",
                    "maximum extreme value" = "MaxExtrVal",
                    "normal" = "Normal", "logistic" = "Logistic")
  if (all(object$x == 0)) {
    shifting <- NULL
  } else {
    shifting <-
      as.basis(
        as.formula(
          paste("~", paste(colnames(object$x), collapse = "+"))
        ), data = data, remove_intercept = TRUE
      )
  }
  mod <- ctm(response = yBasis, shifting = shifting,
             interacting = iBasis, todistr = todistr, data = data)
  coef(mod) <- cfx
  return(mod)
}

#' @importFrom stats model.response model.frame
.get_tramnet_data <- function(object, newdata = NULL) {
  if (is.null(newdata))
    return(cbind(object$model$data, object$x))
  if (is.null(object$process_newdata))
    return(newdata)
  newr <- as.data.frame(model.response(model.frame(
    object$call$formula, data = newdata)))
  colnames(newr) <- variable.names(object$model)
  newx <- object$process_newdata(newdata)
  cbind(newr, newx)
}

.tramnet_sparsity <- function(object, ...) {
  n_coef <- length(coef(object, tol = 0))
  non_zero_coef <- length(coef(object, ...))
  ret <- paste(n_coef, "regression coefficients,", non_zero_coef,
               "of which are non-zero")
  return(ret)
}

.elnet <- function(object) {
  lambda <- object$tuning_parm[1]
  alpha <- object$tuning_parm[2]
  cfx <- coef(object, tol = 0)
  if (is.null(cfx))
    return(0)
  L1 <- sum(abs(cfx))
  L2 <- sqrt(sum(cfx^2))
  ret <- lambda * (0.5 * (1 - alpha) * L2^2 + alpha * L1)
  names(ret) <- NULL
  attr(ret, "L1") <- L1
  attr(ret, "L2") <- L2
  return(ret)
}

.flip_sign <- function(object, parm) {
  stopifnot(inherits(object, "tramnet"))
  wbl <- length(coef(object, tol = 0, with_baseline = TRUE))
  wobl <- length(coef(object, tol = 0))
  idx <- wbl - wobl + seq_len(wobl)
  parm[idx] <- - parm[idx]
  return(parm)
}
