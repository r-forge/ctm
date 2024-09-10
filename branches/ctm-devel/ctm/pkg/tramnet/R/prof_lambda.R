# Profiling tuning parameters

#' Profiling tuning parameters
#'
#' @param model Object of class \code{"tramnet"}.
#' @param min_lambda Minimal value of lambda (default \code{min_lambda = 0}).
#' @param max_lambda Maximal value of lambda (default \code{max_lambda = 15}).
#' @param nprof Number of profiling steps (default \code{nprof = 5}).
#' @param as.lm Return scaled coefficients for class \code{"tramnet_Lm"}.
#'
#' @description Computes the regularization path of all coefficients for a
#'    single tuning parameter, lambda, over a sequence of values.
#'
#'
#' @return Object of class \code{"prof_lambda"} which contains the
#'     regularization path of all coefficients and the log-likelihood over the
#'     penalty parameter lambda
#'
#' @examples
#' \donttest{
#' if (require("survival") & require("penalized")) {
#'   data("nki70", package = "penalized")
#'   nki70$resp <- with(nki70, Surv(time, event))
#'   x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
#'   y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
#'   fit <- tramnet(y, x, lambda = 0, alpha = 1)
#'   pfl <- prof_lambda(fit)
#'   plot_path(pfl)
#' }
#' }
#'
#' @export
prof_lambda <- function(model, min_lambda = 0, max_lambda = 15, nprof = 5,
                        as.lm = FALSE) {
  stopifnot(inherits(model, "tramnet"))
  stopifnot(max_lambda > min_lambda)
  stopifnot(nprof > 0)
  if (min_lambda == 0)
    lambdas <- c(0, 10^seq(-1, log10(max_lambda), length.out = nprof)[-1])
  else
    lambdas <- 10^seq(log10(min_lambda), log10(max_lambda), length.out = nprof)

  cfx <- list()
  lls <- list()

  for (idx in seq_along(lambdas)) {
    message("Step ", idx, "/", length(lambdas),
            " at lambda = ", round(lambdas[idx], 2))
    mod <- try(update(model, lambda = lambdas[idx]))
    if (inherits(mod, "try-error")) {
      cfx[[idx]] <- rep(NA, length(coef(model, tol = 0, as.lm = as.lm)))
      lls[[idx]] <- NA
    } else {
      cfs <- coef(mod, with_baseline = FALSE, tol = 0, as.lm = as.lm)
      cfx[[idx]] <- cfs[names(cfs) != "(Intercept)"]
      lls[[idx]] <- logLik(mod)
    }
  }

  ret <- list(
    lambdas = lambdas,
    cfx = do.call(rbind, cfx),
    lls = do.call(c, lls)
  )

  class(ret) = "prof_lambda"

  return(ret)
}

#' Profiling tuning parameters
#'
#' @description Computes the regularization path of all coefficients for a
#'     single tuning, alpha, parameter over a sequence of values.
#'
#' @param model Object of class \code{"tramnet"}.
#' @param min_alpha Minimal value of alpha (default \code{min_alpha = 0}).
#' @param max_alpha Maximal value of alpha (default \code{max_alpha = 1}).
#' @param nprof Number of profiling steps (default \code{nprof = 5}).
#' @param as.lm Return scaled coefficients for class \code{"tramnet_Lm"}.
#'
#' @return Object of class \code{"prof_alpha"} which contains the regularization
#'     path of all coefficients and the log-likelihood over the mixing parameter
#'     alpha
#'
#' @examples
#' \donttest{
#' if (require("survival") & require("penalized")) {
#'   data("nki70", package = "penalized")
#'   nki70$resp <- with(nki70, Surv(time, event))
#'   x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
#'   y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
#'   fit <- tramnet(y, x, lambda = 1, alpha = 1)
#'   pfa <- prof_alpha(fit)
#'   plot_path(pfa)
#' }
#' }
#'
#' @export
prof_alpha <- function(model, min_alpha = 0, max_alpha = 1, nprof = 5,
                       as.lm = FALSE) {
  stopifnot(inherits(model, "tramnet"))
  stopifnot(.check_bounds(min_alpha, max_alpha))
  stopifnot(nprof > 0)
  alphas <- seq(min_alpha, max_alpha, length.out = nprof)
  cfx <- list()
  lls <- list()

  for (idx in seq_along(alphas)) {
    message("Step ", idx, "/", length(alphas),
            " at alpha = ", round(alphas[idx], 2))
    mod <- try(update(model, alpha = alphas[idx]))
    if (inherits(mod, "try-error")) {
      cfx[[idx]] <- rep(NA, length(coef(model, tol = 0, as.lm = as.lm)))
      lls[[idx]] <- NA
    } else {
      cfs <- coef(mod, with_baseline = FALSE, tol = 0, as.lm = as.lm)
      cfx[[idx]] <- cfs[names(cfs) != "(Intercept)"]
      lls[[idx]] <- logLik(mod)
    }
  }

  ret <- list(
    alphas = alphas,
    cfx = do.call(rbind, cfx),
    lls = do.call(c, lls)
  )

  class(ret) = "prof_alpha"

  return(ret)
}

#' Plot regularization paths
#'
#' @description Plot regularization paths and optionally log-likelihood
#'     trajectories of objects of class \code{"prof_alpha"} and
#'     \code{"prof_lambda"}. Coefficient names are automatically added to the
#'     plot.
#'
#' @param object Object of class \code{"prof_alpha"} or \code{"prof_lambda"}.
#' @param plot_logLik Whether \code{logLik} trajectory should be plotted
#'     (default \code{plot_logLik = FALSE}).
#' @param ... Additional arguments to \code{\link{plot}}
#'
#' @return None.
#'
#' @examples
#' \donttest{
#' if (require("survival") & require("penalized")) {
#'   data("nki70", package = "penalized")
#'   nki70$resp <- with(nki70, Surv(time, event))
#'   x <- scale(model.matrix( ~ 0 + DIAPH3 + NUSAP1 + TSPYL5 + C20orf46, data = nki70))
#'   y <- Coxph(resp ~ 1, data = nki70, order = 10, log_first = TRUE)
#'   fit1 <- tramnet(y, x, lambda = 0, alpha = 1)
#'   pfl <- prof_lambda(fit1)
#'   plot_path(pfl)
#'   fit2 <- tramnet(y, x, lambda = 1, alpha = 1)
#'   pfa <- prof_alpha(fit2)
#'   plot_path(pfa)
#' }
#' }
#'
#' @export
plot_path <- function(object, plot_logLik = FALSE, ...) {
  if (inherits(object, "prof_lambda"))
    .plot_lpath(object, plot_logLik, ...)
  else
    if (inherits(object, "prof_alpha"))
      .plot_apath(object, plot_logLik, ...)
  else
    stop("plot_path() needs an object of class prof_lambda or prof_alpha")
}

.plot_lpath <- function(object, plot_logLik, ...) {
  if (plot_logLik)
    plot(object$lambdas, object$lls, xlab = expression(lambda),
         ylab = "logLik", ...)
  opar <- par("mar")
  on.exit(par(opar))
  par(mar = c(5, 5, 2, 2))
  matplot(object$lambdas, object$cfx, type = "l", xlab = expression(lambda),
          ylab = expression(hat(beta)[j](lambda)), ...)
  text(x = min(object$lambdas) + 0.05 * abs(diff(range(object$lambdas))),
       y = object$cfx[which.min(object$lambdas), ],
       labels = colnames(object$cfx), ...)
}

.plot_apath <- function(object, plot_logLik, ...) {
  if (plot_logLik)
    plot(object$alphas, object$lls, xlab = expression(alpha),
         ylab = "logLik", ...)
  opar <- par("mar")
  par(mar = c(5, 5, 2, 2))
  matplot(object$alphas, object$cfx, type = "l", xlab = expression(alpha),
          ylab = expression(hat(beta)[j](alpha)), ...)
  text(x = min(object$alphas) + 0.05 * abs(diff(range(object$alphas))),
       y = object$cfx[which.min(object$alphas), ],
       labels = colnames(object$cfx), ...)
  par(mar = opar)
}

.check_bounds <- function(mi, ma) {
  if (ma > mi && ma >= 0 && ma <= 1 && mi <= 1 && mi >= 0)
    ret <- TRUE
  else
    ret <- FALSE
  return(ret)
}
