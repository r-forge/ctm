
#' Regularized transformation model classes
#' @rdname model-classes
#' @inheritParams tramnet.formula
#' @param formula Formula specifying the regression. See \code{\link[tram]{tram}}.
#' @return Object of class \code{"tramnet"}.
#' @export
LmNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda, alpha = alpha,
                         tram_fun = "Lm", tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
BoxCoxNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "BoxCox",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
ColrNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "Colr",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
SurvregNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "Survreg",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
CoxphNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "Coxph",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
LehmannNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "Lehmann",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}

#' @rdname model-classes
#' @export
PolrNET <- function(
    formula, data, lambda = 0, alpha = 1, tram_args = NULL, constraints = NULL, ...
) {

  call <- match.call()
  ret <- tramnet.formula(model = formula, data = data, lambda = lambda,
                         alpha = alpha, tram_fun = "Polr",
                         tram_args = tram_args, constraints = constraints,
                         ... = ...)
  ret$call <- call
  ret
}
