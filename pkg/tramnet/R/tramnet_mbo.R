#' Model based optimization for regularized transformation models
#'
#' @description Uses model based optimization to find the optimal tuning
#'     parameter(s) in a regularized transformation model based on
#'     cross-validated log-likelihoods. Here the \code{'tramnet'} package makes
#'     use of the \code{'mlr3mbo'} interface for Bayesian Optimization in
#'     machine learning problems to maximize the cv-logLik as a black-box
#'     function of the tuning parameters alpha and lambda.
#'
#' @param object Object of class \code{"tramnet"}.
#' @param fold Number of cross validation folds.
#' @param n_iter Maximum number of iterations in model-based optimization
#'     routine.
#' @param minlambda Minimum value for lambda (default \code{minlambda = 0}).
#' @param maxlambda Maximum value for lambda (default \code{maxlambda = 16}).
#' @param minalpha Minimum value for alpha (default \code{minalpha = 0}).
#' @param maxalpha Maximum value for alpha (default \code{maxalpha = 1}).
#' @param folds Self specified folds for cross validation (mainly for
#'     reproducibility and comparability purposes).
#' @param noisy indicates whether folds for k-fold cross-validation should
#'     be random for each iteration, leading to a noisy objective function
#'     (default \code{noisy = FALSE}).
#' @param obj_type Objective type, one of \code{"lasso"}, \code{"ridge"} or
#'     \code{"elnet"}.
#' @param verbose Toggle for a verbose output (default \code{verbose = TRUE})
#' @param ... Currently ignored.
#'
#' @return See \code{\link[bbotk]{Optimizer}}'s \code{optimize} function which
#'     returns a \code{data.table::data.table}.
#'
#' @importFrom stats residuals
#' @importFrom bbotk opt trm ObjectiveRFun OptimInstanceSingleCrit
#' @importFrom utils capture.output
#'
#' @export
mbo_tramnet <- function(object, fold = 2, n_iter = 5, minlambda = 0,
                        maxlambda = 16, minalpha = 0, maxalpha = 1,
                        folds = NULL, noisy = FALSE,
                        obj_type = c("lasso", "ridge", "elnet"),
                        verbose = TRUE, ...) {
  if (length(list(...)) > 0) {
    warning("Additional arguments ignored.")
  }
  stopifnot(inherits(object, "tramnet"))
  df <- .get_tramnet_data(object)
  rsp <- variable.names(object$model, "response")
  n <- nrow(df)
  if (!noisy) {
    if (is.null(folds) & !is.null(fold)) {
      folds <- sample(rep(1:fold, ceiling(n/fold)), n)
    } else {
      folds <- round(folds)
      fold <- max(folds)
    }
  } else {
    folds <- NULL
  }
  obj_type <- match.arg(obj_type)
  if (obj_type == "lasso") minalpha <- maxalpha <- 0
  if (obj_type == "ridge") minalpha <- maxalpha <- 1
  inst <- .mbo_obj(object = object, minlambda = minlambda,
                   maxlambda = maxlambda, minalpha = minalpha,
                   maxalpha = maxalpha, folds = folds,
                   fold = fold, n_iter = n_iter)
  sur <- mlr3mbo::default_surrogate(inst)
  acq <- mlr3mbo::acqf("ei")
  acq_opt <- mlr3mbo::acqo(bbotk::opt("nloptr", algorithm = "NLOPT_GN_DIRECT_L"),
                           terminator = bbotk::trm("stagnation", threshold = 1e-8))
  optimizer <- bbotk::opt("mbo", loop_function = mlr3mbo::bayesopt_ego,
                          surrogate = sur, acq_function = acq,
                          acq_optimizer = acq_opt)
  if (verbose)
    ret <- optimizer$optimize(inst)
  else
    utils::capture.output(ret <- optimizer$optimize(inst))
  ret
}

# Elastic net objective function for model based optimization

.mbo_obj <- function(object, minlambda = 0, maxlambda = 16, minalpha = 0,
                      maxalpha = 1, folds, fold, n_iter) {
  fn <- function(xs) {
    pars <- unlist(xs)
    list(y = -2 * cvl_tramnet(object, folds = folds, lambda = pars[1], fold = fold,
                              alpha = pars[2])[["logLik_tab"]][["sum_logLik"]][1])
  }
  dom <- paradox::ps(lambda = paradox::p_dbl(lower = minlambda, upper = maxlambda),
                     alpha = paradox::p_dbl(lower = minalpha, upper = maxalpha))
  cod <- paradox::ps(y = paradox::p_dbl(tags = "minimize"))
  obj <- bbotk::ObjectiveRFun$new(fun = fn, domain = dom, codomain = cod)
  bbotk::OptimInstanceSingleCrit$new(objective = obj, search_space = dom,
                                     terminator = bbotk::trm("evals", n_evals = n_iter))
}

#' Fit recommended regularized tram based on model based optimization output
#'
#' @description Extracts the "optimal" tuning parameters from the output of
#' \code{mbo_tramnet} and fits the corresponding \code{tramnet} model.
#'
#' @param mbo_obj Object returned by \code{\link[tramnet]{mbo_tramnet}}.
#' @param m0 Null model of class \code{"tram"}, see
#'     \code{\link[tramnet]{tramnet.tram}}.
#' @param x Design matrix, see \code{\link[tramnet]{tramnet.tram}}.
#' @param ... Additional arguments to \code{tramnet}.
#'
#' @return Object of class \code{"tramnet"}.
#'
#' @export
mbo_recommended <- function(mbo_obj, m0, x, ...) {
  rec_lambda <- ifelse(is.null(mbo_obj[["x"]][["lmb"]]),
                       0, mbo_obj[["x"]][["lmb"]])
  rec_alpha <- ifelse(is.null(mbo_obj[["x"]][["alp"]]),
                      .get_alpha_from_obj(mbo_obj),
                      mbo_obj[["x"]][["alp"]])
  ret <- tramnet(m0, x = x, lambda = rec_lambda, alpha = rec_alpha, ...)
  return(ret)
}

# Helper Functions

.get_alpha_from_obj <- function(mbo_obj) {
  nm <- attr(mbo_obj$final.opt.state$opt.problem$fun, "name")
  ifelse(nm == "lasso_obj", 1, 0)
}
