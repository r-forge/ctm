#' Cross-validating tramnet models
#'
#' @description
#' k-fold cross validation for \code{"tramnet"} objects over a grid of
#'    the tuning parameters based on out-of-sample log-likelihood.
#'
#' @param object Object of class \code{"tramnet"}.
#' @param fold Number of folds for cross validation.
#' @param lambda Values for lambda to iterate over.
#' @param alpha Values for alpha to iterate over.
#' @param folds Manually specify folds for comparison with other methods.
#' @param fit_opt If \code{TRUE}, returns the full model evaluated at optimal
#'     hyper-parameters
#'
#' @return Returns out-of-sample logLik and coefficient estimates for
#'    corresponding folds and values of the hyper-parameters as an object of
#'    class \code{"cvl_tramnet"}
#'
#' @examples
#' \donttest{
#' set.seed(241068)
#' if (require("survival") & require("TH.data")) {
#'   data("GBSG2", package = "TH.data")
#'   X <- 1 * matrix(GBSG2$horTh == "yes", ncol = 1)
#'   colnames(X) <- "horThyes"
#'   GBSG2$surv <- with(GBSG2, Surv(time, cens))
#'   m <- Coxph(surv ~ 1, data = GBSG2, log_first = TRUE)
#'   mt <- tramnet(model = m, x = X, lambda = 0, alpha = 0)
#'   mc <- Coxph(surv ~ horTh, data = GBSG2)
#'   cvl_tramnet(mt, fold = 2, lambda = c(0, 1), alpha = c(0, 1))
#' }
#' }
#'
#' @export
cvl_tramnet <- function(object, fold = 2, lambda = 0, alpha = 0, folds = NULL,
                        fit_opt = FALSE) {
  if (!inherits(object, "tramnet"))
    stop("Cross validation only for models of class tramnet.")
  df <- .get_tramnet_data(object)
  rsp <- variable.names(object$model, "response")
  n <- nrow(df)
  val_grid <- expand.grid(lambda = lambda, alpha = alpha)
  if (is.null(folds) & !is.null(fold)) {
    folds <- sample(rep(1:fold, ceiling(n/fold)), n)
  } else {
    folds <- round(folds)
    fold <- max(folds)
  }
  out <- .cvl_helper(val_grid = val_grid, df = df, rsp = rsp, folds = folds,
                     object = object, fold = fold, n = n)
  raw_ll <- .get_logLik(out, fold = fold)
  ll_tab <- cbind(val_grid, raw_ll)
  raw_cfx <- .get_cfx(out, fold = fold)
  cfx_tab <- lapply(raw_cfx, function(x) cbind(val_grid, x))
  optimal <- ll_tab[which.max(raw_ll$sum_logLik), ]
  ret <- list(
    logLik_tab = ll_tab,
    optimal = optimal,
    coefficients = cfx_tab,
    folds = folds
  )
  if (fit_opt) {
    optmod <- update(object, lambda = optimal$lambda, alpha = optimal$alpha)
    ret$full_fit <- optmod
  }
  class(ret) <- "cvl_tramnet"
  return(ret)
}

# Helper functions

.cvl_helper <- function(val_grid, df, rsp, folds, object, fold, n) {
  apply(val_grid, 1, function(pars) {
    message("Performing ", fold, "-fold cross validation")
    sapply(seq_len(fold), function(x, lmb = pars[1], alp = pars[2]) {
      message("Fold: ", x)
      idx <- which(x == folds)
      trn <- df[-idx, , drop = FALSE]
      tst <- df[idx, , drop = FALSE]
      xtrn <- object$x[-idx, , drop = FALSE]
      xtst <- object$x[idx, , drop = FALSE]
      fit <- update.default(object$model, data = trn)
      trnt <- try(tramnet(model = fit, x = xtrn, lambda = lmb, alpha = alp,
                          check_dcp = FALSE, solver = "ECOS"))
      ncfx <- length(coef(object, tol = 0))
      if (inherits(trnt, "try-error")) {
        list(ll = NA,
             cfx = rep(NA, ncfx))
      } else {
        list(ll = logLik(trnt, newdata = tst),
             cfx = coef(trnt, tol = 0))
      }
    }, simplify = FALSE)
  })
}

.get_logLik <- function(cvl_helper_out, fold) {
  tmp <- sapply(seq_len(fold), function(x) {
    lapply(cvl_helper_out, function(y) y[[x]]$ll)
  }, simplify = FALSE)
  tmp <- lapply(tmp, function(x) do.call(rbind, x))
  ret <- do.call(cbind, tmp)
  colnames(ret) <- paste("logLik_fold", seq_len(fold), sep = "_")
  ret <- data.frame(ret)
  ret$sum_logLik <- apply(ret, 1, sum, na.rm = FALSE)
  return(ret)
}

.get_cfx <- function(cvl_helper_out, fold) {
  tmp <- sapply(seq_len(fold), function(y) {
    lapply(cvl_helper_out, function(x) x[[y]]$cfx)
  }, simplify = FALSE)
  ret <- lapply(tmp, function(x) do.call(rbind, x))
  names(ret) <- paste("coef_fold", seq_len(fold), sep = "_")
  return(ret)
}

.plot_cvl <- function(object, ...) {
  stopifnot(inherits(object, "cvl_tramnet"))
  ll <- object[["logLik_tab"]]
  cfx <- object[["coefficients"]]
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  par(mar = c(5, 5, 2, 2))
  plot(ll$lambda, ll$sum_logLik, ylab = "CV logLik",
       xlab = expression(lambda))
  lapply(cfx, function(cfxx) {
    matplot(x = cfxx[,1], cfxx[,-(1:2)], type = "l",
            xlab = expression(lambda),
            ylab = expression(hat(beta)[j](lambda)), ...)
  })
}
