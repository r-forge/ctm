robust_score_test <- function(object, ...) {
  UseMethod("robust_score_test")
}

robust_score_test.tram <- function(
    object, parm = names(coef(object)),
    alternative = c("two.sided", "less", "greater"),
    nullvalue = 0, confint = FALSE, level = 0.95,
    ranger_args = NULL, ...) {
  cf <- coef(object)
  stopifnot(all(parm %in% names(cf)))
  alternative <- match.arg(alternative)

  if (length(parm) > 1) {
    ret <- lapply(parm, robust_score_test.tram,
      object = object,
      alternative = alternative,
      nullvalue = nullvalue,
      confint = confint,
      level = level,
      ranger_args = ranger_args,
      ...
    )
    names(ret) <- parm
    class(ret) <- "htests"
    return(ret)
  }

  ### Handle scale terms
  which_residuals <- if (grepl("^scl_", parm)) "scaling" else "shifting"

  hypothesis <- structure(nullvalue, names = parm)
  m <- mlt::as.mlt(object)

  ### Residualize
  rX <- .residualize(object, m, parm, ranger_args)

  ### Compute score statistic under H_0: b = alt
  scfun <- function(alt) {
    names(alt) <- parm
    cf <- stats::coef(stats::update(m, fixed = alt))
    coef(m) <- cf
    rY <- residuals(m, what = which_residuals)
    RR <- c(rY) * c(rX)
    nn <- length(RR)
    R.sq <- RR^2
    meanR <- mean(RR)
    stat <- sqrt(nn) * meanR / sd(RR)
    pval <- switch(alternative,
      "two.sided" = 2 * stats::pnorm(-abs(stat)),
      "greater" = 1 - stats::pnorm(stat),
      "less" = stats::pnorm(stat)
    )
    list(stat = stat, pval = pval, rY = rY)
  }

  h0 <- scfun(hypothesis)
  stat <- h0$stat
  pval <- h0$pval

  parameter <- paste(switch(class(object)[1],
    "Colr" = "Log-odds ratio",
    "Coxph" = "Log-hazard ratio",
    "Lm" = "Standardised difference",
    "Lehmann" = "Lehmann parameter",
    "BoxCox" = "Standardised difference"
  ))
  parameter <- paste(tolower(parameter), "for", parm)
  names(nullvalue) <- parameter

  ci <- NULL
  ci <- if (confint) {
    tryCatch(
      .ci(scfun, level, m, alternative, parm, TRUE),
      error = function(e) {
        warning("Cannot invert score function numerically for paramerer `", parm, "`")
        NULL
      }
    )
  }

  if (!is.null(ci$est)) {
    names(ci$est) <- parameter
  }

  structure(
    list(
      statistic = c("Z" = stat), p.value = pval,
      null.value = nullvalue, alternative = alternative,
      method = paste0("Doubly-robust Transformation Score Test"),
      data.name = paste0(deparse(object$call), collapse = "\n"),
      hypothesis = hypothesis, conf.int = ci$ci, estimate = ci$est,
      rY = h0$rY, rX = rX
    ),
    class = c("tramgcm", "htest")
  )
}

# Helpers -----------------------------------------------------------------

.tmm <- function(object) {
  shi <- model.matrix(object, what = "shifting")
  scl <- tryCatch(model.matrix(object, what = "scaling"), error = function(e) NULL)
  int <- tryCatch(model.matrix(object, what = "interacting"), error = function(e) NULL)
  nshi <- colnames(shi)
  nscl <- colnames(scl)
  scl_rm <- intersect(paste0("scl_", nshi), nscl)
  rm <- numeric(0)
  if (!identical(scl_rm, character(0)) & !is.null(scl_rm)) {
    scl <- scl[, -grep(scl_rm, nscl)]
  }
  cbind(shi, scl, int)
}

.residualize <- function(object, m, parm, ranger_args = NULL) {
  mm <- .tmm(object)
  if (any(grepl("scl_", parm)) & !any(grepl(parm, colnames(mm)))) {
    parm <- substr(parm, start = 5, stop = nchar(parm))
  }
  nx <- mm[, parm, drop = FALSE]
  ox <- mm[, !colnames(mm) %in% parm, drop = FALSE]
  if (length(unique(nx)) == 2) {
    ranger_args$probability <- TRUE
  }

  ### Residualize for GCM
  if (NCOL(ox) == 0) {
    tx <- scale(nx, center = TRUE, scale = FALSE)
  } else {
    if (requireNamespace("ranger")) {
      ranger <- ranger::ranger
      rf <- do.call("ranger", c(list(y = nx, x = ox), ranger_args))
      preds <- stats::predict(rf, data = ox)$predictions
    } else {
      warning("Package 'ranger' required for residualization. Using `lm()` instead.")
      preds <- predict(stats::lm(nx ~ ox))
    }
    if (NCOL(preds) > 1) {
      preds <- preds[, -1]
    }
    tx <- nx - preds
  }
  tx
}

#' @importFrom stats confint qnorm uniroot
.ci <- function(scfun, level, m, alternative, parm, return_est = FALSE) {
  alpha <- 1 - level
  if (alternative == "two.sided") alpha <- alpha / 2
  wci <- try(stats::confint(m, level = 1 - alpha / 5)[parm, ], silent = TRUE)
  if (inherits(wci, "try-error")) {
    wci <- c(-4, 4)
  }
  lwr <- if (alternative == "less") {
    -Inf
  } else {
    stats::uniroot(function(x) scfun(x)$stat - stats::qnorm(alpha),
      lower = wci[1],
      upper = wci[2], extendInt = "yes"
    )$root
  }
  upr <- if (alternative == "greater") {
    Inf
  } else {
    stats::uniroot(function(x) scfun(x)$stat - stats::qnorm(1 - alpha),
      lower = wci[1],
      upper = wci[2], extendInt = "yes"
    )$root
  }
  est <- NULL
  if (return_est) {
    est <- stats::uniroot(function(x) scfun(x)$stat,
      lower = wci[1], upper = wci[2],
      extendInt = "yes"
    )$root
    names(est) <- parm
  }
  ci <- c(lwr, upr)
  attr(ci, "conf.level") <- level
  list(ci = ci, est = est)
}
