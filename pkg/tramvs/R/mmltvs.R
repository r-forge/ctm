#' Select optimal subset based on high dimensional BIC in mmlts
#'
#' @inheritParams abess_mmlt
#' @inheritDotParams abess_mmlt
#' @param supp_max maximum support which to call \code{abess_tram} with.
#' @param verbose show progress bar (default: \code{TRUE})
#' @param parallel toggle for parallel computing via
#'     \code{\link[future.apply]{future_lapply}}
#' @param future_args arguments passed to \code{\link[future]{plan}}; defaults
#'     to a \code{"multisession"} with \code{supp_max} workers
#'
#' @details L0-penalized (i.e., best subset selection) multivariate
#'     transformation models using the abess algorithm.
#'
#' @return object of class \code{"mltvs"}, containing the regularization path
#'     (information criterion \code{SIC} and coefficients \code{coefs}), the
#'      best fit (\code{best_fit}) and all other models (\code{all_fits})
#'
#' @examples
#'
#' @export
mmltVS <- function(mltargs, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, verbose = TRUE, parallel = FALSE,
                   future_args = list(strategy = "multisession",
                                      workers = supp_max), ...) {

  if (is.null(supp_max)) {
    m0 <- do.call("mmlt", mltargs)
    supp_max <- length(ncfs <- coef_mmlt(m0))
  }

  if (verbose & interactive() & !parallel)
    pb <- txtProgressBar(style = 3, width = 50, min = 0, max = supp_max)

  if (parallel) {
    do.call(plan, future_args)
    this.lapply <- \(...) future_lapply(..., future.seed = TRUE)
  } else {
    this.lapply <- lapply
  }

  res <- this.lapply(seq_len(supp_max), \(ts) {
    if (verbose & interactive())
      setTxtProgressBar(pb, ts)
    fit <- abess_mmlt(mltargs, supp = ts, k_max = k_max, thresh = thresh,
                      init = init, m_max = m_max, m0 = m0)
    list(
      fit = fit,
      SIC = -logLik(fit$m) + length(fit$A) * log(length(coef(fit$m))) *
        log(log(nrow(fit$m$data)))
    )
  })

  fits <- lapply(res, \(x) x[[1]])
  SIC <- unlist(lapply(res, \(x) x[[2]]))

  allcf <- structure(rep(0, length(ncfs)), names = names(ncfs))
  traj <- as(do.call("cbind", lapply(fits, \(fit) {
    tnms <- names(tcfx <- coef_mmlt(fit))
    allcf[tnms] <- tcfx
    allcf
  })), "sparseMatrix")
  colnames(traj) <- seq_len(supp_max)

  structure(list(SIC = data.frame(supp = seq_len(supp_max), SIC = SIC),
                 coefs = traj,
                 best_fit = fits[[which.min(SIC)]],
                 all_fits = fits), class = c("mmltvs", "tramvs"))
}

coef_mmlt <- function(obj) {
  all <- coef(obj, type = "all")
  mar <- unlist(unname(coef(obj, type = "marginal")))
  all[setdiff(names(all), names(mar))]
}

#' Optimal subset selection for multivariate transformation models
#'
#' @inheritParams abess_tram
#' @param mltargs Arguments passed to \code{mmlt}
#' @param ... Currently ignored
#'
#' @return List containing the fitted model via \code{mmlt}, active set
#'     \code{A} and inactive set \code{I}.
#'
#' @examples
#'
#' @export
abess_mmlt <- function(mltargs, supp, k_max = supp, thresh = NULL, init = TRUE,
                       m_max = 10, m0 = NULL, ...) {

  if (is.null(k_max))
    k_max <- supp

  stopifnot(k_max <= supp)

  if (is.null(m0))
    m0 <- do.call("mmlt", mltargs)

  ncfs <- names(coef_mmlt(m0))
  p <- length(ncfs)
  n <- nrow(m0$data)

  if (is.null(thresh))
    thresh <- 0.01 * supp * log(p) * log(log(n)) / n

  cors <- cor_init.mmlt(m0, ncfs)

  if (init)
    A0 <- ncfs[.a0_init(cors, supp)]
  else
    A0 <- ncfs[sample.int(ceiling(length(ncfs) / 2), 1)]

  I0 <- setdiff(ncfs, A0)
  fix0 <- numeric(length(I0))
  names(fix0) <- I0
  mltargs$fixed <- fix0
  m0 <- do.call("mmlt", mltargs)

  sm <- s0 <- .splice_mmlt(mltargs, m0, A0, I0, k_max, thresh)

  if (length(s0$A) == length(A0) && all(s0$A == A0)) {
    return(structure(s0, class = c("abess_mmlt", "abess_tram")))
  }

  for (m in seq_len(m_max)) {
    Am <- sm$A
    sm <- .splice_mmlt(mltargs, sm$mod, sm$A, sm$I, k_max, thresh)
    if (length(sm$A) == length(Am) && all(sm$A == Am))
      return(structure(s0, class = c("abess_mmlt", "abess_tram")))
    else
      return(structure(sm, class = c("abess_mmlt", "abess_tram")))
  }
}

#' Method for computing correlations in mmlts
#' @inheritParams cor_init
#' @return Vector of correlation for initializing the active set
#' @exportS3Method cor_init mmlt
cor_init.mmlt <- function(m0, mb) {
  cors <- cor(do.call("cbind", lapply(m0$models$models, residuals)))
  L <- solve(t(chol(cors)))
  LL <- ltMatrices(L[lower.tri(L, diag = TRUE)], diag = TRUE)
  LLL <- invcholD(LL, D = 1 / diagonals(LL))
  structure(c(abs(Lower_tri(LLL))), names = mb)
}

# Helper ------------------------------------------------------------------

.splice_mmlt <- function(args, m, A, I, k_max, thresh, ...) {
  m0 <- m
  A0 <- A
  I0 <- I
  L <- L0 <- - logLik(m) / nrow(m$data)
  cf <- cf0 <- coef_mmlt(m0)
  ncfs <- names(cf)
  cfA <- cf[names(cf) %in% A0]
  cfI <- cf[names(cf) %in% I0]

  bwd_sacrifice <- sapply(seq_along(cfA), \(parm) {
    ncfs <- c(cfA[parm], cfI)
    ncfs[] <- 0
    args$fixed <- ncfs
    m_retrained <- do.call("mmlt", args)
    nll_wo <- - logLik(m_retrained) / nrow(m_retrained$data)
    nll_wo - L
  })

  fwd_sacrifice <- sapply(seq_along(cfI), \(parm) {
    ncfs <- c(cfI, cfA)
    ncfs[names(cfI)] <- 0
    args$fixed <- ncfs[-parm]
    m_retrained <- do.call("mmlt", args)
    nll_wo <- - logLik(m_retrained) / nrow(m$data)
    L - nll_wo
  })

  for (k in seq_len(k_max)) {
    Ak <- ncfs[.ak_compute(bwd_sacrifice, k)]
    Ik <- ncfs[.ik_compute(fwd_sacrifice, k)]

    newA <- sort(union(setdiff(A, Ak), Ik))
    newI <- setdiff(ncfs, newA)

    if (length(newI) == length(I) && all(sort(newI) == sort(I)) |
        length(newI) == length(I0) && all(sort(newI) == sort(I0)) |
        length(newA) > k_max)
      next

    newcfs <- numeric(length(newI))
    names(newcfs) <- newI
    args$fixed <- newcfs
    newm <- do.call("mmlt", args)
    newL <- -logLik(newm) / nrow(m$data)

    if (L > newL) {
      cf <- coef(newm, with_baseline = TRUE)
      cfs <- coef(newm)
      cfb <- cf0[!names(cf) %in% names(cfs)]
      L <- newL
      A <- newA
      I <- newI
    }
  }

  if (L0 - L < thresh)
    ret <- list(mod = m, A = A, I = I)
  else
    ret <- list(mod = m0, A = A0, I = I0)

  ret
}
