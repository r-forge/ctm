#' Optimal subset selection for transformation models
#'
#' @param formula object of class \code{"formula"}.
#' @param data data frame containing the variables in the model.
#' @param modFUN function for fitting a transformation model, e.g., \code{BoxCox()}.
#' @param mandatory formula of mandatory covariates, which will always be included
#'     and estimated in the model. Note that this also changes the intialization
#'     of the active set. The active set is then computed with regards to the
#'     model residuals of \code{modFUN(mandatory, ...)} instead of the unconditional
#'     model.
#' @param supp support size of the coefficient vector
#' @param k_max maximum support size to consider during the splicing algorithm.
#'     Defaults to \code{supp}.
#' @param thresh threshold when to stop splicing. Defaults to
#'     0.01 * \code{supp} * p * log(log(n)) / n$, where p denotes the number of predictors
#'     and n the sample size.
#' @param init initialize active set. Defaults to \code{TRUE} and initializes the
#'     active set with those covariates that are most correlated with score residuals
#'     of an unconditional \code{modFUN(update(formula, . ~ 1))}.
#' @param m_max maximum number of iterating the splicing algorithm.
#' @param m0 Transformation model for initialization
#' @param ... additional arguments supplied to \code{modFUN}.
#'
#' @return List containing the fitted model via \code{modFUN}, active set
#'     \code{A} and inactive set \code{I}.
#'
#' @examples
#' set.seed(24101968)
#' library(tramvs)
#'
#' N <- 1e2
#' P <- 5
#' nz <- 3
#' beta <- rep(c(1, 0), c(nz, P - nz))
#' X <- matrix(rnorm(N * P), nrow = N, ncol = P)
#' Y <- 1 + X %*% beta + rnorm(N)
#'
#' dat <- data.frame(y = Y, x = X)
#'
#' abess_tram(y ~ ., dat, modFUN = Lm, supp = 3)
#'
#' @importFrom stats coef update residuals cor model.matrix
#' @export
abess_tram <- function(formula, data, modFUN, supp, mandatory = NULL, k_max = supp,
                       thresh = NULL, init = TRUE, m_max = 10, m0 = NULL, ...) {
  if (is.null(k_max))
    k_max <- supp

  stopifnot(k_max <= supp)

  if (is.null(m0))
    m0 <- modFUN(formula, data, ... = ...)

  theta_init <- m0$theta

  ncfs <- names(coef(m0))
  p <- length(ncfs)
  n <- nrow(m0$data)

  if (is.null(thresh))
    thresh <- 0.01 * supp * log(p) * log(log(n)) / n

  ### Model for initialization
  if (is.null(mandatory)) {
    fmb <- update(formula, . ~ 1)
  } else {
    fmb <- mandatory
  }
  mb <- modFUN(fmb, data, ... = ...)
  mcfs <- names(coef(mb))

  cors <- cor_init(m0, mb)

  if (init)
    A0 <- ncfs[.a0_init(cors, supp)]
  else
    A0 <- ncfs[sample.int(ceiling(length(ncfs) / 2), 1)]
  if (!is.null(mcfs))
    A0 <- sort(union(mcfs, A0))

  I0 <- setdiff(ncfs, A0)
  fix0 <- numeric(length(I0))
  names(fix0) <- I0
  m0 <- modFUN(formula, data, fixed = fix0,
               theta = theta_init[!names(theta_init) %in% I0],
               ... = ...)

  sm <- s0 <- .splicing(m0, A0, I0, k_max, thresh, modFUN, formula, data,
                        mcfs = mcfs, theta_init = theta_init, ... = ...)

  if (length(s0$A) == length(A0) && all(s0$A == A0)) {
    return(structure(s0, class = "abess_tram"))
  }

  for (m in seq_len(m_max)) {
    Am <- sm$A
    sm <- .splicing(sm$mod, sm$A, sm$I, k_max, thresh, modFUN, formula, data,
                    mcfs = mcfs, theta_init = theta_init, ... = ...)
    if (length(sm$A) == length(Am) && all(sm$A == Am))
      return(structure(s0, class = "abess_tram"))
    else
      return(structure(sm, class = "abess_tram"))
  }
}

#' Compute correlation for initializing the active set
#' @param m0 \code{modFUN(formula, data)}
#' @param mb \code{modFUN(mandatory, data)}
#' @return Vector of correlations for initializing the active set, depends on
#'     type of model (see e.g. \code{\link[tramvs]{cor_init.default}})
#' @export
cor_init <- function(m0, mb) {
  UseMethod("cor_init")
}

#' Default method for computing correlation
#' @inheritParams cor_init
#' @return Vector of correlation for initializing the active set
#' @exportS3Method cor_init default
cor_init.default <- function(m0, mb) {
  res <- residuals(mb)
  mm <- model.matrix(m0)
  cors <- abs(c(cor(res, mm)))
  cors[is.na(cors)] <- 0
  cors
}

#' Shit-scale tram method for computing correlation
#' @inheritParams cor_init
#' @return Vector of correlations for initializing the active set, includes both
#'     shift and scale residuals
#' @exportS3Method cor_init stram
cor_init.stram <- function(m0, mb) {
  mshift <- model.matrix(m0, what = "shifting")
  mscale <- model.matrix(m0, what = "scaling")
  rshift <- residuals(mb, what = "shifting")
  rscale <- residuals(mb, what = "scaling")
  cors <- abs(c(cor(rshift, mshift), cor(rscale, mscale)))
  cors[is.na(cors)] <- 0
  cors
}

#' Select optimal subset based on high dimensional BIC
#'
#' @inheritParams abess_tram
#' @inheritDotParams abess_tram
#' @param supp_max maximum support which to call \code{abess_tram} with.
#' @param verbose show progress bar (default: \code{TRUE})
#' @param parallel toggle for parallel computing via
#'     \code{\link[future.apply]{future_lapply}}
#' @param future_args arguments passed to \code{\link[future]{plan}}; defaults
#'     to a \code{"multisession"} with \code{supp_max} workers
#'
#' @details L0-penalized (i.e., best subset selection) transformation models
#'     using the abess algorithm.
#'
#' @return object of class \code{"tramvs"}, containing the regularization path
#'     (information criterion \code{SIC} and coefficients \code{coefs}), the
#'      best fit (\code{best_fit}) and all other models (\code{all_fits})
#'
#' @examples
#' set.seed(24101968)
#' library("tramvs")
#'
#' N <- 1e2
#' P <- 5
#' nz <- 3
#' beta <- rep(c(1, 0), c(nz, P - nz))
#' X <- matrix(rnorm(N * P), nrow = N, ncol = P)
#' Y <- 1 + X %*% beta + rnorm(N)
#'
#' dat <- data.frame(y = Y, x = X)
#' res <- tramvs(y ~ ., data = dat, modFUN = Lm)
#' plot(res, type = "s")
#' plot(res, which = "path")
#'
#' @importFrom methods as
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @export
tramvs <- function(formula, data, modFUN, mandatory = NULL, supp_max = NULL,
                   k_max = NULL, thresh = NULL, init = TRUE, m_max = 10,
                   m0 = NULL, verbose = TRUE, parallel = FALSE,
                   future_args = list(strategy = "multisession",
                                      workers = supp_max), ...) {
  if (is.null(supp_max)) {
    m0 <- modFUN(formula, data, ... = ...)
    supp_max <- length(coef(m0))
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
    fit <- abess_tram(formula = formula, data = data, modFUN = modFUN,
                      mandatory = mandatory, supp = ts, k_max = k_max,
                      thresh = thresh, init = init, m_max = m_max, m0 = m0,
                      ... = ...)
    list(
      fit = fit,
      SIC = -logLik(fit$m) + length(fit$A) * log(length(coef(fit$m))) *
        log(log(nrow(fit$m$data)))
    )
  })

  fits <- lapply(res, \(x) x[[1]])
  SIC <- unlist(lapply(res, \(x) x[[2]]))

  traj <- as(do.call("cbind", lapply(fits, coef, ... = ...)), "sparseMatrix")
  colnames(traj) <- seq_len(supp_max)

  structure(list(SIC = data.frame(supp = seq_len(supp_max), SIC = SIC),
                 coefs = traj,
                 best_fit = fits[[which.min(SIC)]],
                 all_fits = fits), class = "tramvs")
}

# Helpers

.a0_init <- function(cors, supp) {
  corgrid <- outer(cors, cors, `>=`)
  which(colSums(corgrid) <= supp)
}

#' @importFrom stats logLik coef
.splicing <- function(m, A, I, k_max, thresh, modFUN, formula, data, mcfs,
                      theta_init, ...) {
  m0 <- m
  A0 <- A
  I0 <- I
  L <- L0 <- - logLik(m) / nrow(m$data)
  cf <- cf0 <- coef(m, with_baseline = TRUE)
  cfs <- cfs0 <- coef(m)
  cfb <- cfb0 <- cf0[!names(cf0) %in% names(cfs0)]
  ncfs <- names(cfs0)
  cfA <- cfs[names(cfs) %in% A0]
  cfI <- cfs[names(cfs) %in% I0]

  bwd_sacrifice <- sapply(seq_along(cfA), \(parm) {
    ncfs <- cfs
    ncfs[parm] <- 0
    m_retrained <- modFUN(formula, data, fixed = ncfs,
                          theta = theta_init[!names(theta_init) %in% names(ncfs)],
                          ... = ...)
    nll_wo <- - logLik(m_retrained) / nrow(m$data)
    nll_wo - L
  })

  fwd_sacrifice <- sapply(seq_along(cfI), \(parm) {
    ncfs <- cfs[-parm]
    m_retrained <- modFUN(formula, data, fixed = ncfs,
                          theta = theta_init[!names(theta_init) %in% names(ncfs)],
                          ... = ...)
    nll_wo <- - logLik(m_retrained) / nrow(m$data)
    L - nll_wo
  })

  for (k in seq_len(k_max)) {
    Ak <- ncfs[.ak_compute(bwd_sacrifice, k)]
    Ik <- ncfs[.ik_compute(fwd_sacrifice, k)]

    newA <- sort(union(setdiff(A, Ak), Ik))
    if (!is.null(mcfs))
      newA <- sort(union(newA, mcfs))
    newI <- setdiff(ncfs, newA)

    if (length(newI) == length(I) && all(sort(newI) == sort(I)) |
        length(newI) == length(I0) && all(sort(newI) == sort(I0)) |
        length(newA) > k_max)
      next

    newcfs <- numeric(length(newI))
    names(newcfs) <- newI
    newm <- modFUN(formula, data, fixed = newcfs,
                   theta = theta_init[!names(theta_init) %in% names(newcfs)],
                   ... = ...)
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

.ak_compute <- function(xi, k) {
  which(colSums(outer(xi, xi, `>=`)) <= k)
}

.ik_compute <- function(zeta, k) {
  which(colSums(outer(zeta, zeta, `<=`)) <= k)
}
