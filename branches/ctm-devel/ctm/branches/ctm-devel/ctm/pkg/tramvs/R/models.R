# Aliases for standard trams

#' Optimal subset selection in a BoxCox-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{BoxCox}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram BoxCox
#' @export
BoxCoxVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                     init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, BoxCox, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Colr-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{Colr}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Colr
#' @export
ColrVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Colr, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Coxph-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{Coxph}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Coxph
#' @export
CoxphVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Coxph, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in an Lm-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{Lm}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Lm
#' @export
LmVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Lm, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Lehmann-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Lehmann}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Lehmann
#' @export
LehmannVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                      init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Lehmann, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Polr-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{Polr}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Polr
#' @export
PolrVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Polr, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Survreg model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[tram]{Survreg}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom tram Survreg
#' @export
SurvregVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                      init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, Survreg, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a cotram model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{\link[cotram]{cotram}}
#' @return See \code{\link[tramvs]{tramvs}}
#' @importFrom cotram cotram
#' @export
cotramVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                     init = TRUE, m_max = 10, parallel = FALSE,
                     future_args = list(strategy = "multisession",
                                        workers = supp_max), ...) {
  tramvs(formula, data, cotram, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}
