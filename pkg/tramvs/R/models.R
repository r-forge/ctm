# Aliases for standard trams

#' Optimal subset selection in a BoxCox-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{BoxCox()}
#' @importFrom tram BoxCox
#' @export
BoxCoxVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                     init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, BoxCox, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Colr-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Colr()}
#' @importFrom tram Colr
#' @export
ColrVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Colr, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Coxph-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Coxph()}
#' @importFrom tram Coxph
#' @export
CoxphVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Coxph, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in an Lm-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Lm()}
#' @importFrom tram Lm
#' @export
LmVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Lm, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Lehmann-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Lehmann()}
#' @importFrom tram Lehmann
#' @export
LehmannVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                      init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Lehmann, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Polr-type transformation model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Polr()}
#' @importFrom tram Polr
#' @export
PolrVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                   init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Polr, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}

#' Optimal subset selection in a Survreg model
#' @inheritParams tramvs
#' @param ... Additional arguments supplied to \code{Survreg()}
#' @importFrom tram Survreg
#' @export
SurvregVS <- function(formula, data, supp_max = NULL, k_max = NULL, thresh = NULL,
                      init = TRUE, m_max = 10, ...) {
  tramvs(formula, data, Survreg, supp_max, k_max, thresh, init, m_max,
         ... = ...)
}