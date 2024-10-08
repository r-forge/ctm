% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mmltvs.R
\name{abess_mmlt}
\alias{abess_mmlt}
\title{Optimal subset selection for multivariate transformation models}
\usage{
abess_mmlt(
  mltargs,
  supp,
  k_max = supp,
  thresh = NULL,
  init = TRUE,
  m_max = 10,
  m0 = NULL,
  ...
)
}
\arguments{
\item{mltargs}{Arguments passed to \code{mmlt}}

\item{supp}{support size of the coefficient vector}

\item{k_max}{maximum support size to consider during the splicing algorithm.
Defaults to \code{supp}.}

\item{thresh}{threshold when to stop splicing. Defaults to
0.01 * \code{supp} * p * log(log(n)) / n$, where p denotes the number of predictors
and n the sample size.}

\item{init}{initialize active set. Defaults to \code{TRUE} and initializes the
active set with those covariates that are most correlated with score residuals
of an unconditional \code{modFUN(update(formula, . ~ 1))}.}

\item{m_max}{maximum number of iterating the splicing algorithm.}

\item{m0}{Transformation model for initialization}

\item{...}{Currently ignored}
}
\value{
List containing the fitted model via \code{mmlt}, active set
    \code{A} and inactive set \code{I}.
}
\description{
Optimal subset selection for multivariate transformation models
}
