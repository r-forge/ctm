% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Resp.R
\name{Resp}
\alias{Resp}
\alias{R.Resp}
\alias{print.Resp}
\alias{[.Resp}
\alias{is.na.Resp}
\alias{length.Resp}
\alias{format.Resp}
\title{Response objects}
\usage{
Resp(
  cleft,
  cright,
  tleft,
  tright,
  bounds = c(-Inf, Inf),
  open_lwr_bnd = TRUE,
  tol = sqrt(.Machine$double.eps)
)

\method{R}{Resp}(object, ...)

\method{print}{Resp}(x, ...)

\method{[}{Resp}(x, i, j, drop = FALSE)

\method{is.na}{Resp}(x)

\method{length}{Resp}(x)

\method{format}{Resp}(x, ...)
}
\arguments{
\item{cleft}{A vector of left borders of censoring intervals}

\item{cright}{A vector of right borders of censoring intervals}

\item{tleft}{A vector of left truncation values}

\item{tright}{A vector of right truncation values}

\item{bounds}{An optional numeric vector of two elements (\code{c(a, b)})
that denotes the lower and upper boundaries of the outcome.}

\item{open_lwr_bnd}{Logical; if \code{TRUE}, the lower boundary of the
outcome is open, and we want to enforce this.}

\item{tol}{Tolerance level.}

\item{object}{A \code{Resp} object}

\item{...}{Optional arguments}

\item{x}{A \code{Resp} object}

\item{i}{Row index (typically the only index)}

\item{j}{Column index (typically missing)}

\item{drop}{If \code{TRUE} the result is coerced to the lowest possible dimension}
}
\value{
A \code{Resp} object or a \code{Surv} object
}
\description{
Response objects to represent censored and truncated observations
}
\details{
\code{Resp} extends the functionality of \code{\link[survival]{Surv}} class
  by allowing cases that cannot be defined with it. An example is an
  interval-censored outcome with left truncatation (see Examples).

Censored and exactly observed data can be defined similarly to \code{type =
  "interval2"} objects in \code{\link[survival]{Surv}}. \code{NA} values for
  left or right censoring borders mean left- or right-censored observations,
  respectively. If both borders are \code{NA}, the observation is considered
  \code{NA} by \code{is.na()}.  Truncation times (\code{tleft} and
  \code{tright} arguments) can be omitted or take \code{NA} values, which
  means no truncation. If only the censoring intervals are provided, i.e.,
  no trunction is present, the function returns a \code{Surv} object.

\code{Resp} also provides a limited interface between \code{tramME} and the
  \code{response} class (technically, inherits from it) of \code{mlt} (see
  \code{\link[mlt]{R}}), which uses an internal representation that is not
  compatible with \code{tramME}.

The optional argument \code{open_lwr_bnd} can be used to enforce lower
  boundaries of the outcome. Left boundaries in the \code{Resp} object
  (\code{cleft} and \code{tleft}) that are equal to the first element of
  \code{bounds} will be increased with one \code{tol} value to avoid
  downstream numerical problems in \code{mlt}. This adjustment is recorded
  and reversed when we print the object.
}
\section{Methods (by generic)}{
\itemize{
\item \code{R(Resp)}: Converting \code{Resp} objects to \code{response} (from
\code{mlt}) objects (see \code{\link[mlt]{R}})

\item \code{print(Resp)}: Print method for the \code{Resp} class

\item \code{[}: Subsetting \code{Resp} objects

\item \code{is.na(Resp)}: Missing values

\item \code{length(Resp)}: Length of a \code{Resp} object

\item \code{format(Resp)}: \code{format} method for a \code{Resp} object

}}
\section{Warning}{


This function is experimental and currently limited to continuous outcome
  types. It may be subject to change.
}

\examples{

dat <- data.frame(x1 = 1:10, x2 = c(2:10, NA), x3 = c(NA, 0:8))
dat$r <- with(dat, Resp(x1, x2, x3))

dat$r
dat[1:3, ]$r
dat$r[1:3]

is.na(dat$r)

model.frame(r ~ 1, data = dat, na.action = na.omit)

}
