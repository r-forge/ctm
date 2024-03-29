% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate.tramME}
\alias{simulate.tramME}
\title{Simulate from a \code{tramME} model}
\usage{
\method{simulate}{tramME}(
  object,
  nsim = 1,
  seed = NULL,
  newdata = model.frame(object),
  type = c("ranef", "response", "joint"),
  ...
)
}
\arguments{
\item{object}{A \code{tramME} object.}

\item{nsim}{number of samples to generate}

\item{seed}{optional seed for the random number generator}

\item{newdata}{an optional data frame of observations}

\item{type}{Defaults to \code{"ranef"}. Currently the only avalable option.}

\item{...}{Additional arguments, passed to \code{\link[mlt]{simulate.mlt}}.}
}
\value{
A length \code{nsim} list of draws.
}
\description{
Simulate from a \code{tramME} model
}
\section{Warning}{


This method is under active development and may be subject to change. It is
  currently limited to simulating random effects.
}

\examples{
data("sleepstudy", package = "lme4")
fit <- BoxCoxME(Reaction ~ Days + (Days | Subject), data = sleepstudy)
sim <- simulate(fit, nsim = 10, seed = 123)
}
