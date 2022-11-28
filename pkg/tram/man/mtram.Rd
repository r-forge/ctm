\name{mtram}
\alias{mtram}
\title{
  Transformation Models for Clustered Data
}
\description{
  Marginally interpretable transformation models for clustered data.
}
\usage{
mtram(object, formula, data,
      grd = SparseGrid::createSparseGrid(type = "KPU", 
                dimension = length(rt$cnms[[1]]), k = 10), 
      Hessian = FALSE,  tol = .Machine$double.eps, ...)
}
\arguments{
  \item{object}{
  A \code{tram} object.
}
  \item{formula}{
  A formula specifying the random effects.
}
  \item{data}{
  A data frame. 
}
  \item{grd}{
  A sparse grid used for numerical integration to get the likelihood.
}
  \item{Hessian}{
  A logical, if \code{TRUE}, the hessian is computed and returned.
}
  \item{tol}{numerical tolerance.}
  \item{\dots}{
  Additional argument.
}
}
\details{
  A Gaussian copula with a correlation structure obtained from a random
  intercept or random intercept / random slope model (that is, clustered or
  longitudinal data can by modelled only) is used to capture the
  correlations whereas the marginal distributions are described by a 
  transformation model. The methodology is described in Barbanti and Hothorn
  (2022) and examples are given in the \code{mtram} package vignette.

  This is a proof-of-concept implementation. Only \code{coef()} and 
  \code{logLik()} methods are available at the moment.
}
\value{
  An object of class \code{tram} with \code{coef()} and \code{logLik()}
  methods.
}
\references{

  Luisa Barbanti and Torsten Hothorn (2023). A Transformation Perspective on 
  Marginal and Conditional Models, \emph{Biostatistics}, \doi{10.48550/arXiv.1910.09219}.

}
\examples{

  ### For illustrations see
  ## vignette("mtram", package = "tram")
  ## or
  ## demo("mtram", package = "tram")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
