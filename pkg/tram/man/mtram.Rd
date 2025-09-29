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
      tol = .Machine$double.eps, optim = mltoptim(hessian = TRUE),
      ...)
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
  \item{tol}{numerical tolerance.}
\item{optim}{
a list of optimisers as returned by \code{\link[mlt]{mltoptim}}
}
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

  Only \code{coef()} and \code{logLik()} methods are available at the
  moment, see \code{vignette("mtram", package = "tram")} for worked
  examples.
}
\value{
  An object of class \code{tram} with \code{coef()} and \code{logLik()}
  methods.
}
\seealso{
  \code{vignette("mtram", package = "tram")}
}
\references{

  Luisa Barbanti and Torsten Hothorn (2024). A Transformation Perspective on 
  Marginal and Conditional Models, \emph{Biostatistics}, \bold{25}(2),
  402--428, \doi{10.1093/biostatistics/kxac048}.

}
\examples{

  if (require("lme4")) {
      ### linear mixed model
      sleep_lmer <- lmer(Reaction ~ Days + (Days | Subject), 
                         data = sleepstudy, REML = FALSE)

      ### marginal transformation model
      sleep_LM <- Lm(Reaction ~ Days, data = sleepstudy)
      sleep_LMmer <- mtram(sleep_LM, ~ (Days | Subject), data = sleepstudy)

      ### the same
      logLik(sleep_lmer)
      logLik(sleep_LMmer)

      ### Lm / mtram estimate standardised effects
      sdinv <- 1 / summary(sleep_lmer)$sigma
      fixef(sleep_lmer) * c(-1, 1) * sdinv
      coef(sleep_LMmer)[c("(Intercept)", "Days")]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
