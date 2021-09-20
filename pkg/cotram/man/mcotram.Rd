\name{mcotram}
\alias{mcotram}
\title{
  Multivariate Count Conditional Transformation Models
}
\description{
  A proof-of-concept implementation of multivariate conditional
  transformation models for count data.
}
\usage{
mcotram(..., formula = ~ 1, data, theta = NULL,
        control.outer = list(trace = FALSE), tol = sqrt(.Machine$double.eps), 
        dofit = TRUE)
}
\arguments{
  \item{\dots}{marginal count transformation models, one for each response
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters. The default is set to \code{~ 1} for constant lambdas.
}
  \item{data}{a data.frame
}
  \item{theta}{an optional vector of starting values
}
  \item{control.outer}{a list controlling \code{\link[alabama]{auglag}}
}
  \item{tol}{tolerance}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned
}
}
\details{
  The function implements multivariate count conditional transformation models.
  The response is assumed to be a vector of counts.
}
\value{
An object of class \code{mcotram} and \code{mmlt} with \code{coef} and 
\code{predict} methods.
}
\references{

  Luisa Barbani, Roland Brandl, Torsten Hothorn (2021), Multi-species Count
  Transformation Models. Submitted manuscript, available from the authors.

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2020),
  Multivariate Conditional Transformation Models. \emph{Scandinavian Journal
  of Statistics}, \doi{10.1111/sjos.12501}.

}
\examples{

  op <- options(digits = 2)

  data("spiders", package = "cotram")

  ## fit conditional marginal count transformation models
  m_PF <- cotram(Pardosa_ferruginea ~ Elevation + Canopy_openess, 
                 data = spiders, method = "probit")
  m_HL <- cotram(Harpactea_lepida ~ Elevation + Canopy_openess,
                 data = spiders, method = "probit")
  m_CC <- cotram(Callobius_claustrarius ~ Elevation + Canopy_openess,
                 data = spiders, method = "probit")
  m_CT <- cotram(Coelotes_terrestris ~ Elevation + Canopy_openess,
                 data = spiders, method = "probit")
  m_PL <- cotram(Pardosa_lugubris ~ Elevation + Canopy_openess,
                 data = spiders, method = "probit")
  m_PR <- cotram(Pardosa_riparia ~ Elevation + Canopy_openess,
                 data = spiders, method = "probit")

  ## fit multi-species count transformation model
  ## with constant Cholesky factor of the precision matrix
  m_all_1 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                     formula = ~ 1, data = spiders)

  ## with covariate-dependent Cholesky factor of the precision matrix
  m_all_2 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                     formula = ~ Elevation + Canopy_openess, data = spiders)

  logLik(m_all_1)
  logLik(m_all_2)

  ## lambda defining the Cholesky factor of the precision matrix,
  ## with standard errors
  coef(m_all_1, newdata = spiders[1,], type = "Lambda")
  V <- vcov(m_all_1)[55:69, 55:69]
  (se <- sqrt(diag(V)))

  coef(m_all_2, newdata = spiders[1,], type = "Lambda")
  V <- vcov(m_all_2)[55:99, 55:99]
  (se <- sqrt(diag(V)))

  ## linear correlation, ie Pearson correlation of the models after
  ## transformation to bivariate normality
  (r1 <- coef(m_all_1, newdata = spiders[1,], type = "Corr"))
  (r2 <- coef(m_all_2, newdata = spiders[1,], type = "Corr"))

  options(op)

}
\keyword{models}
