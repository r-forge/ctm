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

  Luisa Barbani, Roland Brandl, Torsten Hothorn (2022), Multi-species Count
  Transformation Models, \doi{10.48550/arXiv.2201.13095}.

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
  ##
  ## define starting values here (this is not necessary but leads
  ## to diffs for ATLAS and OpenBlas)
  cf1 <- c(0.415, 0.921, 0.921, 0.921, 1.257, 1.42, 1.654, 0, 0.001, 0.853,
           1.365, 1.365, 1.365, 1.791, 2.094, 2.428, 0.002, 0, -1.219,
           -1.218, -0.514, -0.403, -0.403, 0.078, 0.863, -0.001, 0, 0.326,
           0.596, 0.91, 0.91, 0.91, 1.481, 2.137, -0.001, 0.001, -1.628,
           -1.188, -1.188, -1.129, -0.434, -0.296, 0.395, 0, -0.001, -6.345,
           -6.015, -5.291, -5.291, -5.291, -4.541, -3.554, -0.006, 0, 0.69,
           -0.345, 0.28, -0.944, 0.787, 0.028, 0.322, -0.559, 0.263, 0.876,
           0.006, -0.207, -0.481, -0.083, 0.161)

  m_all_1 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                     theta = cf1, ## <- not really necessary, please CRAN timings
                     formula = ~ 1, data = spiders)

  ## with covariate-dependent Cholesky factor of the precision matrix
  cf2 <- c(1.649, 2.078, 2.078, 2.078, 2.391, 2.476, 2.694, 0.001, 0, 3.399,
           3.862, 3.862, 3.862, 4.232, 4.593, 4.991, 0.005, -0.001, -3.646, -3.646,
           -2.906, -2.777, -2.777, -2.256, -1.429, -0.003, 0.001, -1.257, -1.01,
           -0.652, -0.652, -0.652, -0.093, 0.601, -0.002, 0.001, 2.294, 2.77, 2.77,
            2.796, 3.593, 3.74, 4.473, 0.006, -0.003, -7.563, -7.173,
           -6.384, -6.384, -6.384, -5.571, -4.509, -0.007, 0.001, -4.788,
           0.007, -0.001, 3.914, -0.005, 0, 1.534, -0.001, 0, 1.098, -0.002,
           -0.001, 2.196, -0.002, 0, -0.072, 0, 0, -0.222, 0.001, -0.001,
           -2.876, 0.003, -0.001, -1.135, 0.002, -0.001, 0.07, 0.001,
           -0.001, 3.279, -0.004, 0, -1.057, 0.001, 0, -0.953, 0.001, 0,
           0.252, 0, 0, 2.443, -0.002, -0.001)

  m_all_2 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                     theta = cf2, ## <- please CRAN timings
                     formula = ~ Elevation + Canopy_openess, data = spiders)

  logLik(m_all_1)
  logLik(m_all_2)

  ## lambda defining the Cholesky factor of the precision matrix
  coef(m_all_1, newdata = spiders[1,], type = "Lambda")
  coef(m_all_2, newdata = spiders[1,], type = "Lambda")

  ## linear correlation, ie Pearson correlation of the models after
  ## transformation to bivariate normality
  (r1 <- coef(m_all_1, newdata = spiders[1,], type = "Corr"))
  (r2 <- coef(m_all_2, newdata = spiders[1,], type = "Corr"))

  options(op)

}
\keyword{models}
