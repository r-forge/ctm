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
mcotram(..., formula = ~ 1, data, conditional = FALSE, theta = NULL,
        fixed = NULL, scale = FALSE, optim = mmltoptim(), 
        M = 1000, dofit = TRUE, domargins = TRUE)
}
\arguments{
  \item{\dots}{marginal count transformation models, one for each response
}
  \item{formula}{a model formula describing a model for the dependency
structure via the lambda parameters. The default is set to \code{~ 1} for constant lambdas.
}
  \item{data}{a data.frame.
}
  \item{conditional}{logical; parameters are defined conditionally (only
possible when all models are probit models). This is the default as
described by Klein et al. (2022). If \code{FALSE}, parameters can be
directly interpreted marginally, this is explained in Section 2.6 by Klein
et al. (2022). Using \code{conditional = FALSE} with probit-only models
gives the same likelihood but different parameter estimates.
}
  \item{theta}{an optional vector of starting values.
}
 \item{fixed}{an optional named numeric vector of predefined parameter values.
}
  \item{scale}{a logical indicating if (internal) scaling shall be applied
          to the model coefficients.
}
  \item{optim}{a list of optimisers as returned by \code{\link[tram]{mmltoptim}}
}
  \item{M}{number of Halton sequences used to approximate the log-likelihood in \code{\link[mvtnorm]{lpmvnorm}}.}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned.
}
  \item{domargins}{logical; all model parameters are fitted by default, 
                   including the parameters of marginal models.
}
}
\details{
  The function implements multivariate count conditional transformation models.
  The response is assumed to be a vector of counts.
}
\value{
An object of class \code{mmlt} with \code{coef} and \code{predict} methods.
}
\references{

  Luisa Barbani, Roland Brandl, Torsten Hothorn (2022), Multi-species Count
  Transformation Models, \doi{10.48550/arXiv.2201.13095}.

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2020),
  Multivariate Conditional Transformation Models. \emph{Scandinavian Journal
  of Statistics}, \doi{10.1111/sjos.12501}.

}
\examples{

library("cotram")
data("spiders", package = "cotram")

### for illustration only
OR <- 1      ### order of transformation function
             ### OR = 1 means log-linear, use OR ~ 6
M <- 100     ### number of Halton sequences, seem sufficient here

## fit conditional marginal count transformation models
## one for each species
## (don't test as it takes too long for CRAN)
\donttest{m_PF <- cotram(Pardosa_ferruginea ~ Elevation + Canopy_openess, 
               data = spiders, method = "probit", order = OR)
m_HL <- cotram(Harpactea_lepida ~ Elevation + Canopy_openess,
               data = spiders, method = "probit", order = OR)
m_CC <- cotram(Callobius_claustrarius ~ Elevation + Canopy_openess,
               data = spiders, method = "probit", order = OR)
m_CT <- cotram(Coelotes_terrestris ~ Elevation + Canopy_openess,
               data = spiders, method = "probit", order = OR)
m_PL <- cotram(Pardosa_lugubris ~ Elevation + Canopy_openess,
               data = spiders, method = "probit", order = OR)
m_PR <- cotram(Pardosa_riparia ~ Elevation + Canopy_openess,
               data = spiders, method = "probit", order = OR)

### fit dependence parameters
mm <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, data = spiders,
              M = M, scale = TRUE)
logLik(mm)

### Kendall's tau: Dependence of species after accounting
### for elevation and canopy openess in marginal models
coef(mm, type = "Kendall")

### regress dependencies on elevation and canopy openess
mmc <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, data = spiders, 
               formula = ~ Elevation + Canopy_openess, M = M, scale = TRUE)
logLik(mmc)

### weak evidence for such effects
pchisq(2 * (logLik(mmc) - logLik(mm)), df = 30, lower.tail = FALSE)

### plot Kendall's tau for different elevations / openess levels
nd <- expand.grid(Elevation = 80:120 * 10, Canopy_openess = 1:10 * 10)
KD <- Lower_tri(coef(mmc, newdata = nd, type = "Kendall"))
f <- factor(rownames(KD))
nd <- cbind(f = rep(f, nrow(nd)), nd[rep(1:nrow(nd), each = nlevels(f)),])
nd$KD <- c(KD)

if (require("lattice")) {
contourplot(KD ~ Elevation + Canopy_openess | f, data = nd, 
            cuts = 18, xlab = "Elevation", ylab = "Canopy openess")
}
}

### for example:
### => constant negative dependence of Pardosa_lugubris and Coelotes_terrestris
### => weak dependence of Harpactea_lepida and Pardosa_ferruginea
###    for low elevations, negative dependence increasing with elevation

}
\keyword{models}
