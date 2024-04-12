\name{mmlt}
\alias{mmlt}
\alias{mmltoptim}
\alias{coef.cmmlt}
\alias{coef.mmmlt}
\alias{predict.mmlt}
\alias{simulate.mmlt}
\title{
  Multivariate Conditional Transformation Models
}
\description{
  Conditional transformation models for multivariate continuous, discrete,
  or a mix of continuous and discrete outcomes
}
\usage{
mmlt(..., formula = ~ 1, data, conditional = FALSE, theta = NULL, fixed = NULL,
     scale = FALSE, optim = mmltoptim(), args = list(seed = 1, M = 1000), 
     dofit = TRUE, domargins = TRUE, sequentialfit = FALSE)
mmltoptim(auglag = list(maxtry = 5), ...)
\method{coef}{cmmlt}(object, newdata, 
     type = c("all", "conditional", "Lambdapar", "Lambda", "Lambdainv", 
              "Precision", "PartialCorr", "Sigma", "Corr", 
              "Spearman", "Kendall"), fixed = FALSE, 
     ...)
\method{coef}{mmmlt}(object, newdata, 
     type = c("all", "marginal", "Lambdapar", "Lambda", "Lambdainv", 
              "Precision", "PartialCorr", "Sigma", "Corr", 
              "Spearman", "Kendall"), fixed = FALSE,
     ...)
\method{predict}{mmlt}(object, newdata, margins = 1:J, 
        type = c("trafo", "distribution", "survivor", "density", "hazard"), 
                 log = FALSE, args = object$args, ...)
\method{simulate}{mmlt}(object, nsim = 1L, seed = NULL, newdata, K = 50, ...)
}
\arguments{
  \item{\dots}{marginal transformation models, one for each response, for
\code{mmlt}. Additional arguments for the methods.
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
  \item{fixed}{an optional named numeric vector of predefined parameter values
or a logical (for \code{coef}) indicating to also return fixed parameters
(only when \code{type = "all"}).
}
  \item{scale}{a logical indicating if (internal) scaling shall be applied
          to the model coefficients.
}
  \item{optim}{a list of optimisers as returned by \code{\link{mmltoptim}}
}
  \item{args}{a list of arguments for \code{\link[mvtnorm]{lpmvnorm}}.}
  \item{dofit}{logical; parameters are fitted by default, otherwise a list
               with log-likelihood and score function is returned.
}
  \item{domargins}{logical; all model parameters are fitted by default, 
                   including the parameters of marginal models.
}
  \item{sequentialfit}{logical; fit parameters sequentially by adding one
       model at a time, keeping all parameters concerning previously added models
       fixed. This does not lead to an optimal estimate but might provide
       attractive starting values.
}
\item{auglag}{a list of arguments to \code{\link[alabama]{auglag}} ensuring
              a numerically determined Hessian is returned.}
  \item{object}{an object of class \code{mmlt}.
}
  \item{newdata}{an optional data.frame coefficients and predictions shall
                 be computed for.
}
  \item{type}{type of coefficient or prediction to be returned.
}
  \item{margins}{indices defining marginal models to be evaluated. Can be
single integers giving the marginal distribution of the corresponding
variable, or multiple integers (currently only \code{1:j} implemented).
}
  \item{log}{logical; return log-probabilities or log-densities of
\code{TRUE}.
}
  \item{nsim}{number of samples to generate.}
  \item{seed}{optional seed for the random number generator.}
  \item{K}{number of grid points to generate.}
}
\details{
  The function implements multivariate conditional transformation models
  as described by Klein et al (2020). 
  Below is a simple example for an unconditional bivariate distribution.
  See \code{demo("undernutrition", package = "tram")} for a conditional
  three-variate example.
}
\value{
An object of class \code{mmlt} with \code{coef} and \code{predict}
methods.
}
\references{

  Nadja Klein, Torsten Hothorn, Luisa Barbanti, Thomas Kneib (2022),
  Multivariate Conditional Transformation Models. \emph{Scandinavian Journal
  of Statistics}, \bold{49}, 116--142, \doi{10.1111/sjos.12501}.

}
\examples{

  data("cars")

  ### fit unconditional bivariate distribution of speed and distance to stop
  ## fit unconditional marginal transformation models
  m_speed <- BoxCox(speed ~ 1, data = cars, support = ss <- c(4, 25), 
                    add = c(-5, 5))
  m_dist <- BoxCox(dist ~ 1, data = cars, support = sd <- c(0, 120), 
                   add = c(-5, 5))

  ## fit multivariate unconditional transformation model
  m_speed_dist <- mmlt(m_speed, m_dist, formula = ~ 1, data = cars)

  ## log-likelihood
  logLik(m_speed_dist)
  sum(predict(m_speed_dist, newdata = cars, type = "density", log = TRUE))

  ## Wald test of independence of speed and dist (the "dist.speed.(Intercept)"
  ## coefficient)
  summary(m_speed_dist)

  ## LR test comparing to independence model
  LR <- 2 * (logLik(m_speed_dist) - logLik(m_speed) - logLik(m_dist))
  pchisq(LR, df = 1, lower.tail = FALSE)

  ## constrain lambda to zero and fit independence model
  ## => log-likelihood is the sum of the marginal log-likelihoods
  mI <- mmlt(m_speed, m_dist, formula = ~1, data = cars, 
             fixed = c("dist.speed.(Intercept)" = 0))
  logLik(m_speed) + logLik(m_dist)
  logLik(mI)

  ## linear correlation, ie Pearson correlation of speed and dist after
  ## transformation to bivariate normality
  (r <- coef(m_speed_dist, type = "Corr"))
  
  ## Spearman's rho (rank correlation) of speed and dist on original scale
  (rs <- coef(m_speed_dist, type = "Spearman"))

  ## evaluate joint and marginal densities (needs to be more user-friendly)
  nd <- expand.grid(c(nd_s <- mkgrid(m_speed, 100), nd_d <- mkgrid(m_dist, 100)))
  nd$d <- predict(m_speed_dist, newdata = nd, type = "density")

  ## compute marginal densities
  nd_s <- as.data.frame(nd_s)
  nd_s$d <- predict(m_speed_dist, newdata = nd_s, margins = 1L,
                    type = "density")
  nd_d <- as.data.frame(nd_d)
  nd_d$d <- predict(m_speed_dist, newdata = nd_d, margins = 2L, 
                    type = "density")

  ## plot bivariate and marginal distribution
  col1 <- rgb(.1, .1, .1, .9)
  col2 <- rgb(.1, .1, .1, .5)
  w <- c(.8, .2)
  layout(matrix(c(2, 1, 4, 3), nrow = 2), width = w, height = rev(w))
  par(mai = c(1, 1, 0, 0) * par("mai"))
  sp <- unique(nd$speed)
  di <- unique(nd$dist)
  d <- matrix(nd$d, nrow = length(sp))
  contour(sp, di, d, xlab = "Speed (in mph)", ylab = "Distance (in ft)", xlim = ss, ylim = sd,
          col = col1)
  points(cars$speed, cars$dist, pch = 19, col = col2)
  mai <- par("mai")
  par(mai = c(0, 1, 0, 1) * mai)
  plot(d ~ speed, data = nd_s, xlim = ss, type = "n", axes = FALSE, 
       xlab = "", ylab = "")
  polygon(nd_s$speed, nd_s$d, col = col2, border = FALSE)
  par(mai = c(1, 0, 1, 0) * mai)
  plot(dist ~ d, data = nd_d, ylim = sd, type = "n", axes = FALSE, 
       xlab = "", ylab = "")
  polygon(nd_d$d, nd_d$dist, col = col2, border = FALSE)

  ### NOTE: marginal densities are NOT normal, nor is the joint
  ### distribution. The non-normal shape comes from the data-driven 
  ### transformation of both variables to joint normality in this model.

}
\keyword{models}
