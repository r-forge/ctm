\name{Compris}
\alias{Compris}
\title{
  Competing Risk Regression
}
\description{
  An alternative approach to competing risk regression via multivariate
  transformation models
}
\usage{
Compris(formula, data, subset, weights, na.action, offset, 
        primary = c("Coxph", "Colr", "BoxCox"), 
        competing = switch(primary, Coxph = "weibull", Colr = "loglogistic", BoxCox = "lognormal"), 
        optim = mmltoptim(), M = 1000, scale = FALSE, tol = 0.001, ...)
}
\arguments{
  \item{formula}{an object of class \code{"formula"}: a symbolic description 
          of the model structure to be
          fitted.  The details of model specification are given under
          Details and in the package vignette.
  }
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}.
  }
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.
  }
  \item{weights}{an optional vector of case weights to be used in the fitting
          process.  Should be \code{NULL} or a numeric vector. If present,
          the weighted log-likelihood is maximised.
  }
  \item{na.action}{a function which indicates what should happen when the data
          contain \code{NA}s.  The default is set to \code{na.omit}.
  }
  \item{offset}{this can be used to specify an _a priori_ known component to
          be included in the linear predictor during fitting.  This
          should be \code{NULL} or a numeric vector of length equal to the
          number of cases.
  }
  \item{primary}{a character defining the marginal model for the primary
  event of interest, that is, the first status level.
}
  \item{competing}{a character defining the marginal models for the
remaining competing events.
}
  \item{optim}{see \code{mmlt}.
}
 \item{args}{a list of arguments for \code{\link[mvtnorm]{lpmvnorm}}.
}
 \item{scale}{logical defining if variables in the linear predictor shall
               be scaled. Scaling is internally used for model estimation,
               rescaled coefficients are reported in model output.
}
  \item{tol}{a tolerance for faking interval censoring.
}
  \item{\dots}{addition arguments.
}
}
\details{
This is a highly experimental approach to an alternative competing risk
regression framework described by Czado and Van Keilegom (2023) and Deresa
and Van Keilegom (2023).
}
\value{
An object of class \code{mmlt}, allowing to derive marginal time-to-event
distributions for the primary event of interest and all competing events.
}
\references{

  Claudia Czado and Ingrid Van Keilegom (2023). Dependent Censoring Based on Parametric Copulas.
  \emph{Biometrika}, \doi{10.1093/biomet/asac067}.

  Negera Wakgari Deresa and Ingrid Van Keilegom (2023).  Copula Based Cox Proportional Hazards 
  Models for Dependent Censoring. \emph{Journal of the American Statistical
  Association}, \doi{10.1080/01621459.2022.2161387}.

}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{survival}
