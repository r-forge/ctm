\name{cyclic_basis}
\alias{cyclic_basis}
\title{
    Cyclic Basis Function
}
\description{
    The cyclic basis function for modelling seasonal effects
}
\usage{
cyclic_basis(var, order = 3, frequency)
}
\arguments{
  \item{var}{a \code{\link[variables]{numeric_var}} object}
  \item{order}{the order of the basis function}
  \item{frequency}{frequency}
}
\details{
  \code{cyclic_basis} returns a set of sin and cosine functions for
  modelling seasonal effects, see Held and Paul (2012), Section 2.2.

  For any choice of coefficients, the function returns the same value
  when evaluated at multiples of \code{frequency}.
}
\references{

  Leonhard Held and Michaela Paul (2012), Modeling Seasonality in Space-time
  Infectious Disease Surveillance Data, \emph{Biometrical Journal},
  \bold{54}(6), 824--843, \doi{10.1002/bimj.201200037}

}
\examples{

  ### set-up basis
  cb <- cyclic_basis(numeric_var("x"), order = 3, frequency = 2 * pi)

  ### generate data + coefficients
  x <- data.frame(x = c(0, pi, 2 * pi))

  ### f(0) = f(2 * pi)
  cb(x) %*% 1:6 / 6
}
