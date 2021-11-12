\name{Bernstein_basis}
\alias{Bernstein_basis}
\title{
    Bernstein Basis Functions
}
\description{
    Basis functions defining a polynomial in Bernstein form
}
\usage{
Bernstein_basis(var, order = 2, ui = c("none", "increasing", "decreasing", 
                                       "cyclic", "zerointegral", "positive",
                                       "negative", "concave", "convex"),
                extrapolate = FALSE, log_first = FALSE)
}
\arguments{
  \item{var}{a \code{\link[variables]{numeric_var}} object}
  \item{order}{the order of the polynomial, one defines a linear function}
  \item{ui}{a character describing possible constraints}
  \item{extrapolate}{logical; if \code{TRUE}, the polynomial is
                     extrapolated linearily outside \code{support(var)}.
                     In particular, the second derivative of the polynomial
                     at \code{support(var)} is constrained to zero.}
  \item{log_first}{logical; the polynomial in Bernstein form is defined on the
                   log-scale if \code{TRUE}. It makes sense to define the
                   support as \code{c(1, q)}$, ie putting the first 
                   basis function of the polynomial on
                   \code{log(1)}.}
}
\details{
  \code{Bernstein_basis} returns a function for the evaluation of
  the basis functions with corresponding \code{model.matrix} and \code{predict}
  methods.
}
\references{

  Rida T. Farouki (2012), The Bernstein Polynomial Basis: A Centennial Retrospective,
  \emph{Computer Aided Geometric Design}, \bold{29}(6), 379--419.
  http://dx.doi.org/10.1016/j.cagd.2012.03.001

}
\examples{

  ### set-up basis
  bb <- Bernstein_basis(numeric_var("x", support = c(0, pi)), 
                        order = 3, ui = "increasing")

  ### generate data + coefficients
  x <- as.data.frame(mkgrid(bb, n = 100))
  cf <- c(1, 2, 2.5, 2.6)

  ### evaluate basis (in two equivalent ways)
  bb(x[1:10,,drop = FALSE])
  model.matrix(bb, data = x[1:10, ,drop = FALSE])

  ### check constraints
  cnstr <- attr(bb(x[1:10,,drop = FALSE]), "constraint")
  all(cnstr$ui \%*\% cf > cnstr$ci)

  ### evaluate and plot Bernstein polynomial defined by
  ### basis and coefficients
  plot(x$x, predict(bb, newdata = x, coef = cf), type = "l")

  ### evaluate and plot first derivative of 
  ### Bernstein polynomial defined by basis and coefficients
  plot(x$x, predict(bb, newdata = x, coef = cf, deriv = c(x = 1)), 
       type = "l")

  ### illustrate constrainted estimation by toy example
  N <- 100
  order <- 10
  x <- seq(from = 0, to = pi, length.out = N)
  y <- rnorm(N, mean = -sin(x) + .5, sd = .5)

  if (require("coneproj")) {
    prnt_est <- function(ui) {
      xv <- numeric_var("x", support = c(0, pi))
      xb <- Bernstein_basis(xv, order = 10, ui = ui)
      X <- model.matrix(xb, data = data.frame(x = x))
      uiM <- as(attr(X, "constraint")$ui, "matrix")
      ci <- attr(X, "constraint")$ci
      if (all(is.finite(ci)))
        parm <- qprog(crossprod(X), crossprod(X, y), 
                      uiM, ci, msg = FALSE)$thetahat
      else
        parm <- coef(lm(y ~ 0 + X))
      plot(x, y, main = ui)
      lines(x, X \%*\% parm, col = col[ui], lwd = 2)
    }
    ui <- eval(formals(Bernstein_basis)$ui)
    col <- 1:length(ui)
    names(col) <- ui
    layout(matrix(1:length(ui), 
                  ncol = ceiling(sqrt(length(ui)))))
    tmp <- sapply(ui, function(x) try(prnt_est(x)))
  }
}
