\name{robust_score_test}
\alias{robust_score_test}
\alias{robust_score_test.tram}
\title{Doubly Robust Transformation Score Test}
\usage{
robust_score_test(object, ...)

\method{robust_score_test}{tram}(
  object,
  parm = names(coef(object)),
  alternative = c("two.sided", "less", "greater"),
  nullvalue = 0,
  confint = FALSE,
  level = 0.95,
  ranger_args = NULL,
  ...
)
}
\arguments{
\item{object}{an object of class '\code{tram}'}

\item{...}{additional arguments, currently ignored.}

\item{parm}{a vector of names of parameters to be tested. These
parameters must be present in \code{object}}

\item{alternative}{a character string specifying the alternative
hypothesis, must be one of \code{"two.sided"} (default),
\code{"greater"}, or \code{"less"}}

\item{nullvalue}{a number specifying an optional parameter used
to form the null hypothesis H_0: parm = nullvalue and defaults
to zero}

\item{confint}{a logical indicating whether to (numerically) invert
the test to obtain a robust score confidence interval}

\item{level}{the confidence level}

\item{ranger_args}{arguments passed to \code{\link[ranger]{ranger}}
for the regression of the column in the design matrix corresponding
to \code{parm} against all others}
}
\value{
An object of class '\code{htest}' or a list thereof.
}
\description{
Doubly robust p-values and confidence intervals for parameters in
(stratified) linear (shift-scale) transformation models obtained
using the tram generalised covariance measure test.
}
\details{
For a (stratified) linear shift (-scale) transformation he tram-GCM test
tests the hypothesis H0: \code{parm} = \code{nullvalue} by re-fitting the
model under the null hypothesis, computing the score residuals (see
\code{\link[tram]{residuals.tram}}), and running an additional regression
of the column in the design matrix corresponding to \code{parm} on the
remaining columns, computing the corresponding residuals, and finally
computing correlation-type test between the score and predictor residuals.
}
\examples{
data("mtcars")
### Linear shift tram
m <- Lm(mpg ~ cyl + disp, data = mtcars)
robust_score_test(m, parm = "cyl")
### Linear shift-scale tram
m2 <- Lm(mpg ~ cyl | disp, data = mtcars)
robust_score_test(m2, parm = "cyl")
robust_score_test(m2, parm = "scl_disp")
### Stratified linear shift tram
m3 <- Lm(mpg | 0 + disp ~ cyl, data = mtcars)
robust_score_test(m3, parm = "cyl")
### Stratified linear shift-scale tram
m4 <- Lm(mpg | 0 + disp ~ cyl | cyl, data = mtcars)
robust_score_test(m4, parm = "cyl")

}
\references{
Kook, L., Saengkyongam, S., Lundborg, A. R., Hothorn, T., &
Peters, J. (2024). Model-based causal feature selection for general response
types. Journal of the American Statistical Association, 1-12.
\doi{10.1080/01621459.2024.2395588}
}
\author{
Lucas Kook, Torsten Hothorn
}
