\name{trafotree}
\alias{trafotree}
\alias{predict.trafotree}
\alias{logLik.trafotree}
\title{Transformation Trees}
\description{
  Partitioned transformation models
}
\usage{
trafotree(object, parm = 1:length(coef(object)), reparm = NULL,
          intercept = c("none", "shift", "scale", "shift-scale"),
          min_update = length(coef(object)) * 2,
          mltargs = list(), ...)
\method{predict}{trafotree}(object, newdata, K = 20, q = NULL,
    type = c("node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    perm = NULL, ...)
\method{logLik}{trafotree}(object, newdata, weights = NULL, perm = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} or \code{\link[mlt]{mlt}} specifying the
                abstract model to be partitioned. For \code{predict} and
                \code{logLik}, \code{object} is an object of class \code{trafotree}.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{reparm}{optional matrix of contrasts for reparameterisation of the scores.
                \code{teststat = "quadratic"} is invariant to this operation
                but \code{teststat = "max"} might be more powerful for
                example when formulating an implicit into an explicit intercept term.}
  \item{intercept}{add optional intercept parameters (constraint to zero) to
                   the model. It may make sense to restrict attention to
                   scores corresponding to those intercept parameters, the additional argument
                   \code{parm = NULL} is needed in this case.}
  \item{min_update}{number of observations necessary to refit the model in a
                    node. If less observations are available, the parameters from the parent
                    node will be reused.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{newdata}{an optional data frame of observations.}
  \item{K}{number of grid points to generate (in the absence of \code{q}).}
  \item{q}{quantiles at which to evaluate the model.}
  \item{type}{type of prediction or plot to generate.}
  \item{weights}{an optional vector of weights.}
  \item{perm}{a vector of integers specifying the variables to be permuted
          prior before splitting (i.e., for computing permutation
          variable importances). The default \code{NULL} doesn't alter the
          data, see \code{\link[partykit]{fitted_node}}.}
  \item{\dots}{arguments to \code{\link[partykit]{ctree}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017). The method can be seen
  in action in Hothorn (2018) and the corresponding code is available as
  \code{demo("BMI")}. \code{demo("applications")} performs transformation
  tree analyses for some standard benchmarking problems.

}
\value{
  An object of class \code{trafotree} with corresponding \code{plot}, \code{logLik} and
  \code{predict} methods.}
\references{

  Torsten Hothorn and Achim Zeileis (2021). Predictive Distribution 
  Modelling Using Transformation Forests. 
  \emph{Journal of Computational and Graphical Statistics},
  \doi{10.1080/10618600.2021.1872581}.

  Torsten Hothorn (2018). Top-Down Transformation Choice. \emph{Statistical
  Modelling}, \bold{3-4}, 274-298. \doi{10.1177/1471082X17748081}

  Natalia Korepanova, Heidi Seibold, Verena Steffen and Torsten Hothorn
  (2019).  Survival Forests under Test: Impact of the Proportional Hazards
  Assumption on Prognostic and Predictive Forests for ALS Survival.
  \doi{10.1177/0962280219862586}.

}
\examples{

### Example: Stratified Medicine Using Partitioned Cox-Models
### A combination of <DOI:10.1515/ijb-2015-0032> and <arXiv:1701.02110>
### based on infrastructure in the mlt R add-on package described in
### https://cran.r-project.org/web/packages/mlt.docreg/vignettes/mlt.pdf

library("trtf")
library("survival")
### German Breast Cancer Study Group 2 data set
data("GBSG2", package = "TH.data")
GBSG2$y <- with(GBSG2, Surv(time, cens))

### set-up Cox model with overall treatment effect in hormonal therapy
cmod <- Coxph(y ~ horTh, data = GBSG2, support = c(100, 2000), order = 5)

### overall log-hazard ratio
coef(cmod)
### roughly the same as 
coef(coxph(y ~ horTh, data = GBSG2))

### partition the model, ie both the baseline hazard function AND the
### treatment effect
(part_cmod <- trafotree(cmod, formula = y ~ horTh | age + menostat + tsize + 
    tgrade + pnodes + progrec + estrec, data = GBSG2))

### compare the log-likelihoods
logLik(cmod)
logLik(part_cmod)

### stronger effects in nodes 2 and 4 and no effect in node 5
coef(part_cmod)[, "horThyes"]

### plot the conditional survivor functions; blue is untreated
### and green is hormonal therapy
nd <- data.frame(horTh = sort(unique(GBSG2$horTh)))
plot(part_cmod, newdata = nd, 
     tp_args = list(type = "survivor", col = c("cadetblue3", "chartreuse4")))

### same model, but with explicit intercept term and max-type statistic
### for _variable_ selection
(part_cmod_max <- trafotree(cmod, formula = y ~ horTh | age + menostat + tsize + 
    tgrade + pnodes + progrec + estrec, data = GBSG2, intercept = "shift",
    control = ctree_control(teststat = "max")))
logLik(part_cmod_max)
coef(part_cmod_max)[, "horThyes"]

### the trees (and log-likelihoods are the same) but the
### p-values are sometimes much smaller in the latter tree
cbind(format.pval(info_node(node_party(part_cmod))$criterion["p.value",]),
      format.pval(info_node(node_party(part_cmod_max))$criterion["p.value",]))

}
\keyword{trees}
