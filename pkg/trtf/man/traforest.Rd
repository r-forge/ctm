\name{traforest}
\alias{traforest}
\alias{predict.traforest}
\alias{logLik.traforest}
\title{Transformation Forests}
\description{
  Partitioned and aggregated transformation models
}
\usage{
traforest(object, parm = 1:length(coef(object)), reparm = NULL,
          intercept = c("none", "shift", "scale", "shift-scale"),
          update = TRUE, min_update = length(coef(object)) * 2,
          mltargs = list(),  ...)
\method{predict}{traforest}(object,  newdata, mnewdata = data.frame(1), K = 20, q = NULL,
    type = c("weights", "node", "coef", "trafo", "distribution", "survivor", "density",
             "logdensity", "hazard", "loghazard", "cumhazard", "quantile"),
    OOB = FALSE, simplify = FALSE, trace = FALSE, updatestart = FALSE, 
    applyfun = NULL, cores = NULL, ...)
\method{logLik}{traforest}(object, newdata, weights = NULL, OOB = FALSE, coef = NULL,  ...)
}
\arguments{
  \item{object}{an object of class \code{\link[mlt]{ctm}} or \code{\link[mlt]{mlt}} specifying the
                abstract model to be partitioned.}
  \item{parm}{parameters of \code{object} those corresponding score is
              used for finding partitions.}
  \item{reparm}{optional matrix of contrasts for reparameterisation of the scores.
                \code{teststat = "quadratic"} is invariant to this operation
                but \code{teststat = "max"} might be more powerful for
                example when formulating an implicit into an explicit intercept term.}
  \item{intercept}{add optional intercept parameters (constraint to zero) to
                   the model.}
  \item{mltargs}{arguments to \code{\link[mlt]{mlt}} for fitting the
                 transformation models.}
  \item{update}{logical, if \code{TRUE}, models and thus scores are updated in
                every node. If \code{FALSE}, the model and scores are
                computed once in the root node. The latter option is faster
                but less accurate.}
  \item{min_update}{number of observations necessary to refit the model in a
                    node. If less observations are available, the parameters from the parent
                    node will be reused.}
  \item{newdata}{an optional data frame of observations for the forest.}
  \item{mnewdata}{an optional data frame of observations for the model.}
  \item{K}{number of grid points to generate (in the absence of \code{q}).}
  \item{q}{quantiles at which to evaluate the model.}
  \item{type}{type of prediction or plot to generate.}
  \item{OOB}{compute out-of-bag predictions.}
  \item{simplify}{simplify predictions (if possible).}
  \item{trace}{a logical indicating if a progress bar shall be printed while
          the predictions are computed.}
  \item{updatestart}{try to be smart about starting values for computing
                     predictions (experimental).}
  \item{applyfun}{an optional \code{lapply}-style function with arguments
          \code{function(X, FUN, ...)} for looping over \code{newdata}.
          The default is to use the
          basic \code{lapply} function unless the \code{cores} argument is
          specified (see below).}
  \item{cores}{numeric. If set to an integer the \code{applyfun} is set to
          \code{mclapply} with the desired number of \code{cores}.}
  \item{weights}{an optional vector of weights.}
  \item{coef}{an optional matrix of precomputed coefficients for
              \code{newdata} (using \code{predict}). Helps to compute the
              coefficients once for later reuse (different weights, for
              example).}
  \item{\dots}{arguments to \code{\link[partykit]{cforest}}, at least
               \code{formula} and \code{data}.}
}
\details{
  Conditional inference trees are used for partitioning likelihood-based transformation
  models as described in Hothorn and Zeileis (2017). The method can be seen
  in action in Hothorn (2018) and the corresponding code is available as
  \code{demo("BMI")}. 

}
\value{
  An object of class \code{traforest} with corresponding \code{logLik} and
  \code{predict} methods.
}
\references{
  Torsten Hothorn and Achim Zeileis (2021). Predictive Distribution 
  Modelling Using Transformation Forests. 
  \emph{Journal of Computational and Graphical Statistics},
  \doi{10.1080/10618600.2021.1872581}.

  Torsten Hothorn (2018). Top-Down Transformation Choice. \emph{Statistical
  Modelling}, \bold{3-4}, 274-298. \doi{10.1177/1471082X17748081}.

  Natalia Korepanova, Heidi Seibold, Verena Steffen and Torsten Hothorn
  (2019).  Survival Forests under Test: Impact of the Proportional Hazards
  Assumption on Prognostic and Predictive Forests for ALS Survival.
  \doi{10.1177/0962280219862586}.


}
\examples{

### Example: Personalised Medicine Using Partitioned and Aggregated Cox-Models
### A combination of <DOI:10.1177/0962280217693034> and <arXiv:1701.02110>
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

\dontrun{

### estimate age-dependent Cox models (here ignoring all other covariates)
ctrl <- ctree_control(minsplit = 50, minbucket = 20, mincriterion = 0)
set.seed(290875)
tf_cmod <- traforest(cmod, formula = y ~ horTh | age, control = ctrl, 
                     ntree = 50, mtry = 1, trace = TRUE, data = GBSG2)

### plot age-dependent treatment effects vs. overall treatment effect
nd <- data.frame(age = 30:70)
cf <- predict(tf_cmod, newdata = nd, type = "coef")
nd$logHR <- sapply(cf, function(x) x["horThyes"])
plot(logHR ~ age, data = nd, pch = 19, xlab = "Age", ylab = "log-Hazard Ratio")
abline(h = coef(cmod <- mlt(m, data = GBSG2))["horThyes"])
### treatment most beneficial in very young patients
### NOTE: scale of log-hazard ratios depends on
### corresponding baseline hazard function which  _differs_
### across age; interpretation of positive / negative treatment effect is,
### however, save.

### mclapply doesn't work in Windows
if (.Platform$OS.type != "windows") {

  ### computing predictions: predicted coefficients
  cf1 <- predict(tf_cmod, newdata = nd, type = "coef")
  ### speedup with plenty of RAM and 4 cores
  cf2 <- predict(tf_cmod, newdata = nd, cores = 4, type = "coef")
  ### memory-efficient with low RAM and _one_ core
  cf3 <- predict(tf_cmod, newdata = nd, cores = 4, applyfun = lapply, type = "coef")
  all.equal(cf1, cf2)
  all.equal(cf1, cf3)

}

}

}
\keyword{trees}
