## check links

library("cotram")
library("survival")

set.seed(29)

### check coefficients (tram vs cotram model)
.check_cf <- function(m, mc) {

  cfc <- coef(as.mlt(mc))
  
  ## tram coefs for negative = TRUE
  cf <- coef(c <- as.mlt(m))
  cf[c$shiftcoef] <- c(-1, 1)[(m$negative & mc$negative) + 1] * cf[c$shiftcoef]
  
  ## check tram vs cotram coefs
  stopifnot(all.equal(cf, cfc, check.attributes = FALSE))
  
  # print(paste("coefficients all equal for tram & cotram, for log_first = ", mc$log_first, " and ", 
  #             mc$model$todistr$name, " inv. link.", sep = "'"))
}


### check log-likelihood (tram vs cotram model)
.check_ll <- function(m, mc) {

  ## simple check wrt to newdata
  stopifnot(logLik(mc) == logLik(mc, newdata = model.frame(mc)))
  
  ## likelihood contributions interval-censored
  L <- predict(m, newdata = data.frame(yi = d$y + as.integer(log_first), x = d$x), type = "distribution") -
    predict(m, newdata = data.frame(yi = d$y + as.integer(log_first) - 1L, x = d$x), type = "distribution")
  
  ## check interval-censored vs cotram log-likelihoods
  stopifnot(all.equal(log(L), mc$logliki(coef(as.mlt(mc)), mc$weights), check.attributes = FALSE))
  
  ## check tram vs cotram log-likelihoods
  stopifnot(all.equal(m$logliki(coef(as.mlt(m)), m$weight), mc$logliki(coef(as.mlt(mc)), mc$weight)))
  
  # print(paste("log-likelihoods all equal for tram & cotram, for log_first = ", mc$log_first, " and ", 
  #             mc$model$todistr$name, " inv. link.", sep = "'"))
}


### run checks for log_first = FALSE / TRUE & all links
## dgp counts
n <- 200
x <- runif(n)
y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))

## trams & cotrams
links <- c("logit", "cloglog", "loglog", "probit")
trams <- c("logit" = "Colr", "cloglog" = "Coxph", "loglog" = "Lehmann", "probit" = "BoxCox")

## run
for (log_first in c(FALSE, TRUE)) {
  
  # print(log_first)
  
  ## plus_one for log_first = TRUE
  plus_one <- as.integer(log_first)
  
  ## counts as interval censored object
  yleft <- y - 1L
  yleft[yleft < 0] <- -Inf
  yi <- Surv(yleft + plus_one, y + plus_one, type = "interval2")
  
  d <- data.frame(y = y, yi = yi, x = x)
  
  for (link in links) {
    # print(link)
    mc <- cotram(as.formula(y ~ x), data = d, method = link, log_first = log_first)
    # print(mc)
    
    ## check model.frame
    stopifnot(all.equal(nd <- model.frame(mc), d[names(nd)], check.attributes = FALSE))
    
    tram <- unname(trams[link])
    m <- do.call(tram, list(formula = yi ~ x, support = mc$support, bounds = mc$bounds, log_first = log_first))
    # print(m)
    
    ## check if same model
    stopifnot(mc$model$todistr$name == m$model$todistr$name)
    
    .check_cf(m, mc)
    .check_ll(m, mc)
  }
}


## additional checks for plus_one
d0 <- data.frame(y = q <- 0L:100L, x = runif(length(q)))

m0 <- cotram(y ~ x, log_first = TRUE, data = d0)

stopifnot(all(mkgrid(m0)$y == q)) ## mkgrid  for log_first

m0 <- cotram(y ~ x, log_first = FALSE, data = d0)

stopifnot(all(mkgrid(m0)$y == q)) ## mkgrid 
