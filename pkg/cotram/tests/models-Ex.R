## check links

library("cotram")
library("survival")

set.seed(29)

### check coefficients (tram vs cotram model)
.check_cf <- function(m, mc, ...) {
  
  cfc <- coef(as.mlt(mc))
  
  ## tram coefs for negative = TRUE
  cf <- coef(c <- as.mlt(m))
  cf[c$shiftcoef] <- c(-1, 1)[(m$negative & mc$negative) + 1] * cf[c$shiftcoef]
  
  ## check tram vs cotram coefs
  print(all.equal(cf, cfc, check.attributes = FALSE, ...))
}


### check log-likelihood (tram vs cotram model)
.check_ll <- function(m, mc, ...) {
  
  ## simple check wrt to newdata
  stopifnot(logLik(mc) == logLik(mc, newdata = model.frame(mc)))
  
  ## likelihood contributions interval-censored
  ### newdata
  nd_m <- d[c("y", "x")]
  colnames(nd_m) <- variable.names(m)
  nd_m1 <- nd_m2 <- nd_m
  nd_m2[, 1] <- nd_m2[, 1] + as.integer(log_first)
  nd_m1[, 1] <- nd_m1[, 1] - 1L + as.integer(log_first)
  
  pm <- predict(m, newdata = nd_m2, type = "distribution")
  pc <- predict(mc, newdata = d, type = "distribution")
  print(all.equal(pm, pc, ...))

  # check tram vs cotram log-likelihood contributions
  print(all.equal(m$logliki(coef(as.mlt(m)), m$weight), mc$logliki(coef(as.mlt(mc)), mc$weight), ...))
}



### run checks for log_first = FALSE / TRUE & all links
## dgp counts
n <- 2000
x <- runif(n, max = 3)
y <- as.integer(rnbinom(n, mu = exp(.5 + .8 * x), size = 10))

## trams & cotrams
links <- c("logit", "cloglog", "loglog", "probit")
trams <- c("logit" = "Colr", "cloglog" = "Coxph", "loglog" = "Lehmann", "probit" = "BoxCox")

## run
for (log_first in c(FALSE, TRUE)) {
  
  set.seed(2)

  LF <- NULL
  if (isTRUE(log_first)) LF <- TRUE  
  
  print(log_first)
  
  ## plus_one for log_first = TRUE
  plus_one <- as.integer(log_first)
  
  ## plus_one for integer y
  yii <- y + plus_one
  
  ## counts as interval censored object
  yleft <- y - 1L
  yleft[yleft < 0] <- -Inf
  yi <- Surv(yleft + plus_one, y + plus_one, type = "interval2")
  
  yi2 <- cotram:::.R.count(as.integer(y), plus_one = log_first)
  print(all.equal(yi, yi2))
  
  d <- data.frame(y = y, yi = yi, x = x)
  
  for (link in links) {
    print(link)
    
    ## fit "cotram"
    mc <- cotram(as.formula(y ~ x), data = d, method = link, log_first = LF)
    
    ## check model.frame
    stopifnot(all.equal(nd <- model.frame(mc), d[names(nd)], check.attributes = FALSE))
    
    ## "tram" with interval-censored response, R() ensures that Surv and
    ## integer use the same preprocessing steps in mlt::pretrafo
    tram <- unname(trams[link])
    m <- do.call(tram, list(formula = R(yi) ~ x, support = mc$support, bounds = mc$bounds, 
                            log_first = LF))
    
    ## check if same model
    stopifnot(mc$model$todistr$name == m$model$todistr$name)
    
    ## check coefs and logLiks
    .check_cf(m, mc)
    .check_ll(m, mc)
    
    ## "tram" with integer response
    m2 <- do.call(tram, list(formula = R(yii) ~ x, support = c(mc$support), bounds = mc$bounds,
      log_first = LF))
    
    ## check if same model
    stopifnot(mc$model$todistr$name == m2$model$todistr$name)
    
    ## check coefs and logLiks
    ### minor differences because integer are treated as interval (y - 1, 1] in "tram" 
    ### NOTE: Differences on Mac seem to be larger
    .check_cf(m2, mc, tol = 1e-3) 
    .check_ll(m2, mc, tol = 1e-3)
  }
}


## NOTE: 
## - Models with integer-y and interval-censored-y (for [y - 1, y]) are not
##   exactly equivalent for "tram": This is because the integer is treated as interval (y - 1, y] 
##   => This also explains the discrepancies between "cotram" and integer-y with "tram"
## - The discrepancy between "cotram" and interval-censored "tram" fit with "prob" given, 
##   comes from the different support (is rounded it "cotram")


## additional checks for plus_one
d0 <- data.frame(y = q <- 0L:100L, x = runif(length(q)))

m0 <- cotram(y ~ x, log_first = TRUE, data = d0)

stopifnot(all(mkgrid(m0)$y == q)) ## mkgrid  for log_first

m0 <- cotram(y ~ x, log_first = FALSE, data = d0)

stopifnot(all(mkgrid(m0)$y == q)) ## mkgrid 
