
library("cotram")
library("survival")

set.seed(25)

### some basic checks for predictions wrt discrete distributions
n <- 50
d <- data.frame(x1 = 0:n, x2 = sample(0:n) + 1, q = sample(0:n))
mc3 <- cotram(q ~ x1 + x2, data = d, log_first = log_first <- FALSE)

di <- d
di$qleft <- with(di, q - 1L)
di$qleft[di$qleft < 0] <- -Inf
di$q <- with(di, Surv(qleft + log_first, q + log_first, type = "interval2"))
m3 <- Colr(q ~ x1 + x2, data = di, support = mc3$support, bounds = mc3$bounds, log_first = log_first)

.chk <- function(x)
  stopifnot(max(abs(x), na.rm = TRUE) < sqrt(.Machine$double.eps))

nd <- d
nd$q <- NULL
q <-0:(n+1)

p <- predict(mc3, newdata = nd, q = q, type = "distribution")
s <- predict(mc3, newdata = nd, q = q, type = "survivor")

## check discrete distribution & survivor function
.chk(predict(mc3, newdata = nd, q = q, type = "distribution", log = TRUE) - log(p))
.chk(predict(mc3, newdata = nd, q = q, type = "distribution", lower.tail = FALSE) - s)
.chk(predict(mc3, newdata = nd, q = q, type = "distribution", 
             lower.tail = FALSE, log = TRUE) - log(s))

## check discrete odds function
o <- predict(mc3, newdata = nd, q = q, type = "odds")
.chk(o - p / s)

## check discrete density function
dd <- predict(mc3, newdata = nd, q = q, type = "density")
.chk(apply(p, 2, function(x) diff(c(0, x))) - dd)

## check discrete (cum)-hazard function
h <- predict(mc3, newdata = nd, q = q , type = "hazard")

.chk(dd / (1 - (p - dd)) - h)

.chk(apply(h, 2, function(x) cumprod(1 - x)) - s)

H <- predict(mc3, newdata = nd, q = q, type = "cumhazard")

.chk(H + log(s))


## simple checks for predict of smooth functions (cotram vs tram)
.check_spr <- function(mc, m) {
  plus_one <- ifelse(type == "quantile", log_first, FALSE)
  .chk(predict(mc, type = type, q = 1:50, newdata = model.frame(mc), smooth = TRUE, prob = 1:9/10) -
         (predict(m, type = type, q = 1:50 + log_first, newdata = model.frame(mc), prob = 1:9/10) - plus_one))
}

for (type in c("trafo", "distribution", "survivor", "density", "logdensity",
               "hazard", "loghazard", "cumhazard", "quantile")) {
  print(type)
.check_spr(mc3, m3)
}
  
### unconditional model for visual checks
y <- rpois(1000, lambda = 2)

## interval-censored response
yleft <- y - 1L
yleft[yleft < 0] <- -Inf

## quantiles & probabilities
q <- 0:10
yy <- seq(from = min(q), to = max(q), length.out = 50)
pr <- 1:9/10

## simple checks for predict wrt to newdata / q (log_first)
.check_prd <- function(mc) {
  ## check q
  stopifnot(all.equal(predict(mc, newdata = data.frame(q), type = "density"), 
                      predict(mc, newdata = data.frame(1), q = q, type = "density")))
  
  ## check q[1] density vs distribution
  stopifnot(all.equal(predict(mc, newdata = data.frame(1), q = 0, type = "density"),
                      predict(mc, newdata = data.frame(1), q = 0, type = "distribution")))
  
  ## rough check discrete density function
  stopifnot(all.equal(predict(mc, newdata = data.frame(1), q = q, type = "density"), 
                      predict(mc, newdata = data.frame(1), q = q, type = "distribution") - 
                        c(0, predict(mc, newdata = data.frame(1), q = q[-1] - 1 , type = "distribution"))))
  
  stopifnot(all.equal(predict(mc, newdata = data.frame(y = q), type = "density"), 
                      predict(mc, newdata = data.frame(y = q), type = "distribution") - 
                        c(0, predict(mc, newdata = data.frame(y = q[-1] - 1), type = "distribution"))))
  
}

## cotram: log_first = FALSE
mc1 <- cotram(y ~ 1, log_first = FALSE, prob = prob <- .99, extrapolate = TRUE)

.check_prd(mc1)

m1 <- Colr(yi ~ 1, data = data.frame(yi = Surv(yleft, y, type = "interval2")), extrapolate = TRUE,
            log_first = FALSE, support = c(0, quantile(y, prob = prob)), bounds = c(-.01, Inf))

## cotram: log_first = TRUE
mc2 <- cotram(y ~ 1, log_first = TRUE, extrapolate = TRUE, prob = .99)

.check_prd(mc2)

## <FIXME> works in package "tram" -- but should it? </FIXME>
try(predict(mc1, type = "quantile", prob = pr, smooth = TRUE, q = 1:10))

## <FIXME> when problems with one newdata configuration, can't predict at all. 
## (This comes from "tram") </FIXME>
try(predict(mc, type = "quantile", prob = pr, smooth = TRUE))

if (FALSE) {
layout(matrix(c(1:4), nrow = 2, byrow = TRUE))

main <- "predict - density: log_first = FALSE"
plot(q, predict(mc1, newdata = data.frame(1), q = q , type = "density"),
     type = "h", main = main)
points(q, predict(mc1, newdata = data.frame(1), q = q, type = "density"), pch = 20)
points(q, dpois(q, lambda = 2), col = "blue")
lines(yy, predict(mc1, newdata = data.frame(1), q = yy, type = "density", smooth = TRUE), col = "grey")
lines(yy, predict(mc1, newdata = data.frame(y = yy), type = "density", smooth = TRUE),
      col = "red")
plot(mc1, newdata = data.frame(1), type = "density", smooth = TRUE,
     col = "red", add = TRUE)
plot(mc1, newdata = data.frame(1), type = "density",
     col = "red", add = TRUE)

main <- "predict - distribution: log_first = FALSE"
plot(q, predict(mc1, newdata = data.frame(1), q = q, type = "distribution"),
     type = "s", ylim = c(0, 1), main = main)
lines(yy, predict(mc1, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE), col = "grey")
points(q, ppois(q, lambda = 2), col = "blue")
points(predict(mc1, newdata = data.frame(1), q = q, prob = pr, type = "quantile",
               smooth = TRUE), pr, col = "darkgreen")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1,
       col = c("blue", "green"))
plot(mc1, newdata = data.frame(1), type = "distribution",
     col = "red", smooth = TRUE, add = TRUE)
plot(mc1, newdata = data.frame(1), q = q, type = "distribution",
     col = "red", add = TRUE)

main <- "predict - density: log_first = TRUE"
plot(q, predict(mc2, newdata = data.frame(1), q = q, type = "density"),
     type = "h", main = main)
points(q, predict(mc2, newdata = data.frame(1), q = q, type = "density"), pch = 20)
points(q, dpois(q, lambda = 2), col = "blue")
lines(yy, predict(mc2, newdata = data.frame(1), type = "density",  q = yy, smooth = TRUE))
plot(mc2, newdata = data.frame(1), type = "density",
     col = "red", smooth = TRUE, add = TRUE, ylim = c(0, .3))
plot(mc2, newdata = data.frame(1), q = q, type = "density",
     col = "red", add = TRUE)

main <- "predict - distribution: log_first = TRUE"
plot(q, predict(mc2, newdata = data.frame(1), q = q, type = "distribution"),
     type = "s", ylim = c(0, 1), main = main)
points(q, ppois(q, lambda = 2), col = "blue")
lines(yy, predict(mc2, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE))
points(predict(mc2, newdata = data.frame(1),
               smooth = TRUE, prob = pr, type = "quantile"), pr, col = "green")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1,
       col = c("blue", "green"))
plot(mc2, newdata = data.frame(1), type = "distribution", q = 0:3,
     col = "red", smooth = TRUE, add = TRUE)
plot(mc2, newdata = data.frame(1), type = "distribution",
     col = "red", add = TRUE)

layout(matrix(c(1:6), nrow = 2, byrow = TRUE))
plot(mc3, newdata = nd,  type = "density", col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "density", col = 2, add = TRUE)
abline(v = 0)

plot(mc3, newdata = nd,  type = "logdensity", ylim = c(-6, 0), col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "logdensity", col = 2, add = TRUE, ylim = c(-6, 0))
abline(v = 0)

plot(mc3, newdata = nd, type = "distribution", col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "distribution", col = 2, add = TRUE)
abline(v = 0)

plot(mc3, newdata = nd, type = "trafo", col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "trafo", col = 2, add = TRUE)   
abline(v = 0)

plot(mc3, newdata = nd, type = "survivor", col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "survivor", col = 2, add = TRUE)   
abline(v = 0)

plot(mc3, newdata = nd, type = "cumhazard", col = rgb(.3, .3, .3, .3))
# plot(m3, newdata = nd, type = "cumhazard", col = 2, add = TRUE)   
abline(v = 0)
}

