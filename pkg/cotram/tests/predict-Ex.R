
library("cotram")

set.seed(25)

y <- rpois(1000, lambda = 2)

## interval-censored response
yleft <- y - 1L
yleft[yleft < 0] <- -Inf

## smooth response
yy <- seq(from = min(q), to = max(q), length.out = 50)

## quantiles & probabilities
q <- 0:10
pr <- 1:9/10

## simple checks for predict
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
mc1 <- cotram(y ~ 1, log_first = FALSE, extrapolate = TRUE, prob = .99)

.check_prd(mc1)

## cotram: log_first = TRUE
mc2 <- cotram(y ~ 1, log_first = TRUE, extrapolate = TRUE, prob = .99)

m2 <- Coxph(yi ~ 1, data = data.frame(yi = Surv(yleft, y, type = "interval2")))

.check_prd(mc2)


### some basic checks for predictions wrt discrete distributions
n <- 50
d <- data.frame(x1 = 1:n, x2 = sample(1:n) + 1, q = sample(1:n))
mod <- cotram(q ~ x1 + x2, data = d)

mmod <- Colr(q ~ x1 + x2, data = d)

.chk <- function(x)
  stopifnot(max(abs(x), na.rm = TRUE) < sqrt(.Machine$double.eps))

nd <- d
nd$q <- NULL
q <-0:(n+1)
p <- predict(mod, newdata = nd, q = q, type = "distribution")
s <- predict(mod, newdata = nd, q = q, type = "survivor")
.chk(predict(mod, newdata = nd, q = q, type = "distribution", log = TRUE) - log(p))
.chk(predict(mod, newdata = nd, q = q, type = "distribution", lower.tail = FALSE) - s)
.chk(predict(mod, newdata = nd, q = q, type = "distribution", 
             lower.tail = FALSE, log = TRUE) - log(s))

o <- predict(mod, newdata = nd, q = q, type = "odds")
.chk(o - p / s)

dd <- predict(mod, newdata = nd, q = q, type = "density")

.chk(apply(p, 2, function(x) diff(c(0, x))) - dd)


h <- predict(mod, newdata = nd, q = q , type = "hazard")

.chk(dd / (1 - (p - dd)) - h)

.chk(apply(h, 2, function(x) cumprod(1 - x)) - s)

H <- predict(mod, newdata = nd, q = q, type = "cumhazard")

.chk(H + log(s))


## predict
stopifnot(all.equal(predict(m, type = "density"), predict(m, newdata = data.frame(q = q, x = x), type = "density")))
stopifnot(all.equal(predict(m, type = "density"), predict(m, newdata = model.frame(m), type = "density")))

ppr <- predict(mm, newdata = nd, type = "density")
pr <- predict(m, newdata = nd, type = "density")

stopifnot(rownames(pr) == rownames(ppr))
stopifnot(!any(abs(ppr - pr) > 0.1))

pprq <- predict(mm, newdata = nd, type = "density", q = 0:10)
prq <- predict(m, newdata = nd, type = "density", q = 0:10)

stopifnot(rownames(prq) == rownames(pprq))
stopifnot(!any(abs(pprq - prq) > 0.1))



## <CHECK>
length(predict(mc2, type = "distribution"))
length(predict(m2, type = "distribution", newdata = model.frame(mc2)))

# mt <- Coxph(yp1 ~ 1, data = data.frame(yp1 = y + 1),
#             support = mc$support,
#             bounds = mc$bounds,
#             log_first = TRUE)
# predict(mc2, type = "distribution") - predict(mt, type = "distribution")
## <\CHECK>

if (FALSE) {
layout(matrix(c(1:4), nrow = 2, byrow = TRUE))

main <- "predict - density: log_first = FALSE"
plot(q, predict(m1, newdata = data.frame(1), q = q , type = "density"),
     type = "h", main = main)
points(q, predict(m1, newdata = data.frame(1), q = q, type = "density"), pch = 20)
points(q, dpois(q, lambda = 2), col = "blue")
lines(yy, predict(m1, newdata = data.frame(1), q = yy, type = "density", smooth = TRUE), col = "grey")
lines(yy, predict(m1, newdata = data.frame(y = yy), type = "density", smooth = TRUE),
      col = "red")
plot(m1, newdata = data.frame(1), type = "density", smooth = TRUE,
     col = "red", add = TRUE)
plot(m1, newdata = data.frame(1), type = "density",
     col = "red", add = TRUE)

main <- "predict - distribution: log_first = FALSE"
plot(q, predict(m1, newdata = data.frame(1), q = q, type = "distribution"),
     type = "s", ylim = c(0, 1), main = main)
lines(yy, predict(m1, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE), col = "grey")
points(q, ppois(q, lambda = 2), col = "blue")
points(predict(m1, newdata = data.frame(1), q = q, prob = pr, type = "quantile",
               smooth = TRUE), pr, col = "darkgreen")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1,
       col = c("blue", "green"))
plot(m1, newdata = data.frame(1), type = "distribution",
     col = "red", smooth = TRUE, add = TRUE)
plot(m1, newdata = data.frame(1), q = q, type = "distribution",
     col = "red", add = TRUE)

main <- "predict - density: log_first = TRUE"
plot(q, predict(m2, newdata = data.frame(1), q = q, type = "density"),
     type = "h", main = main)
points(q, predict(m2, newdata = data.frame(1), q = q, type = "density"), pch = 20)
points(q, dpois(q, lambda = 2), col = "blue")
lines(yy, predict(m2, newdata = data.frame(1), type = "density",  q = yy, smooth = TRUE))
plot(m2, newdata = data.frame(1), type = "density",
     col = "red", smooth = TRUE, add = TRUE, ylim = c(0, .3))
plot(m2, newdata = data.frame(1), q = q, type = "density",
     col = "red", add = TRUE)

main <- "predict - distribution: log_first = TRUE"
plot(q, predict(m2, newdata = data.frame(1), q = q, type = "distribution"),
     type = "s", ylim = c(0, 1), main = main)
points(q, ppois(q, lambda = 2), col = "blue")
lines(yy, predict(m2, newdata = data.frame(1), q = yy, type = "distribution", smooth = TRUE))
points(predict(m2, newdata = data.frame(1),
               smooth = TRUE, prob = pr, type = "quantile"), pr, col = "green")
legend("bottomright", legend = c("ppois", "quantiles"), pch = 1,
       col = c("blue", "green"))
plot(m2, newdata = data.frame(1), type = "distribution", q = 0:3,
     col = "red", smooth = TRUE, add = TRUE)
plot(m2, newdata = data.frame(1), type = "distribution",
     col = "red", add = TRUE)
}

## plot
if (FALSE) {
  plot(m, newdata = nd,  type = "density", col = 1)
  plot(mm, newdata = nd, type = "density", col = 2, add = TRUE)
  abline(v = 0)
  
  plot(m, newdata = nd,  type = "logdensity", col = 1, ylim = c(-6, 0))
  plot(mm, newdata = nd, type = "logdensity", col = 2, add = TRUE, ylim = c(-6, 0))
  abline(v = 0)
  
  plot(m, newdata = nd, type = "distribution", col = 1)
  plot(mm, newdata = nd, type = "distribution", col = 2, add = TRUE)
  abline(v = 0)
  
  plot(m, newdata = nd, type = "trafo", col = 1)
  plot(mm, newdata = nd, type = "trafo", col = 2, add = TRUE)   
  abline(v = 0)
  
  plot(m, newdata = nd, type = "survivor", col = 1)
  plot(mm, newdata = nd, type = "survivor", col = 2, add = TRUE)   
  abline(v = 0)
  
  plot(m, newdata = nd, type = "cumhazard", col = 1)
  plot(mm, newdata = nd, type = "cumhazard", col = 2, add = TRUE)   
  abline(v = 0)
}