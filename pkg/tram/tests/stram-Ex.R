
library("tram")
data("BostonHousing2", package = "mlbench")

m <- Lm(cmedv ~ nox 
        # | rm, 
        , data = BostonHousing2)

d <- predict(as.mlt(m), type = "density")
ll <- sum(log(d))
all.equal(c(logLik(m)), ll, check.attributes = FALSE)
d1 <- predict(as.mlt(m), newdata = BostonHousing2[1:6,], type = "density")
d2 <- diag(predict(as.mlt(m), newdata = BostonHousing2[1:6,-6], q =
             BostonHousing2[1:6,6], type = "density"))
all.equal(d1, d2, check.attributes = FALSE)
prb <- 1:9 / 10
q <- predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], prob = prb,
        type = "quantile")
p1 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,1], 
        type = "distribution")[,1], 2)
p2 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,2], 
        type = "distribution")[,2], 2)
all.equal(p1, prb, check.attributes = FALSE)
all.equal(p2, prb, check.attributes = FALSE)

m <- Lm(cmedv ~ nox 
        | rm, 
        , data = BostonHousing2)

d <- predict(as.mlt(m), type = "density")
ll <- sum(log(d))
all.equal(c(logLik(m)), ll, check.attributes = FALSE)
d1 <- predict(as.mlt(m), newdata = BostonHousing2[1:6,], type = "density")
d2 <- diag(predict(as.mlt(m), newdata = BostonHousing2[1:6,-6], q =
             BostonHousing2[1:6,6], type = "density"))
all.equal(d1, d2, check.attributes = FALSE)
prb <- 1:9 / 10
q <- predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], prob = prb,
        type = "quantile")
p1 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,1], 
        type = "distribution")[,1], 2)
p2 <- round(predict(as.mlt(m), newdata = BostonHousing2[1:2,-6], q = q[,2], 
        type = "distribution")[,2], 2)
all.equal(p1, prb, check.attributes = FALSE)
all.equal(p2, prb, check.attributes = FALSE)

