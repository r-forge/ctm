
library("tram")
library("multcomp")

m <- BoxCox(dist ~ speed, data = cars)
m$negative
coef(m)
nd <- data.frame(speed = as.double(1:5 * 5))

PI(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
PI(m, newdata = nd, reference = lp15)
PI(m, newdata = nd, reference = nd[3,,drop = FALSE])
PI(m, newdata = nd, reference = nd)

F <- m$model$todistr$p
Q <- m$model$todistr$q
f <- function(p, b) F(Q(p) - b)
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])

OVL(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
OVL(m, newdata = nd, reference = lp15)
OVL(m, newdata = nd, reference = nd[3,,drop = FALSE])
OVL(m, newdata = nd, reference = nd)

d <- m$model$todistr$d
Q <- m$model$todistr$q
f <- function(p, b) pmin(1, d(Q(p) - b) / d(Q(p)))
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])

m <- Colr(dist ~ speed, data = cars)
m$negative
coef(m)
nd <- data.frame(speed = as.double(1:5 * 5))

PI(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
PI(m, newdata = nd, reference = lp15)
PI(m, newdata = nd, reference = nd[3,,drop = FALSE])
PI(m, newdata = nd, reference = nd)

F <- m$model$todistr$p
Q <- m$model$todistr$q
f <- function(p, b) F(Q(p) + b)
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])


OVL(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
OVL(m, newdata = nd, reference = lp15)
OVL(m, newdata = nd, reference = nd[3,,drop = FALSE])
OVL(m, newdata = nd, reference = nd)

d <- m$model$todistr$d
Q <- m$model$todistr$q
f <- function(p, b) pmin(1, d(Q(p) + b) / d(Q(p)))
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])


m <- Coxph(dist ~ speed, data = cars)
m$negative
coef(m)
nd <- data.frame(speed = as.double(1:5 * 5))

PI(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
PI(m, newdata = nd, reference = lp15)
PI(m, newdata = nd, reference = nd[3,,drop = FALSE])
PI(m, newdata = nd, reference = nd)

F <- m$model$todistr$p
Q <- m$model$todistr$q
f <- function(p, b) F(Q(p) + b)
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])


OVL(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
OVL(m, newdata = nd, reference = lp15)
OVL(m, newdata = nd, reference = nd[3,,drop = FALSE])
OVL(m, newdata = nd, reference = nd)

d <- m$model$todistr$d
Q <- m$model$todistr$q
f <- function(p, b) pmin(1, d(Q(p) + b) / d(Q(p)))
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])


m <- Lehmann(dist ~ speed, data = cars)
m$negative
coef(m)
nd <- data.frame(speed = as.double(1:5 * 5))

PI(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
PI(m, newdata = nd, reference = lp15)
PI(m, newdata = nd, reference = nd[3,,drop = FALSE])
PI(m, newdata = nd, reference = nd)

F <- m$model$todistr$p
Q <- m$model$todistr$q
f <- function(p, b) F(Q(p) - b)
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])


OVL(m, newdata = nd)
lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
OVL(m, newdata = nd, reference = lp15)
OVL(m, newdata = nd, reference = nd[3,,drop = FALSE])
OVL(m, newdata = nd, reference = nd)

d <- m$model$todistr$d
Q <- m$model$todistr$q
f <- function(p, b) pmin(1, d(Q(p) - b) / d(Q(p)))
lp <- predict(m, newdata = nd, type = "lp")
integrate(f, lower = 0, upper = 1, b = lp[2] - lp[1])

OVL(m, newdata = nd)
OVL(m, newdata = nd, conf.level = .95)

lp15 <- c(predict(m, newdata = data.frame(speed = 15)))
OVL(m, newdata = nd, reference = lp15)
OVL(m, newdata = nd, reference = lp15, conf.level = .95)

OVL(m, newdata = nd[-3,,drop = FALSE], reference = nd[3,,drop = FALSE])
OVL(m, newdata = nd[-3,,drop = FALSE], reference = nd[3,,drop = FALSE], conf.level = .95)

OVL(m, newdata = nd, reference = nd)
OVL(m, newdata = nd, reference = nd, conf.level = .95)

OVL(m, newdata = nd, reference = nd, conf.level = .95, 
    calpha = univariate_calpha())
OVL(m, newdata = nd, reference = nd, conf.level = .95, 
    calpha = adjusted_calpha())
