
set.seed(290875)

ret <- TRUE

if (require("survival") && require("tram") && require("coin")) {

    d <- ovarian
    d$y <- with(d, Surv(futime, fustat))
    d$Ry <- with(d, R(y, as.R.ordered = TRUE))
    d$x <- d$age - 50
    tm <- sort(unique(d$futime[d$fustat == 1]))

    tol <- 1e-4

    mP0 <- Polr(Ry ~ 1, data = d, method = "cloglog")
    s1 <- summary(survfit(y ~ 1, data = d), times = tm)$surv
    s2 <- predict(mP0, newdata = data.frame(1), type = "survivor")
    ret <- all.equal(s1, s2[-length(s2)], check.attributes = FALSE, tol = tol)  &&
    cor(resid(mP0), logrank_trafo(d$y)) > .9999

}

ret
