
set.seed(290875)

ret <- TRUE

if (require("survival") && require("tram") && require("coin")) {

    d <- ovarian
    d$y <- with(d, Surv(futime, fustat))
    d$Ry <- with(d, R(y, as.R.ordered = TRUE))
    d$x <- d$rx
    tm <- sort(unique(d$futime[d$fustat == 1]))

    tol <- 1e-4

    ### nonpara ML, unconditional
    mP0 <- Polr(Ry ~ 1, data = d, method = "cloglog")
    s1 <- summary(survfit(y ~ 1, data = d), times = tm)$surv
    s2 <- predict(mP0, newdata = data.frame(1), type = "survivor")
    ret <- isTRUE(all.equal(s1, s2[-length(s2)], 
                     check.attributes = FALSE, tol = tol))  &&
                  cor(resid(mP0), logrank_trafo(d$y)) > .9999 &&
           isTRUE(all.equal(logLik(mP0, parm = coef(as.mlt(mP0))),
                            logLik(mP0, parm = log(-log(s1)))))

    ### conditional: estimate HR
    mP1 <- Polr(Ry ~ x, data = d, method = "cloglog")
    sf <- survdiff(y ~ x, data = d)
    mc1 <- coxph(y ~ x, data = d)

    ### see Survival Analysis: A Practical Approach, 
    ### Second Edition David Machin, Yin Bun Cheung, Mahesh K.B. Parmar
    ### Page 12
    oe <- sf$obs / sf$exp
    HR <- oe[1] / oe[2]

    ### compare: mP1 is nonpara ML, mc1 partial likelihood
    ### what is HR?
    #exp(coef(mP1))
    #exp(-coef(mc1))
}

ret
