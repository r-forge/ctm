
library("tram")
library("survival")

### check if R(Surv(...), as.R.ordered = TRUE) gives
### results identical to KaplanMeier, also in case
### of ties. 

### Bland & Altman, BMJ 1998;317:1572
prg <- data.frame(time = c(1, 1, 1, 1, 1, 1, 2, 2, 
                           2, 2, 2, 3, 3, 3, 4,
                           4, 4, 6, 6, 9, 9, 9, 10, 13, 16,
                           2, 3, 4, 7, 7, 8, 8, 9, 9, 9, 
                           11, 24, 24),
                  event = rep(c(1, 0), c(25, 13)))
ut <- with(prg, sort(unique(time[event == 1])))

m <- Coxph(R(Surv(time, event), as.R.ordered = TRUE) ~ 1, data = prg)
t1 <- predict(as.mlt(m), newdata = data.frame(time = 1), type = "surv")
t2 <- c(summary(survfit(Surv(time, event) ~ 1, data = prg), time = ut)$surv, 0)
stopifnot(all.equal(t1, t2, check.attributes = FALSE, tol = 1e-4))

### AML
ut <- with(aml, sort(unique(time[status == 1])))

m <- Coxph(R(Surv(time, status), as.R.ordered = TRUE) ~ 1, data = aml)
t1 <- predict(as.mlt(m), newdata = data.frame(1), type = "surv")
t2 <- c(summary(survfit(Surv(time, status) ~ 1, data = aml), time = ut)$surv, 0)

stopifnot(all.equal(t1, t2, check.attributes = FALSE, tol = 1e-4))

