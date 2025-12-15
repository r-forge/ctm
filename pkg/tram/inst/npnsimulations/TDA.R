## ----setup, echo = FALSE, warning = FALSE, message = FALSE--------------------
pkgs <- c("openxlsx", "tram", "mvtnorm", "latticeExtra", "qrng", "survival", "CVXR", 
          "colorspace", "lattice")
req <- sapply(pkgs, require, char = TRUE)
if (!all(req)){
  sapply(pkgs[!req], install.packages)
  req[!req] <- sapply(pkgs[!req], require, char = TRUE)
}
if (!all(req)) 
  stop("cannot load dependencies")

frt <- function(x, digits = 3, width = 3, pval = FALSE) {
    if (pval) {
        x[is.na(x)] <- 1
        idx <- (x < 10^(-width))
        x[idx] <- 10^(-width)
        ret <- formatC(x, digits = digits, width = digits, format = "f")
        return(paste0(ifelse(idx, "<", ""), ret))
    }
    formatC(x, digits = digits, width = digits, format = "f")
}

cols <- diverging_hcl(2, palette = "Blue-Red", alpha = .2)
cols2 <- diverging_hcl(2, palette = "Blue-Red", alpha = .8)


## ----HCCdata, echo = FALSE, cache = TRUE--------------------------------------
### Load data
dat <- read.xlsx("https://datadryad.org/api/v2/files/44697/download", sheet = 1)
d <- with(dat, data.frame(id = 1:nrow(dat),
                          D = HCC_studyGr,
                          AFP = AFP_ng_per_ml,
                          PIV = PIVKA_delete_range,
                          OPN,
                          DKK))
d$D <- factor(d$D, levels = 0:1, labels = c("non-HCC", "HCC"))


## Biomarker columns
bcol <- c("DKK", "OPN", "PIV", "AFP")
J <- length(bcol)
# log all biomarkers
dd <- d
dd[, bcol] <- log(dd[, bcol])
dd$Di <- (0:1)[dd$D]

## censoring
dd$iDKK <- with(dd, R(DKK, as.R.interval = TRUE))
dd$iOPN <- with(dd, R(OPN, as.R.interval = TRUE))
PIVm <- max(dd$PIV)
dd$iPIV <- with(dd, R(Surv(PIV, event = PIV < PIVm), as.R.interval = TRUE))
AFPm <- max(dd$AFP)
dd$iAFP <- with(dd, R(Surv(AFP, event = AFP < AFPm), as.R.interval = TRUE))
dd$out <- with(dd, PIV >= PIVm  | AFP > AFPm)

d2c <- function(x, K = 10) 
    cut(x, breaks = c(-Inf, quantile(x, prob = 1:(K - 1) / K), Inf), ordered_result = TRUE)

dd$cDKK <- with(dd, d2c(DKK))
dd$cOPN <- with(dd, d2c(OPN))
dd$cPIV <- with(dd, d2c(PIV))
dd$cAFP <- with(dd, d2c(AFP))

nd0 <- nd1 <- dd
nd0$D <- sort(unique(dd$D))[1L]
nd0$Di <- (0:1)[nd0$D]
nd1$D <- sort(unique(dd$D))[2L]
nd1$Di <- (0:1)[nd1$D]
Xd0 <- t(model.matrix(~ D, data = nd0))
Xd1 <- t(model.matrix(~ D, data = nd1))

OR <- 6 



## ----HCCplot, echo = FALSE, fig.width = 8, fig.height = 4---------------------
pd <- reshape(dd[, c("D", bcol)], idvar = "id", direction = "long", varying = list(bcol))
pd$Biomarker <- bcol[pd$time]
names(pd)[names(pd) == bcol[1]] <- "obs"
ecdfplot(~ obs | Biomarker, data = pd, groups = D, xlab = "log(Biomarker)", 
         par.settings = list(superpose.line = list(col=cols2)),
         auto.key=list(space="top", lines=TRUE), col = cols2, layout = c(4, 1))


## ----LDA_convex, echo = FALSE, cache = TRUE-----------------------------------
Y <- t(dd[, bcol])
Xd <- t(model.matrix(~ D, data = dd))
L <- Variable(J, J)
B <- Variable(J, 2)
constr <- list(L[upper.tri(L, diag = FALSE)] == 0, diag(L) >= 0)
obj <- 2 * sum(log(diag(L))) - sum_squares(L %*% Y - B %*% Xd) / ncol(Y)
prob <- Problem(Maximize(obj), constr)
result <- solve(prob)#, solver = "SCS")
Lhat <- result$getValue(L)
Lhat <- ltMatrices(Lhat[lower.tri(Lhat, diag = TRUE)], diag = TRUE)
Bhat <- result$getValue(B)
Lhat1 <- ltMatrices(invcholD(Lhat, D = 1 / diagonals(Lhat)), byrow = TRUE)
mu <- solve(Lhat, Bhat)
mu0 <- mu[,1]
mu1 <- rowSums(mu)
mu <- cbind(mu0, mu1)
LL_BB <- ldpmvnorm(obs = Y, invchol = Lhat, mean = mu[,dd$D])
LLR_BB <- .5 * (-colSums((Mult(Lhat, Y) - Bhat %*% Xd0)^2) + colSums((Mult(Lhat, Y) - Bhat %*% Xd1)^2))


## ----CVXRtm, echo = FALSE, cache = TRUE---------------------------------------
itm <- 1:9
tm <- sapply(itm, function(i) system.time(Problem(Maximize(obj), constr))["user.self"])
tm_CVXR <- median(tm)


## ----LDA_mmlt, echo = FALSE, cache = TRUE-------------------------------------
### same normal model, via mmlt
m1 <- Lm(DKK ~ D, data = dd)
m2 <- Lm(OPN ~ D, data = dd)
m3 <- Lm(PIV ~ D, data = dd)
m4 <- Lm(AFP ~ D, data = dd)

m <- Mmlt(m1, m2, m3, m4, formula = ~ 1, data = dd)
LL_LDA <- logLik(m)
# coef(m, type = "Lambdapar")
LLR_LDA <- predict(m, newdata = nd0, type = "density", log = TRUE) -
          predict(m, newdata = nd1, type = "density", log = TRUE)
stopifnot(isTRUE(all.equal(LL_BB, c(LL_LDA), tol = 1e-5)))
stopifnot(isTRUE(all.equal(LLR_BB, LLR_LDA, tol = 1e-4)))


## ----ctab1, echo = FALSE------------------------------------------------------
x <- as.array(Lhat1)[,,1]
ctab <- c(ltMatrices(ltMatrices(x[lower.tri(x)]), byrow = TRUE))
cm <- 1:length(do.call("c", coef(m, type = "marginal")))
cf <- coef(m)
vc <- sqrt(diag(vcov(m)))
ctab <- cbind(ctab, cbind(cf[-cm], vc[-cm]))


## ----LDAtm, echo = FALSE, cache = TRUE----------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3, m4, data = dd, formula = ~ 1))["user.self"])
tm_LDA <- median(tm)


## ----lTDA_mmlt, echo = FALSE, cache = TRUE------------------------------------
### same normal model, via mmlt
m1 <- BoxCox(DKK ~ D, data = dd, order = OR)
m2 <- BoxCox(OPN ~ D, data = dd, order = OR)
m3 <- BoxCox(PIV ~ D, data = dd, order = OR)
m4 <- BoxCox(AFP ~ D, data = dd, order = OR)

m <- Mmlt(m1, m2, m3, m4, formula = ~ 1, data = dd)
LL_lTDA <- logLik(m)
# coef(m, type = "Lambdapar")
LLR_lTDA <- predict(m, newdata = nd0, type = "density", log = TRUE) -
            predict(m, newdata = nd1, type = "density", log = TRUE)


## ----lTDAtm, echo = FALSE, cache = TRUE---------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3, m4, data = dd, formula = ~ 1))["user.self"])
tm_lTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cm <- 1:length(do.call("c", coef(m, type = "marginal")))
cf <- coef(m)
vc <- sqrt(diag(vcov(m)))
ctab <- cbind(ctab, cbind(cf[-cm], vc[-cm]))


## ----lsTDA_mmlt, echo = FALSE, cache = TRUE-----------------------------------
### same normal model, via mmlt
m1 <- BoxCox(DKK ~ D | D, data = dd, order = OR)
m2 <- BoxCox(OPN ~ D | D, data = dd, order = OR)
m3 <- BoxCox(PIV ~ D | D, data = dd, order = OR)
m4 <- BoxCox(AFP ~ D | D, data = dd, order = OR)

m <- Mmlt(m1, m2, m3, m4, formula = ~ 1, data = dd)
LL_lsTDA <- logLik(m)
# coef(m, type = "Lambdapar")
LLR_lsTDA <- predict(m, newdata = nd0, type = "density", log = TRUE) -
          predict(m, newdata = nd1, type = "density", log = TRUE)



## ----lsTDAtm, echo = FALSE, cache = TRUE--------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3, m4, data = dd, formula = ~ 1))["user.self"])
tm_lsTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cm <- 1:length(do.call("c", coef(m, type = "marginal")))
cf <- coef(m)
vc <- sqrt(diag(vcov(m)))
ctab <- cbind(ctab, cbind(cf[-cm], vc[-cm]))


## ----elsTDA_mmlt, echo = FALSE, cache = TRUE----------------------------------
### 2 continuous, 2 interval, shift-scale
m1 <- as.mlt(BoxCox(DKK ~ Di | Di, data = dd, order = OR))
m2 <- as.mlt(BoxCox(OPN ~ Di | Di, data = dd, order = OR))
m3i <- as.mlt(BoxCox(iPIV ~ Di | Di, data = dd, order = OR))
m4i <- as.mlt(BoxCox(iAFP ~ Di | Di, data = dd, order = OR))

M <- 500
args <- list(seed = 1, type = "ghalton", M = M)

mi <- Mmlt(m1, m2, m3i, m4i, formula = ~ 1, data = dd, args = args)
LL_elsTDA <- logLik(mi)
# coef(mi, type = "Lambdapar")

tmp1 <- BoxCox(DKK ~ Di | Di, data = rbind(nd0, nd1), order = OR, dofit = FALSE, LRtest = FALSE, theta = coef(m1))
tmp2 <- BoxCox(OPN ~ Di | Di, data = rbind(nd0, nd1), order = OR, dofit = FALSE, LRtest = FALSE, theta = coef(m2))
tmp3 <- BoxCox(iPIV ~ Di | Di, data = rbind(nd0, nd1), order = OR, dofit = FALSE, LRtest = FALSE, theta = coef(m3i))
tmp4 <- BoxCox(iAFP ~ Di | Di, data = rbind(nd0, nd1), order = OR, dofit = FALSE, LRtest = FALSE, theta = coef(m4i))
tmp <- mmlt(tmp1, tmp2, tmp3, tmp4, formula = ~ 1, data = rbind(nd0, nd1), dofit = FALSE, theta = coef(mi), args = args)
LLR_elsTDA <- matrix(get("ll", env = environment(tmp$ll))(coef(mi)), ncol = 2) %*% c(1, -1)


## ----elsTDAtm, echo = FALSE, cache = TRUE-------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3i, m4i, data = dd, formula = ~ 1, args = args))["user.self"])
tm_elsTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cm <- 1:length(do.call("c", coef(m, type = "marginal")))
cf <- coef(m)
vc <- sqrt(diag(vcov(m)))
ctab <- cbind(ctab, cbind(cf[-cm], vc[-cm]))


## ----table, echo = FALSE, results = "asis"------------------------------------
rn <- gsub(".(Intercept)", "", rownames(ctab), fixed = TRUE)
rn <- gsub(".", ",", rn, fixed = TRUE)
seM <- frt(ctab[nrow(ctab), ncol(ctab)], digits = 3, width = 3)
for (i in 1:nrow(ctab))
    cat(rn[i], " & ", paste(paste("$", frt(ctab[i,], digits = 3, width = 3), "$"), collapse = " &"), "\\\\", "\n")


## ----echo = FALSE-------------------------------------------------------------
ctab <- cf[-cm]


## ----convex-domargins, echo = FALSE-------------------------------------------
mi <- Mmlt(m1, m2, m3i, m4i, formula = ~ 1, data = dd, args = args, domargins = FALSE)
LL_dmlsTDA <- logLik(mi)


## ----dmlsTDAtm, echo = FALSE, cache = TRUE------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3i, m4i, data = dd, formula = ~ 1, args = args, domargins = FALSE))["user.self"])
tm_dmlsTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cf <- coef(mi)
ctab <- cbind(ctab, cf[-cm])


## ----convex-ACS, echo = FALSE, cache = TRUE-----------------------------------
mi <- Mmlt(m1, m2, m3i, m4i, formula = ~ 1, data = dd, args = args, fit = "ACS")
LL_ACSsTDA <- logLik(mi)


## ----ACSsTDAtm, echo = FALSE, cache = TRUE------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3i, m4i, formula = ~ 1, data = dd, args = args, fit = "ACS"))["user.self"])
tm_ACSsTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cf <- coef(mi, fixed = TRUE)
ctab <- cbind(ctab, cf[-cm])


## ----convex-seq, echo = FALSE, cache = TRUE-----------------------------------
mi <- Mmlt(m1, m2, m3i, m4i, formula = ~ 1, data = dd, args = args, fit = "sequential")
LL_sqlsTDA <- logLik(mi)


## ----sqlsTDAtm, echo = FALSE, cache = TRUE------------------------------------
tm <- sapply(itm, function(i) system.time(Mmlt(m1, m2, m3i, m4i, data = dd, formula = ~ 1, 
    args = args, fit = "sequential"))["user.self"])
tm_sqlsTDA <- median(tm)


## ----echo = FALSE-------------------------------------------------------------
cf <- coef(mi, fixed = TRUE)
ctab <- cbind(ctab, cf[-cm])


## ----tableconv, echo = FALSE, results = "asis"--------------------------------
rn <- gsub(".(Intercept)", "", rownames(ctab), fixed = TRUE)
rn <- gsub(".", ",", rn, fixed = TRUE)
for (i in 1:nrow(ctab))
    cat(rn[i], " & ", paste(paste("$", frt(ctab[i,], digits = 3, width = 3), "$"), collapse = " &"), "\\\\", "\n")

