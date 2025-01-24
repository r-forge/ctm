## ----setup, echo = FALSE, message = FALSE-------------------------------------
## Reprocibility Material for
##
##	Adjusted Marginal Inference under Noncollapsibility
##	by Susanne Dandl & Torsten Hothorn, UZH
##

## packages
pkgs <- c("tram", "TH.data", "multcomp", "survival")

ip <- rownames(installed.packages())
if (any(!pkgs %in% ip))
    install.packages(pkgs[!pkgs %in% ip], repos = "https://stat.ethz.ch/CRAN/")

OK <- sapply(pkgs, require, character.only = TRUE)
if (!all(OK)) 
    stop("package(s) ", paste(pkgs[!OK], collapse = ", "), " not available")

## some helper functions
frmt1 <- function(x, digits = 1) {
    if (!is.numeric(x)) return(x)
    formatC(round(x, digits = digits), digits = digits, format = "f") 
}
frmt3 <- function(x) 
    frmt1(x, digits = 3)
frmtci <- function(x, digits = 3) {
    if (!is.numeric(x)) return(x)
    if (length(x) != 2) stop("not a confidence interval")
    return(paste("(", frmt1(x[1], digits = digits), 
                 ",", frmt1(x[2], digits = digits), ")"))
}

## model checks
CHECK <- interactive()

## discrete nonparanormal likelihood relies on Monte Carlo, set seed
set.seed(221224L)



## ----Chloramine, echo = FALSE, message = FALSE, results = "hide"--------------
## data taken from package SiTuR (not on CRAN, therefore dumped)
# dataset: create baseline var as time point = week 1 
# data("immun", package = "SiTuR")
# immun <- do.call(rbind, lapply(split(immun, immun$Anino), 
#                                function(x) {
#                                  cbind(x, baseline = with(x, x$weight[x$time == 1]))
#                                }
# ))
# immun <- immun[immun$time != 1,]
# y <- "weight"
# w <- "dose"
# x <- "baseline"
# 
# generate_data <- function(dose, time) {
#     dt <- immun[immun$dose %in% c(0, dose) & immun$time == time,]
#     dt <- dt[c(y, w, x)] 
#     names(dt) <- c("y", "w", "x")
#     dt
# }
#
# d <- generate_data(100, 29)
# d$w <- d$w[, drop = TRUE]
d <-
structure(list(y = c(21.699999999999999, 23.899999999999999, 
22.699999999999999, 23.399999999999999, 26.800000000000001, 24.800000000000001, 
23.399999999999999, 25.100000000000001, 24.100000000000001, 23.300000000000001, 
25.399999999999999, 25.100000000000001, 23.800000000000001, 23.100000000000001, 
24, 24.199999999999999, 27.399999999999999, 23.300000000000001, 
22.600000000000001, 23.300000000000001, 24.800000000000001, 23.899999999999999, 
22.199999999999999, 21.399999999999999, 22.800000000000001, 22.300000000000001, 
22.399999999999999, 30.399999999999999, 30.600000000000001, 21.699999999999999, 
24.800000000000001, 25.600000000000001, 21.399999999999999, 24.300000000000001, 
23.5, 25.800000000000001, 21.600000000000001, 22.899999999999999, 
23.800000000000001, 22.600000000000001, 24.199999999999999, 24.300000000000001, 
25.699999999999999, 23.199999999999999, 24.600000000000001, 24.5, 
22.699999999999999, 26.300000000000001, 27.199999999999999, 27.100000000000001, 
22.699999999999999, 24.600000000000001, 23, 23.199999999999999, 
23.899999999999999, 23, 20.800000000000001, 23.399999999999999, 
24.300000000000001, 24.399999999999999, 22.600000000000001, 22.100000000000001, 
22.199999999999999, 24.100000000000001, 28.100000000000001, 23.399999999999999, 
26.800000000000001, 24, 25.899999999999999, 24.699999999999999, 
24.100000000000001, 26.899999999999999, 23.899999999999999, 24.399999999999999, 
25.199999999999999, 22.899999999999999, 25.699999999999999, 24.300000000000001, 
25.199999999999999, 24.100000000000001), w = structure(c(1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), levels = c("0", 
"100"), class = "factor"), x = c(19, 22, 19.899999999999999, 
20.800000000000001, 22.899999999999999, 21, 21.800000000000001, 
22.600000000000001, 21.5, 22.600000000000001, 21.5, 22, 20.699999999999999, 
21.699999999999999, 20.800000000000001, 22.300000000000001, 21.800000000000001, 
19.800000000000001, 20.699999999999999, 20.600000000000001, 19.899999999999999, 
20.199999999999999, 20.199999999999999, 19, 19.600000000000001, 
18.699999999999999, 18.699999999999999, 21.699999999999999, 21.300000000000001, 
19.199999999999999, 21, 21.100000000000001, 18.699999999999999, 
20.300000000000001, 21.600000000000001, 21.899999999999999, 18.100000000000001, 
19.199999999999999, 20, 18.899999999999999, 20.399999999999999, 
20.5, 21.300000000000001, 19.199999999999999, 19.399999999999999, 
20.199999999999999, 18.199999999999999, 21, 21.5, 22.100000000000001, 
19.600000000000001, 21.399999999999999, 19.199999999999999, 21, 
20.899999999999999, 19.600000000000001, 19, 21.300000000000001, 
22.199999999999999, 21.699999999999999, 20.800000000000001, 20.600000000000001, 
19.199999999999999, 20.800000000000001, 23.600000000000001, 21.199999999999999, 
22.600000000000001, 21.100000000000001, 23.300000000000001, 22.300000000000001, 
21.600000000000001, 22.600000000000001, 21.199999999999999, 23.899999999999999, 
23.800000000000001, 20.699999999999999, 22.199999999999999, 22.5, 
22.199999999999999, 22.699999999999999)), row.names = c("1.5", 
"2.10", "3.15", "4.20", "5.25", "6.30", "7.35", "8.40", "33.165", 
"34.170", "35.175", "36.180", "37.185", "38.190", "39.195", "40.200", 
"57.285", "58.290", "59.295", "60.300", "61.305", "62.310", "63.315", 
"64.320", "89.445", "90.450", "91.455", "92.460", "93.465", "94.470", 
"95.475", "96.480", "121.605", "122.610", "123.615", "124.620", 
"125.625", "126.630", "127.635", "128.640", "153.765", "154.770", 
"155.775", "156.780", "157.785", "158.790", "159.795", "160.800", 
"177.885", "178.890", "179.895", "180.900", "181.905", "182.910", 
"183.915", "184.920", "209.1045", "210.1050", "211.1055", "212.1060", 
"213.1065", "214.1070", "215.1075", "216.1080", "233.1165", "234.1170", 
"235.1175", "236.1180", "237.1185", "238.1190", "239.1195", "240.1200", 
"265.1325", "266.1330", "267.1335", "268.1340", "269.1345", "270.1350", 
"271.1355", "272.1360"), class = "data.frame")

## marginal outcome normal model feat. Cohen's d
m0 <- Lm(y ~ w, data = d)
## marginal model for baseline weight
m1 <- BoxCox(x ~ 1, data = d)
## multivariate transformation model
m <- mmlt(m0, m1, formula = ~ 1, data = d)

## marginal model
(cf0 <- coef(m0))
(ci0 <- confint(m0))
## multivariate model
(cf1 <- coef(m)["y.w100"])
(ci1 <- confint(m)["y.w100",])
## latent correlation
(r <- c(unclass(coef(m, type = "Corr"))))

## observed and 
(se0 <- sqrt(vcov(m0)))
(se1 <- sqrt(diag(vcov(m))["y.w100"]))
## theoretical standard errors
(se0t <- sqrt(2/80 * (cf0^2/4 + 2)))
(lambda <- coef(m, type = "Lambdapar"))
(se1t <- sqrt(2/80 * ((1 + lambda^2)*cf1^2 + (2 + lambda^2)*4)/(2*lambda^4 + 6*lambda^2+4)))

## check model assumptions
## fit conditional model (4) and check if h_y is linear, h_baseline
## monotone
if (CHECK) {
    m <- BoxCoxME(y ~ w + s(x), data = d)
    plot(m) ## transformation function linear? Not really.
    plot(smooth_terms(m)) ## s(x) monotone? Yes!
}

## fit more flexible model
m0 <- BoxCox(y ~ w, data = d)
m1 <- BoxCox(x ~ 1, data = d)
m <- mmlt(m0, m1, formula = ~ 1, data = d)

(cf1f <- coef(m)["y.w100"])
(ci1f <- confint(m)["y.w100",])


## ----CAOdata, echo = FALSE, message = FALSE-----------------------------------
load(system.file("rda", "Primary_endpoint_data.rda", package = "TH.data"))
## outcome
CAOsurv$ypT0ypN0 <- factor(CAOsurv$path_stad == "ypT0ypN0")
CAOsurv$strata <- with(CAOsurv, interaction(strat_n, strat_t)) 
rt <- table(CAOsurv$randarm)


## ----CAO-glm, echo = FALSE, message = FALSE, results = "hide"-----------------
## estimate marginal log-odds ratio
mg_w <- glm(ypT0ypN0 ~ randarm,
            data = CAOsurv, family = binomial())
(lOR <- coef(mg_w)["randarm5-FU + Oxaliplatin"])
## Wald confidence interval
(ci <- confint(glht(mg_w), calpha = univariate_calpha())$confint[2,-1])


## ----CAO-mmlt, echo = FALSE, cache = TRUE-------------------------------------
## marginal models
## for outcome (equivalent to binary logistic GLM)
mpCR <- Polr(ypT0ypN0 ~ randarm, data = CAOsurv, na.action = na.pass, method = "logistic")
## for clinical T: binary probit GLM
mT <- Polr(strat_t ~ 1, data = CAOsurv, na.action = na.pass, method = "probit")
## for clinical N: binary probit GLM
mN <- Polr(strat_n ~ 1, data = CAOsurv, na.action = na.pass, method = "probit")
## for distance to anal verge (dichotomised, binary probit GLM)
mentf <- Polr(bentf ~ 1, data = CAOsurv, na.action = na.pass, method = "probit")
## for age (continuous)
mage <- BoxCox(age ~ 1, data = CAOsurv, na.action = na.pass)
## for sex: binary probit GLM
msex <- Polr(geschlecht ~ 1, data = CAOsurv, na.action = na.pass, method = "probit")
## for ECOG status (ordered probit model) 
CAOsurv$ecog_o <- as.ordered(CAOsurv$ecog_b)
mecog <- Polr(ecog_o ~ 1, data = CAOsurv, na.action = na.pass, method = "probit")
## multivariate transformation model: uses mixed continuous-discrete
## likelihood (arXiv:2408.17346), integrates out variables with missings
m <- Mmlt(mT, mN, mentf, mage, msex, mecog, mpCR, 
          data = CAOsurv, args = list(type = "ghalton", M = 250))


## ----CAO-output, echo = FALSE, results = "hide"-------------------------------
## adjusted analysis
prm <- "ypT0ypN0.randarm5-FU + Oxaliplatin"
(lOR <- coef(m)[prm])
ci <- confint(glht(m, coef. = function(...) coef(..., fixed = FALSE)), calpha = univariate_calpha())$confint
(ci <- ci[prm,-1])


## ----CAO-correlation, echo = FALSE--------------------------------------------
## correlation with outcome
mr <- as.array(coef(m, type = "Cor"))["ypT0ypN0",,1]
i <- which.max(abs(mr[-length(mr)]))
ni <- names(mr)[i]
mr <- mr[i]


## ----flies, echo = FALSE, results = "hide"------------------------------------
library("Stat2Data")
data("FruitFlies", package = "Stat2Data")

## 2 groups only
flies <- FruitFlies 
flies <- flies[flies$Treatment %in% c("8 virgin", "8 pregnant"),]
flies$Treatment <- flies$Treatment[, drop = TRUE]
flies$Longevity <- as.double(flies$Longevity)

## marginal model for thorax length
xmod <- BoxCox(Thorax ~ 1, data = flies)
## marginal Cox model for survival
flies$survival <- Surv(flies$Longevity)
coxph_w <- Coxph(survival ~ Treatment, data = flies)
## multivariate transformation model
m <- mmlt(xmod, coxph_w, data = flies, formula = ~ 1)

## marginal log-hazard ratio + Wald CI
(cf0 <- coef(coxph_w))
(ci0 <- confint(coxph_w))

## adjusted marginal log-hazard ratio + Wald CI
(cf1 <- coef(m)["survival.Treatment8 virgin"])
(ci1 <- confint(m)["survival.Treatment8 virgin",])

## conditional log-hazard ratio
(cf2 <- coef(Coxph(survival ~ Treatment + Thorax, data = flies))["Treatment8 virgin"])

## check marginal proportional hazards assumption
if (CHECK) {
    ## marginal survivor functions parallel on cloglog scale? Yes.
    plot(survfit(survival ~ Treatment, data = flies), fun = "cloglog")
    ## h_1(Thorax) monotone, in both groups? This is conditional model
    ## (10).
    m <- CoxphME(survival ~ s(Thorax, k = 5), data = flies, 
                 subset = Treatment == levels(Treatment)[1])
    plot(smooth_terms(m)) ## monotone? Yes.
    m <- CoxphME(survival ~ s(Thorax, k = 5), data = flies, 
                 subset = Treatment == levels(Treatment)[2])
    plot(smooth_terms(m)) ## monotone? Yes.
}

