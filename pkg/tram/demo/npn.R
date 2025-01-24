
###
### Fit models for HCC data presented in
### 
###	On Nonparanormal Likelihoods
###	by Torsten Hothorn, UZH
###

set.seed(290875)

pkgs <- c("openxlsx", "tram", "survival")

ip <- rownames(installed.packages())
if (any(!pkgs %in% ip))
    install.packages(pkgs[!pkgs %in% ip], repos = "https://stat.ethz.ch/CRAN/")

OK <- sapply(pkgs, require, character.only = TRUE)
if (!all(OK)) 
    stop("package(s) ", paste(pkgs[!OK], collapse = ", "), " not available")

### Load data
dat <- read.xlsx("https://datadryad.org/api/v2/files/44697/download", sheet = 1)
HCC <- with(dat, data.frame(id = 1:nrow(dat),
                            x = factor(HCC_studyGr),
                            AFP = log(AFP_ng_per_ml),
                            PIV = log(PIVKA_delete_range),
                            OPN = log(OPN),
                            DKK = log(DKK)))
### limits of detection
PIVm <- max(HCC$PIV)
AFPm <- max(HCC$AFP)

### marginal location-scale models
mDKK <- BoxCox(
    DKK ~                               ### probit, h(DKK) via Bernstein
    x                                   ### location non-HCC / HCC
    | x,                                ### scale non-HCC / HCC
    data = HCC)
mOPN <- BoxCox(OPN ~ x | x, data = HCC)
HCC$PIVi <- with(HCC, R(       
    Surv(PIV, event = PIV < PIVm),      ### right censoring
    as.R.interval = TRUE))              ### empirical likelihood
mPIV <- BoxCox(PIVi ~
    x | x,                              ### location-scale
    data = HCC)
HCC$AFPi <- with(HCC, R(       
    Surv(AFP, event = AFP < AFPm),      ### right censoring
    as.R.interval = TRUE))              ### empirical likelihood
mAFP <- BoxCox(AFPi ~ x | x,  data = HCC)

### joint estimation of marginal and Gaussian copula parameters, s = 2
### location-scale transformation discriminant analysis
m <- Mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC)
logLik(m)
### marginal parameters
coef(m, type = "marginal")
### copula parameter: Lambda
coef(m, type = "Lambdapar")
### standard errors for all parameters
sqrt(diag(vcov(m)))

### convex approximations
## pseudo
mp <- Mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC, fit = "pseudo")
logLik(mp)
## sequential
ms <- Mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC, fit = "sequential")
logLik(ms)
## ACS
ma <- Mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC, fit = "ACS", ACSiter = 1)
logLik(ma)
