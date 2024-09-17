
###
### Fit models for HCC data presented in
### 
###	On Nonparanormal Likelihoods
###	by Torsten Hothorn
###

pkgs <- c("openxlsx", "tram", "survival")
req <- sapply(pkgs, require, char = TRUE)
if (!all(req)){
  sapply(pkgs[!req], install.packages)
  req[!req] <- sapply(pkgs[!req], require, char = TRUE)
}
if (!all(req)) 
  stop("cannot load dependencies")

### tram version 1.0-6 only
mmlt <- Mmlt

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
m <- mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC)
### marginal parameters
coef(m, type = "marginal")
### copula parameter: Lambda
coef(m, type = "Lambdapar")
### standard errors for all parameters
sqrt(diag(vcov(m)))

### convex approximations
## pseudo
mm <- mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC, domargins = FALSE)
## sequential
ms <- mmlt(mDKK, mOPN, mPIV, mAFP, data = HCC, sequentialfit = TRUE)
