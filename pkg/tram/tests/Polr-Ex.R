
library("MASS")
library("tram")

### Windows diffs...
options(digits = 3)

tol <- .Machine$double.eps^(1/4)

cmp <- function(x, y)
    stopifnot(isTRUE(all.equal(x, y, tolerance = tol)))

(house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing))
summary(house.plr)

(house.plr2 <- Polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing))
summary(house.plr2)

cmp(coef(house.plr), coef(house.plr2))
ll <- logLik(house.plr)
attr(ll, "nobs") <- NULL
cmp(ll, logLik(house.plr2))

if (require("TH.data")) {

    ### blood loss data
    load(system.file("rda", "bloodloss.rda", package = "TH.data"))
    sMBL <- sort(unique(blood$MBL))
    blood$MBLc <- cut(blood$MBL, breaks = c(-Inf, sMBL), ordered_result = TRUE)

    op <- mltoptim()[c("auglag", "nlminb")]
    m <- lapply(1:length(op), function(i) 
        Polr(MBLc ~ 1, data = blood, method = "probit", optim = op[i])
    )
    stopifnot(all(diff(sapply(m, logLik)) < 1e-5))

    m <- lapply(1:length(op), function(i) 
        Polr(MBLc ~ IOL + DAUER.ap + FET.GEW, data = blood, method = "probit", optim = op[i])
    )
    stopifnot(all(diff(sapply(m, logLik)) < 1e-5))
}
