# Demo ordinal

set.seed(24101968)
library(tramvs)

data("wine", package = "ordinal")
wine$noise <- rnorm(nrow(wine))

# Fixed support size
abess_tram(rating ~ temp + contact + noise, data = wine, modFUN = Polr,
           supp = 2)

# Estimate support size via HBIC
res <- tramvs(rating ~ temp + contact + noise, data = wine, modFUN = Polr)
plot(res, type = "b")
plot(res, which = "path")

coef(res)
coef(res, with_baseline = TRUE)

# Active set
support(res)
