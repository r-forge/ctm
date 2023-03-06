# Demo ordinal

# Windows diffs...
old <- options(digits = 3)

set.seed(24101968)
library("tramvs")

if (require("ordinal")) {
  data("wine", package = "ordinal")
  wine$noise <- rnorm(nrow(wine))

  # Estimate support size via HBIC
  res <- tramvs(rating ~ temp + contact + noise, data = wine, modFUN = Polr)
  plot(res, type = "b")
  plot(res, which = "path")

  coef(res)
  coef(res, with_baseline = TRUE)

  # Active set
  support(res)

  options(old)
}
