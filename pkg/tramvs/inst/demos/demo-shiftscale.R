# Demo sstram

library(tramvs)
library(sstram)

data("BostonHousing2", package = "mlbench")

abess_tram(cmedv ~ nox + rm + age + lat | rm + tax,
           data = BostonHousing2,
           modFUN = BoxCoxS, supp = 3, stabilizer.sqrt = TRUE, scale_shift = TRUE)

res <- tramvs(cmedv ~ nox + rm + age + lat | rm + tax,
            data = BostonHousing2,
            modFUN = BoxCoxS, supp_max = 5)

plot(res, type = "b")
plot(res, which = "path")

# Active set
support(res)
