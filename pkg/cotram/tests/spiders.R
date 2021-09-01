## example with new dataset
data("spiders")

## fit conditional marginal count transformation models
m_PF <- cotram(Pardosa_ferruginea ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")
m_HL <- cotram(Harpactea_lepida ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")
m_CC <- cotram(Callobius_claustrarius ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")
m_CT <- cotram(Coelotes_terrestris ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")
m_PL <- cotram(Pardosa_lugubris ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")
m_PR <- cotram(Pardosa_riparia ~ Elevation + Canopy_openess, data = spiders,
               method = "probit")

## fit multivariate count conditional transformation model
m_all_1 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                   formula = ~ 1, data = spiders)

m_all_2 <- mcotram(m_PF, m_HL, m_CC, m_CT, m_PL, m_PR, 
                   formula = ~ Elevation + Canopy_openess, data = spiders)

logLik(m_all_1)
logLik(m_all_2)

## lambda defining the Cholesky factor of the precision matrix,
## with standard error
coef(m_all_1, newdata = spiders[1,], type = "Lambda")
V <- vcov(m_all_1)[55:69, 55:69]
(se <- sqrt(diag(V)))

coef(m_all_2, newdata = spiders[1,], type = "Lambda")
V <- vcov(m_all_2)[55:99, 55:99]
(se <- sqrt(diag(V)))

## linear correlation, ie Pearson correlation of the models after
## transformation to bivariate normality
(r1 <- coef(m_all_1, newdata = spiders[1,], type = "Corr"))
(r2 <- coef(m_all_2, newdata = spiders[1,], type = "Corr"))