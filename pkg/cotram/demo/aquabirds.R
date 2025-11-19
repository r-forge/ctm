
############################################################
#
#    Reproducibility Material for
#
#    Multi-species count transformation models
#
#    by Lukas Graz, Luisa Barbanti, Roland Brandl and Torsten Hothorn
#
#############################################################

if (!all(sapply(c("TH.data", "tram", "viridis"), require, character.only = TRUE)))
  stop("Please install the required packages: TH.data, tram, viridis")

### Load and prep data --------------------
###########################################
load(system.file("rda", "aquabirds.rda", package = "TH.data"))

# add day of Year
aquabirds$dayY <- as.numeric(format(aquabirds$Date, "%j"))
# add factor Year
aquabirds$Year <- as.factor(format(aquabirds$Date, "%Y"))
# Add cyclic basis for day of Year
cb <- cyclic_basis(numeric_var("dayY"), order = 4, frequency = 365)
XTag <- as.data.frame(model.matrix(cb, data = aquabirds))
colnames(XTag) <- paste0("t", 1:ncol(XTag))
tvars <- paste0("t", 1:ncol(XTag))
toy <- paste0(tvars, collapse = "+") #time of Year

aquabirds <- cbind(aquabirds, XTag)
rm(XTag, cb);

### Model Functions   --------------------
###########################################

# turns out that the marginal coefs can be extracted via:
#   m_m$models$parm(coef(m_m, fixed = TRUE))
get_Mcoefs <- function(m_m){
  Mmodel_coefs <- list()
  Lambda_coefs <- c()
  Model_coefs_ind <- list()
  for(mod in m_m$models$models){
    response <- mod$model$response
    Model_coefs_ind[[response]] <- match(
      paste0(response, ".", names(mod$theta)), 
      names(m_m$theta))
    tmp <- m_m$theta[Model_coefs_ind[[response]]]
    names(tmp) <- gsub(paste0("^", response, "\\."), "", names(tmp))
    Mmodel_coefs[[response]] <- tmp
  }
  Lambda_coefs <- m_m$theta[-unlist(Model_coefs_ind)]
  return(list(
    Mmodels = Mmodel_coefs,
    Lambdas = Lambda_coefs))
}


get_Mmodels <- function(m_m){
  CoefMMs <- get_Mcoefs(m_m)$Mmodels
  Mmodels <- list()
  for(i in seq_along(m_m$models$models)){
    mod <- m_m$models$models[[i]]
    call.  <- mod$call
    call.$theta <- CoefMMs[[i]]
    call.$dofit <- FALSE # keep marginally fitted coefs
    Mmodels[[mod$model$response]] <- eval(call.)
  }
  Mmodels
}


### Plotting Setups    --------------------
###########################################

cor.cols <- c("#1E88E5", "#DCA707", "#DE1760")

## newdata
cb <- cyclic_basis(numeric_var("dayY"), order = 4, frequency = 365)
# XTag <- as.data.frame(model.matrix(cb, data = aquabirds))
# colnames(XTag) <- paste0("t", 1:ncol(XTag))
nd <- data.frame(dayY = 1:365)
X <- model.matrix(cb, data = nd)
colnames(X) <- paste0("t", 1:ncol(X))
nd <- cbind(nd, as.data.frame(X))


# Dates for x-axis
ats <- c(min(aquabirds$Date), "2003-01-01", "2004-01-01",
         "2005-01-01", "2006-01-01", "2007-01-01", "2008-01-01",
         "2009-01-01", "2010-01-01", "2011-01-01", "2012-01-01",
         "2013-01-01", "2014-01-01", "2015-01-01", "2016-01-01",
         max(aquabirds$Date))
labs <- c("05.2002", " ", "2004", " ", "2006", " ", "2008", " ",
          "2010", " ", "2012", " ", "2014", " ", "2016", "11.2016")
xtick <- cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31))
xlabs <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
          "Oct", "Nov", "Dec", "Jan")
# colors with alpha:
qs <- c(seq(5, 95, by = 5))/100
qs <- c(0.01, qs, 0.99)
col1 <- viridis::plasma((length(qs) + 1)/2)
cols_quantiles <- c(col1, rev(col1[-length(col1)]))
alphas <- c(qs * I(qs <= 0.5) + (1-qs) * I(qs > 0.5))
alphas <- alphas + 0.08

s_rho <- function(x) {
  return(6 * asin(x / 2) / pi)
}

plot_corr.df <- function(df, newplot=TRUE, main=NULL, legend = TRUE,
              lwd=NULL, ylab = expression(paste("Spearman\'s  ", "",rho)), ...){
  if (is.null(lwd)) lwd <- 2.5

  df <- as.data.frame(df)
  stopifnot(c("dayY", "GreatCormorant.GreatCrestedGrebe", 
   "Goosander.GreatCrestedGrebe", "Goosander.GreatCormorant") %in% colnames(df))
  if(newplot){
    plot(df$dayY, df$GreatCormorant.GreatCrestedGrebe, ylim = c(-1, 1), 
      type = "l", xlab = "", ylab = ylab, las = 1, col = cor.cols[1], main=main,
      xaxt = "n", lwd = lwd, ...)
    axis(side = 1, at = xtick, labels = xlabs)
  } else {
    lines(df$dayY, df$GreatCormorant.GreatCrestedGrebe, col = cor.cols[1], 
      lwd = lwd, ...)
  }
  lines(df$dayY, df$Goosander.GreatCrestedGrebe, col = cor.cols[2], 
    lwd = lwd, ...)
  lines(df$dayY, df$Goosander.GreatCormorant, col = cor.cols[3], lwd = lwd, ...)
  abline(h = 0, lty = 2, col = "grey")
  if (newplot && legend){
    legend("bottomleft", col = cor.cols, lty = 1, lwd = lwd, 
        legend = c("Great Crested Grebe - Great Cormorant", 
        "Great Crested Grebe - Goosander", "Great Cormorant - Goosander"), 
        bty = "n")
  }
}


plot_corr.Mmlt <- function(Mmlt_obj, main = NULL, ...){
  stopifnot(inherits(Mmlt_obj, "Mmlt"))
  tmp <- cbind(nd, t(unclass(coef(Mmlt_obj, newdata = nd, type = "Spearman"))))

  # remove R. coding and ".." -> "." 
  if (any(grepl("^R\\.", colnames(tmp)))){
    colnames(tmp) <- gsub("..", ".", colnames(tmp), fixed=TRUE)
    colnames(tmp) <- gsub("R.", "", colnames(tmp), fixed=TRUE)
    colnames(tmp) <- gsub("\\.$", "", colnames(tmp))
  }
  
  # potentially Komoran is called Komoran1p. If so, rename collumns that have 
  # combinations (same for GreatCrestedGrebe1p and Goosander1p)
  colnames(tmp) <- gsub("GreatCormorant1p", "GreatCormorant", colnames(tmp))
  colnames(tmp) <- gsub("GreatCrestedGrebe1p","GreatCrestedGrebe",colnames(tmp))
  colnames(tmp) <- gsub("Goosander1p", "Goosander", colnames(tmp))

  plot_corr.df(tmp, main = main, ...)
}


#' Plot aquabirds density predictions from tram model
#' @param gather_years Logical, gather data from all years? reasonable if tram 
#'                      model has no year effect 
plot_aquabirds_denisty.tram <- function(model, gather_years = FALSE, ...){
  stopifnot(inherits(model, "tram"))
  # ylims <- c(0, 5.5) 

  # qs <- 1:19/20
  pred <- predict(as.mlt(model), type = "quantile", smooth = TRUE, prob = qs)
  # pred <- predict(as.mlt(model), type = "quantile", smooth = TRUE, prob=0.001)
  # if censored information, pred is a different object, 
  # convert back again to quantiles:
  if(inherits(pred, "response"))
    pred <- matrix(pred$approx, nrow = length(qs), 
                    dimnames = list(prob=as.character(qs), NULL))
  # pred <- pred$approxy
  if (!is.null(model$log_first) && model$log_first)
    pred <- pred - 1

  x_axis <- if (gather_years) aquabirds$dayY else aquabirds$Date

  # discrete count quantiles
  pred <- floor(pred)

  # data points
  plot(x_axis, model$response$approxy + 1, 
    # data = aquabirds, 
    log = "y",
    ylab = "counts + 1",
    xaxt = "n",
    xlab = " ",
    col = rgb(.1, .1, .1, .35),
    las = 1
    , ... )
  if (gather_years){
    axis(side = 1, at = xtick, labels = xlabs)
  } else {
    axis(side = 1, at = ats, labels = labs)
  }

  # quantile lines
  ind <- if(gather_years) 
    aquabirds$Year == "2003" else rep(TRUE, nrow(aquabirds))
  for(i in seq_along(qs)) {
    lines(x_axis[ind], pred[i, ind] + 1, 
      col = adjustcolor(cols_quantiles[i], alpha = alphas[i]*2.5), lwd = 1.3)
  }
}



plot_triple <- function(h, k, s, gather_years = FALSE){
  # get marginal models first, if Mmlt provided
  if (inherits(h, "mmlt")){
    m <- h
    MMs <- get_Mmodels(h)
    return(do.call(plot_triple, c(unname(MMs), gather_years = gather_years)))
  }

  # Set up layout with 4 rows, where the legend (4th row) is smaller
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 1)
  layout(layout_matrix, heights = c(1, 1, 1, 0.36))# Legend gets 0.36 rel-height
  
  # Set margins for the main plots
  mymai <- c(1, 1, 1, 1)
  # par(mai = par("mai") * mymai)
  par(mai = c(1.02, 0.82, 0.82, 0.42))
  par(mai = c(.3, .6, .3, .3))

  # Plot the three main plots
  plot_aquabirds_denisty.tram(h, main = "Great Crested Grebe", 
    gather_years=gather_years)
  plot_aquabirds_denisty.tram(k, main = "Great Cormorant", 
    gather_years=gather_years)
  plot_aquabirds_denisty.tram(s, main = "Goosander", 
    gather_years=gather_years)
  
  # Add the legend with smaller margins
  par(mai = c(0, 0.6, 0, 0.6))  # Smaller top/bottom margins for legend

  add_leg()
  
  # Reset layout
  layout(1)
}

# quantile-legend
add_leg <- function(x = TRUE) {
  xl <- 1
  yb <- 0.1
  xr <- 1.2
  yt <- 1.9
  
  n. <- length(qs)

  plot(NA, type = "n", ann = FALSE, 
       xlim = c(0, 2), ylim = c(0, 2), 
       xaxt = "n", yaxt = "n", bty = "n")
  rect(
    head(seq(yb, yt, (yt - yb)/n.), -1), xl,
    tail(seq(yb, yt, (yt - yb)/n.), -1), xr,
    col = diag(sapply(alphas, function(alpha) 
      adjustcolor(col = cols_quantiles, alpha = alpha*2.5)))
  )
  
  mtext(qs, side = 1, at = tail(seq(yb, yt, (yt - yb)/n.), -1) - 0.05, las = 1, 
        cex = 0.5, line = -2.5)
  mtext("Quantiles", side = 3, line = -2, cex = .7)
}


plot_aquabirds_data <- function(...){
  op <- par(mfrow = c(3, 1), mai = c(0.35, 0.6, 0.3, 0.1))
  plot(I(GreatCrestedGrebe + 1) ~ Date#, col = adjustcolor("#009E73", alpha = .3)
    , main = "Great Crested Grebe", xaxt = "n", ...)
  axis(side = 1, at = ats, labels = labs)
  abline(v = ats, lty=2, col="lightgray")
  plot(I(GreatCormorant + 1) ~ Date#, col = adjustcolor("#0072B2", alpha = .3)
    , main = "Great Cormorant", xaxt = "n", ...)
  axis(side = 1, at = ats, labels = labs)
  abline(v = ats, lty=2, col="lightgray")
  plot(I(Goosander + 1) ~ Date#, col = adjustcolor("#56B4E9", alpha = .3)
    , main = "Goosander", xaxt = "n", ...)
  axis(side = 1, at = ats, labels = labs)
  abline(v = ats, lty=2, col="lightgray")
  par(op)
}


# function to add confidence intervals to existing correlation plot
add_cor_CI <- function(samples, alpha=0.25, cols=cor.cols){
  stopifnot(
    class(samples) == "list",
    length(samples) == 3,
    all(sapply(samples, is.matrix)),
    all(sapply(samples, nrow) == 365)
  )
  
  for (i in 1:3){
    rho <- s_rho(samples[[i]])
    q_low <- apply(rho, MARGIN = 1, FUN = function(x){quantile(x, prob =0.025)})
    q_upp <- apply(rho, MARGIN = 1, FUN = function(x){quantile(x, prob =0.975)})
    polygon(c(1:365, 365:1), c(q_low, rev(q_upp)),
            col = adjustcolor(cor.cols[i], alpha = alpha), border = NA)
  }
}


# Wald/asymptotic parametric bootstrap confidence bands for correlations

### sampling nsamp values from the asymptotic distribution of the parameters
sample.parboot <- function(mod, nsamp = 200) {
  set.seed(1)
  H <- Hessian(mod)
  H <- (H + t(H))/2 # ensure symmetry

  V <- solve(H)
  P <- rmvnorm(nsamp, mean = coef(mod), sigma = (V + t(V))/2)
  m_tmp <- mod
  CR <- vector(mode = "list", length = nrow(P))
  
  # setup temp-model and retrieve corr via coef()
  for (i in 1:nsamp) {
    m_tmp$par[] <- P[i, ]
    CR[[i]] <- t(unclass(coef(m_tmp, newdata = nd, type = "Corr")))
  }
  
  # save results to matricies
  r21 <- r31 <- r32 <- matrix(NA, nrow = nrow(nd), ncol = nsamp)
  for(l in 1:nsamp) {
    r21[, l] <- CR[[l]][, 1]
    r31[, l] <- CR[[l]][, 2]
    r32[, l] <- CR[[l]][, 3]
  }
  
  return(list(r21, r31, r32
    # ,r21s = s_rho(r21), r31s = s_rho(r31), r32s = s_rho(r32)
    ))
}


get_kendall_df <- function(window_length){
  d <- 1:max(aquabirds$dayY)
  wd <- data.frame(dayY = d, GreatCormorant.GreatCrestedGrebe = NA, 
            Goosander.GreatCrestedGrebe = NA, Goosander.GreatCormorant = NA)

  for (i in d) {
      # tmp <- subset(aquabirds, dayY %in% i:(i + wl))
      tmp <- subset(aquabirds, 
        pmin((i - aquabirds$dayY) %% 365, 
        (aquabirds$dayY - i)  %% 365) <= window_length/2)
      CR <- suppressWarnings(cor(tmp[, c("GreatCrestedGrebe", "GreatCormorant", 
        "Goosander")], meth = "spearman", use = "pairwise.complete.obs"))
      wd[i,c("GreatCormorant.GreatCrestedGrebe", "Goosander.GreatCrestedGrebe", 
        "Goosander.GreatCormorant")] <- CR[lower.tri(CR)]
  }
  wd
}


## ----spearman-plot, echo = FALSE, warning = FALSE, message = FALSE, fig.width = 6, fig.height = 2.3----
op <- par(mai = c(1.02/2, 0.82, 0.82/5, 0.42/2))
get_kendall_df(window_length = 21) |> 
  plot_corr.df(ylab = expression(paste("Spearman\'s  ", rho)))
par(op)

## HMSC computation -----------------------------------------------------
message("to compute HMSC correlations, set 'executeHMSC <- TRUE'")
if (exists("executeHMSC") && executeHMSC){
  library("Hmsc")
  set.seed(270161)

  X <- model.matrix(as.formula(paste("~ Year + ", toy)), data = aquabirds)
  Y <- as.matrix(
    aquabirds[, c("GreatCrestedGrebe", "GreatCormorant", "Goosander")])

  ##########################################
  ## Hmsc with covariate-dependent Omega  
  ## extension by Tikhonov et al. (2017) is needed
  set.seed(270161)
  studyDesign <- data.frame(sample = as.factor(1:nrow(X)), 
                            dayY = as.factor(aquabirds$dayY))
  # rL1 <- HmscRandomLevel(units = levels(studyDesign$dayY))
  rL2 <- HmscRandomLevel(xData = data.frame(intercept = rep(1, nrow(aquabirds)),
                                            t1 = aquabirds$t1, 
                                            t2 = aquabirds$t2,
                                            t3 = aquabirds$t3,
                                            t4 = aquabirds$t4,
                                            t5 = aquabirds$t5,
                                            t6 = aquabirds$t6,
                                            t7 = aquabirds$t7,
                                            t8 = aquabirds$t8))
  m0 <- Hmsc(Y = Y, X = X, 
                distr = "poisson", studyDesign = studyDesign, 
                ranLevels = list(sample = rL2))#, dayY = rL1))
  # takes a long time
  m_var <- sampleMcmc(m0, thin = 1, samples = 2000, transient = 500,
                nChains = 8, verbose = 1, nParallel = 8)
  ## computeAssociations for complex model
  ## --> HMSC book: equation 7.12
  OmegaCor <- vector("list", m_var$nr)
  start <- thin <- 1
  postList <- poolMcmcChains(m_var$postList, start = start, thin = thin)
  # get lambda array from postlist
  Lambdas <- array(NA, dim = c(
                        dim(postList[[1]]$Lambda[[1]])[1],
                        dim(postList[[1]]$Lambda[[1]])[2],
                        dim(postList[[1]]$Lambda[[1]])[3], # = #covariates
                        length(postList)), # = #MCMC-samples
                        dimnames = list(
                          NULL,NULL,
                          covariate = names(rL2$x),
                          sample = paste0("MCMCsample", 1:length(postList))
  ))
  # Fill the array
  for (i in seq_along(postList)) {
    Lambdas[, , , i] <- postList[[i]]$Lambda[[1]]
  }
  # get an array of dim [1:2, 1:3, 1:nrow(xData), 1:8] 
  #   where we multiply each Lambda[i,j,,k] with its xData[,m]
  ind <- match(1:365, aquabirds$dayY)
  tmp <- as.matrix(rL2$x)[ind,] # design matrix for dayY
  LamX <- apply(Lambdas, c(1,2,4), function(a) c(sample.nr = tmp %*% a)) |> 
    aperm(c(2,3,1,4))
  Omega <- apply(LamX, 3:4, crossprod) # first dim (length 9) is 3x3 matrix
  Omega_cor <- apply(Omega, 2:3, \(x) cov2cor(matrix(x,3,3)))
  # average over MCMC samples
  Omega_cor_avg <- apply(Omega_cor, 1:2, mean) |> 
    t() |> _[, c(2,3,6)] # only unique pairs
  dimnames(Omega_cor_avg) <- list(
    dayY = 1:nrow(Omega_cor_avg),
    pair = c("GreatCormorant.GreatCrestedGrebe", 
            "Goosander.GreatCrestedGrebe", 
            "Goosander.GreatCormorant")
  )
  Omega_cor_avg <- cbind(Omega_cor_avg, dayY = 1:nrow(Omega_cor_avg))
  stopifnot(c("dayY", "GreatCormorant.GreatCrestedGrebe", 
    "Goosander.GreatCrestedGrebe", "Goosander.GreatCormorant") 
    %in% colnames(Omega_cor_avg))
  # saveRDS(Omega_cor_avg, file = "R/cache-manual/hmsc_omega_cor_avg.rds")
  subset_ind <- seq(1, dim(Omega_cor)[3], length.out=2000) |> round() |>unique()
  # get sample replicates for CIs
  HMSC_samples_for_CIs <- list(
    HK = Omega_cor[2,,subset_ind], # GreatCormorant - GreatCrestedGrebe
    HS = Omega_cor[3,,subset_ind], # Goosander - GreatCrestedGrebe
    KS = Omega_cor[6,,subset_ind]  # Goosander - GreatCormorant
  )
  # saveRDS(HMSC_samples_for_CIs,file="R/cache-manual/hmsc_samples_for_CIs.rds")
} else {
  message("Skipping HMSC computations")
}

## ----model_fits, cache = TRUE-------------------------------------------------
# first fit three models
m_h <- BoxCox(R(GreatCrestedGrebe)
         ~ t1+t2+t3+t4+t5+t6+t7+t8 | t1+t2+t3+t4+t5+t6+t7+t8, 
         data = aquabirds, LRtest = FALSE, na.action = na.pass)
m_k <- BoxCox(R(GreatCormorant)   
         ~ t1+t2+t3+t4+t5+t6+t7+t8 | t1+t2+t3+t4+t5+t6+t7+t8, 
         data = aquabirds, LRtest = FALSE, na.action = na.pass)
m_s <- BoxCox(R(Goosander)        
         ~ t1+t2+t3+t4+t5+t6+t7+t8 | t1+t2+t3+t4+t5+t6+t7+t8, 
         data = aquabirds, LRtest = FALSE, na.action = na.pass)

# then fit a multi- speciesmodel
m_m <- Mmlt(m_h, m_k, m_s, formula = ~t1+t2+t3+t4+t5+t6+t7+t8, data = aquabirds
  , args = list(M = 250, seed = 1, type = "ghalton")
)

Mods <- list(m_h = m_h, m_k = m_k, m_s = m_s, m_m = m_m)


## ----timeseries, echo = FALSE, cache = FALSE, message = FALSE, warning = FALSE, fig.width = 6, fig.height=6----
plot_aquabirds_data(data = aquabirds, log = "y", ylab = "Counts + 1", xlab = ""
   , pch = 1, cex = .6, las = 1, col=rgb(0, 0, 0, .25))


## ----marginal_dist, eval = TRUE, fig.width = 6, fig.height = 6----------------
plot_triple(Mods[["m_m"]], gather_years = TRUE) # plot mariganls here


## ----triple_cor_plot, eval = TRUE, fig.width = 6, fig.height = 6--------------
op <- par(mfrow = c(3,1), adj = 0, mai = c(1.02/3, 0.82/1.5, 0.82/3, 0.42/2))

df <- get_kendall_df(window_length = 21)
plot_corr.df(df, ylab = expression(paste("Spearman\'s  ", rho)), 
  main = "(a) Empirical trajectories")

plot_corr.Mmlt(Mods[["m_m"]], main="(b) MSCTM based trajectories", legend=FALSE)
add_cor_CI(sample.parboot(Mods[["m_m"]]), alpha = 0.15)
plot_corr.Mmlt(Mods[["m_m"]], newplot=FALSE)

## HMSC
## RDS generated from R/competitors/Hmsc.R

if (exists("executeHMSC") && executeHMSC){
  # Omega_cor_avg <- readRDS("R/cache-manual/hmsc_omega_cor_avg.rds"); message("read hmsc_omega_cor_avg.rds")
  # HMSC_samples_for_CIs <- readRDS("R/cache-manual/hmsc_samples_for_CIs.rds")
  plot_corr.df(Omega_cor_avg, ylab = "Latent residual correlation", 
    main = "(c) HMSC based trajectories", legend=FALSE)
  add_cor_CI(HMSC_samples_for_CIs)
  plot_corr.df(Omega_cor_avg, newplot = FALSE)
} else {
  plot.new()
  text(0.5, 0.5, "HMSC computations skipped", cex = 2)
}

par(op)
