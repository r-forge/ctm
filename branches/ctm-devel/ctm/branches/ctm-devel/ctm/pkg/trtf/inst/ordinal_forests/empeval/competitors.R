
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages --------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library("trtf")
library("party")
library("partykit")
library("tram")
library("ranger")
library("ordinalForest")

library("psych")
# https://personality-project.org/r/psych/
# https://personality-project.org/r/psych-manual.pdf


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define global variables ----------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

set.seed(12345)

nsim <- c(100)
ntree <- c(250, 2000)

sim_para <- expand.grid(nsim = nsim, ntree = ntree)
rownames(sim_para) <- paste("nsim", sim_para$nsim, "_",
                            "ntree", sim_para$ntree, sep = "")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kullback-Leibler-Divergence ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

KL_div <- function(p,      # orig probability density function of DGP
                   y,      # orig y realisations of DGP test data set
                   q,      # estimated probability density function of random forest algorithm
                   yhat){  # predicted y observations based on random forest algorithm
  
  # Approach I: KL(p(y_k|x) || q(y_k|x))
  # estimated probability density function (with original y observations)
  q_y <- q[cbind(1:nrow(q), unclass(y))]
  KL_I <- mean( (-1) * log(q_y) )
  # library("MASS")
  # truehist((-1) * log(q_y))
  
  # Approach II: KL(q(\hat{y_k}|x) || p(\hat{y_k}|x))
  # model based predicted y observations
  p_yhat <- p[cbind(1:nrow(p), unclass(yhat))]
  KL_II <- mean( (-1) * log(p_yhat) )
  # truehist((-1) * log(p_yhat))
  
  # Approach III: KL(p(#|x) || q(#|x))
  # comparing whole probability density functions
  KL_III <- mean( rowSums(p * log(p / q)) )
  # i = 2
  # j = 3
  # p[i, j] * log(p[i, j] / q[i, j])
  
  return(list(KL_I = KL_I,
              KL_II = KL_II,
              KL_III = KL_III))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal Forest - perffunction = "equal" ------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ordforest_equal <- function(train, test, test_prb, mtry) {
  
  rf <- ordfor(depvar = "y", data = train,
               nsets = 1000,       # number of considered score sets
               ntreeperdiv = 100,  # number of trees considered per tried score
               ntreefinal = ntree, # number of trees for the forest
               perffunction = "equal",
               min.node.size = max(minbucket, nrow(train) / 2^nodedepth), 
               mtry = mtry, 
               keep.inbag = TRUE)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test)$classfreqtree
  rf_y <- factor(predict(rf, newdata = test)$ypred,
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  of_eq_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  # Kullback-Leibler Divergence
  of_eq_KL_div <- KL_div(p = test_prb, y = test$y,
                         q = rf_prb, yhat = rf_y)
  
  return(list(of_eq_ll = of_eq_ll,
              of_eq_KL_I = of_eq_KL_div$KL_I,
              of_eq_KL_II = of_eq_KL_div$KL_II,
              of_eq_KL_III = of_eq_KL_div$KL_III))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal Forest - perffunction = "proportional" -----------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ordforest_prop <- function(train, test, test_prb, mtry) {
  
  rf <- ordfor(depvar = "y", data = train,
               nsets = 1000,       # number of considered score sets
               ntreeperdiv = 100,  # number of trees considered per tried score
               ntreefinal = ntree, # number of trees for the forest
               perffunction = "proportional",
               min.node.size = max(minbucket, nrow(train) / 2^nodedepth), 
               mtry = mtry, 
               keep.inbag = TRUE)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test)$classfreqtree
  rf_y <- factor(predict(rf, newdata = test)$ypred,
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  of_prop_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  # Kullback-Leibler Divergence
  of_prop_KL_div <- KL_div(p = test_prb, y = test$y,
                           q = rf_prb, yhat = rf_y)
  
  return(list(of_prop_ll = of_prop_ll,
              of_prop_KL_I = of_prop_KL_div$KL_I,
              of_prop_KL_II = of_prop_KL_div$KL_II,
              of_prop_KL_III = of_prop_KL_div$KL_III))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal transfromation forest with Bernstein basis and general score -------
# Bs(theta) --> non-proportional odds deviations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

traforest_theta <- function(train, test, test_prb, mtry) {
  
  m0 <- as.mlt(Polr(y ~ 1, data = train))
  rf <- traforest(m0, 
                  formula = y ~ ., 
                  data = train,
                  ntree = ntree, 
                  trace = FALSE, 
                  mtry = mtry,
                  control = ctrl_partykit)
  
  # model-based probability density function and y predictions
  rf_prb <- do.call("rbind", lapply(predict(rf, newdata = test, type = "density"), c))
  rf_prb <- data.frame(rf_prb)
  colnames(rf_prb) <- levels(test$y)
  rf_y <- factor(predict.cforest(rf, newdata = test, type = "response"),
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  names(rf_y) <- NULL
  
  # log-likelihood with “nearest neighbour” forest weights prediction
  cf <- predict(rf, newdata = test, type = "coef")
  # Warning message:
  #   In model.frame.default(object$predictf, data = newdata,
  #   na.action = na.pass,  : variable 'y' is not a factor
  tf_theta_ll <- logLik(rf, newdata = test, coef = cf)
  
  # Kullback-Leibler Divergence
  tf_theta_KL_div <- KL_div(p = test_prb, y = test$y,
                            q = rf_prb, yhat = rf_y)
  
  return(list(tf_theta_ll = tf_theta_ll,
              tf_theta_KL_I = tf_theta_KL_div$KL_I,
              tf_theta_KL_II = tf_theta_KL_div$KL_II,
              tf_theta_KL_III = tf_theta_KL_div$KL_III))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Ordinal transformation forest with Bernstein basis and log-rank score ------
# Bs(alpha) --> proportional odds deviations from the model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

traforest_alpha <- function(train, test, test_prb, mtry) {
  
  ### There is no intercept in Bernstein-basis, therefor to perform
  ### log-rank splitting constant variable is added to the data implicitly
  train$alpha <- 1
  test$alpha <- 1
  m0 <- as.mlt(Polr(y ~ alpha, data = train, fixed = c("alpha" = 0)))
  ### split wrt to intercept ("log-rank scores") only
  rf <- traforest(m0, 
                  formula = y | alpha ~ ., 
                  data = train,
                  parm = "alpha",
                  mltargs = list(fixed = c("alpha" = 0)),
                  ntree = ntree, 
                  trace = FALSE, 
                  mtry = mtry,
                  control = ctrl_partykit)
  
  # model-based probability density function and y predictions
  rf_prb <- do.call("rbind", lapply(predict(rf, newdata = test,
                                            mnewdata = data.frame(alpha = 1),
                                            type = "density"), c))
  rf_prb <- data.frame(rf_prb)
  colnames(rf_prb) <- levels(test$y)
  rf_y <- factor(predict.cforest(rf, newdata = test, type = "response")[,"y"],
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood with “nearest neighbour” forest weights prediction
  cf <- predict(rf, newdata = test, type = "coef")
  tf_alpha_ll <- logLik(rf, newdata = test, coef = cf)
  
  # Kullback-Leibler Divergence
  tf_alpha_KL_div <- KL_div(p = test_prb, y = test$y,
                            q = rf_prb, yhat = rf_y)
  
  return(list(tf_alpha_ll = tf_alpha_ll,
              tf_alpha_KL_I = tf_alpha_KL_div$KL_I,
              tf_alpha_KL_II = tf_alpha_KL_div$KL_II,
              tf_alpha_KL_III = tf_alpha_KL_div$KL_III))
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Conditional Inference Trees --> Forests ------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mycforest <- function(train, test, test_prb, mtry) {
  
  ctrl_party <- party:::cforest_unbiased(ntree = ntree, 
                                         maxdepth = nodedepth, 
                                         mtry = mtry,
                                         minbucket = minbucket)
  
  ### log-rank splitting
  rf <- party::cforest(y ~ ., data = train, 
                       control = ctrl_party)
  
  # model-based probability density function and y predictions
  rf_prb <- predict(rf, newdata = test, type = "prob")
  rf_prb <- do.call("rbind", rf_prb)
  rf_y <- factor(predict(rf, newdata = test, type = "response"),
                 levels = levels(train$y),
                 labels = levels(train$y),
                 ordered = TRUE)
  
  # log-likelihood
  prb <- rf_prb[cbind(1:nrow(test), unclass(test$y))]
  cf_ll <- sum(log(pmax(sqrt(.Machine$double.eps), prb)))
  
  # Kullback-Leibler Divergence
  cf_KL_div <- KL_div(p = test_prb, y = test$y,
                      q = rf_prb, yhat = rf_y)
  
  return(list(cf_ll = cf_ll,
              cf_KL_I = cf_KL_div$KL_I,
              cf_KL_II = cf_KL_div$KL_II,
              cf_KL_III = cf_KL_div$KL_III))
  
}

