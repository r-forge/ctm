
### catch constraint violations here
.log <- function(x) {
  ret <- log(pmax(.Machine$double.eps, x))
  dim(ret) <- dim(x)
  ret
}

### compute more sensible starting values
.start <- function(by, bx = NULL, data, ...) {
  
  J <- length(by)
  
  mctm <- vector(mode = "list", length = J)
  mmlt <- vector(mode = "list", length = J)
  mctm[[1]] <- by[[1]]$model
  mmlt[[1]] <- mlt(mctm[[1]], data = data, ...)
  pdat <- data
  htotal <- "~ 1"
  
  for (j in 2:J) {
    hhat <- paste("hhat", j - 1, sep = "_")
    htotal <- c(htotal, hhat)
    data[[hhat]] <- predict(mmlt[[j - 1]], newdata = pdat, 
                            type = "trafo")
    pdat[[hhat]] <- 0
    bhi <- as.basis(as.formula(paste(htotal, collapse = "+")), 
                    data = data, remove_intercept = TRUE)
    if (!is.null(bx)) {
      shift <- b(bh = bhi, bx = bx)
      if (!is.null(by[[j]]$model$bases$shifting))
        shift <- c(shift = by[[j]]$model$bases$shifting, bhbx = b(bh = bhi, bx = bx))
      mctm[[j]] <- ctm(by[[j]]$model$bases$response, 
                       interacting = by[[j]]$model$bases$interacting,
                       shifting = shift,
                       todistr = "Normal")
    } else {
      shift <- bhi
      if (!is.null(by[[j]]$model$bases$shifting))
        shift <- c(shift = by[[j]]$model$bases$shifting, bhbx = b(bh = bhi, bx = bx))
      mctm[[j]] <- ctm(by[[j]]$model$bases$response, 
                       interacting = by[[j]]$model$bases$interacting,
                       shifting = shift,
                       todistr = "Normal")
    }
    ### set todistr
    mctm[[j]]$todistr <- by[[j]]$todistr
    ### get marginal parameters as starting values
    theta <- coef(mctm[[j]])
    theta[] <- 0
    theta[names(coef(by[[j]]))] <- coef(by[[j]])
    mmlt[[j]] <- mlt(mctm[[j]], data = data, theta = theta, ...)
  }
  
  ### postprocess parameters
  p <- ncol(model.matrix(bx, data = data))
  cf <- lapply(mmlt, coef)
  mpar <- c()
  for (i in 1:length(cf))
    mpar <- c(mpar, cf[[i]][names(coef(by[[i]]))])
  cpar <- c()
  j <- 1
  for (i in 2:length(cf)) {
    cp <- cf[[i]][grep("hhat", names(cf[[i]]))]
    cpar <- rbind(cpar, matrix(cp, ncol = p))
  }
  list(mpar = mpar, cpar = cpar)
}


mmlt <- function(..., formula = ~ 1, data, 
                 theta = NULL, control.outer = list(trace = FALSE), scale = FALSE,
                 dofit = TRUE) {
  
  call <- match.call()
  
  m <- lapply(list(...), function(x) as.mlt(x))
  nm <- abbreviate(sapply(m, function(x) x$model$response), 4)
  J <- length(m)
  
  ### weights are not yet allowed
  w <- unique(do.call("c", lapply(m, weights)))
  stopifnot(isTRUE(all.equal(w, 1)))
  
  ### check if data is continuous and branch to discrete version here
  
  lu <- lapply(m, function(mod) {
    eY <- get("eY", environment(mod$parm))
    if (is.null(eY)) 
      stop("Only continuous outcomes without censoring implemented. For count outcomes
           consider using the mcotram function in the cotram package.")
    fixed <- get("fixed", environment(mod$parm))
    offset <- get("offset", environment(mod$parm))
    tmp <- attr(eY$Y, "constraint")
    wf <- !colnames(eY$Y) %in% names(fixed)
    eY$Y <- eY$Y[, wf,drop = FALSE]
    attr(eY$Y, "constraint") <- tmp
    list(exact = eY$Y, prime = eY$Yprime)
  })
  
  Jp <- J * (J - 1) / 2
  GAUSSIAN <- all(unlist(lapply(m, function(x) { x$todistr$name == "normal"})))

  bx <- formula
  if (inherits(formula, "formula"))
    bx <- as.basis(formula, data)
  lX <- model.matrix(bx, data = data)
  
  N <- nrow(lX)
  nobs <- sapply(lu, function(m) nrow(m$exact))
  stopifnot(length(unique(nobs)) == 1L)
  
  Y <- do.call("bdiag", lapply(lu, function(m) m$exact))
  Yprime <- do.call("bdiag", lapply(lu, function(m) m$prime))
  
  cnstr <- do.call("bdiag", 
                   lapply(lu, function(m) attr(m$exact, "constraint")$ui))
  ui <- bdiag(cnstr, Diagonal(Jp * ncol(lX)))
  ci <- do.call("c", lapply(lu, function(m) attr(m$exact, "constraint")$ci))
  ci <- c(ci, rep(-Inf, Jp * ncol(lX)))
  ui <- ui[is.finite(ci),]
  ci <- ci[is.finite(ci)]
  ui <- as(ui, "matrix")
  
  
  ll <- function(par) {
      
    mpar <- par[1:ncol(Y)]
    cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
      
    Yp <- matrix(Y %*% mpar, nrow = N)
    Yprimep <- matrix(Yprime %*% mpar, nrow = N)
    Xp <- ltmatrices(lX %*% cpar, byrow = TRUE, diag = FALSE, names = nm)

      
    if (GAUSSIAN) {
      C <- .mult(Xp, Yp)
      ret <- sum(.log(Yprimep))
      ret <- ret + sum(dnorm(C, log = TRUE))
    } else { ## !GAUSSIAN
        
      Sigmas2 <- .tcrossprod.ltmatrices(solve(Xp), diag_only = TRUE)
      Sigmas <- sqrt(Sigmas2)
        
      F_Zj_Yp <- Phi_01_inv <- Phi_Sigmas_inv <- Yp
        
      for (j in 1:J) {
        if (m[[j]]$todistr$name != "normal") {
          F_Zj_Yp[, j] <- m[[j]]$todistr$p(Yp[, j], log.p = TRUE)
          Phi_01_inv[, j] <- qnorm(F_Zj_Yp[, j], log.p = TRUE)
          Phi_Sigmas_inv[, j] <- Sigmas[, j] * Phi_01_inv[, j]
        }
      }
      C <- .mult(Xp, Phi_Sigmas_inv)
        
      ret <- sum(.log(Yprimep))
      for (j in 1:J) {
        ret <- ret + sum(dnorm(C[, j], log = TRUE))
        if (m[[j]]$todistr$name != "normal") {
          ret <- ret + 
            sum(m[[j]]$todistr$d(Yp[, j], log = TRUE)) +
            sum(.log(Sigmas[, j])) +
            0.5 * sum(Phi_01_inv[, j]^2)
        } 
      }
      ret <- ret - J * N * log(1 / sqrt(2 * pi))
    }
    return(-ret)
  }
   
  sc <- function(par) {
      
    mpar <- par[1:ncol(Y)] 
    cpar <- matrix(par[-(1:ncol(Y))], nrow = ncol(lX))
      
    Yp <- matrix(Y %*% mpar, nrow = N)
    Yprimep <- matrix(Yprime %*% mpar, nrow = N)
    Xp <- ltmatrices(lX %*% cpar, byrow = TRUE, diag = FALSE, names = nm)

    L <- diag(0, J)
    if (attr(Xp, "byrow")) {
        L[upper.tri(L)] <- 1:Jp
        L <- t(L)
    } else {
        L[lower.tri(L)] <- 1:Jp
    }
      
    if (GAUSSIAN) {
        
      C <- .mult(Xp, Yp)
      C1 <- -C
      B <- C[, -1L, drop = FALSE]
        
      mret <- vector(length = J, mode = "list")
      for (k in 1:J) {
        Lk <- L[,k]
        D <- cbind(matrix(rep(0, (k-1)*N), nrow = N), 1, unclass(Xp)[,Lk[Lk > 0]])
        mret[[k]] <- colSums(rowSums(C1 * D) * lu[[k]]$exact) +
                     colSums(lu[[k]]$prime / Yprimep[,k])
      }
        
      cret <- vector(length = J - 1, mode = "list")
      for (k in 1:(J - 1)) {  # go over rows
        B1 <- matrix(rep(B[,k], k), ncol = k)
        tmp <- -B1 * Yp[,1:k]
        ret <- c()
        for (i in 1:k) {
          tmp1 <- matrix(rep(tmp[,i], ncol(lX)), ncol = ncol(lX))
          ret <- c(ret, colSums(tmp1 * lX))
        }
        cret[[k]] <- ret
      }
    } else { ## !GAUSSIAN
        
      Sigmas2 <- .tcrossprod.ltmatrices(solve(Xp), diag_only = TRUE)
      Sigmas <- sqrt(Sigmas2)
        
      F_Zj_Yp <- Phi_01_inv <- Phi_Sigmas_inv <- Yp
        
      for (j in 1:J) {
        if (m[[j]]$todistr$name != "normal") {
          F_Zj_Yp[, j] <- m[[j]]$todistr$p(Yp[, j], log.p = TRUE)
          Phi_01_inv[, j] <- qnorm(F_Zj_Yp[, j], log.p = TRUE)
          Phi_Sigmas_inv[, j] <- Sigmas[, j] * Phi_01_inv[, j]
        }
      }
        
      C <- .mult(Xp, Phi_Sigmas_inv)
      C1 <- -C
      B <- C[, -1L, drop = FALSE]
        
      mret <- vector(length = J, mode = "list")
      for (k in 1:J) {
        Lk <- L[,k]
        D <- cbind(matrix(rep(0, (k-1)*N), nrow = N), 1, unclass(Xp)[,Lk[Lk > 0]])
          
        if(m[[k]]$todistr$name == "normal") {
          mret[[k]] <- colSums(rowSums(C1 * D) * lu[[k]]$exact) +
                       colSums(lu[[k]]$prime / Yprimep[,k])
        } else { ## k-th model not BoxCox
          f_k <- m[[k]]$todistr$d
          omega_k <- m[[k]]$todistr$dd2d
          mret[[k]] <- colSums(rowSums(C1 * D) * Sigmas[, k] * 
                                   f_k(Yp[, k]) * lu[[k]]$exact / dnorm(Phi_01_inv[, k])) +
              colSums(Phi_01_inv[, k] * f_k(Yp[, k]) * lu[[k]]$exact / dnorm(Phi_01_inv[, k])) +
              colSums(omega_k(Yp[, k]) * lu[[k]]$exact) +
              colSums(lu[[k]]$prime / Yprimep[, k])
        }
      }
        
      cret <- vector(length = J - 1, mode = "list")
      for (k in 1:(J - 1)) {  # go over rows
        B1 <- matrix(rep(B[,k], k), ncol = k)
        tmp1 <- - B1 * Phi_Sigmas_inv[,1:k]
        tmp2 <- - B1 * Phi_01_inv[,k+1]
        ret <- c()

        if(m[[k+1]]$todistr$name == "normal") {
          for (i in 1:k) {
            tmp3 <- matrix(rep(tmp1[,i], ncol(lX)), ncol = ncol(lX))
            ret <- c(ret, colSums(tmp3 * lX))
          }
        } else {
          Lk <- L[k+1, ]
          lambda_ktk <- unclass(Xp)[, Lk[Lk > 0]]
          tmp4 <- tmp1 + 
              (lambda_ktk / Sigmas[, k+1]) * tmp2  +
              lambda_ktk / Sigmas2[, k+1]
          for (i in 1:k) {
            tmp3 <- matrix(rep(tmp4[,i], ncol(lX)), ncol = ncol(lX))
            ret <- c(ret, colSums(tmp3 * lX))
          }
        }
        cret[[k]] <- ret
      }
    }
    mret <- -do.call("c", mret)
    cret <- -do.call("c", cret)
    return(c(mret, cret))
  }

  ### user-defined starting parameters for optimization
  if(!is.null(theta)) {
    start <- unname(theta)
  }
  else {
    if(inherits(formula, "formula") && formula == ~1) {
      ### don't bother with .start(), simply use the marginal coefficients
      ### and zero for the lambda parameters
      start <- do.call("c", lapply(m, function(mod) coef(as.mlt(mod))))
      start <- c(start, rep(0, Jp * ncol(lX)))
    }
    else { # formula != ~ 1
      start <- .start(m, bx = bx, data = data)
      start <- c(start$mpar, c(t(start$cpar)))
    }
  }
  
  if (scale) {
    Ytmp <- cbind(do.call("cbind", lapply(lu, function(m) m$exact)), 
                  kronecker(matrix(1, ncol = Jp), lX))
    Ytmp[!is.finite(Ytmp)] <- NA
    scl <- apply(abs(Ytmp), 2, max, na.rm = TRUE)
    lt1 <- scl < 1.1
    gt1 <- scl >= 1.1
    scl[gt1] <- 1 / scl[gt1]
    scl[lt1] <- 1
    start <- start / scl
    if (!is.null(ui))
      ui <- t(t(ui) * scl)
    f <- function(par) ll(scl * par)
    g <- function(par) sc(scl * par) * scl
  } else {
    f <- function(par) ll(par)
    g <- sc
  }
  
  if (!dofit)
    return(list(ll = ll, sc = sc))
  
  opt <- alabama::auglag(par = start, fn = f, gr = g,
                         hin = function(par) ui %*% par - ci,
                         hin.jac = function(par) ui,
                         control.outer = control.outer)[c("par",
                                                          "value",
                                                          "gradient",
                                                          "hessian")]
  
  if (scale) opt$par <- opt$par * scl
  
  mpar <- opt$par[1:(sum(sapply(lu, function(m) ncol(m$exact))))]
  
  mlist <- split(mpar, sf <- rep(factor(1:J), sapply(lu, function(m) ncol(m$exact))))
  mmod <- vector(mode = "list", length = J)
  for (j in 1:J) {
    mmod[[j]] <- as.mlt(m[[j]])
    coef(mmod[[j]]) <- mlist[[j]]
  }
  cpar <- matrix(opt$par[-(1:length(mpar))], ncol = Jp)
  tmp <- ltmatrices(cpar, byrow = TRUE, diag = FALSE, names = nm)
  args <- expand.grid(colnames(lX), colnames(unclass(tmp)))[,2:1]
  colnames(cpar) <- colnames(unclass(tmp))
  rownames(cpar) <- colnames(lX)
  args$sep <- "."
  names(opt$par) <- c(sapply(1:J, function(j) 
                             paste(nm[j], names(coef(mmod[[j]])), sep = ".")),
                      do.call("paste", args))
  
  ret <- list(marginals = mmod, formula = formula, bx = bx, data = data,
              call = call, diag = FALSE,
              gaussian = GAUSSIAN, 
              pars = list(mpar = mpar, cpar = cpar),
              par = opt$par, ll = ll, sc = sc, logLik = -opt$value,
              hessian = opt$hessian, names = nm)
  class(ret) <- "mmlt"
  ret
}

predict.mmlt <- function(object, newdata, margins = 1:J, 
                         type = c("trafo", "distribution", "density"), log = FALSE, ...) {

  type <- match.arg(type)

  if (type == "trafo")
    if (log) warning("argument log ignored")

  J <- length(object$marginals)
  yvar <- sapply(object$marginals, function(mg) mg$model$response)
  dx <- rep(1, J)
  names(dx) <- yvar

  link <- sapply(object$marginals, function(m) m$todistr$name)

  if (length(margins) == 1L) {

    ### Section 2.6: tilde{h} are already marginals for F_Z != Phi
    if (link[margins] != "normal")
        return(predict(object$marginals[[margins]], newdata = newdata, type = type, ...))

    ### lists currently not allowed
    stopifnot(is.data.frame(newdata)) 

    ### F_Z = Phi: need to rescale
    tr <- predict(object$marginals[[margins]], newdata = newdata, type = "trafo", ...)
    if (type == "trafo") return(tr)
    Vx <- coef(object, newdata = newdata, type = "Sigma")
    sdg <- matrix(sqrt(diagonals(Vx))[, margins], nrow = NROW(tr), ncol = NCOL(tr), byrow = TRUE)
    if (type == "distribution")
      return(pnorm(tr / sdg, log.p = log))
    trp <- predict(object$marginals[[margins]], newdata = newdata, type = "trafo", deriv = dx[margins], ...)
    ret <- dnorm(tr / sdg, log = TRUE) - .log(sdg) + .log(trp)
    if (log) return(ret)
    return(exp(ret))
  }

  tr <- do.call("cbind", lapply(margins, function(i)
                c(predict(object$marginals[[i]], newdata = newdata, type = "trafo"))))
  if (type == "trafo")
    return(tr)

  ret <- numeric(nrow(newdata))

  gm <- sapply(margins, function(i) object$marginals[[i]]$todistr$name == "normal")
  if (any(!gm)) {
    logF <- do.call("cbind", lapply(margins[!gm], function(i)
            c(predict(object$marginals[[i]], newdata = newdata, type = "distribution", log = TRUE))))
    Vx <- coef(object, newdata = newdata, type = "Sigma")
    sdg <- sqrt(diagonals(Vx))[, margins[!gm]]
    h <- tr
    h[,!gm] <- (Z <- qnorm(logF, log.p = TRUE)) * sdg
  } else {
    h <- tr
  }

  if (type == "distribution") {
    Smat <- as.array(coef(object, newdata = newdata, type = "Sigma")[, margins])
    for (i in 1:nrow(newdata))
      ret[i] <- pmvnorm(lower = rep(-Inf, length(margins)), upper = h[i,], sigma = Smat[,, i])
    if (log) return(.log(ret))
    return(ret)
  } else {
    if (1 %in% margins && all(diff(margins) == 1L)) {
      Lmat <- coef(object, newdata = newdata, type = "Lambda")[, margins]
    } else {
      if (length(margins) == 1L) {
        Vx <- coef(object, newdata = newdata, type = "Sigma")
        Lmat <- ltmatrices(sqrt(diagonals(Vx))[, margins, drop = FALSE], diag = TRUE)
      } else {
        stop("cannot evaluate density for selected margins; reorder and refit such that margins = 1:j")
      }
    }
    if (any(gm)) {
      trp <- do.call("cbind", lapply(margins[gm], function(i)
                     c(predict(object$marginals[[i]], newdata = newdata, type = "trafo", deriv = dx[i]))))
    } else {
      trp <- 1
    }
    ret <- rowSums(dnorm(.mult(Lmat, h), log = TRUE) + .log(trp))
      
    if (any(!gm)) {
      ld <- do.call("cbind", lapply(margins[!gm], function(i)
                    c(predict(object$marginals[[i]], newdata = newdata, 
                              type = "logdensity"))))
      ret <- ret + rowSums(ld) + rowSums(.log(sdg)) + .5 * rowSums(Z^2)
      ret <- ret - J * log(1 / sqrt(2 * pi))
    }

    if (log) return(ret)
    return(exp(ret))
  }
}

logLik.mmlt <- function(object, parm = coef(object), ...) {

  args <- list(...)
    if (length(args) > 0) 
      warning("Arguments ", names(args), " are ignored")

  ret <- -object$ll(parm)
  attr(ret, "df") <- length(object$par)
  class(ret) <- "logLik"
  ret
}

coef.mmlt <- function(object, newdata, 
                      type = c("all", "marginal", "Lambda", "Lambdainv", "Sigma", "Corr", "Spearman"), 
                      ...)
{
  
  type <- match.arg(type)
  if (type == "all") return(object$par)
  if (type == "marginal") return(lapply(object$marginals, coef))
  
  if (missing(newdata)) {
      if (nrow(object$pars$cpar) > 1L)
          stop("newdata not specified")
      ret <- ltmatrices(object$pars$cpar, byrow = TRUE, diag = FALSE, names = object$names)
  } else {
      X <- model.matrix(object$bx, data = newdata)
      ret <- ltmatrices(X %*% object$pars$cpar, byrow = TRUE, diag = FALSE, names = object$names)
  }

  if (type == "Spearman")
    return(6 * asin(coef(object, newdata = newdata, type = "Cor") / 2) / pi)

  ret <- switch(type, "Lambda" = ret,
                      "Lambdainv" = solve(ret),
                      "Sigma" = .tcrossprod.ltmatrices(solve(ret)),
                      "Corr" = {
                        inv <- solve(ret)
                        ret <- .tcrossprod.ltmatrices(inv)
                        isd <- 1 / sqrt(.tcrossprod.ltmatrices(inv, diag_only = TRUE))
                        J <- length(object$marginals)

                        if (attr(ret, "diag")) {
                            ### remove diagonal elements from ret
                            L <- diag(0, J)
                            L[upper.tri(L, diag = TRUE)] <- 1:ncol(unclass(ret))
                            L <- t(L)
                            ret <- unclass(ret)[, -diag(L), drop = FALSE]
                        } else {
                            ret <- unclass(ret)
                        }

                        L1 <- matrix(1:J, nrow = J, ncol = J)
                        L2 <- matrix(1:J, nrow = J, ncol = J, byrow = TRUE)
                        tmp <- ltmatrices(isd[, L2[lower.tri(L2)], drop = FALSE] * 
                                          isd[, L1[lower.tri(L1)], drop = FALSE], byrow = FALSE, diag = FALSE)
                        ret <- ret * unclass(.reorder(tmp, byrow = TRUE))
                        ret <- ltmatrices(ret, byrow = TRUE, diag = FALSE, names = object$names)
                        class(ret)[1L] <- "symatrices"
                        ret
                      })
  return(ret)
}

vcov.mmlt <- function(object, ...) {
  step <- 0
  lam <- 1e-6
  H <- object$hessian
  while((step <- step + 1) <= 3) {
        ret <- try(solve(H + (step - 1) * lam * diag(nrow(H))))
        if (!inherits(ret, "try-error")) break
  }
  if (inherits(ret, "try-error"))
      stop("Hessian is not invertible")
  if (step > 1)
      warning("Hessian is not invertible, an approximation is used")
  rownames(ret) <- colnames(ret) <- names(coef(object))
  ret
}

summary.mmlt <- function(object, ...) {
  ret <- list(call = object$call,
              #                tram = object$tram,
              test = cftest(object, parm = names(coef(object, with_baseline = FALSE))),
              ll = logLik(object))
  class(ret) <- "summary.mmlt"
  ret
}

print.summary.mmlt <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\n", "Multivariate conditional transformation model", "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  pq <- x$test$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  colnames(mtests) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  sig <- .Machine$double.eps
  printCoefmat(mtests, digits = digits, has.Pvalue = TRUE, 
               P.values = TRUE, eps.Pvalue = sig)
  cat("\nLog-Likelihood:\n ", x$ll, " (df = ", attr(x$ll, "df"), ")", sep = "")
  cat("\n\n")
  invisible(x)
}

print.mmlt <- function(x, ...) {
  cat("\n", "Multivariate conditional transformation model", "\n")
  cat("\nCall:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(coef(x))
  if (x$diag) {
     cat("\nDiagonal:\n", "elements are estimated.\n")
  } else { 
    cat("\nDiagonal:\n", "elements are constrained to 1.\n")
    }
  invisible(x)
}
