
### N (nrow(x)) lower triangular J x J matrices (ncol(x) = J * (J - 1) / 2 + diag * J)
### with unit (!diag) or non-unit (diag) diagonal in column order (!byrow) or row order (byrow)
### rcnames is a J-dim unique char vector
ltmatrices <- function(x, diag = FALSE, byrow = FALSE, names = TRUE) {

    if (!is.matrix(x)) x <- matrix(x, nrow = 1L)

    J <- floor((1 + sqrt(1 + 4 * 2 * ncol(x))) / 2 - diag)
    stopifnot(ncol(x) == J * (J - 1) / 2 + diag * J)

    if (!isTRUE(names)) {
        stopifnot(is.character(names) &&
                  length(unique(names)) == J)
    } else {
        names <- as.character(1:J)
    }

    L1 <- matrix(names, nrow = J, ncol = J)
    L2 <- matrix(names, nrow = J, ncol = J, byrow = TRUE)
    L <- matrix(paste(L1, L2, sep = "."), nrow = J, ncol = J)
    if (byrow)
        colnames(x) <- t(L)[upper.tri(L, diag = diag)]
    else
        colnames(x) <- L[lower.tri(L, diag = diag)]
         
    attr(x, "diag") <- diag
    attr(x, "byrow") <- byrow
    attr(x, "rcnames") <- names
    class(x) <- c("ltmatrices", class(x))
    x
}

### dimensions: N x J x J
dim.ltmatrices <- function(x) {
    J <- length(attr(x, "rcnames"))
    class(x) <- class(x)[-1L]
    return(c(nrow(x), J, J))
} 

### convert ltmatrices to J x J x N array for pretty printing
### symmetric returns the symmetric matrix where lower = upper
as.array.ltmatrices <- function(x, symmetric = FALSE, ...) {

    diag <- attr(x, "diag")
    byrow <- attr(x, "byrow")
    rcnames <- attr(x, "rcnames")
    class(x) <- class(x)[-1L]
    J <- length(rcnames)

    L <- matrix(1L, nrow = J, ncol = J)
    diag(L) <- 2L
    if (byrow) {
        L[upper.tri(L, diag = diag)] <- floor(2L + 1:(J * (J - 1) / 2L + diag * J))
        L <- t(L)
    } else {
        L[lower.tri(L, diag = diag)] <- floor(2L + 1:(J * (J - 1) / 2L + diag * J))
    }
    if (symmetric) {
        L[upper.tri(L)] <- 0L
        dg <- diag(L)
        L <- L + t(L)
        diag(L) <- dg
    }
    ret <- t(cbind(0, 1, x)[, c(L), drop = FALSE])
    class(ret) <- "array"
    dim(ret) <- c(J, J, nrow(x))
    dimnames(ret) <- list(rcnames, rcnames, rownames(x))
    return(ret)
}

as.array.symatrices <- function(x, ...)
    return(as.array.ltmatrices(x, symmetric = TRUE))

print.ltmatrices <- function(x, ...)
    print(as.array(x))

print.symatrices <- function(x, ...)
    print(as.array(x))

### change storage from column to row order or vice versa
.reorder <- function(x, byrow = FALSE) {
    
    stopifnot(inherits(x, "ltmatrices"))
    if (attr(x, "byrow") == byrow) return(x)

    diag <- attr(x, "diag")
    rcnames <- attr(x, "rcnames")
    J <- length(rcnames)
    class(x) <- class(x)[-1L]

    rL <- cL <- diag(0, nrow = J)
    rL[lower.tri(rL, diag = diag)] <- cL[upper.tri(cL, diag = diag)] <- 1:ncol(x)
    cL <- t(cL)
    if (attr(x, "byrow")) ### row -> col order
        return(ltmatrices(x[, cL[lower.tri(cL, diag = diag)], drop = FALSE], 
                          diag = diag, byrow = FALSE, names = rcnames))
    ### col -> row order
    return(ltmatrices(x[, t(rL)[upper.tri(rL, diag = diag)], drop = FALSE], 
                      diag = diag, byrow = TRUE, names = rcnames))
}

### subset: i selects rows out of 1:N and j columns out of 1:J
### returns ltmatrices in the same storage order as x
"[.ltmatrices" <- function(x, i, j, ..., drop = FALSE) {

    if (drop) warning("argument drop is ignored")
    if (missing(i) && missing(j)) return(x)
    diag <- attr(x, "diag")
    byrow <- attr(x, "byrow")
    rcnames <- attr(x, "rcnames")
    class(x) <- class(x)[-1L]
    J <- length(rcnames)
    if (!missing(j)) {
        if (length(j) == 1L && !diag)
            return(ltmatrices(matrix(1, nrow = nrow(x), ncol = 1), diag = TRUE, 
                              names = rcnames[j]))
        L <- diag(0L, nrow = J)
        if (byrow) {
            L[upper.tri(L, diag = diag)] <- 1:ncol(x)
            L <- L[j, j, drop = FALSE]
            L <- L[upper.tri(L, diag = diag)]
        } else {
            L[lower.tri(L, diag = diag)] <- 1:ncol(x)
            L <- L[j, j, drop = FALSE]
            L <- L[lower.tri(L, diag = diag)]
        }
        if (missing(i))
            return(ltmatrices(x[, c(L), drop = FALSE], diag = diag, 
                              byrow = byrow, names = rcnames[j]))
        return(ltmatrices(x[i, c(L), drop = FALSE], diag = diag, 
                          byrow = byrow, names = rcnames[j]))
    }
    return(ltmatrices(x[i, , drop = FALSE], diag = diag, 
                      byrow = byrow, names = rcnames))
}

"[.symatrices" <- function(x, i, j, ..., drop = FALSE) {
    class(x)[1L] <- "ltmatrices"
    ret <- x[i, j, ..., drop = drop]
    class(ret)[1L] <- "symatrices"
    return(ret)
}

### inverse of ltmatrices
### returns inverse matrices as ltmatrices in same storage order (missing b)
### or mult(solve(a), b)
solve.ltmatrices <- function(a, b, ...) {

    byrow_orig <- attr(a, "byrow")

    x <- .reorder(a, byrow = FALSE)
    diag <- attr(x, "diag")
    byrow <- attr(x, "byrow")
    rcnames <- attr(x, "rcnames")
    class(x) <- class(x)[-1L]
    xdim <- dim(x)
    J <- length(rcnames)
  
    if (diag) {
        L <- diag(0, nrow = J)
        L[!upper.tri(L)] <- 1:ncol(x)  ## column-wise
        idx_d <- diag(L)
        idx_l <- L[lower.tri(L)] ## column-wise
    
        x_orig <- x
        x_diag <- x[, idx_d, drop = FALSE]
        x <- x[, -idx_d, drop = FALSE]
    
        idx_f <- rep(1, J - 1L)
        if(J > 2) {
            for (j in 2:J)
                idx_f <- c(idx_f, rep(j, J - j))
        }
    
        x <- x / x_diag[, idx_f, drop = FALSE]
        if(J == 2) x <- t(x)
    }
  
    xij <- function(x = NULL, i, j) {
        if (i == j) return(1)
        if (j == 1) {
            ret <- i - 1
        } else {
            idx <- J - (1:(J - 1L))
            ret <- sum(idx[1:(j - 1)]) + (i - (J - idx[j]))
        }
        if (is.null(x))
            return(ret)
        return(x[, ret, drop = FALSE])
    }

    ret <- matrix(0, nrow = nrow(x), ncol(x))
    for (i in 2:J) {
        for (j in 1:(i - 1L)) {
            s <- 0L
            for (k in j:(i - 1L))
                s <- s + xij(x, i, k) * xij(ret, k, j)
            ret[, xij(NULL, i, j)] <- -s
        }
    }
    if (!diag) {
        ret <- .reorder(ltmatrices(ret, diag = diag, byrow = byrow, 
                                   names = rcnames), 
                        byrow = byrow_orig)
        if (missing(b)) return(ret)
        return(.mult(ret, b))
    }

    x_norm <- ret
    
    idx_f1 <- 2:J
    if(J > 2) {
        for (j in 3:J)
          idx_f1 <- c(idx_f1, j:J)
    }
    
    ret <- matrix(0, nrow = xdim[1L], ncol = xdim[2L])
    ret[, idx_l] <- x_norm / x_diag[, idx_f1]
    ret[, idx_d] <- 1 / x_diag
    
    ret <- .reorder(ltmatrices(ret, diag = diag, byrow = byrow, 
                               names = rcnames), 
                     byrow = byrow_orig)
    if (missing(b)) return(ret)
    return(.mult(ret, b))
} 

### L %*% t(L)
### NOTE: this returns symatrices (symmetric)
.tcrossprod.ltmatrices <- function(x, y = NULL, diag_only = FALSE) {

    stopifnot(is.null(y))
    byrow_orig <- attr(x, "byrow")
    rcnames <- attr(x, "rcnames")
    diag <- attr(x, "diag")
    J <- length(rcnames)

    x <- .reorder(x, byrow = FALSE)
    byrow <- attr(x, "byrow")
    class(x) <- class(x)[-1L]
    xdim <- dim(x)
    N <- xdim[1L]

    L <- diag(0, J)
    if (diag) {
        L[lower.tri(L, diag = diag)] <- 1:ncol(x)
    } else {
        L[lower.tri(L, diag = diag)] <- 1L + 1:ncol(x)
        diag(L) <- 1L
        x <- cbind(1, x)
    }
    tL <- t(L)

    if (diag_only) {
        ret <- matrix(0, nrow = nrow(x), ncol = J)
        colnames(ret) <- rcnames
        rownames(ret) <- rownames(x)
        k <- 0
        for (i in 1:J) {
            idx <- L[i,] * tL[,i] > 0
            k <- k + 1
            ret[, k] <- rowSums(x[, L[i,idx], drop = FALSE]^2)
        }
        return(ret)
    }

    ret <- matrix(0, nrow = nrow(x), ncol = J * (J - 1) / 2 + J)
    k <- 0
    for (i in 1:J) {
        for (j in 1:i) {
            idx <- L[i,] * tL[,j] > 0
            k <- k + 1
            ret[, k] <- rowSums(x[, L[i,idx], drop = FALSE] * x[,tL[idx,j], drop = FALSE])
        }
    }

    ret <- .reorder(ltmatrices(ret, diag = TRUE, byrow = TRUE, names = rcnames), 
                    byrow = byrow_orig)
    class(ret)[1L] <- "symatrices"
    ret
}

diagonals <- function(x, ...)
    UseMethod("diagonals")

diagonals.ltmatrices <- function(x, ...) {

    rcnames <- attr(x, "rcnames")
    diag <- attr(x, "diag")
    byrow <- attr(x, "byrow")
    diag <- attr(x, "diag")
    J <- length(rcnames)
    class(x) <- class(x)[-1L]

    if (!diag) {
        ret <- matrix(1, nrow = nrow(x), ncol = J)
        rownames(ret) <- rownames(x)
        colnames(ret) <- rcnames
        return(ret)
    } else {
        L <- diag(0, J)
        if (byrow) {
            L[upper.tri(L, diag = TRUE)] <- 1:ncol(x)
            L <- t(L)
            idx <- diag(L)
        } else {
            L[lower.tri(L, diag = TRUE)] <- 1:ncol(x)
            idx <- diag(L)
        }
        ret <- x[, idx, drop = FALSE]
        colnames(ret) <- rcnames
        return(ret)
    }
}

diagonals.symatrices <- diagonals.ltmatrices

### L %*% y
.mult <- function(x, y) {

    stopifnot(inherits(x, "ltmatrices"))

    rcnames <- attr(x, "rcnames")
    diag <- attr(x, "diag")
    idx <- attr(x, "idx")
    J <- length(rcnames)
    mx <- ifelse(J > 10, Matrix, matrix)
    x <- .reorder(x, byrow = TRUE)
    class(x) <- class(x)[-1L]

    if (!diag) {
        idx <- 1
        S <- 1
        if (J > 2) {
            S <- mx(rep(rep(1:0, (J - 1)), c(rbind(1:(J - 1), ncol(x)))), nrow = ncol(x))[, -J,drop = FALSE]
            idx <- unlist(lapply(colSums(S), seq_len))
        }
    } else {
        S <- mx(rep(rep(1:0, J),
                    c(rbind(1:J, ncol(x)))), nrow = ncol(x))[, -(J + 1), drop = FALSE]
        idx <- unlist(lapply(colSums(S), seq_len))
    }

    if (!diag) {
        A <- y[, idx] * x
        B <- A %*% S + y[, -1L, drop = FALSE]
        ret <- cbind(y[, 1L, drop = FALSE], as(B, "matrix"))
    } else {
	A <- y[, idx] * x
        ret <- as(A %*% S, "matrix")
    }
    colnames(ret) <- rcnames
    rownames(ret) <- rownames(x)
    return(ret)
}
