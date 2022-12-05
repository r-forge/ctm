
library("tram")

chk <- function(...) stopifnot(all.equal(...))

### some checks
## dimensions
set.seed(290875)
N <- 100
J <- 10
Jn <- J * (J - 1) / 2
## data
xn <- matrix(runif(N * Jn), nrow = N, byrow = TRUE)
xd <- matrix(runif(N * (Jn + J)), nrow = N, byrow = TRUE)

## constructor + .reorder + as.array
a <- as.array(ltmatrices(xn, byrow = TRUE))
b <- as.array(tram:::.reorder(ltmatrices(xn, byrow = TRUE), byrow = FALSE))
chk(a, b)

a <- as.array(ltmatrices(xn, byrow = FALSE))
b <- as.array(tram:::.reorder(ltmatrices(xn, byrow = FALSE), byrow = TRUE))
chk(a, b)

a <- as.array(ltmatrices(xd, byrow = TRUE, diag = TRUE))
b <- as.array(tram:::.reorder(ltmatrices(xd, byrow = TRUE, diag = TRUE), byrow = FALSE))
chk(a, b)

a <- as.array(ltmatrices(xd, byrow = FALSE, diag = TRUE))
b <- as.array(tram:::.reorder(ltmatrices(xd, byrow = FALSE, diag = TRUE), byrow = TRUE))
chk(a, b)

## subset
a <- as.array(ltmatrices(xn, byrow = FALSE)[1:2, 2:4])
b <- as.array(ltmatrices(xn, byrow = FALSE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltmatrices(xn, byrow = TRUE)[1:2, 2:4])
b <- as.array(ltmatrices(xn, byrow = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltmatrices(xd, byrow = FALSE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltmatrices(xd, byrow = FALSE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltmatrices(xd, byrow = TRUE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltmatrices(xd, byrow = TRUE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

## solve
A <- as.array(lxn <- ltmatrices(xn, byrow = FALSE))
a <- solve(lxn)
a <- as.array(a)
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

A <- as.array(lxd <- ltmatrices(xd, byrow = FALSE, diag = TRUE))
a <- as.array(solve(lxd))
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

## tcrossprod
a <- as.array(tram:::.tcrossprod.ltmatrices(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) tcrossprod(x), simplify = TRUE), dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- tram:::.tcrossprod.ltmatrices(lxn, diag_only = TRUE)
chk(d, t(apply(a, 3, diag)))
chk(d, diagonals(tram:::.tcrossprod.ltmatrices(lxn)))

a <- as.array(tram:::.tcrossprod.ltmatrices(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) tcrossprod(x), simplify = TRUE), dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- tram:::.tcrossprod.ltmatrices(lxd, diag_only = TRUE)
chk(d, t(apply(a, 3, diag)))
chk(d, diagonals(tram:::.tcrossprod.ltmatrices(lxd)))

## multiplication
y <- matrix(runif(N * J), nrow = N)
a <- tram:::.mult(lxn, y)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:nrow(y), function(i) t(A[,,i] %*% t(y[i,,drop = FALSE]))))
chk(a, b)

a <- tram:::.mult(lxd, y)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:nrow(y), function(i) t(A[,,i] %*% t(y[i,,drop = FALSE]))))
chk(a, b)

### tcrossprod as multiplication
i <- sample(1:N)[1]
M <- t(as.array(lxn)[,,i])
a <- sapply(1:J, function(j) tram:::.mult(lxn[i,], t(M[,j,drop = FALSE])))
rownames(a) <- colnames(a) <- attr(lxn, "rcnames")
b <- as.array(tram:::.tcrossprod.ltmatrices(lxn[i,]))[,,1]
chk(a, b)

head(unclass(ltmatrices(xn, byrow = TRUE)))
head(unclass(ltmatrices(xn, byrow = FALSE)))

head(unclass(ltmatrices(xd, byrow = TRUE, diag = TRUE)))
head(unclass(ltmatrices(xd, byrow = FALSE, diag = TRUE)))

