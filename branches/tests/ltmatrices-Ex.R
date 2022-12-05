
library("mvtnorm")

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
a <- as.array(ltMatrices(xn, byrow = TRUE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = TRUE), byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = FALSE))
b <- as.array(ltMatrices(ltMatrices(xn, byrow = FALSE), byrow = TRUE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = TRUE, diag = TRUE), byrow = FALSE))
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))
b <- as.array(ltMatrices(ltMatrices(xd, byrow = FALSE, diag = TRUE), byrow = TRUE))
chk(a, b)

## subset
a <- as.array(ltMatrices(xn, byrow = FALSE)[1:2, 2:4])
b <- as.array(ltMatrices(xn, byrow = FALSE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xn, byrow = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xn, byrow = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xd, byrow = FALSE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

a <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE)[1:2, 2:4])
b <- as.array(ltMatrices(xd, byrow = TRUE, diag = TRUE))[2:4, 2:4, 1:2]
chk(a, b)

## solve
A <- as.array(lxn <- ltMatrices(xn, byrow = FALSE))
a <- solve(lxn)
a <- as.array(a)
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

A <- as.array(lxd <- ltMatrices(xd, byrow = FALSE, diag = TRUE))
a <- as.array(solve(lxd))
b <- array(apply(A, 3L, function(x) solve(x), simplify = TRUE), dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

## tcrossprod
a <- as.array(Tcrossprod(lxn))
b <- array(apply(as.array(lxn), 3L, function(x) tcrossprod(x), simplify = TRUE), dim = rev(dim(lxn)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxn, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxn)))

a <- as.array(Tcrossprod(lxd))
b <- array(apply(as.array(lxd), 3L, function(x) tcrossprod(x), simplify = TRUE), dim = rev(dim(lxd)))
chk(a, b, check.attributes = FALSE)

# diagonal elements only
d <- Tcrossprod(lxd, diag_only = TRUE)
chk(d, apply(a, 3, diag))
chk(d, diagonals(Tcrossprod(lxd)))

## multiplication
y <- matrix(runif(N * J), nrow = J)
a <- Mult(lxn, y)
A <- as.array(lxn)
b <- do.call("rbind", lapply(1:ncol(y), function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b))

a <- Mult(lxd, y)
A <- as.array(lxd)
b <- do.call("rbind", lapply(1:ncol(y), function(i) t(A[,,i] %*% y[,i,drop = FALSE])))
chk(a, t(b))

### tcrossprod as multiplication
i <- sample(1:N)[1]
M <- t(as.array(lxn)[,,i])
a <- sapply(1:J, function(j) Mult(lxn[i,], M[,j,drop = FALSE]))
rownames(a) <- colnames(a) <- dimnames(lxn)[[2L]]
b <- as.array(Tcrossprod(lxn[i,]))[,,1]
chk(a, b, check.attributes = FALSE)

head(unclass(ltMatrices(xn, byrow = TRUE)))
head(unclass(ltMatrices(xn, byrow = FALSE)))

head(unclass(ltMatrices(xd, byrow = TRUE, diag = TRUE)))
head(unclass(ltMatrices(xd, byrow = FALSE, diag = TRUE)))

