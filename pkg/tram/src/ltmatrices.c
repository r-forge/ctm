
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h> /* for dtptri */

/*
    Solves A x = b 
        or
    computes A^{-1}
    
    for a set of lower triangular J x J matrices 
    (possibly with non-unit diagonals)
    
    A:	lower triangular J x J matrices in packed column-wise storage
        N x p (p = J * (J - 1) / 2 + diag * J)
    b:  missing or J x N
    N:  number of rows in A / columns of b
    J:  number of columns of A / rows of b
    diag: 1 if diag(A) != 1
    
    Returns A^{-1} WITH diagonal elements (missing b) 
        or 
    solutions x (J x N)
*/

SEXP R_ltmatrices_solve (SEXP A, SEXP b, SEXP N, SEXP J, SEXP diag)
{

    SEXP ans, ansx;
    double *dA, *dans, *dansx, *ptrA, *ptrb, *db;
    int n, p, info, k, len, j, jj, idx, ONE = 1;
    int iJ = INTEGER(J)[0];

    /* number of matrices */
    int iN = INTEGER(N)[0];
    /* p = J * (J - 1) / 2 + diag * J */
    p = LENGTH(A) / iN;
    dA = REAL(A);
    Rboolean Rdiag = asLogical(diag);
    char di, lo = 'L', tr = 'N';
    
    if (Rdiag) {
        /* non-unit diagonal elements */
        di = 'N';
    } else {
        /* unit diagonal elements */
        di = 'U';
    }

    /* return object: include unit diagnonal elements if Rdiag == 0 */
    len = iN * (p + (1 - Rdiag) * iJ);
    PROTECT(ans = allocVector(REALSXP, len));
    dans = REAL(ans);

    if (b != R_NilValue) {
        PROTECT(ansx = allocVector(REALSXP, iN * iJ));
        dansx = REAL(ansx);
        db = REAL(b);
    }
    
    /* loop over matrices, ie rows of x */    
    for (int n = 0; n < iN; n++) {

        /* point to subvector of ans */    
        ptrA = dans + n * (p + (1 - Rdiag) * iJ);

        /* copy data and insert unit diagonal elements when necessary */
        jj = 0;
        k = 0;
        idx = 0;
        j = 0;
        while(j < p) {
            if (!Rdiag && (jj == idx)) {
                ptrA[jj] = 1.0;
                idx = idx + (iJ - k);
                k++;
            } else {
                ptrA[jj] = dA[j * iN + n];
                j++;
            }
            jj++;
        }
        if (!Rdiag) ptrA[idx] = 1.0;

        if (b == R_NilValue) {
            /* compute inverse */
            F77_CALL(dtptri)(&lo, &di, &iJ, ptrA, &info FCONE FCONE);
            if (info != 0)
                error("Cannot solve ltmatices");
        } else {
            /* solve linear system */
            ptrb = dansx + n * iJ;
            for (j = 0; j < iJ; j++)
                ptrb[j] = db[n * iJ + j];
            F77_CALL(dtpsv)(&lo, &tr, &di, &iJ, ptrA, ptrb, &ONE FCONE FCONE);
        }
    }
    
    if (b == R_NilValue) {
        UNPROTECT(1);
        /* note: ans always includes diagonal elements */
        return(ans);
    } else {
        UNPROTECT(2);
        return(ansx);
    }
}

#define IDX(i, j, n, d) ((i) >= (j) ? (n) * ((j) - 1) - ((j) - 2) * ((j) - 1)/2 + (i) - (j) - (!d) * (j) : 0)

SEXP R_ltmatrices_tcrossprod (SEXP A, SEXP N, SEXP J, SEXP diag, SEXP diag_only) {

    SEXP ans;
    int iJ = INTEGER(J)[0];
    int iN = INTEGER(N)[0];
    int i, j, k, ix;
    double *dans, *dA;
    Rboolean Rdiag_only = asLogical(diag_only);
    Rboolean Rdiag = asLogical(diag);

    dA = REAL(A);
    
    if (Rdiag_only) {
        PROTECT(ans = allocVector(REALSXP, iN * iJ));
        dans = REAL(ans);
        for (int n = 0; n < iN; n++) {
            dans[n] = 1.0;
            if (Rdiag) dans[n] = pow(dA[n], 2);
            for (i = 1; i < iJ; i++) {
                ix = i * iN + n;
                dans[ix] = 0.0;
                for (k = 0; k < i; k++)
                    dans[ix] += pow(dA[IDX(i + 1, k + 1, iJ, Rdiag) * iN + n], 2);
                if (Rdiag) {
                    dans[ix] += pow(dA[IDX(i + 1, i + 1, iJ, Rdiag) * iN + n], 2);
                } else {
                    dans[ix] += 1.0;
                }
            }
        }
    } else {
        PROTECT(ans = allocVector(REALSXP, iN * iJ * (iJ + 1) / 2));
        dans = REAL(ans);
        for (int n = 0; n < INTEGER(N)[0]; n++) {
            dans[n] = 1.0;
            if (Rdiag) dans[n] = pow(dA[n], 2);
            for (i = 1; i < iJ; i++) {
                for (j = 0; j <= i; j++) {
                    ix =  IDX(i + 1, j + 1, iJ, 1) * iN + n;
                    dans[ix] = 0.0;
                    for (k = 0; k < j; k++)
                        dans[ix] += 
                            dA[IDX(i + 1, k + 1, iJ, Rdiag) * iN + n] *
                            dA[IDX(j + 1, k + 1, iJ, Rdiag) * iN + n];
                    if (Rdiag) {
                        dans[ix] += 
                            dA[IDX(i + 1, j + 1, iJ, Rdiag) * iN + n] *
                            dA[IDX(j + 1, j + 1, iJ, Rdiag) * iN + n];
                    } else {
                        if (j < i) {
                            dans[ix] += dA[IDX(i + 1, j + 1, iJ, Rdiag) * iN + n];
                        } else {
                            dans[ix] += 1.0;
                        }
                    }
                }
            }
        }
    }
    UNPROTECT(1);
    return(ans);
}
