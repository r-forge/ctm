
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

    /* number of matrices */
    n = INTEGER(N)[0];
    /* p = J * (J - 1) / 2 + diag * J */
    p = LENGTH(A) / n;
    dA = REAL(A);
    int idiag = LOGICAL(diag)[0];
    char di, lo = 'L', tr = 'N';
    
    if (idiag) {
        /* non-unit diagonal elements */
        di = 'N';
    } else {
        /* unit diagonal elements */
        di = 'U';
    }

    /* return object: include unit diagnonal elements if idiag == 0 */
    len = n * (p + (1 - idiag) * INTEGER(J)[0]);
    PROTECT(ans = allocVector(REALSXP, len));
    dans = REAL(ans);

    if (b != R_NilValue) {
        PROTECT(ansx = allocVector(REALSXP, n * INTEGER(J)[0]));
        dansx = REAL(ansx);
        db = REAL(b);
    }
    
    /* loop over matrices, ie rows of x */    
    for (int i = 0; i < n; i++) {

        /* point to subvector of ans */    
        ptrA = dans + i * (p + (1 - idiag) * INTEGER(J)[0]);

        /* copy data and insert unit diagonal elements when necessary */
        jj = 0;
        k = 0;
        idx = 0;
        j = 0;
        while(j < p) {
            if (idiag == 0 && (jj == idx)) {
                ptrA[jj] = 1.0;
                idx = idx + (INTEGER(J)[0] - k);
                k++;
            } else {
                ptrA[jj] = dA[j * n + i];
                j++;
            }
            jj++;
        }
        if (idiag == 0) ptrA[idx] = 1.0;

        if (b == R_NilValue) {
            /* compute inverse */
            F77_CALL(dtptri)(&lo, &di, INTEGER(J), ptrA, &info FCONE FCONE);
            if (info != 0)
                error("Cannot solve ltmatices");
        } else {
            /* solve linear system */
            ptrb = dansx + i * INTEGER(J)[0];
            for (j = 0; j < INTEGER(J)[0]; j++)
                ptrb[j] = db[i * INTEGER(J)[0] + j];
            F77_CALL(dtpsv)(&lo, &tr, &di, INTEGER(J), ptrA, ptrb, &ONE FCONE FCONE);
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
