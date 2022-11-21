
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h> /* for dtptri */

SEXP R_ltmatrices_solve (SEXP x, SEXP N, SEXP J, SEXP diag)
{

    SEXP ans;
    double *dx, *dans, *ptr;
    int n, p, info, k, len, j, jj, idx;

    n = INTEGER(N)[0];
    p = LENGTH(x) / n;
    dx = REAL(x);
    int idiag = LOGICAL(diag)[0];
    char di, lo = 'L';
    
    if (idiag) {
        di = 'N';
    } else {
        di = 'U';
    }

    len = n * (p + (1 - idiag) * INTEGER(J)[0]);
    PROTECT(ans = allocVector(REALSXP, len));
    dans = REAL(ans);
    
    for (int i = 0; i < n; i++) {
    
        ptr = dans + i * (p + (1 - idiag) * INTEGER(J)[0]);
        jj = 0;
        k = 0;
        idx = 0;
        j = 0;
        while(j < p) {
            if (idiag == 0 && (jj == idx)) {
                ptr[jj] = 1.0;
                idx = idx + (INTEGER(J)[0] - k);
                k++;
            } else {
                ptr[jj] = dx[j * n + i];
                j++;
            }
            jj++;
        }
        if (idiag == 0) ptr[idx] = 1.0;

        F77_CALL(dtptri)(&lo, &di, INTEGER(J), ptr, &info FCONE FCONE);

        if (info != 0)
            error("Cannot solve ltmatices");
    }
    UNPROTECT(1);
    return(ans);
}
