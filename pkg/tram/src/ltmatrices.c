
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

    /* number of matrices */
    n = INTEGER(N)[0];
    /* if (diag)
           p = J * (J + 1) / 2
       else 
           p = J * (J - 1) / 2
    */
    p = LENGTH(x) / n;
    dx = REAL(x);
    int idiag = LOGICAL(diag)[0];
    char di, lo = 'L';
    
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

    /* loop over matrices, ie rows of x */    
    for (int i = 0; i < n; i++) {

        /* point to subvector of ans */    
        ptr = dans + i * (p + (1 - idiag) * INTEGER(J)[0]);

        /* copy data and insert unit diagonal elements when necessary */
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

        /* compute inverse */
        F77_CALL(dtptri)(&lo, &di, INTEGER(J), ptr, &info FCONE FCONE);

        if (info != 0)
            error("Cannot solve ltmatices");
    }
    UNPROTECT(1);

    /* note: ans always includes diagonal elements */
    return(ans);
}
