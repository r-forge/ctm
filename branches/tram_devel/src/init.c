
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rconfig.h>
#include <R_ext/Lapack.h> /* for dspev */
#include <R_ext/Visibility.h>

/* .Call calls */
extern SEXP R_pnormMRS(SEXP);
extern SEXP R_inner(SEXP, SEXP);
extern SEXP R_ltmatrices_solve(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP R_ltmatrices_tcrossprod(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_pnormMRS", (DL_FUNC) &R_pnormMRS, 1},
    {"R_inner", (DL_FUNC) &R_inner, 2},
    {"R_ltmatrices_solve", (DL_FUNC) &R_ltmatrices_solve, 5},
    {"R_ltmatrices_tcrossprod", (DL_FUNC) &R_ltmatrices_tcrossprod, 5},
    {NULL, NULL, 0}
};

void attribute_visible R_init_tram(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
