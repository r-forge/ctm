
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

static const R_CallMethodDef CallEntries[] = {
    {"R_pnormMRS", (DL_FUNC) &R_pnormMRS, 1},
    {"R_inner", (DL_FUNC) &R_inner, 2},
    {NULL, NULL, 0}
};

void attribute_visible R_init_tram(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
