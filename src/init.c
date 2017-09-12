#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP adfica(SEXP, SEXP, SEXP);
extern SEXP dgbn(SEXP, SEXP);
extern SEXP dgln(SEXP, SEXP);
extern SEXP dgpow3(SEXP);
extern SEXP dgrn(SEXP, SEXP);
extern SEXP ficadef(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ficasym(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gbn(SEXP, SEXP);
extern SEXP Gbn(SEXP, SEXP);
extern SEXP gln(SEXP, SEXP);
extern SEXP Gln(SEXP, SEXP);
extern SEXP gpow3(SEXP);
extern SEXP Gpow3(SEXP);
extern SEXP grn(SEXP, SEXP);
extern SEXP Grn(SEXP, SEXP);
extern SEXP relfica(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"adfica",  (DL_FUNC) &adfica,  3},
    {"dgbn",    (DL_FUNC) &dgbn,    2},
    {"dgln",    (DL_FUNC) &dgln,    2},
    {"dgpow3",  (DL_FUNC) &dgpow3,  1},
    {"dgrn",    (DL_FUNC) &dgrn,    2},
    {"ficadef", (DL_FUNC) &ficadef, 5},
    {"ficasym", (DL_FUNC) &ficasym, 5},
    {"gbn",     (DL_FUNC) &gbn,     2},
    {"Gbn",     (DL_FUNC) &Gbn,     2},
    {"gln",     (DL_FUNC) &gln,     2},
    {"Gln",     (DL_FUNC) &Gln,     2},
    {"gpow3",   (DL_FUNC) &gpow3,   1},
    {"Gpow3",   (DL_FUNC) &Gpow3,   1},
    {"grn",     (DL_FUNC) &grn,     2},
    {"Grn",     (DL_FUNC) &Grn,     2},
    {"relfica", (DL_FUNC) &relfica, 4},
    {NULL, NULL, 0}
};

void R_init_fICA(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
