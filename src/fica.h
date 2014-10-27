#include <RcppArmadillo.h>

RcppExport SEXP adfica(SEXP x, SEXP epsilon, SEXP maxiter);
RcppExport SEXP relfica(SEXP x, SEXP g, SEXP epsilon, SEXP maxiter);
RcppExport SEXP ficadef(SEXP x, SEXP g, SEXP w0, SEXP epsilon, SEXP maxiter);
RcppExport SEXP ficasym(SEXP x, SEXP g, SEXP w0, SEXP epsilon, SEXP maxiter);
RcppExport SEXP gpow3(SEXP x);
RcppExport SEXP dgpow3(SEXP x);
RcppExport SEXP Gpow3(SEXP x);
RcppExport SEXP grn(SEXP x, SEXP a1);
RcppExport SEXP dgrn(SEXP x, SEXP a1);
RcppExport SEXP Grn(SEXP x, SEXP a1);
RcppExport SEXP gln(SEXP x, SEXP a1);
RcppExport SEXP dgln(SEXP x, SEXP a1);
RcppExport SEXP Gln(SEXP x, SEXP a1);
RcppExport SEXP gbn(SEXP x, SEXP a1);
RcppExport SEXP dgbn(SEXP x, SEXP a1);
RcppExport SEXP Gbn(SEXP x, SEXP a1);
