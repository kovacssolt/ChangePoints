#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>

void F77_NAME(intervalF)(int *n2, double *dec, int *dep2, double *ilen, int *nint2, int *nsum2, int *bou);

void F77_NAME(findcptsF)(double *x, int *n2, double *dec, int *dep2, double *ilen, int *nint2, int *nsum2, int *par, int *stats2, int *pen2, double *thr, int *t, int *t2, double *tot, double *bou, int *bou2, double *cpts, int *cpts2);

void F77_NAME(narrowcptsF)(double *x, int *n2, double *dec, int *dep2, double *ilen, int *nint2, int *nsum2, int *par, double *thr, int *t, double *bou, int *bou2, int *cpts2);


extern SEXP intervalC(SEXP n, SEXP dec, SEXP dep, SEXP ilen, SEXP nint, SEXP nsum){
  SEXP bou;
  int *n2 = INTEGER(n);
  int *dep2 = INTEGER(dep);
  int *nint2 = INTEGER(nint);
  int *nsum2 = INTEGER(nsum);
  PROTECT(bou = allocMatrix(INTSXP, *nsum2, 2));
  F77_CALL(intervalF)(n2, REAL(dec), dep2, REAL(ilen), nint2, nsum2, INTEGER(bou));
  UNPROTECT(1);
  return(bou);
}

extern SEXP findcptsC(SEXP x, SEXP n, SEXP dec, SEXP dep, SEXP ilen, SEXP nint, SEXP nsum, SEXP par, SEXP stats, SEXP pen, SEXP thr){
  SEXP bou;
  SEXP bou2;
  SEXP cpts;
  SEXP cpts2;
  SEXP t;
  SEXP t2;
  SEXP tot;
  int *n2 = INTEGER(n);
  int *nsum2 = INTEGER(nsum);
  /* allocate and populate list */
  SEXP bou3 = PROTECT(allocVector(VECSXP, 7));
  PROTECT(bou = allocVector(REALSXP, *nsum2));
  PROTECT(bou2 = allocMatrix(INTSXP, *nsum2, 3));
  PROTECT(cpts = allocMatrix(REALSXP, *n2, 3));
  PROTECT(cpts2 = allocMatrix(INTSXP, *n2, 3));
  PROTECT(t = allocVector(INTSXP, 1));
  PROTECT(t2 = allocVector(INTSXP, 1));
  PROTECT(tot = allocVector(REALSXP, 1));
  F77_CALL(findcptsF)(REAL(x), n2, REAL(dec), INTEGER(dep), REAL(ilen), INTEGER(nint), nsum2, INTEGER(par), INTEGER(stats), INTEGER(pen), REAL(thr), INTEGER(t), INTEGER(t2), REAL(tot), REAL(bou), INTEGER(bou2), REAL(cpts), INTEGER(cpts2));
  SET_VECTOR_ELT(bou3, 0, bou);
  SET_VECTOR_ELT(bou3, 1, bou2);;
  SET_VECTOR_ELT(bou3, 2, cpts);
  SET_VECTOR_ELT(bou3, 3, cpts2);
  SET_VECTOR_ELT(bou3, 4, t);
  SET_VECTOR_ELT(bou3, 5, t2);
  SET_VECTOR_ELT(bou3, 6, tot);
  /* create names */
  SEXP nms = PROTECT(allocVector(STRSXP, 7));
  SET_STRING_ELT(nms, 0, mkChar("value"));
  SET_STRING_ELT(nms, 1, mkChar("int"));
  SET_STRING_ELT(nms, 2, mkChar("Cpts"));
  SET_STRING_ELT(nms, 3, mkChar("Cpts2"));
  SET_STRING_ELT(nms, 4, mkChar("t"));
  SET_STRING_ELT(nms, 5, mkChar("t2"));
  SET_STRING_ELT(nms, 6, mkChar("tot"));
  /* assign names to list */
  setAttrib(bou3, R_NamesSymbol, nms);
  UNPROTECT(9);
  return(bou3);
}


extern SEXP narrowcptsC(SEXP x, SEXP n, SEXP dec, SEXP dep, SEXP ilen, SEXP nint, SEXP nsum, SEXP par, SEXP thr){
  SEXP bou;
  SEXP bou2;
  SEXP cpts2;
  SEXP t;
  int *n2 = INTEGER(n);
  int *nsum2 = INTEGER(nsum);
  /* allocate and populate list */
  SEXP bou3 = PROTECT(allocVector(VECSXP, 4));
  PROTECT(bou = allocVector(REALSXP, *nsum2));
  PROTECT(bou2 = allocMatrix(INTSXP, *nsum2, 3));
  PROTECT(cpts2 = allocVector(INTSXP, *n2));
  PROTECT(t = allocVector(INTSXP, 1));
  F77_CALL(narrowcptsF)(REAL(x), n2, REAL(dec), INTEGER(dep), REAL(ilen), INTEGER(nint), nsum2, INTEGER(par), REAL(thr), INTEGER(t), REAL(bou), INTEGER(bou2), INTEGER(cpts2));
  SET_VECTOR_ELT(bou3, 0, bou);
  SET_VECTOR_ELT(bou3, 1, bou2);;
  SET_VECTOR_ELT(bou3, 2, cpts2);
  SET_VECTOR_ELT(bou3, 3, t);
  /* create names */
  SEXP nms = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, mkChar("value"));
  SET_STRING_ELT(nms, 1, mkChar("int"));
  SET_STRING_ELT(nms, 2, mkChar("Cpts2"));
  SET_STRING_ELT(nms, 3, mkChar("t"));
  /* assign names to list */
  setAttrib(bou3, R_NamesSymbol, nms);
  UNPROTECT(6);
  return(bou3);
}

static const R_CallMethodDef CallEntries[] = {
  {"intervalC",   (DL_FUNC) &intervalC,   6},
  {"findcptsC",   (DL_FUNC) &findcptsC,   11},
  {"narrowcptsC",   (DL_FUNC) &narrowcptsC,   9},
  {NULL, NULL, 0}
};


void R_init_ChangePoints(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}


