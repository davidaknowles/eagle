#pragma once
#ifndef _VBGLM_H
#define _VBGLM_H
#include <Rcpp.h>

using namespace Rcpp ;

RcppExport SEXP runvb(SEXP alt_sexp, SEXP n_sexp, SEXP x_sexp, SEXP settings); 

RcppExport SEXP fitnull(SEXP alt_sexp, SEXP n_sexp, SEXP max_its, SEXP tol, SEXP rev, SEXP debug); 

#endif
