/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * carrqr.h
 *
 * Code generation for function 'carrqr'
 *
 */

#ifndef CARRQR_H
#define CARRQR_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "carrqr_types.h"

/* Function Declarations */
extern void carrqr(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax,
                   const real_T A[100], int32_T lda, const int32_T jpvt[10],
                   const real_T tau[10], int32_T info);

#endif

/* End of code generation (carrqr.h) */
