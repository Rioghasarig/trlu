/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mcarrqr.h
 *
 * Code generation for function 'mcarrqr'
 *
 */

#ifndef MCARRQR_H
#define MCARRQR_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "mcarrqr_types.h"

/* Function Declarations */
extern void mcarrqr(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax,
                    real_T A[100], int32_T lda, int32_T jpvt[10], real_T tau[10],
                    int32_T info);

#endif

/* End of code generation (mcarrqr.h) */
