/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mcarrqr.c
 *
 * Code generation for function 'mcarrqr'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mcarrqr.h"

/* Custom Source Code */
#include "ccarrqr.h"

/* Function Definitions */
void mcarrqr(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax, real_T
             A[100], int32_T lda, int32_T jpvt[10], real_T tau[10], int32_T info)
{
  (void)sp;
  ccarrqr(m, n, blkMax, A, lda, jpvt, tau, info);
}

/* End of code generation (mcarrqr.c) */
