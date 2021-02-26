/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * carrqr.c
 *
 * Code generation for function 'carrqr'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "carrqr.h"

/* Custom Source Code */
#include "ccarrqr.h"

/* Function Definitions */
void carrqr(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax, const
            real_T A[100], int32_T lda, const int32_T jpvt[10], const real_T
            tau[10], int32_T info)
{
  int32_T b_m;
  int32_T b_n;
  int32_T b_blkMax;
  real_T b_A[100];
  int32_T b_lda;
  int32_T b_jpvt[10];
  real_T b_tau[10];
  int32_T b_info;
  (void)sp;
  (void)m;
  (void)n;
  (void)blkMax;
  (void)A;
  (void)lda;
  (void)jpvt;
  (void)tau;
  (void)info;
  carrqr(&b_m, &b_n, &b_blkMax, b_A, &b_lda, b_jpvt, b_tau, &b_info);
}

/* End of code generation (carrqr.c) */
