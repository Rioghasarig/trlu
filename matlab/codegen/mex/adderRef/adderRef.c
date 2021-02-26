/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * adderRef.c
 *
 * Code generation for function 'adderRef'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "adderRef.h"

/* Custom Source Code */
#include "cAdd.h"

/* Function Definitions */
void adderRef(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax, real_T
              in1[100], int32_T lda, int32_T jpvt[10], real_T in2[100], int32_T
              info)
{
  (void)sp;

  /*  the input numel(in1) is converted to integer type  */
  /*  to match the cAdd function signature */
  cAdd(m, n, blkMax, in1, lda, jpvt, in2, info);
}

/* End of code generation (adderRef.c) */
