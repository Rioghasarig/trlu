/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * adderRef.h
 *
 * Code generation for function 'adderRef'
 *
 */

#ifndef ADDERREF_H
#define ADDERREF_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "adderRef_types.h"

/* Function Declarations */
extern void adderRef(const emlrtStack *sp, int32_T m, int32_T n, int32_T blkMax,
                     real_T in1[100], int32_T lda, int32_T jpvt[10], real_T in2
                     [100], int32_T info);

#endif

/* End of code generation (adderRef.h) */
