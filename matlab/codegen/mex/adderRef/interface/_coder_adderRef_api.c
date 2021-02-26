/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_adderRef_api.c
 *
 * Code generation for function '_coder_adderRef_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "adderRef.h"
#include "_coder_adderRef_api.h"
#include "adderRef_data.h"

/* Function Declarations */
static int32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *in1,
  const char_T *identifier))[100];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[100];
static int32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *jpvt,
  const char_T *identifier))[10];
static int32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *m, const
  char_T *identifier);
static void emlrt_marshallOut(const real_T u[100], const mxArray *y);
static int32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[10];
static int32_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100];
static int32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10];

/* Function Definitions */
static int32_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  int32_T y;
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *in1,
  const char_T *identifier))[100]
{
  real_T (*y)[100];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(in1), &thisId);
  emlrtDestroyArray(&in1);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[100]
{
  real_T (*y)[100];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static int32_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *jpvt,
  const char_T *identifier))[10]
{
  int32_T (*y)[10];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(jpvt), &thisId);
  emlrtDestroyArray(&jpvt);
  return y;
}
  static int32_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *m, const
  char_T *identifier)
{
  int32_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(m), &thisId);
  emlrtDestroyArray(&m);
  return y;
}

static void emlrt_marshallOut(const real_T u[100], const mxArray *y)
{
  static const int32_T iv1[2] = { 10, 10 };

  emlrtMxSetData((mxArray *)y, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)y, iv1, 2);
}

static int32_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[10]
{
  int32_T (*y)[10];
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static int32_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  int32_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "int32", false, 0U, &dims);
  ret = *(int32_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[100]
{
  real_T (*ret)[100];
  static const int32_T dims[2] = { 10, 10 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[100])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static int32_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[10]
{
  int32_T (*ret)[10];
  static const int32_T dims[1] = { 10 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "int32", false, 1U, dims);
  ret = (int32_T (*)[10])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void adderRef_api(const mxArray *prhs[8], const mxArray *plhs[1])
{
  int32_T m;
  int32_T n;
  int32_T blkMax;
  real_T (*in1)[100];
  int32_T lda;
  int32_T (*jpvt)[10];
  real_T (*in2)[100];
  int32_T info;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  prhs[3] = emlrtProtectR2012b(prhs[3], 3, true, -1);
  prhs[5] = emlrtProtectR2012b(prhs[5], 5, false, -1);
  prhs[6] = emlrtProtectR2012b(prhs[6], 6, false, -1);

  /* Marshall function inputs */
  m = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "m");
  n = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "n");
  blkMax = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "blkMax");
  in1 = c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "in1");
  lda = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "lda");
  jpvt = e_emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "jpvt");
  in2 = c_emlrt_marshallIn(&st, emlrtAlias(prhs[6]), "in2");
  info = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "info");

  /* Invoke the target function */
  adderRef(&st, m, n, blkMax, *in1, lda, *jpvt, *in2, info);

  /* Marshall function outputs */
  emlrt_marshallOut(*in1, prhs[3]);
  plhs[0] = prhs[3];
}

/* End of code generation (_coder_adderRef_api.c) */
