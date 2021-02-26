/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_carrqr_mex.c
 *
 * Code generation for function '_coder_carrqr_mex'
 *
 */

/* Include files */
#include "carrqr.h"
#include "_coder_carrqr_mex.h"
#include "carrqr_terminate.h"
#include "_coder_carrqr_api.h"
#include "carrqr_initialize.h"
#include "carrqr_data.h"

/* Function Declarations */
static void carrqr_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs[8]);

/* Function Definitions */
static void carrqr_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs[8])
{
  int32_T n;
  const mxArray *inputs[8];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4, 6,
                        "carrqr");
  }

  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 6,
                        "carrqr");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  carrqr_api(inputs);

  /* Module termination. */
  carrqr_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  (void)plhs;
  mexAtExit(carrqr_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  carrqr_initialize();

  /* Dispatch the entry-point. */
  carrqr_mexFunction(nlhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_carrqr_mex.c) */
