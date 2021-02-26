/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_mcarrqr_mex.c
 *
 * Code generation for function '_coder_mcarrqr_mex'
 *
 */

/* Include files */
#include "mcarrqr.h"
#include "_coder_mcarrqr_mex.h"
#include "mcarrqr_terminate.h"
#include "_coder_mcarrqr_api.h"
#include "mcarrqr_initialize.h"
#include "mcarrqr_data.h"

/* Function Declarations */
static void mcarrqr_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs
  [8]);

/* Function Definitions */
static void mcarrqr_mexFunction(int32_T nlhs, int32_T nrhs, const mxArray *prhs
  [8])
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
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4, 7,
                        "mcarrqr");
  }

  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "mcarrqr");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  mcarrqr_api(inputs);

  /* Module termination. */
  mcarrqr_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  (void)plhs;
  mexAtExit(mcarrqr_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  mcarrqr_initialize();

  /* Dispatch the entry-point. */
  mcarrqr_mexFunction(nlhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_mcarrqr_mex.c) */
