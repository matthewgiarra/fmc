/*
 * transformationResiduals_terminate.c
 *
 * Code generation for function 'transformationResiduals_terminate'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "transformationResiduals_terminate.h"

/* Function Definitions */
void transformationResiduals_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void transformationResiduals_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (transformationResiduals_terminate.c) */
