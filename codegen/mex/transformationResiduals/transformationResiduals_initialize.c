/*
 * transformationResiduals_initialize.c
 *
 * Code generation for function 'transformationResiduals_initialize'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "transformationResiduals_initialize.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar;

/* Function Definitions */
void transformationResiduals_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (transformationResiduals_initialize.c) */
