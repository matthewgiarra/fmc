/*
 * testCompile_initialize.c
 *
 * Code generation for function 'testCompile_initialize'
 *
 * C source code generated on: Wed May 14 18:06:02 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "testCompile.h"
#include "testCompile_initialize.h"

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar;

/* Function Definitions */
void testCompile_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (testCompile_initialize.c) */
