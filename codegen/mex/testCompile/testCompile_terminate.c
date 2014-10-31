/*
 * testCompile_terminate.c
 *
 * Code generation for function 'testCompile_terminate'
 *
 * C source code generated on: Wed May 14 18:06:02 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "testCompile.h"
#include "testCompile_terminate.h"

/* Function Definitions */
void testCompile_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void testCompile_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (testCompile_terminate.c) */
