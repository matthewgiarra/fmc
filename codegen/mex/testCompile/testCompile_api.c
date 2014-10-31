/*
 * testCompile_api.c
 *
 * Code generation for function 'testCompile_api'
 *
 * C source code generated on: Wed May 14 18:06:02 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "testCompile.h"
#include "testCompile_api.h"

/* Function Declarations */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId);
static real_T c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId);
static real_T emlrt_marshallIn(const mxArray *a, const char_T *identifier);
static const mxArray *emlrt_marshallOut(real_T u);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId)
{
  real_T ret;
  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 0U, 0);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const mxArray *a, const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = b_emlrt_marshallIn(emlrtAlias(a), &thisId);
  emlrtDestroyArray(&a);
  return y;
}

static const mxArray *emlrt_marshallOut(real_T u)
{
  const mxArray *y;
  const mxArray *m1;
  y = NULL;
  m1 = mxCreateDoubleScalar(u);
  emlrtAssign(&y, m1);
  return y;
}

const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  static const int32_T iv0[2] = { 0, 1 };

  const mxArray *m0;
  nameCaptureInfo = NULL;
  m0 = mxCreateNumericArray(2, (int32_T *)&iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtAssign(&nameCaptureInfo, m0);
  return nameCaptureInfo;
}

void testCompile_api(const mxArray * const prhs[2], const mxArray *plhs[1])
{
  real_T a;
  real_T b;

  /* Marshall function inputs */
  a = emlrt_marshallIn(emlrtAliasP(prhs[0]), "a");
  b = emlrt_marshallIn(emlrtAliasP(prhs[1]), "b");

  /* Invoke the target function */
  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(testCompile(a, b));
}

/* End of code generation (testCompile_api.c) */
