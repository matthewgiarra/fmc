/*
 * inv.c
 *
 * Code generation for function 'inv'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "inv.h"
#include "interp2.h"
#include "norm.h"
#include "transformationResiduals_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo c_emlrtRSI = { 27, "inv",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtRSInfo d_emlrtRSI = { 40, "inv",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtRSInfo e_emlrtRSI = { 44, "inv",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/matfun/inv.m" };

static emlrtRSInfo f_emlrtRSI = { 16, "eml_warning",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtRSInfo g_emlrtRSI = { 29, "eml_flt2str",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtMCInfo emlrtMCI = { 16, 13, "eml_warning",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtMCInfo b_emlrtMCI = { 16, 5, "eml_warning",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtMCInfo c_emlrtMCI = { 29, 23, "eml_flt2str",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtMCInfo d_emlrtMCI = { 29, 15, "eml_flt2str",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

/* Function Declarations */
static void b_eml_warning(const char_T varargin_2[14]);
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, char_T y[14]);
static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location);
static const mxArray *b_sprintf(const mxArray *b, const mxArray *c, const
  mxArray *d, emlrtMCInfo *location);
static const mxArray *c_sprintf(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location);
static void eml_warning(void);
static void emlrt_marshallIn(const mxArray *d_sprintf, const char_T *identifier,
  char_T y[14]);
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, char_T ret[14]);
static void warning(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void b_eml_warning(const char_T varargin_2[14])
{
  const mxArray *y;
  static const int32_T iv3[2] = { 1, 33 };

  const mxArray *m2;
  char_T cv4[33];
  int32_T i;
  static const char_T cv5[33] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'l', 'l', 'C', 'o', 'n', 'd', 'i', 't', 'i', 'o',
    'n', 'e', 'd', 'M', 'a', 't', 'r', 'i', 'x' };

  const mxArray *b_y;
  static const int32_T iv4[2] = { 1, 14 };

  char_T b_varargin_2[14];
  emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
  y = NULL;
  m2 = mxCreateCharArray(2, iv3);
  for (i = 0; i < 33; i++) {
    cv4[i] = cv5[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 33, m2, cv4);
  emlrtAssign(&y, m2);
  b_y = NULL;
  m2 = mxCreateCharArray(2, iv4);
  for (i = 0; i < 14; i++) {
    b_varargin_2[i] = varargin_2[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 14, m2, b_varargin_2);
  emlrtAssign(&b_y, m2);
  warning(b_message(y, b_y, &emlrtMCI), &b_emlrtMCI);
  emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
}

static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, char_T y[14])
{
  g_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_message(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location)
{
  const mxArray *pArrays[2];
  const mxArray *m9;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m9, 2, pArrays, "message",
    TRUE, location);
}

static const mxArray *b_sprintf(const mxArray *b, const mxArray *c, const
  mxArray *d, emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  const mxArray *m7;
  pArrays[0] = b;
  pArrays[1] = c;
  pArrays[2] = d;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m7, 3, pArrays, "sprintf",
    TRUE, location);
}

static const mxArray *c_sprintf(const mxArray *b, const mxArray *c, emlrtMCInfo *
  location)
{
  const mxArray *pArrays[2];
  const mxArray *m8;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 1, &m8, 2, pArrays, "sprintf",
    TRUE, location);
}

static void eml_warning(void)
{
  const mxArray *y;
  static const int32_T iv2[2] = { 1, 27 };

  const mxArray *m1;
  char_T cv2[27];
  int32_T i;
  static const char_T cv3[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a', 't',
    'r', 'i', 'x' };

  emlrtPushRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
  y = NULL;
  m1 = mxCreateCharArray(2, iv2);
  for (i = 0; i < 27; i++) {
    cv2[i] = cv3[i];
  }

  emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 27, m1, cv2);
  emlrtAssign(&y, m1);
  warning(message(y, &emlrtMCI), &b_emlrtMCI);
  emlrtPopRtStackR2012b(&f_emlrtRSI, emlrtRootTLSGlobal);
}

static void emlrt_marshallIn(const mxArray *d_sprintf, const char_T *identifier,
  char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(emlrtAlias(d_sprintf), &thisId, y);
  emlrtDestroyArray(&d_sprintf);
}

static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, char_T ret[14])
{
  int32_T iv9[2];
  int32_T i4;
  for (i4 = 0; i4 < 2; i4++) {
    iv9[i4] = 1 + 13 * i4;
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "char", FALSE, 2U, iv9);
  emlrtImportCharArray(src, ret, 14);
  emlrtDestroyArray(&src);
}

static void warning(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "warning", TRUE,
                        location);
}

void inv(const real_T x[9], real_T y[9])
{
  real_T b_x[9];
  int32_T p1;
  int32_T p2;
  int32_T p3;
  real_T absx11;
  real_T absx21;
  real_T absx31;
  int32_T itmp;
  const mxArray *b_y;
  static const int32_T iv1[2] = { 1, 8 };

  const mxArray *m0;
  char_T cv0[8];
  static const char_T cv1[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  char_T str[14];
  memcpy(&b_x[0], &x[0], 9U * sizeof(real_T));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = muDoubleScalarAbs(x[0]);
  absx21 = muDoubleScalarAbs(x[1]);
  absx31 = muDoubleScalarAbs(x[2]);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    b_x[0] = x[1];
    b_x[1] = x[0];
    b_x[3] = x[4];
    b_x[4] = x[3];
    b_x[6] = x[7];
    b_x[7] = x[6];
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      b_x[0] = x[2];
      b_x[2] = x[0];
      b_x[3] = x[5];
      b_x[5] = x[3];
      b_x[6] = x[8];
      b_x[8] = x[6];
    }
  }

  absx11 = b_x[1] / b_x[0];
  b_x[1] /= b_x[0];
  absx21 = b_x[2] / b_x[0];
  b_x[2] /= b_x[0];
  b_x[4] -= absx11 * b_x[3];
  b_x[5] -= absx21 * b_x[3];
  b_x[7] -= absx11 * b_x[6];
  b_x[8] -= absx21 * b_x[6];
  if (muDoubleScalarAbs(b_x[5]) > muDoubleScalarAbs(b_x[4])) {
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    b_x[1] = absx21;
    b_x[2] = absx11;
    absx11 = b_x[4];
    b_x[4] = b_x[5];
    b_x[5] = absx11;
    absx11 = b_x[7];
    b_x[7] = b_x[8];
    b_x[8] = absx11;
  }

  absx11 = b_x[5] / b_x[4];
  b_x[5] /= b_x[4];
  b_x[8] -= absx11 * b_x[7];
  absx11 = (b_x[5] * b_x[1] - b_x[2]) / b_x[8];
  absx21 = -(b_x[1] + b_x[7] * absx11) / b_x[4];
  y[p1] = ((1.0 - b_x[3] * absx21) - b_x[6] * absx11) / b_x[0];
  y[p1 + 1] = absx21;
  y[p1 + 2] = absx11;
  absx11 = -b_x[5] / b_x[8];
  absx21 = (1.0 - b_x[7] * absx11) / b_x[4];
  y[p2] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p2 + 1] = absx21;
  y[p2 + 2] = absx11;
  absx11 = 1.0 / b_x[8];
  absx21 = -b_x[7] * absx11 / b_x[4];
  y[p3] = -(b_x[3] * absx21 + b_x[6] * absx11) / b_x[0];
  y[p3 + 1] = absx21;
  y[p3 + 2] = absx11;
  emlrtPushRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
  absx11 = norm(x);
  absx21 = norm(y);
  absx31 = 1.0 / (absx11 * absx21);
  if ((absx11 == 0.0) || (absx21 == 0.0) || (absx31 == 0.0)) {
    emlrtPushRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
    eml_warning();
    emlrtPopRtStackR2012b(&d_emlrtRSI, emlrtRootTLSGlobal);
  } else {
    if (muDoubleScalarIsNaN(absx31) || (absx31 < 2.2204460492503131E-16)) {
      emlrtPushRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      b_y = NULL;
      m0 = mxCreateCharArray(2, iv1);
      for (p1 = 0; p1 < 8; p1++) {
        cv0[p1] = cv1[p1];
      }

      emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 8, m0, cv0);
      emlrtAssign(&b_y, m0);
      c_y = NULL;
      m0 = mxCreateDoubleScalar(14.0);
      emlrtAssign(&c_y, m0);
      d_y = NULL;
      m0 = mxCreateDoubleScalar(6.0);
      emlrtAssign(&d_y, m0);
      e_y = NULL;
      m0 = mxCreateDoubleScalar(absx31);
      emlrtAssign(&e_y, m0);
      emlrt_marshallIn(c_sprintf(b_sprintf(b_y, c_y, d_y, &c_emlrtMCI), e_y,
        &d_emlrtMCI), "sprintf", str);
      emlrtPopRtStackR2012b(&g_emlrtRSI, emlrtRootTLSGlobal);
      b_eml_warning(str);
      emlrtPopRtStackR2012b(&e_emlrtRSI, emlrtRootTLSGlobal);
    }
  }

  emlrtPopRtStackR2012b(&c_emlrtRSI, emlrtRootTLSGlobal);
}

/* End of code generation (inv.c) */
