/*
 * interp2.c
 *
 * Code generation for function 'interp2'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "interp2.h"
#include "transformationResiduals_mexutil.h"

/* Variable Definitions */
static emlrtRSInfo y_emlrtRSI = { 81, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

static emlrtRSInfo eb_emlrtRSI = { 92, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

static emlrtRSInfo fb_emlrtRSI = { 101, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

static emlrtRSInfo gb_emlrtRSI = { 102, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

static emlrtRSInfo hb_emlrtRSI = { 31, "eml_bsearch",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/eml_bsearch.m" };

static emlrtMCInfo g_emlrtMCI = { 93, 5, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

static emlrtMCInfo h_emlrtMCI = { 92, 15, "interp2",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/polyfun/interp2.m" };

/* Function Declarations */
static void error(const mxArray *b, emlrtMCInfo *location);

/* Function Definitions */
static void error(const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(emlrtRootTLSGlobal, 0, NULL, 1, &pArray, "error", TRUE,
                        location);
}

void interp2(const real_T x[16384], const real_T y[16384], const real_T z[16384],
             const real_T xi[16384], const real_T yi[16384], real_T zi[16384])
{
  int32_T k;
  int32_T exitg2;
  boolean_T p;
  boolean_T guard1 = FALSE;
  int32_T exitg1;
  const mxArray *b_y;
  static const int32_T iv5[2] = { 1, 37 };

  const mxArray *m3;
  char_T cv6[37];
  int32_T i;
  static const char_T cv7[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'i', 'n', 't', 'e', 'r', 'p', '2', '_', 'n', 'o',
    'n', 'i', 'n', 'c', 'r', 'e', 'a', 's', 'i', 'n', 'g', 'X', 'Y' };

  int32_T low_i;
  int32_T low_ip1;
  int32_T high_i;
  int32_T b_low_i;
  real_T rx;
  real_T qx1;
  real_T ry;
  emlrtPushRtStackR2012b(&y_emlrtRSI, emlrtRootTLSGlobal);
  k = 0;
  do {
    exitg2 = 0;
    if (k < 127) {
      if (!(x[k << 7] < x[(k + 1) << 7])) {
        p = FALSE;
        exitg2 = 1;
      } else {
        k++;
      }
    } else {
      p = TRUE;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  guard1 = FALSE;
  if (p) {
    k = 0;
    do {
      exitg1 = 0;
      if (k < 127) {
        if (!(y[k] < y[1 + k])) {
          p = FALSE;
          exitg1 = 1;
        } else {
          k++;
        }
      } else {
        p = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    if (p) {
      p = TRUE;
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    p = FALSE;
  }

  if (p) {
  } else {
    emlrtPushRtStackR2012b(&eb_emlrtRSI, emlrtRootTLSGlobal);
    b_y = NULL;
    m3 = mxCreateCharArray(2, iv5);
    for (i = 0; i < 37; i++) {
      cv6[i] = cv7[i];
    }

    emlrtInitCharArrayR2013a(emlrtRootTLSGlobal, 37, m3, cv6);
    emlrtAssign(&b_y, m3);
    error(message(b_y, &g_emlrtMCI), &h_emlrtMCI);
    emlrtPopRtStackR2012b(&eb_emlrtRSI, emlrtRootTLSGlobal);
  }

  for (k = 0; k < 16384; k++) {
    if ((xi[k] >= x[0]) && (xi[k] <= x[16256]) && (yi[k] >= y[0]) && (yi[k] <=
         y[127])) {
      emlrtPushRtStackR2012b(&fb_emlrtRSI, emlrtRootTLSGlobal);
      low_i = 0;
      low_ip1 = 2;
      high_i = 128;
      while (high_i > low_ip1) {
        emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
        i = (low_i + high_i) + 1;
        if (i >= 0) {
          i = (int32_T)((uint32_T)i >> 1);
        } else {
          i = ~(int32_T)((uint32_T)~i >> 1);
        }

        emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
        if (xi[k] >= x[(i - 1) << 7]) {
          low_i = i - 1;
          low_ip1 = i + 1;
        } else {
          high_i = i;
        }
      }

      emlrtPopRtStackR2012b(&fb_emlrtRSI, emlrtRootTLSGlobal);
      emlrtPushRtStackR2012b(&gb_emlrtRSI, emlrtRootTLSGlobal);
      b_low_i = 1;
      low_ip1 = 2;
      high_i = 128;
      while (high_i > low_ip1) {
        emlrtPushRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
        i = b_low_i + high_i;
        if (i >= 0) {
          i = (int32_T)((uint32_T)i >> 1);
        } else {
          i = ~(int32_T)((uint32_T)~i >> 1);
        }

        emlrtPopRtStackR2012b(&hb_emlrtRSI, emlrtRootTLSGlobal);
        if (yi[k] >= y[i - 1]) {
          b_low_i = i;
          low_ip1 = i + 1;
        } else {
          high_i = i;
        }
      }

      emlrtPopRtStackR2012b(&gb_emlrtRSI, emlrtRootTLSGlobal);
      rx = (xi[k] - x[low_i << 7]) / (x[(low_i + 1) << 7] - x[low_i << 7]);
      if (z[(b_low_i + (low_i << 7)) - 1] == z[(b_low_i + ((low_i + 1) << 7)) -
          1]) {
        qx1 = z[(b_low_i + (low_i << 7)) - 1];
      } else {
        qx1 = (1.0 - rx) * z[(b_low_i + (low_i << 7)) - 1] + rx * z[(b_low_i +
          ((low_i + 1) << 7)) - 1];
      }

      if (z[b_low_i + (low_i << 7)] == z[b_low_i + ((low_i + 1) << 7)]) {
        rx = z[b_low_i + (low_i << 7)];
      } else {
        rx = (1.0 - rx) * z[b_low_i + (low_i << 7)] + rx * z[b_low_i + ((low_i +
          1) << 7)];
      }

      if (qx1 == rx) {
      } else {
        ry = (yi[k] - y[b_low_i - 1]) / (y[b_low_i] - y[b_low_i - 1]);
        qx1 = (1.0 - ry) * qx1 + ry * rx;
      }

      zi[k] = qx1;
    } else {
      zi[k] = 0.0;
    }
  }

  emlrtPopRtStackR2012b(&y_emlrtRSI, emlrtRootTLSGlobal);
}

/* End of code generation (interp2.c) */
