/*
 * transformationResiduals.c
 *
 * Code generation for function 'transformationResiduals'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "transformImageCoordinates.h"
#include "interp2.h"
#include "inv.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 32, "transformationResiduals",
  "/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m"
};

static emlrtRSInfo b_emlrtRSI = { 34, "transformationResiduals",
  "/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m"
};

/* Function Definitions */
void transformationResiduals(c_transformationResidualsStackD *SD, const real_T
  MATRIXELEMENTS[6], const real_T IMAGE1[16384], const real_T IMAGE2[16384],
  const real_T SPATIALWINDOW[16384], const real_T XGRID[16384], const real_T
  YGRID[16384], real_T RESIDUAL[16384])
{
  real_T dv0[9];
  real_T dv1[9];
  real_T b_MATRIXELEMENTS[9];
  int32_T i0;
  static const int8_T iv0[3] = { 0, 0, 1 };

  real_T dv2[9];
  real_T dv3[9];
  int32_T i1;
  int32_T i2;

  /*  transformationResiduals(IMAGE1, MATRIXELEMENTS) calculates the residual intensity between images  */
  /*  INPUTS */
  /*    IMAGE1 = The original first image that is to be transformed (before windowing) */
  /*  X-scaling */
  /*  Rotation angle */
  /*  X-shearing */
  /*  Y-shearing */
  /*  X-translation */
  /*  Y-translation */
  /*  % Matrix elements */
  /*  Isotropic scaling matrix */
  /*  Rotation matrix */
  /*  Shearing matrix */
  /*  Translation matrix */
  /*  Affine transformation matrix  */
  /*  These are the vertical and horizontal coordinates of the inverse transformation. These */
  /*  are the non-integer image coordinates whose (interpolated) image values will */
  /*  make up the integer image coordinates in the output image. */
  emlrtPushRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
  dv1[0] = 1.0;
  dv1[3] = 0.0;
  dv1[6] = MATRIXELEMENTS[4];
  dv1[1] = 0.0;
  dv1[4] = 1.0;
  dv1[7] = MATRIXELEMENTS[5];
  b_MATRIXELEMENTS[0] = MATRIXELEMENTS[0];
  b_MATRIXELEMENTS[3] = 0.0;
  b_MATRIXELEMENTS[6] = 0.0;
  b_MATRIXELEMENTS[1] = 0.0;
  b_MATRIXELEMENTS[4] = MATRIXELEMENTS[0];
  b_MATRIXELEMENTS[7] = 0.0;
  for (i0 = 0; i0 < 3; i0++) {
    dv1[2 + 3 * i0] = iv0[i0];
    b_MATRIXELEMENTS[2 + 3 * i0] = iv0[i0];
  }

  dv3[0] = muDoubleScalarCos(MATRIXELEMENTS[1]);
  dv3[3] = -muDoubleScalarSin(MATRIXELEMENTS[1]);
  dv3[6] = 0.0;
  dv3[1] = muDoubleScalarSin(MATRIXELEMENTS[1]);
  dv3[4] = muDoubleScalarCos(MATRIXELEMENTS[1]);
  dv3[7] = 0.0;
  for (i0 = 0; i0 < 3; i0++) {
    for (i1 = 0; i1 < 3; i1++) {
      dv2[i0 + 3 * i1] = 0.0;
      for (i2 = 0; i2 < 3; i2++) {
        dv2[i0 + 3 * i1] += dv1[i0 + 3 * i2] * b_MATRIXELEMENTS[i2 + 3 * i1];
      }
    }

    dv3[2 + 3 * i0] = iv0[i0];
  }

  b_MATRIXELEMENTS[0] = 1.0;
  b_MATRIXELEMENTS[3] = MATRIXELEMENTS[2];
  b_MATRIXELEMENTS[6] = 0.0;
  b_MATRIXELEMENTS[1] = MATRIXELEMENTS[3];
  b_MATRIXELEMENTS[4] = 1.0;
  b_MATRIXELEMENTS[7] = 0.0;
  for (i0 = 0; i0 < 3; i0++) {
    for (i1 = 0; i1 < 3; i1++) {
      dv1[i0 + 3 * i1] = 0.0;
      for (i2 = 0; i2 < 3; i2++) {
        dv1[i0 + 3 * i1] += dv2[i0 + 3 * i2] * dv3[i2 + 3 * i1];
      }
    }

    b_MATRIXELEMENTS[2 + 3 * i0] = iv0[i0];
  }

  for (i0 = 0; i0 < 3; i0++) {
    for (i1 = 0; i1 < 3; i1++) {
      dv0[i0 + 3 * i1] = 0.0;
      for (i2 = 0; i2 < 3; i2++) {
        dv0[i0 + 3 * i1] += dv1[i0 + 3 * i2] * b_MATRIXELEMENTS[i2 + 3 * i1];
      }
    }
  }

  inv(dv0, dv1);
  transformImageCoordinates(SD, dv1, XGRID, YGRID, SD->f1.transformedFirstImage,
    SD->f1.xInv);
  emlrtPopRtStackR2012b(&emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);
  memcpy(&SD->f1.b_transformedFirstImage[0], &SD->f1.transformedFirstImage[0],
         sizeof(real_T) << 14);
  interp2(XGRID, YGRID, IMAGE1, SD->f1.xInv, SD->f1.b_transformedFirstImage,
          SD->f1.transformedFirstImage);
  for (i0 = 0; i0 < 16384; i0++) {
    SD->f1.transformedFirstImage[i0] *= SPATIALWINDOW[i0];
  }

  emlrtPopRtStackR2012b(&b_emlrtRSI, emlrtRootTLSGlobal);

  /*  Calculate residuals */
  for (i0 = 0; i0 < 16384; i0++) {
    RESIDUAL[i0] = SD->f1.transformedFirstImage[i0] - IMAGE2[i0];
  }
}

/* End of code generation (transformationResiduals.c) */
