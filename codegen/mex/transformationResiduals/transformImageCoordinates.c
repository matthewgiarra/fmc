/*
 * transformImageCoordinates.c
 *
 * Code generation for function 'transformImageCoordinates'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "transformImageCoordinates.h"

/* Variable Definitions */
static emlrtRSInfo h_emlrtRSI = { 31, "transformImageCoordinates",
  "/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m" };

static emlrtRSInfo i_emlrtRSI = { 55, "mtimes",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo j_emlrtRSI = { 54, "eml_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m" };

static emlrtRSInfo k_emlrtRSI = { 32, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo l_emlrtRSI = { 110, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo m_emlrtRSI = { 111, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo n_emlrtRSI = { 112, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo o_emlrtRSI = { 113, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo p_emlrtRSI = { 114, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo q_emlrtRSI = { 115, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo r_emlrtRSI = { 119, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo s_emlrtRSI = { 122, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo t_emlrtRSI = { 125, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo u_emlrtRSI = { 128, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo v_emlrtRSI = { 131, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo w_emlrtRSI = { 134, "eml_blas_xgemm",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"
};

static emlrtRSInfo x_emlrtRSI = { 14, "eml_c_cast",
  "/usr/local/MATLAB/R2013a/toolbox/eml/lib/matlab/eml/blas/external/eml_c_cast.m"
};

/* Function Definitions */
void transformImageCoordinates(c_transformationResidualsStackD *SD, const real_T
  TRANSFORM[9], const real_T XGRID[16384], const real_T YGRID[16384], real_T
  YOUT[16384], real_T XOUT[16384])
{
  real_T yPointsVect[16384];
  int32_T k;
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;

  /*  TRANSFORM = affine matrix */
  /*  IMAGESIZE = size of image (pixels); IMAGESIZE = [nRows nCols] */
  /*  Default to not cropping the image. */
  /*  Count rows and columns in the original image  */
  /*  Determine the center of the transformation. Default to center of grid */
  /*  Shift grid origin to the center of grid */
  /*  Reshape coordinate matrices into vectors */
  for (k = 0; k < 16384; k++) {
    SD->f0.xPointsVect[k] = XGRID[k] - 64.5;
    yPointsVect[k] = YGRID[k] - 64.5;
  }

  /*  Apply the transformation matrix to the shifted points */
  emlrtPushRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);
  for (k = 0; k < 16384; k++) {
    SD->f0.b[3 * k] = SD->f0.xPointsVect[k];
    SD->f0.b[1 + 3 * k] = yPointsVect[k];
    SD->f0.b[2 + 3 * k] = 1.0;
  }

  emlrtPushRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
  alpha1 = 1.0;
  beta1 = 0.0;
  TRANSB = 'N';
  TRANSA = 'N';
  memset(&SD->f0.transformedPoints[0], 0, 49152U * sizeof(real_T));
  emlrtPushRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  m_t = (ptrdiff_t)(3);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&l_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  n_t = (ptrdiff_t)(16384);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&m_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  k_t = (ptrdiff_t)(3);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&n_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  lda_t = (ptrdiff_t)(3);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&o_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&p_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  ldb_t = (ptrdiff_t)(3);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&p_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&q_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  ldc_t = (ptrdiff_t)(3);
  emlrtPopRtStackR2012b(&x_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&q_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
  alpha1_t = (double *)(&alpha1);
  emlrtPopRtStackR2012b(&r_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
  Aia0_t = (double *)(&TRANSFORM[0]);
  emlrtPopRtStackR2012b(&s_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&t_emlrtRSI, emlrtRootTLSGlobal);
  Bib0_t = (double *)(&SD->f0.b[0]);
  emlrtPopRtStackR2012b(&t_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&u_emlrtRSI, emlrtRootTLSGlobal);
  beta1_t = (double *)(&beta1);
  emlrtPopRtStackR2012b(&u_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&v_emlrtRSI, emlrtRootTLSGlobal);
  Cic0_t = (double *)(&SD->f0.transformedPoints[0]);
  emlrtPopRtStackR2012b(&v_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPushRtStackR2012b(&w_emlrtRSI, emlrtRootTLSGlobal);
  dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t, Bib0_t,
        &ldb_t, beta1_t, Cic0_t, &ldc_t);
  emlrtPopRtStackR2012b(&w_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&k_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&j_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&i_emlrtRSI, emlrtRootTLSGlobal);
  emlrtPopRtStackR2012b(&h_emlrtRSI, emlrtRootTLSGlobal);

  /*  Extract the X- and Y coordinates of the transformed points. */
  /*  Turn transformed X- and Y- coordinates back into matrices and shift the origin back to (1, 1) */
  for (k = 0; k < 16384; k++) {
    /*  Crop the image if specified. */
    /*  Don't crop the image if cropping wasn't specified. */
    XOUT[k] = SD->f0.transformedPoints[3 * k] + 64.5;
    YOUT[k] = SD->f0.transformedPoints[1 + 3 * k] + 64.5;
  }
}

/* End of code generation (transformImageCoordinates.c) */
