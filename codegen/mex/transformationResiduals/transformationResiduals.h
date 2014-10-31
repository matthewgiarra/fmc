/*
 * transformationResiduals.h
 *
 * Code generation for function 'transformationResiduals'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

#ifndef __TRANSFORMATIONRESIDUALS_H__
#define __TRANSFORMATIONRESIDUALS_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "transformationResiduals_types.h"

/* Function Declarations */
extern void transformationResiduals(c_transformationResidualsStackD *SD, const real_T MATRIXELEMENTS[6], const real_T IMAGE1[16384], const real_T IMAGE2[16384], const real_T SPATIALWINDOW[16384], const real_T XGRID[16384], const real_T YGRID[16384], real_T RESIDUAL[16384]);
#endif
/* End of code generation (transformationResiduals.h) */
