/*
 * transformImageCoordinates.h
 *
 * Code generation for function 'transformImageCoordinates'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

#ifndef __TRANSFORMIMAGECOORDINATES_H__
#define __TRANSFORMIMAGECOORDINATES_H__
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
extern void transformImageCoordinates(c_transformationResidualsStackD *SD, const real_T TRANSFORM[9], const real_T XGRID[16384], const real_T YGRID[16384], real_T YOUT[16384], real_T XOUT[16384]);
#endif
/* End of code generation (transformImageCoordinates.h) */
