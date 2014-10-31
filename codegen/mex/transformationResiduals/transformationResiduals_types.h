/*
 * transformationResiduals_types.h
 *
 * Code generation for function 'transformationResiduals'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

#ifndef __TRANSFORMATIONRESIDUALS_TYPES_H__
#define __TRANSFORMATIONRESIDUALS_TYPES_H__

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_ResolvedFunctionInfo
#define typedef_ResolvedFunctionInfo
typedef struct
{
    const char * context;
    const char * name;
    const char * dominantType;
    const char * resolved;
    uint32_T fileTimeLo;
    uint32_T fileTimeHi;
    uint32_T mFileTimeLo;
    uint32_T mFileTimeHi;
} ResolvedFunctionInfo;
#endif /*typedef_ResolvedFunctionInfo*/
#ifndef typedef_c_transformationResidualsStackD
#define typedef_c_transformationResidualsStackD
typedef struct
{
    struct
    {
        real_T b[49152];
        real_T transformedPoints[49152];
        real_T xPointsVect[16384];
    } f0;
    struct
    {
        real_T xInv[16384];
        real_T transformedFirstImage[16384];
        real_T b_transformedFirstImage[16384];
    } f1;
} c_transformationResidualsStackD;
#endif /*typedef_c_transformationResidualsStackD*/

#endif
/* End of code generation (transformationResiduals_types.h) */
