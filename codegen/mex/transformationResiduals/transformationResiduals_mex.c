/*
 * transformationResiduals_mex.c
 *
 * Code generation for function 'transformationResiduals'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "mex.h"
#include "transformationResiduals_api.h"
#include "transformationResiduals_initialize.h"
#include "transformationResiduals_terminate.h"

/* Function Declarations */
static void transformationResiduals_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
MEXFUNCTION_LINKAGE mxArray *emlrtMexFcnProperties(void);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "transformationResiduals", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
emlrtCTX emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void transformationResiduals_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *outputs[1];
  mxArray *inputs[6];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  c_transformationResidualsStackD* d_transformationResidualsStackD = (c_transformationResidualsStackD*)mxCalloc(1,sizeof(c_transformationResidualsStackD));
  /* Module initialization. */
  transformationResiduals_initialize(&emlrtContextGlobal);
  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 6, mxCHAR_CLASS, 23, "transformationResiduals");
  } else if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 23, "transformationResiduals");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = (mxArray *)prhs[n];
  }
  /* Call the function. */
  transformationResiduals_api(d_transformationResidualsStackD, (const mxArray**)inputs, (const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  transformationResiduals_terminate();
  mxFree(d_transformationResidualsStackD);
}

void transformationResiduals_atexit_wrapper(void)
{
   transformationResiduals_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(transformationResiduals_atexit_wrapper);
  /* Dispatch the entry-point. */
  transformationResiduals_mexFunction(nlhs, plhs, nrhs, prhs);
}

mxArray *emlrtMexFcnProperties(void)
{
  const char *mexProperties[] = {
    "Version",
    "ResolvedFunctions",
    "EntryPoints"};
  const char *epProperties[] = {
    "Name",
    "NumberOfInputs",
    "NumberOfOutputs",
    "ConstantInputs"};
  mxArray *xResult = mxCreateStructMatrix(1,1,3,mexProperties);
  mxArray *xEntryPoints = mxCreateStructMatrix(1,1,4,epProperties);
  mxArray *xInputs = NULL;
  xInputs = mxCreateLogicalMatrix(1, 6);
  mxSetFieldByNumber(xEntryPoints, 0, 0, mxCreateString("transformationResiduals"));
  mxSetFieldByNumber(xEntryPoints, 0, 1, mxCreateDoubleScalar(6));
  mxSetFieldByNumber(xEntryPoints, 0, 2, mxCreateDoubleScalar(1));
  mxSetFieldByNumber(xEntryPoints, 0, 3, xInputs);
  mxSetFieldByNumber(xResult, 0, 0, mxCreateString("8.1.0.604 (R2013a)"));
  mxSetFieldByNumber(xResult, 0, 1, (mxArray*)emlrtMexFcnResolvedFunctionsInfo());
  mxSetFieldByNumber(xResult, 0, 2, xEntryPoints);

  return xResult;
}
/* End of code generation (transformationResiduals_mex.c) */
