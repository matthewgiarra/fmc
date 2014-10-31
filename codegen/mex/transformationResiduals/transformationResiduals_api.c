/*
 * transformationResiduals_api.c
 *
 * Code generation for function 'transformationResiduals_api'
 *
 * C source code generated on: Thu May  1 18:35:45 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "transformationResiduals.h"
#include "transformationResiduals_api.h"

/* Function Declarations */
static const mxArray *b_emlrt_marshallOut(real_T u[16384]);
static void b_info_helper(ResolvedFunctionInfo info[108]);
static real_T (*c_emlrt_marshallIn(const mxArray *MATRIXELEMENTS, const char_T
  *identifier))[6];
static real_T (*d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId))[6];
static real_T (*e_emlrt_marshallIn(const mxArray *IMAGE1, const char_T
  *identifier))[16384];
static const mxArray *emlrt_marshallOut(ResolvedFunctionInfo u[108]);
static real_T (*f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId))[16384];
static real_T (*h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId))[6];
static real_T (*i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId))[16384];
static void info_helper(ResolvedFunctionInfo info[108]);

/* Function Definitions */
static const mxArray *b_emlrt_marshallOut(real_T u[16384])
{
  const mxArray *y;
  static const int32_T iv7[1] = { 0 };

  const mxArray *m5;
  static const int32_T iv8[1] = { 16384 };

  y = NULL;
  m5 = mxCreateNumericArray(1, (int32_T *)&iv7, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m5, (void *)u);
  mxSetDimensions((mxArray *)m5, iv8, 1);
  emlrtAssign(&y, m5);
  return y;
}

static void b_info_helper(ResolvedFunctionInfo info[108])
{
  info[64].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/strcmpi.m";
  info[64].name = "min";
  info[64].dominantType = "double";
  info[64].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  info[64].fileTimeLo = 1311280518U;
  info[64].fileTimeHi = 0U;
  info[64].mFileTimeLo = 0U;
  info[64].mFileTimeHi = 0U;
  info[65].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/datafun/min.m";
  info[65].name = "eml_min_or_max";
  info[65].dominantType = "char";
  info[65].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m";
  info[65].fileTimeLo = 1334096690U;
  info[65].fileTimeHi = 0U;
  info[65].mFileTimeLo = 0U;
  info[65].mFileTimeHi = 0U;
  info[66].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  info[66].name = "eml_scalar_eg";
  info[66].dominantType = "double";
  info[66].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[66].fileTimeLo = 1286843996U;
  info[66].fileTimeHi = 0U;
  info[66].mFileTimeLo = 0U;
  info[66].mFileTimeHi = 0U;
  info[67].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  info[67].name = "eml_scalexp_alloc";
  info[67].dominantType = "double";
  info[67].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m";
  info[67].fileTimeLo = 1352446460U;
  info[67].fileTimeHi = 0U;
  info[67].mFileTimeLo = 0U;
  info[67].mFileTimeHi = 0U;
  info[68].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_bin_extremum";
  info[68].name = "eml_index_class";
  info[68].dominantType = "";
  info[68].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[68].fileTimeLo = 1323192178U;
  info[68].fileTimeHi = 0U;
  info[68].mFileTimeLo = 0U;
  info[68].mFileTimeHi = 0U;
  info[69].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_min_or_max.m!eml_scalar_bin_extremum";
  info[69].name = "eml_scalar_eg";
  info[69].dominantType = "double";
  info[69].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[69].fileTimeLo = 1286843996U;
  info[69].fileTimeHi = 0U;
  info[69].mFileTimeLo = 0U;
  info[69].mFileTimeHi = 0U;
  info[70].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/strcmpi.m";
  info[70].name = "lower";
  info[70].dominantType = "char";
  info[70].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/lower.m";
  info[70].fileTimeLo = 1327440710U;
  info[70].fileTimeHi = 0U;
  info[70].mFileTimeLo = 0U;
  info[70].mFileTimeHi = 0U;
  info[71].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/lower.m";
  info[71].name = "eml_string_transform";
  info[71].dominantType = "char";
  info[71].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_string_transform.m";
  info[71].fileTimeLo = 1327440710U;
  info[71].fileTimeHi = 0U;
  info[71].mFileTimeLo = 0U;
  info[71].mFileTimeHi = 0U;
  info[72].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_string_transform.m";
  info[72].name = "eml_assert_supported_string";
  info[72].dominantType = "char";
  info[72].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_assert_supported_string.m";
  info[72].fileTimeLo = 1327440710U;
  info[72].fileTimeHi = 0U;
  info[72].mFileTimeLo = 0U;
  info[72].mFileTimeHi = 0U;
  info[73].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_assert_supported_string.m";
  info[73].name = "eml_charmax";
  info[73].dominantType = "";
  info[73].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_charmax.m";
  info[73].fileTimeLo = 1327440710U;
  info[73].fileTimeHi = 0U;
  info[73].mFileTimeLo = 0U;
  info[73].mFileTimeHi = 0U;
  info[74].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_string_transform.m";
  info[74].name = "eml_charmax";
  info[74].dominantType = "";
  info[74].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_charmax.m";
  info[74].fileTimeLo = 1327440710U;
  info[74].fileTimeHi = 0U;
  info[74].mFileTimeLo = 0U;
  info[74].mFileTimeHi = 0U;
  info[75].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_string_transform.m";
  info[75].name = "colon";
  info[75].dominantType = "int8";
  info[75].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  info[75].fileTimeLo = 1348217128U;
  info[75].fileTimeHi = 0U;
  info[75].mFileTimeLo = 0U;
  info[75].mFileTimeHi = 0U;
  info[76].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  info[76].name = "colon";
  info[76].dominantType = "int8";
  info[76].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  info[76].fileTimeLo = 1348217128U;
  info[76].fileTimeHi = 0U;
  info[76].mFileTimeLo = 0U;
  info[76].mFileTimeHi = 0U;
  info[77].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m";
  info[77].name = "floor";
  info[77].dominantType = "double";
  info[77].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  info[77].fileTimeLo = 1343855580U;
  info[77].fileTimeHi = 0U;
  info[77].mFileTimeLo = 0U;
  info[77].mFileTimeHi = 0U;
  info[78].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m";
  info[78].name = "eml_scalar_floor";
  info[78].dominantType = "double";
  info[78].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m";
  info[78].fileTimeLo = 1286843926U;
  info[78].fileTimeHi = 0U;
  info[78].mFileTimeLo = 0U;
  info[78].mFileTimeHi = 0U;
  info[79].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  info[79].name = "intmin";
  info[79].dominantType = "char";
  info[79].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  info[79].fileTimeLo = 1311280518U;
  info[79].fileTimeHi = 0U;
  info[79].mFileTimeLo = 0U;
  info[79].mFileTimeHi = 0U;
  info[80].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!checkrange";
  info[80].name = "intmax";
  info[80].dominantType = "char";
  info[80].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[80].fileTimeLo = 1311280516U;
  info[80].fileTimeHi = 0U;
  info[80].mFileTimeLo = 0U;
  info[80].mFileTimeHi = 0U;
  info[81].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  info[81].name = "intmin";
  info[81].dominantType = "char";
  info[81].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmin.m";
  info[81].fileTimeLo = 1311280518U;
  info[81].fileTimeHi = 0U;
  info[81].mFileTimeLo = 0U;
  info[81].mFileTimeHi = 0U;
  info[82].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  info[82].name = "intmax";
  info[82].dominantType = "char";
  info[82].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[82].fileTimeLo = 1311280516U;
  info[82].fileTimeHi = 0U;
  info[82].mFileTimeLo = 0U;
  info[82].mFileTimeHi = 0U;
  info[83].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_integer_colon_dispatcher";
  info[83].name = "eml_isa_uint";
  info[83].dominantType = "int8";
  info[83].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  info[83].fileTimeLo = 1286843984U;
  info[83].fileTimeHi = 0U;
  info[83].mFileTimeLo = 0U;
  info[83].mFileTimeHi = 0U;
  info[84].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  info[84].name = "eml_unsigned_class";
  info[84].dominantType = "char";
  info[84].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  info[84].fileTimeLo = 1323192180U;
  info[84].fileTimeHi = 0U;
  info[84].mFileTimeLo = 0U;
  info[84].mFileTimeHi = 0U;
  info[85].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  info[85].name = "eml_index_class";
  info[85].dominantType = "";
  info[85].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[85].fileTimeLo = 1323192178U;
  info[85].fileTimeHi = 0U;
  info[85].mFileTimeLo = 0U;
  info[85].mFileTimeHi = 0U;
  info[86].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  info[86].name = "intmax";
  info[86].dominantType = "char";
  info[86].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[86].fileTimeLo = 1311280516U;
  info[86].fileTimeHi = 0U;
  info[86].mFileTimeLo = 0U;
  info[86].mFileTimeHi = 0U;
  info[87].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  info[87].name = "eml_isa_uint";
  info[87].dominantType = "int8";
  info[87].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isa_uint.m";
  info[87].fileTimeLo = 1286843984U;
  info[87].fileTimeHi = 0U;
  info[87].mFileTimeLo = 0U;
  info[87].mFileTimeHi = 0U;
  info[88].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!integer_colon_length_nonnegd";
  info[88].name = "eml_index_plus";
  info[88].dominantType = "double";
  info[88].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[88].fileTimeLo = 1286843978U;
  info[88].fileTimeHi = 0U;
  info[88].mFileTimeLo = 0U;
  info[88].mFileTimeHi = 0U;
  info[89].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/colon.m!eml_signed_integer_colon";
  info[89].name = "eml_int_forloop_overflow_check";
  info[89].dominantType = "";
  info[89].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[89].fileTimeLo = 1346535540U;
  info[89].fileTimeHi = 0U;
  info[89].mFileTimeLo = 0U;
  info[89].mFileTimeHi = 0U;
  info[90].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_string_transform.m";
  info[90].name = "char";
  info[90].dominantType = "int8";
  info[90].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m";
  info[90].fileTimeLo = 1319755168U;
  info[90].fileTimeHi = 0U;
  info[90].mFileTimeLo = 0U;
  info[90].mFileTimeHi = 0U;
  info[91].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[91].name = "isequal";
  info[91].dominantType = "double";
  info[91].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  info[91].fileTimeLo = 1286843958U;
  info[91].fileTimeHi = 0U;
  info[91].mFileTimeLo = 0U;
  info[91].mFileTimeHi = 0U;
  info[92].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isequal.m";
  info[92].name = "eml_isequal_core";
  info[92].dominantType = "double";
  info[92].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m";
  info[92].fileTimeLo = 1286843986U;
  info[92].fileTimeHi = 0U;
  info[92].mFileTimeLo = 0U;
  info[92].mFileTimeHi = 0U;
  info[93].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_isequal_core.m!isequal_scalar";
  info[93].name = "isnan";
  info[93].dominantType = "double";
  info[93].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[93].fileTimeLo = 1286843960U;
  info[93].fileTimeHi = 0U;
  info[93].mFileTimeLo = 0U;
  info[93].mFileTimeHi = 0U;
  info[94].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[94].name = "eml_scalar_eg";
  info[94].dominantType = "double";
  info[94].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[94].fileTimeLo = 1286843996U;
  info[94].fileTimeHi = 0U;
  info[94].mFileTimeLo = 0U;
  info[94].mFileTimeHi = 0U;
  info[95].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[95].name = "length";
  info[95].dominantType = "double";
  info[95].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  info[95].fileTimeLo = 1303171406U;
  info[95].fileTimeHi = 0U;
  info[95].mFileTimeLo = 0U;
  info[95].mFileTimeHi = 0U;
  info[96].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m!intlength";
  info[96].name = "eml_index_class";
  info[96].dominantType = "";
  info[96].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[96].fileTimeLo = 1323192178U;
  info[96].fileTimeHi = 0U;
  info[96].mFileTimeLo = 0U;
  info[96].mFileTimeHi = 0U;
  info[97].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[97].name = "ismatrix";
  info[97].dominantType = "double";
  info[97].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m";
  info[97].fileTimeLo = 1331326458U;
  info[97].fileTimeHi = 0U;
  info[97].mFileTimeLo = 0U;
  info[97].mFileTimeHi = 0U;
  info[98].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!isincreasing";
  info[98].name = "length";
  info[98].dominantType = "double";
  info[98].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  info[98].fileTimeLo = 1303171406U;
  info[98].fileTimeHi = 0U;
  info[98].mFileTimeLo = 0U;
  info[98].mFileTimeHi = 0U;
  info[99].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!interp2_linear_or_nearest";
  info[99].name = "eml_int_forloop_overflow_check";
  info[99].dominantType = "";
  info[99].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[99].fileTimeLo = 1346535540U;
  info[99].fileTimeHi = 0U;
  info[99].mFileTimeLo = 0U;
  info[99].mFileTimeHi = 0U;
  info[100].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!interp2_linear_or_nearest";
  info[100].name = "eml_bsearch";
  info[100].dominantType = "double";
  info[100].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_bsearch.m";
  info[100].fileTimeLo = 1286843896U;
  info[100].fileTimeHi = 0U;
  info[100].mFileTimeLo = 0U;
  info[100].mFileTimeHi = 0U;
  info[101].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_bsearch.m";
  info[101].name = "eml_index_class";
  info[101].dominantType = "";
  info[101].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[101].fileTimeLo = 1323192178U;
  info[101].fileTimeHi = 0U;
  info[101].mFileTimeLo = 0U;
  info[101].mFileTimeHi = 0U;
  info[102].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_bsearch.m";
  info[102].name = "eml_unsigned_class";
  info[102].dominantType = "char";
  info[102].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  info[102].fileTimeLo = 1323192180U;
  info[102].fileTimeHi = 0U;
  info[102].mFileTimeLo = 0U;
  info[102].mFileTimeHi = 0U;
  info[103].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_unsigned_class.m";
  info[103].name = "eml_index_class";
  info[103].dominantType = "";
  info[103].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[103].fileTimeLo = 1323192178U;
  info[103].fileTimeHi = 0U;
  info[103].mFileTimeLo = 0U;
  info[103].mFileTimeHi = 0U;
  info[104].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_bsearch.m";
  info[104].name = "intmax";
  info[104].dominantType = "char";
  info[104].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[104].fileTimeLo = 1311280516U;
  info[104].fileTimeHi = 0U;
  info[104].mFileTimeLo = 0U;
  info[104].mFileTimeHi = 0U;
  info[105].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!scalar_bilinear_interp";
  info[105].name = "eml_scalar_eg";
  info[105].dominantType = "double";
  info[105].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[105].fileTimeLo = 1286843996U;
  info[105].fileTimeHi = 0U;
  info[105].mFileTimeLo = 0U;
  info[105].mFileTimeHi = 0U;
  info[106].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!scalar_bilinear_interp";
  info[106].name = "mrdivide";
  info[106].dominantType = "double";
  info[106].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[106].fileTimeLo = 1357973148U;
  info[106].fileTimeHi = 0U;
  info[106].mFileTimeLo = 1319755166U;
  info[106].mFileTimeHi = 0U;
  info[107].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m!scalar_bilinear_interp";
  info[107].name = "mtimes";
  info[107].dominantType = "double";
  info[107].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[107].fileTimeLo = 1289541292U;
  info[107].fileTimeHi = 0U;
  info[107].mFileTimeLo = 0U;
  info[107].mFileTimeHi = 0U;
}

static real_T (*c_emlrt_marshallIn(const mxArray *MATRIXELEMENTS, const char_T
  *identifier))[6]
{
  real_T (*y)[6];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = d_emlrt_marshallIn(emlrtAlias(MATRIXELEMENTS), &thisId);
  emlrtDestroyArray(&MATRIXELEMENTS);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *
  parentId))[6]
{
  real_T (*y)[6];
  y = h_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const mxArray *IMAGE1, const char_T
  *identifier))[16384]
{
  real_T (*y)[16384];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(emlrtAlias(IMAGE1), &thisId);
  emlrtDestroyArray(&IMAGE1);
  return y;
}
  static const mxArray *emlrt_marshallOut(ResolvedFunctionInfo u[108])
{
  const mxArray *y;
  int32_T iv6[1];
  int32_T i3;
  ResolvedFunctionInfo *r0;
  const char * b_u;
  const mxArray *b_y;
  const mxArray *m4;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  uint32_T c_u;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  y = NULL;
  iv6[0] = 108;
  emlrtAssign(&y, mxCreateStructArray(1, iv6, 0, NULL));
  for (i3 = 0; i3 < 108; i3++) {
    r0 = &u[i3];
    b_u = r0->context;
    b_y = NULL;
    m4 = mxCreateString(b_u);
    emlrtAssign(&b_y, m4);
    emlrtAddField(y, b_y, "context", i3);
    b_u = r0->name;
    c_y = NULL;
    m4 = mxCreateString(b_u);
    emlrtAssign(&c_y, m4);
    emlrtAddField(y, c_y, "name", i3);
    b_u = r0->dominantType;
    d_y = NULL;
    m4 = mxCreateString(b_u);
    emlrtAssign(&d_y, m4);
    emlrtAddField(y, d_y, "dominantType", i3);
    b_u = r0->resolved;
    e_y = NULL;
    m4 = mxCreateString(b_u);
    emlrtAssign(&e_y, m4);
    emlrtAddField(y, e_y, "resolved", i3);
    c_u = r0->fileTimeLo;
    f_y = NULL;
    m4 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m4) = c_u;
    emlrtAssign(&f_y, m4);
    emlrtAddField(y, f_y, "fileTimeLo", i3);
    c_u = r0->fileTimeHi;
    g_y = NULL;
    m4 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m4) = c_u;
    emlrtAssign(&g_y, m4);
    emlrtAddField(y, g_y, "fileTimeHi", i3);
    c_u = r0->mFileTimeLo;
    h_y = NULL;
    m4 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m4) = c_u;
    emlrtAssign(&h_y, m4);
    emlrtAddField(y, h_y, "mFileTimeLo", i3);
    c_u = r0->mFileTimeHi;
    i_y = NULL;
    m4 = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    *(uint32_T *)mxGetData(m4) = c_u;
    emlrtAssign(&i_y, m4);
    emlrtAddField(y, i_y, "mFileTimeHi", i3);
  }

  return y;
}

static real_T (*f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId))[16384]
{
  real_T (*y)[16384];
  y = i_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*h_emlrt_marshallIn(const mxArray *src, const
  emlrtMsgIdentifier *msgId))[6]
{
  real_T (*ret)[6];
  int32_T iv10[2];
  int32_T i5;
  for (i5 = 0; i5 < 2; i5++) {
    iv10[i5] = 1 + 5 * i5;
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv10);
  ret = (real_T (*)[6])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T (*i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId))[16384]
{
  real_T (*ret)[16384];
  int32_T iv11[2];
  int32_T i;
  for (i = 0; i < 2; i++) {
    iv11[i] = 128;
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv11);
  ret = (real_T (*)[16384])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static void info_helper(ResolvedFunctionInfo info[108])
{
  info[0].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[0].name = "cos";
  info[0].dominantType = "double";
  info[0].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  info[0].fileTimeLo = 1343855572U;
  info[0].fileTimeHi = 0U;
  info[0].mFileTimeLo = 0U;
  info[0].mFileTimeHi = 0U;
  info[1].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m";
  info[1].name = "eml_scalar_cos";
  info[1].dominantType = "double";
  info[1].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m";
  info[1].fileTimeLo = 1286843922U;
  info[1].fileTimeHi = 0U;
  info[1].mFileTimeLo = 0U;
  info[1].mFileTimeHi = 0U;
  info[2].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[2].name = "sin";
  info[2].dominantType = "double";
  info[2].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  info[2].fileTimeLo = 1343855586U;
  info[2].fileTimeHi = 0U;
  info[2].mFileTimeLo = 0U;
  info[2].mFileTimeHi = 0U;
  info[3].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m";
  info[3].name = "eml_scalar_sin";
  info[3].dominantType = "double";
  info[3].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m";
  info[3].fileTimeLo = 1286843936U;
  info[3].fileTimeHi = 0U;
  info[3].mFileTimeLo = 0U;
  info[3].mFileTimeHi = 0U;
  info[4].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[4].name = "mtimes";
  info[4].dominantType = "double";
  info[4].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[4].fileTimeLo = 1289541292U;
  info[4].fileTimeHi = 0U;
  info[4].mFileTimeLo = 0U;
  info[4].mFileTimeHi = 0U;
  info[5].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[5].name = "eml_index_class";
  info[5].dominantType = "";
  info[5].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[5].fileTimeLo = 1323192178U;
  info[5].fileTimeHi = 0U;
  info[5].mFileTimeLo = 0U;
  info[5].mFileTimeHi = 0U;
  info[6].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[6].name = "eml_scalar_eg";
  info[6].dominantType = "double";
  info[6].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[6].fileTimeLo = 1286843996U;
  info[6].fileTimeHi = 0U;
  info[6].mFileTimeLo = 0U;
  info[6].mFileTimeHi = 0U;
  info[7].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[7].name = "eml_xgemm";
  info[7].dominantType = "char";
  info[7].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  info[7].fileTimeLo = 1299098372U;
  info[7].fileTimeHi = 0U;
  info[7].mFileTimeLo = 0U;
  info[7].mFileTimeHi = 0U;
  info[8].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m";
  info[8].name = "eml_blas_inline";
  info[8].dominantType = "";
  info[8].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m";
  info[8].fileTimeLo = 1299098368U;
  info[8].fileTimeHi = 0U;
  info[8].mFileTimeLo = 0U;
  info[8].mFileTimeHi = 0U;
  info[9].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold";
  info[9].name = "mtimes";
  info[9].dominantType = "double";
  info[9].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[9].fileTimeLo = 1289541292U;
  info[9].fileTimeHi = 0U;
  info[9].mFileTimeLo = 0U;
  info[9].mFileTimeHi = 0U;
  info[10].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  info[10].name = "eml_index_class";
  info[10].dominantType = "";
  info[10].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[10].fileTimeLo = 1323192178U;
  info[10].fileTimeHi = 0U;
  info[10].mFileTimeLo = 0U;
  info[10].mFileTimeHi = 0U;
  info[11].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  info[11].name = "eml_scalar_eg";
  info[11].dominantType = "double";
  info[11].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[11].fileTimeLo = 1286843996U;
  info[11].fileTimeHi = 0U;
  info[11].mFileTimeLo = 0U;
  info[11].mFileTimeHi = 0U;
  info[12].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m";
  info[12].name = "eml_refblas_xgemm";
  info[12].dominantType = "char";
  info[12].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m";
  info[12].fileTimeLo = 1299098374U;
  info[12].fileTimeHi = 0U;
  info[12].mFileTimeLo = 0U;
  info[12].mFileTimeHi = 0U;
  info[13].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[13].name = "inv";
  info[13].dominantType = "double";
  info[13].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m";
  info[13].fileTimeLo = 1305343200U;
  info[13].fileTimeHi = 0U;
  info[13].mFileTimeLo = 0U;
  info[13].mFileTimeHi = 0U;
  info[14].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3";
  info[14].name = "eml_index_class";
  info[14].dominantType = "";
  info[14].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[14].fileTimeLo = 1323192178U;
  info[14].fileTimeHi = 0U;
  info[14].mFileTimeLo = 0U;
  info[14].mFileTimeHi = 0U;
  info[15].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3";
  info[15].name = "abs";
  info[15].dominantType = "double";
  info[15].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[15].fileTimeLo = 1343855566U;
  info[15].fileTimeHi = 0U;
  info[15].mFileTimeLo = 0U;
  info[15].mFileTimeHi = 0U;
  info[16].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[16].name = "eml_scalar_abs";
  info[16].dominantType = "double";
  info[16].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m";
  info[16].fileTimeLo = 1286843912U;
  info[16].fileTimeHi = 0U;
  info[16].mFileTimeLo = 0U;
  info[16].mFileTimeHi = 0U;
  info[17].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3";
  info[17].name = "eml_div";
  info[17].dominantType = "double";
  info[17].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  info[17].fileTimeLo = 1313373010U;
  info[17].fileTimeHi = 0U;
  info[17].mFileTimeLo = 0U;
  info[17].mFileTimeHi = 0U;
  info[18].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3";
  info[18].name = "mtimes";
  info[18].dominantType = "double";
  info[18].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[18].fileTimeLo = 1289541292U;
  info[18].fileTimeHi = 0U;
  info[18].mFileTimeLo = 0U;
  info[18].mFileTimeHi = 0U;
  info[19].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!inv3x3";
  info[19].name = "eml_index_plus";
  info[19].dominantType = "double";
  info[19].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[19].fileTimeLo = 1286843978U;
  info[19].fileTimeHi = 0U;
  info[19].mFileTimeLo = 0U;
  info[19].mFileTimeHi = 0U;
  info[20].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m";
  info[20].name = "eml_index_class";
  info[20].dominantType = "";
  info[20].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[20].fileTimeLo = 1323192178U;
  info[20].fileTimeHi = 0U;
  info[20].mFileTimeLo = 0U;
  info[20].mFileTimeHi = 0U;
  info[21].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[21].name = "norm";
  info[21].dominantType = "double";
  info[21].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m";
  info[21].fileTimeLo = 1336547294U;
  info[21].fileTimeHi = 0U;
  info[21].mFileTimeLo = 0U;
  info[21].mFileTimeHi = 0U;
  info[22].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm";
  info[22].name = "abs";
  info[22].dominantType = "double";
  info[22].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m";
  info[22].fileTimeLo = 1343855566U;
  info[22].fileTimeHi = 0U;
  info[22].mFileTimeLo = 0U;
  info[22].mFileTimeHi = 0U;
  info[23].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm";
  info[23].name = "isnan";
  info[23].dominantType = "double";
  info[23].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[23].fileTimeLo = 1286843960U;
  info[23].fileTimeHi = 0U;
  info[23].mFileTimeLo = 0U;
  info[23].mFileTimeHi = 0U;
  info[24].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/norm.m!mat1norm";
  info[24].name = "eml_guarded_nan";
  info[24].dominantType = "char";
  info[24].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m";
  info[24].fileTimeLo = 1286843976U;
  info[24].fileTimeHi = 0U;
  info[24].mFileTimeLo = 0U;
  info[24].mFileTimeHi = 0U;
  info[25].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_guarded_nan.m";
  info[25].name = "eml_is_float_class";
  info[25].dominantType = "char";
  info[25].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  info[25].fileTimeLo = 1286843982U;
  info[25].fileTimeHi = 0U;
  info[25].mFileTimeLo = 0U;
  info[25].mFileTimeHi = 0U;
  info[26].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[26].name = "mtimes";
  info[26].dominantType = "double";
  info[26].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[26].fileTimeLo = 1289541292U;
  info[26].fileTimeHi = 0U;
  info[26].mFileTimeLo = 0U;
  info[26].mFileTimeHi = 0U;
  info[27].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[27].name = "eml_warning";
  info[27].dominantType = "char";
  info[27].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_warning.m";
  info[27].fileTimeLo = 1286844002U;
  info[27].fileTimeHi = 0U;
  info[27].mFileTimeLo = 0U;
  info[27].mFileTimeHi = 0U;
  info[28].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[28].name = "isnan";
  info[28].dominantType = "double";
  info[28].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isnan.m";
  info[28].fileTimeLo = 1286843960U;
  info[28].fileTimeHi = 0U;
  info[28].mFileTimeLo = 0U;
  info[28].mFileTimeHi = 0U;
  info[29].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[29].name = "eps";
  info[29].dominantType = "char";
  info[29].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  info[29].fileTimeLo = 1326749596U;
  info[29].fileTimeHi = 0U;
  info[29].mFileTimeLo = 0U;
  info[29].mFileTimeHi = 0U;
  info[30].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  info[30].name = "eml_is_float_class";
  info[30].dominantType = "char";
  info[30].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_is_float_class.m";
  info[30].fileTimeLo = 1286843982U;
  info[30].fileTimeHi = 0U;
  info[30].mFileTimeLo = 0U;
  info[30].mFileTimeHi = 0U;
  info[31].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/eps.m";
  info[31].name = "eml_eps";
  info[31].dominantType = "char";
  info[31].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  info[31].fileTimeLo = 1326749596U;
  info[31].fileTimeHi = 0U;
  info[31].mFileTimeLo = 0U;
  info[31].mFileTimeHi = 0U;
  info[32].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_eps.m";
  info[32].name = "eml_float_model";
  info[32].dominantType = "char";
  info[32].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_float_model.m";
  info[32].fileTimeLo = 1326749596U;
  info[32].fileTimeHi = 0U;
  info[32].mFileTimeLo = 0U;
  info[32].mFileTimeHi = 0U;
  info[33].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/matfun/inv.m!checkcond";
  info[33].name = "eml_flt2str";
  info[33].dominantType = "double";
  info[33].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m";
  info[33].fileTimeLo = 1309476396U;
  info[33].fileTimeHi = 0U;
  info[33].mFileTimeLo = 0U;
  info[33].mFileTimeHi = 0U;
  info[34].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_flt2str.m";
  info[34].name = "char";
  info[34].dominantType = "double";
  info[34].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/char.m";
  info[34].fileTimeLo = 1319755168U;
  info[34].fileTimeHi = 0U;
  info[34].mFileTimeLo = 0U;
  info[34].mFileTimeHi = 0U;
  info[35].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[35].name = "transformImageCoordinates";
  info[35].dominantType = "double";
  info[35].resolved =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m";
  info[35].fileTimeLo = 1398983717U;
  info[35].fileTimeHi = 0U;
  info[35].mFileTimeLo = 0U;
  info[35].mFileTimeHi = 0U;
  info[36].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m";
  info[36].name = "mrdivide";
  info[36].dominantType = "double";
  info[36].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[36].fileTimeLo = 1357973148U;
  info[36].fileTimeHi = 0U;
  info[36].mFileTimeLo = 1319755166U;
  info[36].mFileTimeHi = 0U;
  info[37].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p";
  info[37].name = "rdivide";
  info[37].dominantType = "double";
  info[37].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[37].fileTimeLo = 1346535588U;
  info[37].fileTimeHi = 0U;
  info[37].mFileTimeLo = 0U;
  info[37].mFileTimeHi = 0U;
  info[38].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[38].name = "eml_scalexp_compatible";
  info[38].dominantType = "double";
  info[38].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m";
  info[38].fileTimeLo = 1286843996U;
  info[38].fileTimeHi = 0U;
  info[38].mFileTimeLo = 0U;
  info[38].mFileTimeHi = 0U;
  info[39].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m";
  info[39].name = "eml_div";
  info[39].dominantType = "double";
  info[39].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m";
  info[39].fileTimeLo = 1313373010U;
  info[39].fileTimeHi = 0U;
  info[39].mFileTimeLo = 0U;
  info[39].mFileTimeHi = 0U;
  info[40].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m";
  info[40].name = "reshape";
  info[40].dominantType = "double";
  info[40].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m";
  info[40].fileTimeLo = 1286843968U;
  info[40].fileTimeHi = 0U;
  info[40].mFileTimeLo = 0U;
  info[40].mFileTimeHi = 0U;
  info[41].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m";
  info[41].name = "eml_index_class";
  info[41].dominantType = "";
  info[41].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[41].fileTimeLo = 1323192178U;
  info[41].fileTimeHi = 0U;
  info[41].mFileTimeLo = 0U;
  info[41].mFileTimeHi = 0U;
  info[42].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size";
  info[42].name = "eml_index_class";
  info[42].dominantType = "";
  info[42].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[42].fileTimeLo = 1323192178U;
  info[42].fileTimeHi = 0U;
  info[42].mFileTimeLo = 0U;
  info[42].mFileTimeHi = 0U;
  info[43].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!varargin_nempty";
  info[43].name = "eml_index_class";
  info[43].dominantType = "";
  info[43].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[43].fileTimeLo = 1323192178U;
  info[43].fileTimeHi = 0U;
  info[43].mFileTimeLo = 0U;
  info[43].mFileTimeHi = 0U;
  info[44].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size";
  info[44].name = "eml_assert_valid_size_arg";
  info[44].dominantType = "double";
  info[44].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  info[44].fileTimeLo = 1286843894U;
  info[44].fileTimeHi = 0U;
  info[44].mFileTimeLo = 0U;
  info[44].mFileTimeHi = 0U;
  info[45].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!isintegral";
  info[45].name = "isinf";
  info[45].dominantType = "double";
  info[45].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/isinf.m";
  info[45].fileTimeLo = 1286843960U;
  info[45].fileTimeHi = 0U;
  info[45].mFileTimeLo = 0U;
  info[45].mFileTimeHi = 0U;
  info[46].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m!numel_for_size";
  info[46].name = "mtimes";
  info[46].dominantType = "double";
  info[46].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[46].fileTimeLo = 1289541292U;
  info[46].fileTimeHi = 0U;
  info[46].mFileTimeLo = 0U;
  info[46].mFileTimeHi = 0U;
  info[47].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  info[47].name = "eml_index_class";
  info[47].dominantType = "";
  info[47].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[47].fileTimeLo = 1323192178U;
  info[47].fileTimeHi = 0U;
  info[47].mFileTimeLo = 0U;
  info[47].mFileTimeHi = 0U;
  info[48].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m";
  info[48].name = "intmax";
  info[48].dominantType = "char";
  info[48].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[48].fileTimeLo = 1311280516U;
  info[48].fileTimeHi = 0U;
  info[48].mFileTimeLo = 0U;
  info[48].mFileTimeHi = 0U;
  info[49].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m!reshape_varargin_to_size";
  info[49].name = "eml_index_prod";
  info[49].dominantType = "coder.internal.indexInt";
  info[49].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m";
  info[49].fileTimeLo = 1286843980U;
  info[49].fileTimeHi = 0U;
  info[49].mFileTimeLo = 0U;
  info[49].mFileTimeHi = 0U;
  info[50].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m";
  info[50].name = "eml_index_class";
  info[50].dominantType = "";
  info[50].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[50].fileTimeLo = 1323192178U;
  info[50].fileTimeHi = 0U;
  info[50].mFileTimeLo = 0U;
  info[50].mFileTimeHi = 0U;
  info[51].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m";
  info[51].name = "eml_int_forloop_overflow_check";
  info[51].dominantType = "";
  info[51].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[51].fileTimeLo = 1346535540U;
  info[51].fileTimeHi = 0U;
  info[51].mFileTimeLo = 0U;
  info[51].mFileTimeHi = 0U;
  info[52].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper";
  info[52].name = "intmax";
  info[52].dominantType = "char";
  info[52].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[52].fileTimeLo = 1311280516U;
  info[52].fileTimeHi = 0U;
  info[52].mFileTimeLo = 0U;
  info[52].mFileTimeHi = 0U;
  info[53].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_prod.m";
  info[53].name = "eml_index_times";
  info[53].dominantType = "coder.internal.indexInt";
  info[53].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  info[53].fileTimeLo = 1286843980U;
  info[53].fileTimeHi = 0U;
  info[53].mFileTimeLo = 0U;
  info[53].mFileTimeHi = 0U;
  info[54].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m";
  info[54].name = "eml_index_class";
  info[54].dominantType = "";
  info[54].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m";
  info[54].fileTimeLo = 1323192178U;
  info[54].fileTimeHi = 0U;
  info[54].mFileTimeLo = 0U;
  info[54].mFileTimeHi = 0U;
  info[55].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m";
  info[55].name = "eml_scalar_eg";
  info[55].dominantType = "double";
  info[55].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m";
  info[55].fileTimeLo = 1286843996U;
  info[55].fileTimeHi = 0U;
  info[55].mFileTimeLo = 0U;
  info[55].mFileTimeHi = 0U;
  info[56].context = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/reshape.m";
  info[56].name = "eml_int_forloop_overflow_check";
  info[56].dominantType = "";
  info[56].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m";
  info[56].fileTimeLo = 1346535540U;
  info[56].fileTimeHi = 0U;
  info[56].mFileTimeLo = 0U;
  info[56].mFileTimeHi = 0U;
  info[57].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m";
  info[57].name = "length";
  info[57].dominantType = "double";
  info[57].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/length.m";
  info[57].fileTimeLo = 1303171406U;
  info[57].fileTimeHi = 0U;
  info[57].mFileTimeLo = 0U;
  info[57].mFileTimeHi = 0U;
  info[58].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/transformImageCoordinates.m";
  info[58].name = "mtimes";
  info[58].dominantType = "double";
  info[58].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m";
  info[58].fileTimeLo = 1289541292U;
  info[58].fileTimeHi = 0U;
  info[58].mFileTimeLo = 0U;
  info[58].mFileTimeHi = 0U;
  info[59].context =
    "[E]/home/voodoo/Desktop/branch_fastLogPolar/testCodes/transformationResiduals.m";
  info[59].name = "interp2";
  info[59].dominantType = "double";
  info[59].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[59].fileTimeLo = 1332190274U;
  info[59].fileTimeHi = 0U;
  info[59].mFileTimeLo = 0U;
  info[59].mFileTimeHi = 0U;
  info[60].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/polyfun/interp2.m";
  info[60].name = "strcmpi";
  info[60].dominantType = "char";
  info[60].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/strcmpi.m";
  info[60].fileTimeLo = 1327440710U;
  info[60].fileTimeHi = 0U;
  info[60].mFileTimeLo = 0U;
  info[60].mFileTimeHi = 0U;
  info[61].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/strcmpi.m";
  info[61].name = "eml_assert_supported_string";
  info[61].dominantType = "char";
  info[61].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_assert_supported_string.m";
  info[61].fileTimeLo = 1327440710U;
  info[61].fileTimeHi = 0U;
  info[61].mFileTimeLo = 0U;
  info[61].mFileTimeHi = 0U;
  info[62].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_assert_supported_string.m!inrange";
  info[62].name = "eml_charmax";
  info[62].dominantType = "";
  info[62].resolved =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_charmax.m";
  info[62].fileTimeLo = 1327440710U;
  info[62].fileTimeHi = 0U;
  info[62].mFileTimeLo = 0U;
  info[62].mFileTimeHi = 0U;
  info[63].context =
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/strfun/eml_charmax.m";
  info[63].name = "intmax";
  info[63].dominantType = "char";
  info[63].resolved = "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m";
  info[63].fileTimeLo = 1311280516U;
  info[63].fileTimeHi = 0U;
  info[63].mFileTimeLo = 0U;
  info[63].mFileTimeHi = 0U;
}

const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  ResolvedFunctionInfo info[108];
  nameCaptureInfo = NULL;
  info_helper(info);
  b_info_helper(info);
  emlrtAssign(&nameCaptureInfo, emlrt_marshallOut(info));
  emlrtNameCapturePostProcessR2012a(emlrtAlias(nameCaptureInfo));
  return nameCaptureInfo;
}

void transformationResiduals_api(c_transformationResidualsStackD *SD, const
  mxArray * const prhs[6], const mxArray *plhs[1])
{
  real_T (*RESIDUAL)[16384];
  real_T (*MATRIXELEMENTS)[6];
  real_T (*IMAGE1)[16384];
  real_T (*IMAGE2)[16384];
  real_T (*SPATIALWINDOW)[16384];
  real_T (*XGRID)[16384];
  real_T (*YGRID)[16384];
  RESIDUAL = (real_T (*)[16384])mxMalloc(sizeof(real_T [16384]));

  /* Marshall function inputs */
  MATRIXELEMENTS = c_emlrt_marshallIn(emlrtAlias(prhs[0]), "MATRIXELEMENTS");
  IMAGE1 = e_emlrt_marshallIn(emlrtAlias(prhs[1]), "IMAGE1");
  IMAGE2 = e_emlrt_marshallIn(emlrtAlias(prhs[2]), "IMAGE2");
  SPATIALWINDOW = e_emlrt_marshallIn(emlrtAlias(prhs[3]), "SPATIALWINDOW");
  XGRID = e_emlrt_marshallIn(emlrtAlias(prhs[4]), "XGRID");
  YGRID = e_emlrt_marshallIn(emlrtAlias(prhs[5]), "YGRID");

  /* Invoke the target function */
  transformationResiduals(SD, *MATRIXELEMENTS, *IMAGE1, *IMAGE2, *SPATIALWINDOW,
    *XGRID, *YGRID, *RESIDUAL);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(*RESIDUAL);
}

/* End of code generation (transformationResiduals_api.c) */
