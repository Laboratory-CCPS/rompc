#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "mx_helper.h"
#include "algorithm/alg.h"
#include "algorithm/linalg.h"

static const char* const MEX_FILE_ID = "ROMPC:multiplier_line_search";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly two input arguments expected.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most one output argument expected.");
    }

    const mat_t nu = cview_assert_nonempty_double_vector(prhs[0], 0, MEX_FILE_ID);
    const mat_t Delta_nu = cview_assert_nonempty_double_vector(prhs[1], 1, MEX_FILE_ID);

    assert_mex(nu.m == Delta_nu.m, "nu and Delta_nu expected to have same length.", MEX_FILE_ID);

    mxArray* mxAlpha = mxCreateDoubleMatrix(1, 1, mxREAL);

    *mxGetPr(mxAlpha) = multiplier_line_search(&nu, &Delta_nu);

    plhs[0] = mxAlpha;

    return;
}
