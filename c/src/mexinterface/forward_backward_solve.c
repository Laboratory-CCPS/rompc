#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:forward_backward_solve";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly two input arguments expected.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most one output argument expected.");
    }

    const mat_t G = cview_assert_nonempty_square_double_matrix(prhs[0], 0, MEX_FILE_ID);
    const mat_t b = cview_assert_nonempty_double_vector(prhs[1], 1, MEX_FILE_ID);

    assert_mex(b.m == G.n, "length of be expected to be equal to the number of columns of G",
               MEX_FILE_ID);

    mxArray* mxX = mxCreateDoubleMatrix(G.n, 1, mxREAL);

    mat_t bx = mview_matrix_from_mx(mxX);
    mat_copy(&b, &bx);

    forward_backward_solve(&G, &bx);

    plhs[0] = mxX;

    return;
}
