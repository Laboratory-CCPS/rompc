#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:get_condensed_hessian";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 5) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly five input arguments expected.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most one output argument expected.");
    }

    const mat_t W = cview_assert_nonempty_double_matrix(prhs[0], 0, MEX_FILE_ID);
    const mat_t Q = cview_assert_nonempty_square_double_matrix(prhs[1], 1, MEX_FILE_ID);
    const mat_t R = cview_assert_nonempty_square_double_matrix(prhs[2], 2, MEX_FILE_ID);
    const mat_t P = cview_assert_nonempty_square_double_matrix(prhs[3], 3, MEX_FILE_ID);
    const lapackint_t N = get_assert_positive_integer(prhs[4], 4, MEX_FILE_ID);

    assert_mex(Q.m == P.m, "Q and P expected to have same dimension.", MEX_FILE_ID);
    assert_mex(W.m == N * (Q.n + R.n),
               "W expected to have N * (nb. of columns of Q + nb. of columns of R) rows.",
               MEX_FILE_ID);

    mxArray* mxH = mxCreateDoubleMatrix(W.n, W.n, mxREAL);
    mat_t H = mview_matrix_from_mx(mxH);

    sizedbuffer_t work = sizedbuffer_create(W.n * ((Q.m > R.m) ? Q.m : R.m));

    get_condensed_hessian(&W, &Q, &R, &P, N, &H, &work);

    sizedbuffer_free(&work);

    plhs[0] = mxH;

    return;
}