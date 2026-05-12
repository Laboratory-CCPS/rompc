#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:get_condensed_quadratic_constraints";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 7) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly seven input arguments expected.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most one output argument expected.");
    }

    const mat_t W = cview_assert_nonempty_double_matrix(prhs[0], 0, MEX_FILE_ID);
    const mat_t Gx = cview_assert_nonempty_double_matrix(prhs[1], 1, MEX_FILE_ID);
    const mat_t Gu = cview_assert_nonempty_double_matrix(prhs[2], 2, MEX_FILE_ID);
    const mat_t Gf = cview_assert_nonempty_double_matrix(prhs[3], 3, MEX_FILE_ID);
    const mat_t y = cview_assert_nonempty_double_vector(prhs[4], 4, MEX_FILE_ID);
    const mat_t lambda = cview_assert_nonempty_double_vector(prhs[5], 5, MEX_FILE_ID);
    const lapackint_t N = get_assert_positive_integer(prhs[6], 6, MEX_FILE_ID);

    assert_mex(Gf.n == Gx.n, "Gf and Gx expected to have the same number of columns.", MEX_FILE_ID);
    assert_mex(W.m == N * (Gx.n + Gu.n),
               "W expected to have N * (nb. of columns of Gx + nb. of columns of Gu) rows.",
               MEX_FILE_ID);
    assert_mex(y.m == lambda.m, "y and lambda expected to have same length.", MEX_FILE_ID);

    mxArray* mxD = mxCreateDoubleMatrix(W.n, W.n, mxREAL);
    mat_t D = mview_matrix_from_mx(mxD);

    lapackint_t buf_dim = (Gx.m > Gu.m) ? Gx.m : Gu.m;
    buf_dim = (buf_dim > Gf.m) ? buf_dim : Gf.m;
    buf_dim = buf_dim * W.n * 2;

    sizedbuffer_t work = sizedbuffer_create(buf_dim);

    get_condensed_quadratic_constraints(&W, &Gx, &Gu, &Gf, &y, &lambda, N, &D, &work);

    sizedbuffer_free(&work);

    plhs[0] = mxD;

    return;
}
