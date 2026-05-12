#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:get_condensed_primal_feasibility";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 9) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly nine input arguments expected.");
    }

    if (nlhs > 1) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most one output argument expected.");
    }

    const mat_t w = cview_assert_nonempty_double_vector(prhs[0], 0, MEX_FILE_ID);
    const mat_t W = cview_assert_nonempty_double_matrix(prhs[1], 1, MEX_FILE_ID);
    const mat_t Gx = cview_assert_nonempty_double_matrix(prhs[2], 2, MEX_FILE_ID);
    const mat_t Gu = cview_assert_nonempty_double_matrix(prhs[3], 3, MEX_FILE_ID);
    const mat_t Gf = cview_assert_nonempty_double_matrix(prhs[4], 4, MEX_FILE_ID);
    const mat_t fx = cview_assert_nonempty_double_vector(prhs[5], 5, MEX_FILE_ID);
    const mat_t fu = cview_assert_nonempty_double_vector(prhs[6], 6, MEX_FILE_ID);
    const mat_t ff = cview_assert_nonempty_double_vector(prhs[7], 7, MEX_FILE_ID);
    const lapackint_t N = get_assert_positive_integer(prhs[8], 8, MEX_FILE_ID);

    const lapackint_t q = N * Gu.m + (N - 1) * Gx.m + Gf.m;

    assert_mex(w.m == W.m, "length of w expected to be equal to the number of rows of W.",
               MEX_FILE_ID);
    assert_mex(fx.m == Gx.m, "length of fx expected to be equal to the number of rows of Gx.",
               MEX_FILE_ID);
    assert_mex(fu.m == Gu.m, "length of fu expected to be equal to the number of rows of Gu.",
               MEX_FILE_ID);
    assert_mex(ff.m == Gf.m, "length of ff expected to be equal to the number of rows of Gf.",
               MEX_FILE_ID);

    assert_mex(W.m == N * (Gx.n + Gu.n),
               "W expected to have N * (nb. of columns of Gx + nb. of columns of Gu) rows.",
               MEX_FILE_ID);

    mxArray* mxDelta = mxCreateDoubleMatrix(q, 1, mxREAL);
    mat_t delta = mview_matrix_from_mx(mxDelta);

    get_condensed_primal_feasibility(&w, &W, &Gx, &Gu, &Gf, &fx, &fu, &ff, N, &delta);

    plhs[0] = mxDelta;

    return;
}
