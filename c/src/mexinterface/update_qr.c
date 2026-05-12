#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:update_qr";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs != 2) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "exactly two input arguments expected.");
    }

    if (nlhs > 2) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most two output arguments expected.");
    }

    const mat_t Hcheck = cview_assert_nonempty_double_matrix(prhs[0], 0, MEX_FILE_ID);
    const mat_t Qcheck2 = cview_assert_nonempty_double_matrix(prhs[1], 1, MEX_FILE_ID);

    assert_mex(Hcheck.m == Hcheck.n + 1,
               "Hcheck expected to have exactly one row more than columns.", MEX_FILE_ID);
    assert_mex(Hcheck.n == Qcheck2.n,
               "Hcheck and Qcheck2 expected to have the same number of columns.", MEX_FILE_ID);

    mxArray* mxRt = mxCreateDoubleMatrix(Hcheck.m - 1, Hcheck.n, mxREAL);
    mxArray* mxQhat2 = mxCreateDoubleMatrix(Qcheck2.m + 1, Qcheck2.n + 1, mxREAL);

    mat_t Rt = mview_matrix_from_mx(mxRt);
    mat_t Qhat2 = mview_matrix_from_mx(mxQhat2);

    sizedbuffer_t work = sizedbuffer_create(Hcheck.m * Hcheck.n);

    update_qr(&Hcheck, &Qcheck2, &Rt, &Qhat2, &work);

    sizedbuffer_free(&work);

    plhs[0] = mxRt;

    if (nlhs >= 2) {
        plhs[1] = mxQhat2;
    } else {
        mxDestroyArray(mxQhat2);
    }

    return;
}
