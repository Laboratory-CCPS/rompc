#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"
#include "rompc_mex_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:condensed_rompc_step";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 7 || nrhs > 8) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "seven or eight input arguments expected.");
    }

    if (nlhs > 4) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most four output arguments expected.");
    }

    const mat_t xt = cview_assert_nonempty_double_vector(prhs[0], 0, MEX_FILE_ID);
    const mat_t zt = cview_assert_nonempty_double_vector(prhs[1], 1, MEX_FILE_ID);
    const mat_t Rcheck = cview_assert_nonempty_square_double_matrix(prhs[2], 2, MEX_FILE_ID);
    const mat_t Qcheck = cview_assert_nonempty_square_double_matrix(prhs[3], 3, MEX_FILE_ID);

    const mat_t A =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[4], "A", 4, MEX_FILE_ID);
    const mat_t B =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[4], "B", 4, MEX_FILE_ID);
    const mat_t Qlq =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[4], "Q", 4, MEX_FILE_ID);
    const mat_t Rlq =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[4], "R", 4, MEX_FILE_ID);
    const mat_t Gx =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[4], "Gx", 4, MEX_FILE_ID);
    const mat_t fx =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[4], "fx", 4, MEX_FILE_ID);
    const mat_t Gu =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[4], "Gu", 4, MEX_FILE_ID);
    const mat_t fu =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[4], "fu", 4, MEX_FILE_ID);
    const sys_t sys = {
        .A = A, .B = B, .Q = Qlq, .R = Rlq, .Gx = Gx, .fx = fx, .Gu = Gu, .fu = fu
    };

    const lapackint_t N =
        get_assert_positive_integer_from_struct_field(prhs[5], "N", 5, MEX_FILE_ID);
    const mat_t Plq =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[5], "P", 5, MEX_FILE_ID);
    const mat_t K =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[5], "K", 5, MEX_FILE_ID);
    const mat_t Gf =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[5], "Gf", 5, MEX_FILE_ID);
    const mat_t ff =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[5], "ff", 5, MEX_FILE_ID);
    const rompc_t rompc = { .N = N, .P = Plq, .K = K, .Gf = Gf, .ff = ff };

    const mat_t U =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[6], "U", 6, MEX_FILE_ID);
    const mat_t xi =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[6], "xi", 6, MEX_FILE_ID);

    const lapackint_t n = sys.A.n;
    const lapackint_t m = sys.B.n;
    const lapackint_t r = U.n;

    assert_mex(xt.m == sys.A.n, "length of xt expected to be equal the number of columns of sys.A",
               MEX_FILE_ID);
    assert_mex(zt.m == N * m, "length of zt expected to be equal to N * m", MEX_FILE_ID);
    assert_mex(Rcheck.m == N * (n + m), "Rcheck expected to have dimension (N*(n+m) x N*(n+m))",
               MEX_FILE_ID);
    assert_mex(Qcheck.m == N * (n + m) + r,
               "Rcheck expected to have dimension (N*(n+m)+r x N*(n+m)+r)", MEX_FILE_ID);
    assert_sys(&sys, MEX_FILE_ID);
    assert_rompc(&rompc, sys.A.n, sys.B.n, MEX_FILE_ID);

    assert_mex(U.m == N * m, "number of rows of U expected to be equal to N * m", MEX_FILE_ID);
    assert_mex(xi.m == N * m, "length of xi expected to be equal to N * m, MEX_FILE_ID",
               MEX_FILE_ID);

    rompc_options_t opts = { .tol = 1e-6, .max_iter = 50 };

    if (nrhs > 7 && mxIsStruct(prhs[nrhs - 1])) {
        opts.tol = get_assert_positive_double_from_struct_field(prhs[nrhs - 1], "tol", nrhs - 1,
                                                                MEX_FILE_ID);

        opts.max_iter = get_assert_positive_integer_from_struct_field(prhs[nrhs - 1], "max_iter",
                                                                      nrhs - 1, MEX_FILE_ID);
    }

    mxArray* mxUtt = mxCreateDoubleMatrix(m, 1, mxREAL);
    mxArray* mxZtt = mxCreateDoubleMatrix(N * m, 1, mxREAL);

    mat_t ut_ret = mview_matrix_from_mx(mxUtt);
    mat_t zt_ret = mview_matrix_from_mx(mxZtt);

    const rompc_stats_t stats = condensed_rompc_step(&xt, &zt, &Rcheck, &Qcheck, &sys, &rompc, &U,
                                                     &xi, &ut_ret, &zt_ret, opts);

    plhs[0] = mxUtt;

    if (nlhs >= 2) {
        plhs[1] = mxZtt;
    } else {
        mxDestroyArray(mxZtt);
    }

    if (nlhs >= 3) {
        mxArray* mxFlag = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetPr(mxFlag) = (double)stats.flag;
        plhs[2] = mxFlag;
    }

    if (nlhs >= 4) {
        plhs[3] = create_rompc_stats(stats);
    }

    return;
}
