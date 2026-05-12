#include <stdio.h>
#include <stdlib.h>

#include "mex.h"

#include "algorithm/alg.h"
#include "algorithm/linalg.h"

#include "mx_helper.h"
#include "rompc_mex_helper.h"

static const char* const MEX_FILE_ID = "ROMPC:condensed_predictor_corrector";

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    if (nrhs < 5 || nrhs > 6) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "five or six input arguments expected.");
    }

    if (nlhs > 4) {
        mexErrMsgIdAndTxt(MEX_FILE_ID, "at most four output arguments expected.");
    }

    const mat_t x0 = cview_assert_nonempty_double_vector(prhs[0], 0, MEX_FILE_ID);
    const mat_t zt = cview_assert_nonempty_double_vector(prhs[1], 1, MEX_FILE_ID);

    const mat_t A =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[2], "A", 2, MEX_FILE_ID);
    const mat_t B =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[2], "B", 2, MEX_FILE_ID);
    const mat_t Q =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[2], "Q", 2, MEX_FILE_ID);
    const mat_t R =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[2], "R", 2, MEX_FILE_ID);
    const mat_t Gx =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[2], "Gx", 2, MEX_FILE_ID);
    const mat_t fx =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[2], "fx", 2, MEX_FILE_ID);
    const mat_t Gu =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[2], "Gu", 2, MEX_FILE_ID);
    const mat_t fu =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[2], "fu", 2, MEX_FILE_ID);

    const lapackint_t N =
        get_assert_positive_integer_from_struct_field(prhs[3], "N", 3, MEX_FILE_ID);
    const mat_t P =
        cview_assert_nonempty_square_double_matrix_from_struct_field(prhs[3], "P", 3, MEX_FILE_ID);
    const mat_t K =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[3], "K", 3, MEX_FILE_ID);
    const mat_t Gf =
        cview_assert_nonempty_double_matrix_from_struct_field(prhs[3], "Gf", 3, MEX_FILE_ID);
    const mat_t ff =
        cview_assert_nonempty_double_vector_from_struct_field(prhs[3], "ff", 3, MEX_FILE_ID);

    const mat_t eqQ1 =
        cview_assert_double_matrix_from_struct_struct_field(prhs[3], "eq", "Q1", 3, MEX_FILE_ID);
    const mat_t eqQ2 =
        cview_assert_double_matrix_from_struct_struct_field(prhs[3], "eq", "Q2", 3, MEX_FILE_ID);
    const mat_t eqR =
        cview_assert_double_matrix_from_struct_struct_field(prhs[3], "eq", "R", 3, MEX_FILE_ID);

    const mat_t p_init = cview_assert_nonempty_double_vector(prhs[4], 4, MEX_FILE_ID);

    const sys_t sys = { .A = A, .B = B, .Q = Q, .R = R, .Gx = Gx, .fx = fx, .Gu = Gu, .fu = fu };
    const rompc_t rompc = { .N = N, .P = P, .K = K, .Gf = Gf, .ff = ff };
    const eqc_t eqc = { .Q1 = eqQ1, .Q2 = eqQ2, .R = eqR };

    const lapackint_t s = p_init.m;
    const lapackint_t Nnm = N * (Q.m + R.m);

    assert_sys(&sys, MEX_FILE_ID);
    assert_rompc(&rompc, sys.A.n, sys.B.n, MEX_FILE_ID);
    assert_eqc(&eqc, Nnm, s, MEX_FILE_ID);

    rompc_options_t opts = { .tol = 1e-6, .max_iter = 50 };
    if (nrhs == 6) {
        assert_mex(mxIsStruct(prhs[5]), "options must be a struct", MEX_FILE_ID);

        opts.tol = get_assert_positive_double_from_struct_field(prhs[5], "tol", 5, MEX_FILE_ID);

        opts.max_iter =
            get_assert_positive_integer_from_struct_field(prhs[5], "max_iter", 5, MEX_FILE_ID);
    }

    mxArray* mxW_opt = mxCreateDoubleMatrix(Nnm, 1, mxREAL);
    mxArray* mxBeta_opt = mxCreateDoubleMatrix(s, 1, mxREAL);

    mat_t w_opt = mview_matrix_from_mx(mxW_opt);
    mat_t beta_opt = mview_matrix_from_mx(mxBeta_opt);

    const rompc_stats_t stats = condensed_predictor_corrector(&x0, &zt, &sys, &rompc, &eqc, &p_init,
                                                              &w_opt, &beta_opt, opts);

    plhs[0] = mxW_opt;

    if (nlhs >= 2) {
        plhs[1] = mxBeta_opt;
    } else {
        mxDestroyArray(mxBeta_opt);
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
