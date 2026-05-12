#pragma once

#include "mx_helper.h"
#include "algorithm/alg.h"

static void assert_sys(const sys_t* sys, const char* const mex_file_id) {
    const lapackint_t nx = sys->A.m;
    const lapackint_t nu = sys->B.n;

    assert_mex(sys->B.m == nx, "sys.A and sys.B expected to have the same number of rows.",
               mex_file_id);
    assert_mex(sys->Q.n == nx, "sys.Q expected to have same number of columns as sys.A.",
               mex_file_id);
    assert_mex(sys->R.n == nu, "sys.R expected to have same number of columns as sys.B.",
               mex_file_id);

    assert_mex(sys->Gx.m == sys->fx.m, "sys.Gx and sys.fx expected to have same number of rows.",
               mex_file_id);
    assert_mex(sys->Gx.n == nx, "sys.Gx expected to have the same number of columns as sys.A",
               mex_file_id);

    assert_mex(sys->Gu.m == sys->fu.m, "sys.Gu and sys.fu expected to have same number of rows.",
               mex_file_id);
    assert_mex(sys->Gu.n == nu, "sys.Gu expected to have the same number of columns as sys.B",
               mex_file_id);
}

static void assert_rompc(const rompc_t* rompc, lapackint_t nx, lapackint_t nu,
                         const char* const mex_file_id) {

    assert_mex(rompc->P.n == nx, "rompc.P expected to have the same number of columns as sys.A",
               mex_file_id);
    assert_mex(rompc->Gf.m == rompc->ff.m,
               "rompc.Gf and rompc.ff expected to have same number of rows.", mex_file_id);
    assert_mex(rompc->Gf.n == nx, "rompc.Gf expected to have the same number of columns as sys.A",
               mex_file_id);
    assert_mex(rompc->K.n == nx, "rompc.K expected to have the same number of columns as sys.A",
               mex_file_id);
    assert_mex(rompc->K.m == nu,
               "number of rows of rompc.K expected to be equal to the number of columns of sys.B",
               mex_file_id);
}

static void assert_eqc(const eqc_t* eqc, lapackint_t Nnm, lapackint_t s,
                       const char* const mex_file_id) {
    assert_mex(eqc->R.m == eqc->R.n, "rompc.eq.R expected to be square", mex_file_id);
    assert_mex(eqc->R.m == Nnm, "number of rows of rompc.eq.R expected to be N * (n + m)",
               mex_file_id);
    assert_mex(eqc->Q1.m == eqc->Q2.m,
               "rompc.eq.Q1 and rompc.eq.Q2 expected to have the same number of rows", mex_file_id);
    assert_mex(eqc->Q1.n == Nnm, "number of columns of rompc.eq.Q1 expected to be N * (n + m)",
               mex_file_id);
    assert_mex(eqc->Q2.n == s,
               "number of columns of rompc.eq.Q2 expected to be equal the size of p_init (s)",
               mex_file_id);
}

#define NB_STAT_FIELDS 7
static const char* stat_fields[NB_STAT_FIELDS] = { "flag",   "msg",         "iters",       "gam_rp",
                                                   "gam_rd", "duration_qr", "duration_cpc" };

static mxArray* create_rompc_stats(rompc_stats_t stats) {
    mxArray* mxFlag = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxFlag) = (double)stats.flag;

    const char* msg;

    if (stats.flag == ROMPC_FLAG_CONVERGED) {
        msg = "successfully converged";
    } else if (stats.flag == ROMPC_FLAG_ITER_LIMIT) {
        msg = "not converged: iteration limit reached";
    } else if (stats.flag == ROMPC_FLAG_CHOLESKY_FAILED) {
        msg = "not converged: cholesky factorization failed";
    } else {
        msg = "unexpected error";
    }

    mxArray* mxMsg = mxCreateString(msg);

    mxArray* mxIters = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxIters) = (double)stats.iters;

    mxArray* mxGamRp = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxGamRp) = stats.gam_rp;

    mxArray* mxGamRd = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxGamRd) = stats.gam_rd;

    mxArray* mxDurationQr = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxDurationQr) = stats.duration_qr;

    mxArray* mxDurationCpc = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(mxDurationCpc) = stats.duration_cpc;

    mxArray* mx = mxCreateStructMatrix(1, 1, NB_STAT_FIELDS, stat_fields);

    mxSetFieldByNumber(mx, 0, 0, mxFlag);
    mxSetFieldByNumber(mx, 0, 1, mxMsg);
    mxSetFieldByNumber(mx, 0, 2, mxIters);
    mxSetFieldByNumber(mx, 0, 3, mxGamRp);
    mxSetFieldByNumber(mx, 0, 4, mxGamRd);
    mxSetFieldByNumber(mx, 0, 5, mxDurationQr);
    mxSetFieldByNumber(mx, 0, 6, mxDurationCpc);

    return mx;
}
