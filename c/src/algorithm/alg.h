#pragma once

#include "linalg.h"

// The function signatures follow the convention of BLAS/lapack in the sense, that the inputs are
// given first and the outputs last, such that the ordering is:
//      1. input arguments (that are not changed by the function)
//      2. input/output arguments
//      3. output arguments
//      4. working buffers

// work: buffer with a capacity of at least  double values
void update_qr(const mat_t* Hcheck, const mat_t* Qcheck2, mat_t* Rt, mat_t* Qhat2,
               sizedbuffer_t* work);

void forward_backward_solve(const mat_t* G, mat_t* bx);

double multiplier_line_search(const mat_t* nu, const mat_t* Delta_nu);

// work: buffer with a capacity of at least max(n, m) double values
void get_condensed_adjoint(const mat_t* lambda, const mat_t* W, const mat_t* Gx, const mat_t* Gu,
                           const mat_t* Gf, const lapackint_t N, mat_t* theta, sizedbuffer_t* work);

// work: buffer with a capacity of at least max(n, m) double values
void get_condensed_gradient(const mat_t* w, const mat_t* W, const mat_t* Q, const mat_t* R,
                            const mat_t* P, const lapackint_t N, mat_t* eta, sizedbuffer_t* work);

// work: buffer with a capacity of at least max(n, m) * s double values
void get_condensed_hessian(const mat_t* W, const mat_t* Q, const mat_t* R, const mat_t* P,
                           const lapackint_t N, mat_t* H, sizedbuffer_t* work);

void get_condensed_primal_feasibility(const mat_t* w, const mat_t* W, const mat_t* Gx,
                                      const mat_t* Gu, const mat_t* Gf, const mat_t* fx,
                                      const mat_t* fu, const mat_t* ff, const lapackint_t N,
                                      mat_t* delta);

// work: buffer with a capacity of at least max(qx, qu, qf) * 2 * s double values
void get_condensed_quadratic_constraints(const mat_t* W, const mat_t* Gx, const mat_t* Gu,
                                         const mat_t* Gf, const mat_t* y, const mat_t* lambda,
                                         const lapackint_t N, mat_t* D, sizedbuffer_t* work);

typedef struct {
    mat_t A;
    mat_t B;
    mat_t Q;
    mat_t R;
    mat_t Gx;
    mat_t fx;
    mat_t Gu;
    mat_t fu;
} sys_t;

typedef struct {
    lapackint_t N;
    mat_t P;
    mat_t K;
    mat_t Gf;
    mat_t ff;
} rompc_t;

typedef struct {
    mat_t Q1;
    mat_t Q2;
    mat_t R;
} eqc_t;

typedef struct {
    double tol;
    lapackint_t max_iter;
} rompc_options_t;

typedef enum {
    ROMPC_FLAG_UNEXPECTED_ERROR = -10,
    ROMPC_FLAG_CHOLESKY_FAILED = -2,
    ROMPC_FLAG_ITER_LIMIT = -1,
    ROMPC_FLAG_CONVERGED = 1,
} rompc_flag_t;

typedef struct {
    rompc_flag_t flag;
    lapackint_t iters;
    double gam_rp;
    double gam_rd;
    double duration_qr;
    double duration_cpc;
} rompc_stats_t;

rompc_stats_t condensed_predictor_corrector(const mat_t* x0, const mat_t* zt, const sys_t* sys,
                                            const rompc_t* rompc, const eqc_t* eqc,
                                            const mat_t* p_init, mat_t* w_opt, mat_t* beta_opt,
                                            rompc_options_t opts);

rompc_stats_t condensed_rompc_step(const mat_t* xt, const mat_t* zt, const mat_t* Rcheck,
                                   const mat_t* Qcheck, const sys_t* sys, const rompc_t* rompc,
                                   const mat_t* U, const mat_t* xi, mat_t* ut_ret, mat_t* zt_ret,
                                   rompc_options_t opts);
