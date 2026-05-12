#include "alg.h"

#include "utils.h"

void update_qr(const mat_t* Hcheck, const mat_t* Qcheck2, mat_t* Rt_ret, mat_t* Qhat2,
               sizedbuffer_t* work) {

    ASSERT(Hcheck->n == Qcheck2->n);

    const lapackint_t d = Hcheck->n;

    mat_t Rt = mat_from_sizedbuffer(work, Hcheck->m, Hcheck->n);
    mat_copy(Hcheck, &Rt);

    mat_set_zero(Qhat2);
    MAT_ENTRY(Qhat2, Qhat2->m - 1, 0) = 1.0;

    mat_copy_2(Qcheck2, mview_north_east_dim(Qhat2, Qcheck2->m, Qcheck2->n));

    for (lapackint_t k = 0; k < d; ++k) {
        const givens_t ck_sk = givens(MAT_ENTRY(&Rt, k, k), MAT_ENTRY(&Rt, k + 1, k));

        // for (lapackint_t j = k; j < d; ++j) {
        //     const double tau1 = MAT_ENTRY(&Rt, k, j);
        //     const double tau2 = MAT_ENTRY(&Rt, k + 1, j);
        //     MAT_ENTRY(&Rt, k, j) = ck_sk.c * tau1 - ck_sk.s * tau2;
        //     MAT_ENTRY(&Rt, k + 1, j) = ck_sk.s * tau1 + ck_sk.c * tau2;
        // }
        const lapackint_t r = d - k;
        drot(&r, MAT_ENTRY_PTR(&Rt, k + 1, k), &Rt.ld, MAT_ENTRY_PTR(&Rt, k, k), &Rt.ld, &ck_sk.c,
             &ck_sk.s);

        // for (lapackint_t j = 0; j < Qhat2->m; ++j) {
        //     const double t1 = MAT_ENTRY(Qhat2, j, k);
        //     const double t2 = MAT_ENTRY(Qhat2, j, k + 1);
        //     MAT_ENTRY(Qhat2, j, k) = ck_sk.c * t1 - ck_sk.s * t2;
        //     MAT_ENTRY(Qhat2, j, k + 1) = ck_sk.s * t1 + ck_sk.c * t2;
        // }
        drot(&Qhat2->m, MAT_ENTRY_PTR(Qhat2, 0, k + 1), &ione, MAT_ENTRY_PTR(Qhat2, 0, k), &ione,
             &ck_sk.c, &ck_sk.s);
    }

    mat_copy_1(cview_north(&Rt, Rt_ret->m), Rt_ret);
}

void forward_backward_solve(const mat_t* G, mat_t* bx) {
    ASSERT(G->m == G->n);
    ASSERT(G->n == bx->m);
    ASSERT(bx->n == 1);

    dtrsv("u", "t", "n", &G->m, G->p, &G->ld, bx->p, &ione);
    dtrsv("u", "n", "n", &G->m, G->p, &G->ld, bx->p, &ione);
}

double multiplier_line_search(const mat_t* nu, const mat_t* Delta_nu) {
    ASSERT(nu->m == Delta_nu->m);
    ASSERT(nu->n == 1);
    ASSERT(Delta_nu->n == 1);

    double alpha = 1.0;

    const double* pnu = nu->p;
    const double* pDnu = Delta_nu->p;
    const double* const pend = pnu + nu->m;

    while (pnu < pend) {
        if (*pDnu < 0) {
            const double mu = -*pnu / *pDnu;

            if (mu < alpha) {
                alpha = mu;
            }
        }

        ++pnu;
        ++pDnu;
    }

    return alpha;
}

double multiplier_line_search_tau(double tau, const mat_t* nu, const mat_t* Delta_nu) {
    ASSERT(nu->m == Delta_nu->m);
    ASSERT(nu->n == 1);
    ASSERT(Delta_nu->n == 1);

    double alpha = 1.0;

    const double* pnu = nu->p;
    const double* pDnu = Delta_nu->p;
    const double* const pend = pnu + nu->m;

    while (pnu < pend) {
        if (*pDnu < 0) {
            const double mu = -tau * *pnu / *pDnu;

            if (mu < alpha) {
                alpha = mu;
            }
        }

        ++pnu;
        ++pDnu;
    }

    return alpha;
}

void get_condensed_adjoint(const mat_t* lambda, const mat_t* W, const mat_t* Gx, const mat_t* Gu,
                           const mat_t* Gf, lapackint_t N, mat_t* theta, sizedbuffer_t* work) {

    const lapackint_t qx = Gx->m;
    const lapackint_t n = Gx->n;
    const lapackint_t qu = Gu->m;
    const lapackint_t m = Gu->n;
    const lapackint_t qf = Gf->m;

    mat_t tmpx = mat_from_sizedbuffer(work, n, 1);
    mat_t tmpu = mat_from_sizedbuffer(work, m, 1);

    mat_set_zero(theta);

    for (lapackint_t k = 0; k < N; ++k) {

        const mat_t Wu_k = cview_rows_idx_dim(W, k * (m + n), m);
        const mat_t lambda_u_k = cview_rows_idx_dim(lambda, k * (qu + qx), qu);

        matT_vec_mult(Gu, &lambda_u_k, &tmpu);
        matT_vec_mult_add(&Wu_k, &tmpu, theta);

        if (k == N - 1) {
            const mat_t Wx_k = cview_rows_idx_dim(W, k * (m + n) + m, n);
            const mat_t lambda_x_k = cview_rows_idx_dim(lambda, k * (qu + qx) + qu, qf);

            matT_vec_mult(Gf, &lambda_x_k, &tmpx);
            matT_vec_mult_add(&Wx_k, &tmpx, theta);
        } else {
            const mat_t Wx_k = cview_rows_idx_dim(W, k * (m + n) + m, n);
            const mat_t lambda_x_k = cview_rows_idx_dim(lambda, k * (qu + qx) + qu, qx);

            matT_vec_mult(Gx, &lambda_x_k, &tmpx);
            matT_vec_mult_add(&Wx_k, &tmpx, theta);
        }
    }
}

void get_condensed_gradient(const mat_t* w, const mat_t* W, const mat_t* Q, const mat_t* R,
                            const mat_t* P, lapackint_t N, mat_t* eta, sizedbuffer_t* work) {
    const lapackint_t n = Q->m;
    const lapackint_t m = R->m;
    const lapackint_t s = W->n;

    mat_t tmpx = mat_from_sizedbuffer(work, n, 1);
    mat_t tmpu = mat_from_sizedbuffer(work, m, 1);

    mat_set_zero(eta);

    for (lapackint_t k = 0; k < N; ++k) {

        const mat_t Wu_k = cview_rows_idx_dim(W, k * (m + n), m);
        const mat_t w_u_k = cview_rows_idx_dim(w, k * (m + n), m);

        matT_vec_mult(R, &w_u_k, &tmpu);
        matT_vec_mult_add(&Wu_k, &tmpu, eta);

        const mat_t Wx_k = cview_rows_idx_dim(W, k * (m + n) + m, n);
        const mat_t w_x_k = cview_rows_idx_dim(w, k * (m + n) + m, n);

        if (k == N - 1) {
            matT_vec_mult(P, &w_x_k, &tmpx);
        } else {
            matT_vec_mult(Q, &w_x_k, &tmpx);
        }

        matT_vec_mult_add(&Wx_k, &tmpx, eta);
    }
}

void get_condensed_hessian(const mat_t* W, const mat_t* Q, const mat_t* R, const mat_t* P,
                           lapackint_t N, mat_t* H, sizedbuffer_t* work) {
    const lapackint_t n = Q->m;
    const lapackint_t m = R->m;
    const lapackint_t s = W->n;

    mat_t tmpx = mat_from_sizedbuffer(work, n, s);
    mat_t tmpu = mat_from_sizedbuffer(work, m, s);

    mat_set_zero(H);

    for (lapackint_t k = 0; k < N; ++k) {
        const mat_t Wu_k = cview_rows_idx_dim(W, k * (m + n), m);

        mat_mult(R, &Wu_k, &tmpu);
        matT_mat_mult_add(&Wu_k, &tmpu, H);

        const mat_t Wx_k = cview_rows_idx_dim(W, k * (m + n) + m, n);

        if (k == N - 1) {
            mat_mult(P, &Wx_k, &tmpx);
        } else {
            mat_mult(Q, &Wx_k, &tmpx);
        }

        matT_mat_mult_add(&Wx_k, &tmpx, H);
    }
}

void get_condensed_primal_feasibility(const mat_t* w, const mat_t* W, const mat_t* Gx,
                                      const mat_t* Gu, const mat_t* Gf, const mat_t* fx,
                                      const mat_t* fu, const mat_t* ff, lapackint_t N,
                                      mat_t* delta) {
    const lapackint_t qx = Gx->m;
    const lapackint_t n = Gx->n;
    const lapackint_t qu = Gu->m;
    const lapackint_t m = Gu->n;
    const lapackint_t qf = Gf->m;

    for (lapackint_t k = 0; k < N; ++k) {

        const mat_t w_u_k = cview_rows_idx_dim(w, k * (m + n), m);
        mat_t delta_u_k = mview_rows_idx_dim(delta, k * (qu + qx), qu);

        mat_copy(fu, &delta_u_k);
        mat_vec_mult_minus(Gu, &w_u_k, &delta_u_k);

        if (k == N - 1) {
            const mat_t w_x_k = cview_rows_idx_dim(w, k * (m + n) + m, n);
            mat_t delta_f_k = mview_rows_idx_dim(delta, k * (qu + qx) + qu, qf);

            mat_copy(ff, &delta_f_k);
            mat_vec_mult_minus(Gf, &w_x_k, &delta_f_k);
        } else {
            const mat_t w_x_k = cview_rows_idx_dim(w, k * (m + n) + m, n);
            mat_t delta_x_k = mview_rows_idx_dim(delta, k * (qu + qx) + qu, qx);

            mat_copy(fx, &delta_x_k);
            mat_vec_mult_minus(Gx, &w_x_k, &delta_x_k);
        }
    }
}

void get_condensed_quadratic_constraints(const mat_t* W, const mat_t* Gx, const mat_t* Gu,
                                         const mat_t* Gf, const mat_t* y, const mat_t* lambda,
                                         lapackint_t N, mat_t* D, sizedbuffer_t* work) {
    const lapackint_t qx = Gx->m;
    const lapackint_t n = Gx->n;
    const lapackint_t qu = Gu->m;
    const lapackint_t m = Gu->n;
    const lapackint_t qf = Gf->m;
    const lapackint_t s = W->n;

    lapackint_t dim = (qx > qu) ? qx : qu;
    dim = (dim > qf) ? dim : qf;

    ASSERT(work->n >= dim * s * 2);

    mat_t tmpGfWxL = mat_from_raw_buffer(work->p, qf, s);
    mat_t tmpGfWxR = mat_from_raw_buffer(work->p + dim * s, qf, s);

    mat_t tmpGxWxL = mat_from_raw_buffer(work->p, qx, s);
    mat_t tmpGxWxR = mat_from_raw_buffer(work->p + dim * s, qx, s);

    mat_t tmpGuWuL = mat_from_raw_buffer(work->p, qu, s);
    mat_t tmpGuWuR = mat_from_raw_buffer(work->p + dim * s, qu, s);

    mat_set_zero(D);

    for (lapackint_t k = 0; k < N; ++k) {
        const mat_t Wu_k = cview_rows_idx_dim(W, k * (m + n), m);
        const mat_t lambda_u_k = cview_rows_idx_dim(lambda, k * (qx + qu), qu);
        const mat_t y_u_k = cview_rows_idx_dim(y, k * (qx + qu), qu);

        mat_mult(Gu, &Wu_k, &tmpGuWuL);
        mat_copy(&tmpGuWuL, &tmpGuWuR);

        mat_scale_recip_diag(&y_u_k, &tmpGuWuL);
        mat_scale_diag(&lambda_u_k, &tmpGuWuR);

        matT_mat_mult_add(&tmpGuWuL, &tmpGuWuR, D);

        const mat_t Wx_k = cview_rows_idx_dim(W, k * (m + n) + m, n);

        if (k == N - 1) {
            const mat_t lambda_x_k = cview_rows_idx_dim(lambda, k * (qx + qu) + qu, qf);
            const mat_t y_x_k = cview_rows_idx_dim(y, k * (qx + qu) + qu, qf);

            mat_mult(Gf, &Wx_k, &tmpGfWxL);
            mat_copy(&tmpGfWxL, &tmpGfWxR);

            mat_scale_recip_diag(&y_x_k, &tmpGfWxL);
            mat_scale_diag(&lambda_x_k, &tmpGfWxR);

            matT_mat_mult_add(&tmpGfWxL, &tmpGfWxR, D);
        } else {
            const mat_t lambda_x_k = cview_rows_idx_dim(lambda, k * (qx + qu) + qu, qx);
            const mat_t y_x_k = cview_rows_idx_dim(y, k * (qx + qu) + qu, qx);

            mat_mult(Gx, &Wx_k, &tmpGxWxL);
            mat_copy(&tmpGxWxL, &tmpGxWxR);

            mat_scale_recip_diag(&y_x_k, &tmpGxWxL);
            mat_scale_diag(&lambda_x_k, &tmpGxWxR);

            matT_mat_mult_add(&tmpGxWxL, &tmpGxWxR, D);
        }
    }
}

rompc_stats_t condensed_predictor_corrector(const mat_t* x0, const mat_t* zt, const sys_t* sys,
                                            const rompc_t* rompc, const eqc_t* eqc,
                                            const mat_t* p_init, mat_t* w_opt, mat_t* beta_opt,
                                            rompc_options_t opts) {

    rompc_stats_t stats = { .flag = ROMPC_FLAG_UNEXPECTED_ERROR,
                            .iters = 0,
                            .gam_rp = NAN,
                            .gam_rd = NAN,
                            .duration_qr = NAN,
                            .duration_cpc = NAN };

    const double tol = opts.tol;
    const int max_iter = opts.max_iter;

    const lapackint_t n = sys->Q.m;
    const lapackint_t m = sys->R.n;
    const lapackint_t N = rompc->N;

    const lapackint_t q = N * sys->Gu.m + (N - 1) * sys->Gx.m + rompc->Gf.m;

    const lapackint_t qx = sys->Gx.m;
    const lapackint_t qu = sys->Gu.m;
    const lapackint_t qf = rompc->Gf.m;

    const lapackint_t s = eqc->Q2.n;

    lapackint_t buf_dim = (qx > qu) ? qx : qu;
    buf_dim = (buf_dim > qf) ? buf_dim : qf;
    buf_dim = buf_dim * s * 2;
    buf_dim = (buf_dim > m * s) ? buf_dim : m * s;
    buf_dim = (buf_dim > n * s) ? buf_dim : n * s;
    buf_dim = (buf_dim > N * (n + m)) ? buf_dim : N * (n + m);

    sizedbuffer_t work = sizedbuffer_create(buf_dim);

    const mat_t Z = cview_north(&eqc->Q1, N * (n + m));
    const mat_t W = cview_north(&eqc->Q2, N * (n + m));

    mat_t g = mat_create_zeros(N * (n + m), 1);

    mat_t Ax0 = mat_from_sizedbuffer(&work, n, 1);
    mat_vec_mult(&sys->A, x0, &Ax0);
    mat_copy_2(&Ax0, mview_north(&g, n));

    mat_t Kx0 = mat_from_sizedbuffer(&work, m, 1);
    mat_vec_neg_mult(&rompc->K, x0, &Kx0);
    mat_copy_2(&Kx0, mview_idx_dim(&g, N * n, 0, m, 1));
    vec_minus_2(zt, mview_rows_idx_dim(&g, N * n, N * m));

    g.p[0] /= eqc->R.p[0];

    for (lapackint_t i = 1; i < N * (n + m); ++i) {
        const double v = ddot(&i, MAT_ENTRY_PTR(&eqc->R, 0, i), &ione, g.p, &ione);

        g.p[i] = (g.p[i] - v) / MAT_ENTRY(&eqc->R, i, i);
    }

    mat_t H = mat_create_zeros(s, s);
    get_condensed_hessian(&W, &sys->Q, &sys->R, &rompc->P, N, &H, &work);

    mat_t p = mat_create_copy(p_init);

    mat_t y_lambda = mat_create_ones(2 * q, 1);
    mat_t y = mview_rows_idx_dim(&y_lambda, 0, q);
    mat_t lambda = mview_rows_idx_dim(&y_lambda, q, q);

    // w is just an alias for w_opt
    mat_t* w = w_opt;

    mat_t eta = mat_create(s, 1);
    mat_t delta = mat_create(q, 1);

    mat_t rd = mat_create(s, 1);
    mat_t rp = mat_create(q, 1);

    mat_t L = mat_create(s, s);

    mat_t phat = mat_create(s, 1);
    mat_t Dp_aff = mat_create(s, 1);
    mat_t w_aff = mat_create(N * (n + m), 1);

    mat_t Dyl_aff = mat_create(2 * q, 1);
    mat_t Dy_aff = mview_rows_idx_dim(&Dyl_aff, 0, q);
    mat_t Dlambda_aff = mview_rows_idx_dim(&Dyl_aff, q, q);

    mat_t kappa = mat_create(q, 1);

    mat_t phat_tmp = mat_create(s, 1);

    lapackint_t k = 0;

    while (1) {
        mat_vec_mult(&Z, &g, w);
        mat_vec_mult_add(&W, &p, w);

        get_condensed_gradient(w, &W, &sys->Q, &sys->R, &rompc->P, rompc->N, &eta, &work);
        get_condensed_primal_feasibility(w, &W, &sys->Gx, &sys->Gu, &rompc->Gf, &sys->fx, &sys->fu,
                                         &rompc->ff, rompc->N, &delta);

        mat_set_zero(&rd);
        get_condensed_adjoint(&lambda, &W, &sys->Gx, &sys->Gu, &rompc->Gf, N, &rd, &work);
        vec_add(&eta, &rd);

        mat_copy(&delta, &rp);
        vec_add(&y, &rp);
        vec_scalar_mult(-1.0, &rp);

        const double gam_rd = vec_inf_norm(&rd);
        const double gam_rp = vec_inf_norm(&rp);

        stats.iters = k;
        stats.gam_rd = gam_rd;
        stats.gam_rp = gam_rp;

        if ((gam_rd < tol) && (gam_rp < tol)) {
            stats.flag = ROMPC_FLAG_CONVERGED;
            break;
        } else if (k == max_iter) {
            stats.flag = ROMPC_FLAG_ITER_LIMIT;
            break;
        }

        get_condensed_quadratic_constraints(&W, &sys->Gx, &sys->Gu, &rompc->Gf, &y, &lambda, N, &L,
                                            &work);
        mat_add_cont(&H, &L);

        const lapackint_t chol_flag = chol(&L);

        if (chol_flag != 0) {
            stats.flag = ROMPC_FLAG_CHOLESKY_FAILED;
            break;
        }

        for (lapackint_t i = 0; i < delta.m; ++i) {
            delta.p[i] = (lambda.p[i] * delta.p[i]) / y.p[i] + lambda.p[i];
        }

        mat_set_zero(&phat);
        get_condensed_adjoint(&delta, &W, &sys->Gx, &sys->Gu, &rompc->Gf, N, &phat, &work);
        vec_add(&eta, &phat);
        vec_scalar_mult(-1.0, &phat);

        mat_copy(&phat, &Dp_aff);
        forward_backward_solve(&L, &Dp_aff);

        mat_copy(w, &w_aff);
        mat_vec_mult_add(&W, &Dp_aff, &w_aff);

        mat_set_zero(&Dlambda_aff);
        get_condensed_primal_feasibility(&w_aff, &W, &sys->Gx, &sys->Gu, &rompc->Gf, &sys->fx,
                                         &sys->fu, &rompc->ff, N, &Dlambda_aff);

        mat_copy(&Dlambda_aff, &Dy_aff);
        vec_add(&y, &Dy_aff);
        vec_scalar_mult(-1.0, &Dy_aff);

        for (lapackint_t i = 0; i < Dlambda_aff.m; ++i) {
            Dlambda_aff.p[i] = (lambda.p[i] * Dlambda_aff.p[i]) / y.p[i];
        }

        const double mu = vec_dot(&y, &lambda) / q;
        const double alpha_aff = multiplier_line_search(&y_lambda, &Dyl_aff);

        double mu_aff = 0;

        for (lapackint_t i = 0; i < Dy_aff.m; ++i) {
            mu_aff +=
                (y.p[i] + alpha_aff * Dy_aff.p[i]) * (lambda.p[i] + alpha_aff * Dlambda_aff.p[i]);
        }
        mu_aff /= q;

        double sigma = (mu_aff / mu);
        sigma = sigma * sigma * sigma;

        const double sigma_mu = sigma * mu;

        for (lapackint_t i = 0; i < kappa.m; ++i) {
            kappa.p[i] = (Dlambda_aff.p[i] * Dy_aff.p[i] - sigma_mu) / y.p[i];
        }

        get_condensed_adjoint(&kappa, &W, &sys->Gx, &sys->Gu, &rompc->Gf, N, &Dp_aff, &work);
        vec_add(&phat, &Dp_aff);

        forward_backward_solve(&L, &Dp_aff);

        mat_copy(w, &w_aff);
        mat_vec_mult_add(&W, &Dp_aff, &w_aff);

        mat_set_zero(&Dlambda_aff);
        get_condensed_primal_feasibility(&w_aff, &W, &sys->Gx, &sys->Gu, &rompc->Gf, &sys->fx,
                                         &sys->fu, &rompc->ff, N, &Dlambda_aff);

        mat_copy(&Dlambda_aff, &Dy_aff);
        vec_add(&y, &Dy_aff);
        vec_scalar_mult(-1.0, &Dy_aff);

        for (lapackint_t i = 0; i < Dlambda_aff.m; ++i) {
            Dlambda_aff.p[i] = (lambda.p[i] * Dlambda_aff.p[i]) / y.p[i] - kappa.p[i];
        }

        const double tau = 0.5;
        const double alpha_primal = multiplier_line_search_tau(tau, &y, &Dy_aff);
        const double alpha_dual = multiplier_line_search_tau(tau, &lambda, &Dlambda_aff);
        const double alpha = (alpha_primal < alpha_dual) ? alpha_primal : alpha_dual;

        vec_add_factor(alpha, &Dp_aff, &p);
        vec_add_factor(alpha, &Dy_aff, &y);
        vec_add_factor(alpha, &Dlambda_aff, &lambda);

        k += 1;
    }

    sizedbuffer_free(&work);
    // printf("niters: %d\n", k);
    // w_opt, beta_opt

    mat_vec_mult_1(cview_south_idx(&eqc->Q1, N * (n + m)), &g, beta_opt);
    mat_vec_mult_add_1(cview_south_idx(&eqc->Q2, N * (n + m)), &p, beta_opt);

    mat_free(&phat_tmp);
    mat_free(&Dyl_aff);
    mat_free(&w_aff);
    mat_free(&Dp_aff);
    mat_free(&phat);
    mat_free(&L);
    mat_free(&rd);
    mat_free(&rp);
    mat_free(&delta);
    mat_free(&eta);
    mat_free(&y_lambda);
    mat_free(&p);
    mat_free(&H);
    mat_free(&g);

    return stats;
}

rompc_stats_t condensed_rompc_step(const mat_t* xt, const mat_t* zt, const mat_t* Rcheck,
                                   const mat_t* Qcheck, const sys_t* sys, const rompc_t* rompc,
                                   const mat_t* U, const mat_t* xi, mat_t* ut_ret, mat_t* zt_ret,
                                   rompc_options_t opts) {

    const lapackint_t N = rompc->N;
    const lapackint_t n = sys->Q.m;
    const lapackint_t m = sys->R.n;
    const lapackint_t Nnm = N * (n + m);
    const lapackint_t r = U->n;
    const lapackint_t s = r + 1;

    mat_t z0 = mat_create_copy(zt);
    vec_add_factor(-1.0, xi, &z0);
    vec_scalar_mult(-1.0, &z0);

    mat_t Hcheck = mat_create(N * m + 1, N * m);
    mat_t z0T = mat_from_raw_buffer(z0.p, 1, z0.m); // basically a transposed view
    mat_copy_2(&z0T, mview_north(&Hcheck, 1));

    mat_copy_3(cview_south_east_idx(Rcheck, N * n, N * n), mview_south_idx(&Hcheck, 1));

    mat_t Q = mat_create_zeros(Nnm + s, Nnm + s);

    mat_copy_3(cview_west(Qcheck, N * n), mview_idx_dim(&Q, 0, 0, Q.m - 1, N * n));
    mat_copy_3(cview_east_idx(Qcheck, Nnm), mview_idx_dim(&Q, 0, Nnm + 1, Q.m - 1, r));

    mat_t Qhat2 = mview_columns_idx_dim(&Q, N * n, N * m + 1);

    mat_t R = mat_create_zeros(Nnm, Nnm);
    mat_copy_3(cview_north(Rcheck, N * n), mview_north(&R, N * n));

    mat_t Rtilde = mview_south_east_dim(&R, N * m, N * m);

    const mat_t Qcheck2 = cview_columns_idx_dim(Qcheck, N * n, N * m);

    sizedbuffer_t work = sizedbuffer_create(Hcheck.m * Hcheck.n);

    struct timespec ts_started, ts_finished;
    clock_gettime(CLOCK_MONOTONIC, &ts_started);
    update_qr(&Hcheck, &Qcheck2, &Rtilde, &Qhat2, &work);
    clock_gettime(CLOCK_MONOTONIC, &ts_finished);

    const double duration_qr = get_duration(ts_started, ts_finished);

    sizedbuffer_free(&work);

    const mat_t Q1 = cview_columns_idx_dim(&Q, 0, Nnm);
    const mat_t Q2 = cview_columns_idx_dim(&Q, Nnm, Q.n - Nnm);

    // mat_dbg("Q", Q);
    // mat_dbg("R", R);
    // mat_dbg("Q1", Q1);
    // mat_dbg("Q2", Q2);

    const eqc_t eqc = { .Q1 = Q1, .Q2 = Q2, .R = R };

    mat_t p_init = mat_create_zeros(s, 1);
    mat_t w_opt = mat_create_zeros(Nnm, 1);
    mat_t beta_opt = mat_create_zeros(s, 1);


    clock_gettime(CLOCK_MONOTONIC, &ts_started);
    rompc_stats_t stats =
        condensed_predictor_corrector(xt, zt, sys, rompc, &eqc, &p_init, &w_opt, &beta_opt, opts);
    clock_gettime(CLOCK_MONOTONIC, &ts_finished);

    const double duration_cpc = get_duration(ts_started, ts_finished);
    stats.duration_qr = duration_qr;
    stats.duration_cpc = duration_cpc;

    mat_copy_1(cview_north(&w_opt, m), ut_ret);

    mat_mult_2(U, cview_north(&beta_opt, U->n), zt_ret);
    vec_add_factor(MAT_ENTRY(&beta_opt, beta_opt.m - 1, 0), &z0, zt_ret);
    vec_add(zt, zt_ret);

    memmove(zt_ret->p, zt_ret->p + m, sizeof(double) * (zt_ret->m - m));
    memset(zt_ret->p + zt_ret->m - m, 0, sizeof(double) * m);

    mat_free(&beta_opt);
    mat_free(&w_opt);
    mat_free(&p_init);
    mat_free(&R);
    mat_free(&Q);
    mat_free(&Hcheck);
    mat_free(&z0);

    return stats;
}