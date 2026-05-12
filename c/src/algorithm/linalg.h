#pragma once

#include <math.h>
#include <string.h>

#include "blas.h"
#include "lapack.h"

#include "assertions.h"

#define INLINE inline static
// #define INLINE __forceinline

typedef ptrdiff_t lapackint_t;

#ifdef COMPILE_MEX_FUNCTION
#include "mex.h"
#define mallocX(a) mxMalloc(a)
#define callocX(a, b) mxCalloc(a, b)
#define freeX(a) mxFree(a)
#else
#define mallocX(a) malloc(a)
#define callocX(a, b) calloc(a, b)
#define freeX(a) free(a)
#endif

// By default, all matrices (mat_t) are passed by-pointer.
// Sometimes it is more convenient for the functions to accept a mat_t value instead of a pointer
// to mat_t, to avoid introducing a local variable.
// This is the case in a situation where one would like to write something like
//      mat_copy(&cview_north(A, 5), B)
// which is invalid as the address operator cannot be applied to unnamed objects.
// For this cases there are variants of the functions marked with an underscore and a number that
// accept mat_t values direclty, for example
//      mat_copy_1(cview(A, 5), B)
// The number encodes the number of the arguments that are passed by-value in a flag-like manner,
// e.g. _1 means the first argument, _2 the second, _4 the third, and _3 the first and the second.
// These variants are not implemented for all possible combinations, but added only for these cases
// where they were needed.

// A function wit "vec" in its name assumes that the corresponding operand is a column vector, that
// is, a matrix with one column. (This is asserted in debug builds)


// We are using "plain" BLAS/lapack (and not lapacke). Therefore we need to pass also scalar values
// as pointers. The commonly used values are defined here for convenience.
// (This might not be optimal, as they reside in the "static" memory area of the program, and
// potentially must be loaded from outside the cache. It might be more efficient to define the
// values needed locally within the functions.)
static const double fone = 1;
static const double fzero = 0;
static const double fmone = -1;
static const lapackint_t ione = 1;

typedef struct {
    double* p;
    lapackint_t n;
} sizedbuffer_t;

INLINE sizedbuffer_t sizedbuffer_create(lapackint_t n) {
    sizedbuffer_t buf = { .p = mallocX(sizeof(double) * n), .n = n };
    return buf;
}

INLINE void sizedbuffer_free(sizedbuffer_t* buf) {
    freeX(buf->p);
    buf->p = NULL;
    buf->n = 0;
}

typedef struct {
    double* p;
    lapackint_t m;
    lapackint_t n;
    lapackint_t ld;
} mat_t;

#define MAT_ENTRY_PTR(A, i, j) ((A)->p + (i) + (A)->ld * (j))
#define MAT_ENTRY(A, i, j) (*MAT_ENTRY_PTR(A, i, j))

void mat_print_dim(const mat_t* A);
void mat_print(const mat_t* A);
void mat_dbg(const char* name, const mat_t* A);

INLINE mat_t mat_create(lapackint_t m, lapackint_t n) {
    mat_t A = { .p = (double*)mallocX(sizeof(double) * m * n), .m = m, .n = n, .ld = m };

    return A;
}

INLINE mat_t mat_create_zeros(lapackint_t m, lapackint_t n) {
    mat_t A = { .p = (double*)callocX(sizeof(double), m * n), .m = m, .n = n, .ld = m };

    return A;
}

INLINE mat_t mat_create_ones(lapackint_t m, lapackint_t n) {
    mat_t A = { .p = (double*)mallocX(sizeof(double) * m * n), .m = m, .n = n, .ld = m };

    for (lapackint_t i = 0; i < m * n; ++i) {
        A.p[i] = 1.0;
    }

    return A;
}

INLINE void mat_copy(const mat_t* A, mat_t* B) {
    ASSERT(A->m == B->m);
    ASSERT(A->n == B->n);

    dlacpy("x", &A->m, &A->n, A->p, &A->ld, B->p, &B->ld);
}

INLINE void mat_copy_1(const mat_t A, mat_t* B) { mat_copy(&A, B); }

INLINE void mat_copy_2(const mat_t* A, mat_t B) { mat_copy(A, &B); }

INLINE void mat_copy_3(const mat_t A, mat_t B) { mat_copy(&A, &B); }

INLINE mat_t mat_create_copy(const mat_t* A) {
    mat_t B = {
        .p = (double*)mallocX(sizeof(double) * A->m * A->n), .m = A->m, .n = A->n, .ld = A->m
    };

    if (A->ld == A->m) {
        memcpy(B.p, A->p, A->m * A->n * sizeof(double));
    } else {
        mat_copy(A, &B);
    }

    return B;
}

INLINE mat_t mat_create_copy_1(const mat_t A) { return mat_create_copy(&A); }

INLINE void mat_free(mat_t* A) {
    freeX(A->p);
    A->p = NULL;
}

INLINE mat_t mat_from_raw_buffer(double* p, lapackint_t m, lapackint_t n) {
    mat_t A = { .p = p, .m = m, .n = n, .ld = m };
    return A;
}

INLINE mat_t mat_from_sizedbuffer(sizedbuffer_t* work, lapackint_t m, lapackint_t n) {
    ASSERT(work->n >= m * n);

    mat_t A = { .p = work->p, .m = m, .n = n, .ld = m };
    return A;
}

INLINE void mat_set_zero(mat_t* A) {
    if (A->m == A->ld) {
        memset(A->p, 0, sizeof(double) * A->m * A->n);
    } else {
        for (lapackint_t i = 0; i < A->n; ++i) {
            memset(A->p + i * A->ld, 0, sizeof(double) * A->n);
        }
    }
}

// ==================================================
// Functions to create "views" into an existing matrix.
// A "view" is actually just another mat_t object that describes (parts of) another matrix.
// These functions help to specify these submatrices in a readable fashion.
//
// mview_... vs cview_...
//      Every function is provided in two variants:
//      - mview... (mutable view): Takes a pointer to non-const mat_t and returns a mat_t
//      - cview... (const view): Takes a pointer to const mat_t and returns a const mat_t
//      The intention is to use the corresponding function to make ones intent clear. From a
//      language point of view the cview-version would work in each case without any warnings.
//      (A pointer to non-const can always be casted into a pointer to const without warnings,
//      and the const qualifier on a returned value type hasn't any significance, i.e. there will
//      be not even a warning if the returned value is assigned to a non-const variable.)
//
// _idx_ vs _dim_
//      Often the specification of the submatrix can be given in two ways:
//      - idx: zero-based indices (if this specifies an interval end, it is exclusive)
//      - dim: size
//      The used variant is indicated in the function name using "_idx_" and "_dim_".
//      Not all possible variants are implemented.
//
//      For example:
//      - mview_south_east_idx(A, m0, n0): Returns a view of the lower right part of the matrix A
//                                         starting at the position (m0, n0).
//                                         (Matlab: A(m0 + 1:end, n0 + 1:end))
//      - mview_south_east_dim(A, m, n)  : Returns a view of the lower right part of the matrix A
//                                         consisting of m rows and n columns.
//                                         (Matlab: A(end - m + 1:end, end - n + 1:end))
//
//      - mview_idx_dim(A, m0, n0, m, n): Returns a view the submatrix von A that starts at the
//                                        index (m0, n0) and has m rows and n columns.
//                                        (Matlab: A(m0 + 1:m0 + m, n0 + 1: n0 + n))
//
// In debug builds it is asserted that the submatrix specified lies completely within the "parent"
// matrix.

INLINE mat_t mview_idx_dim(mat_t* A, lapackint_t m0, lapackint_t n0, lapackint_t m, lapackint_t n) {
    ASSERT(A->m >= m0 + m);
    ASSERT(A->n >= n0 + n);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, m0, n0), .m = m, .n = n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_idx_dim(const mat_t* A, lapackint_t m0, lapackint_t n0, lapackint_t m,
                                 lapackint_t n) {

    return mview_idx_dim((mat_t*)A, m0, n0, m, n);
}

INLINE mat_t mview_columns_idx_dim(mat_t* A, lapackint_t n0, lapackint_t n) {
    ASSERT(A->n >= n0 + n);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, 0, n0), .m = A->m, .n = n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_columns_idx_dim(const mat_t* A, lapackint_t n0, lapackint_t n) {
    return mview_columns_idx_dim((mat_t*)A, n0, n);
}

INLINE mat_t mview_rows_idx_dim(mat_t* A, lapackint_t m0, lapackint_t m) {
    ASSERT(A->m >= m0 + m);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, m0, 0), .m = m, .n = A->n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_rows_idx_dim(const mat_t* A, lapackint_t m0, lapackint_t m) {
    return mview_rows_idx_dim((mat_t*)A, m0, m);
}

INLINE mat_t mview_north(mat_t* A, lapackint_t m) {
    ASSERT(A->m >= m);

    mat_t blk = { .p = A->p, .m = m, .n = A->n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_north(const mat_t* A, lapackint_t m) { return mview_north((mat_t*)A, m); }

INLINE mat_t mview_south_idx(mat_t* A, lapackint_t m0) {
    ASSERT(A->m > m0);

    mat_t blk = { .p = A->p + m0, .m = A->m - m0, .n = A->n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_south_idx(const mat_t* A, lapackint_t m0) {
    return mview_south_idx((mat_t*)A, m0);
}

INLINE mat_t mview_west(mat_t* A, lapackint_t n) {
    ASSERT(A->n > n);

    mat_t blk = { .p = A->p, .m = A->m, .n = n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_west(const mat_t* A, lapackint_t n) { return mview_west((mat_t*)A, n); }

INLINE mat_t mview_east_idx(mat_t* A, lapackint_t n0) {
    ASSERT(A->n > n0);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, 0, n0), .m = A->m, .n = A->n - n0, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_east_idx(const mat_t* A, lapackint_t n0) {
    return mview_east_idx((mat_t*)A, n0);
}

INLINE mat_t mview_north_east_dim(mat_t* A, lapackint_t m, lapackint_t n) {
    ASSERT(A->m >= m);
    ASSERT(A->n >= n);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, 0, A->n - n), .m = m, .n = n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_north_east_dim(const mat_t* A, lapackint_t m, lapackint_t n) {
    return mview_north_east_dim((mat_t*)A, m, n);
}

INLINE mat_t mview_south_east_idx(mat_t* A, lapackint_t m0, lapackint_t n0) {
    ASSERT(A->m > m0);
    ASSERT(A->n > n0);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, m0, n0), .m = A->m - m0, .n = A->n - n0, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_south_east_idx(const mat_t* A, lapackint_t m0, lapackint_t n0) {
    return mview_south_east_idx((mat_t*)A, m0, n0);
}

INLINE mat_t mview_south_east_dim(mat_t* A, lapackint_t m, lapackint_t n) {
    ASSERT(A->m >= m);
    ASSERT(A->n >= n);

    mat_t blk = { .p = MAT_ENTRY_PTR(A, A->m - m, A->n - n), .m = m, .n = n, .ld = A->ld };

    return blk;
}

INLINE const mat_t cview_south_east_dim(const mat_t* A, lapackint_t m, lapackint_t n) {
    return mview_south_east_dim((mat_t*)A, m, n);
}

// ==================================================
// Operations

// C <- A * B
INLINE void mat_mult(const mat_t* A, const mat_t* B, mat_t* C) {
    ASSERT(A->n == B->m);
    ASSERT(A->m == C->m);
    ASSERT(B->n == C->n);

    dgemm("n", "n", &A->m, &B->n, &A->n, &fone, A->p, &A->ld, B->p, &B->ld, &fzero, C->p, &C->ld);
}

INLINE void mat_mult_2(const mat_t* A, const mat_t B, mat_t* C) { mat_mult(A, &B, C); }

// C <- A^T * B + C
INLINE void matT_mat_mult_add(const mat_t* A, const mat_t* B, mat_t* C) {
    ASSERT(A->m == B->m);
    ASSERT(A->n == C->m);
    ASSERT(B->n == C->n);

    dgemm("t", "n", &A->n, &B->n, &A->m, &fone, A->p, &A->ld, B->p, &B->ld, &fone, C->p, &C->ld);
}

// c <- A * b - c
INLINE void mat_vec_mult_minus(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->n == b->m);
    ASSERT(c->m == A->m);
    ASSERT(c->n == 1);
    ASSERT(b->n == 1);

    dgemv("n", &A->m, &A->n, &fone, A->p, &A->ld, b->p, &ione, &fmone, c->p, &ione);
}

// c <- A^T * b
INLINE void matT_vec_mult(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->m == b->m);
    ASSERT(c->m == A->n);
    ASSERT(b->n == 1);
    ASSERT(c->n == 1);

    dgemv("t", &A->m, &A->n, &fone, A->p, &A->ld, b->p, &ione, &fzero, c->p, &ione);
}

// c <- A^T * b + c
INLINE void matT_vec_mult_add(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->m == b->m);
    ASSERT(c->m == A->n);
    ASSERT(b->n == 1);
    ASSERT(c->n == 1);

    dgemv("t", &A->m, &A->n, &fone, A->p, &A->ld, b->p, &ione, &fone, c->p, &ione);
}

// c <- A * b
INLINE void mat_vec_mult(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->n == b->m);
    ASSERT(c->m == A->m);
    ASSERT(b->n == 1);
    ASSERT(c->n == 1);

    dgemv("n", &A->m, &A->n, &fone, A->p, &A->ld, b->p, &ione, &fzero, c->p, &ione);
}

INLINE void mat_vec_mult_1(const mat_t A, const mat_t* b, mat_t* c) { mat_vec_mult(&A, b, c); }

// c <- A * b + c
INLINE void mat_vec_mult_add(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->n == b->m);
    ASSERT(c->m == A->m);
    ASSERT(b->n == 1);
    ASSERT(c->n == 1);

    dgemv("n", &A->m, &A->n, &fone, A->p, &A->ld, b->p, &ione, &fone, c->p, &ione);
}

INLINE void mat_vec_mult_add_1(const mat_t A, const mat_t* b, mat_t* c) {
    mat_vec_mult_add(&A, b, c);
}

// c <- -A * b
INLINE void mat_vec_neg_mult(const mat_t* A, const mat_t* b, mat_t* c) {
    ASSERT(A->n == b->m);
    ASSERT(c->m == A->m);
    ASSERT(b->n == 1);
    ASSERT(c->n == 1);

    dgemv("n", &A->m, &A->n, &fmone, A->p, &A->ld, b->p, &ione, &fzero, c->p, &ione);
}

// ret <- A * B
INLINE mat_t mat_mult_create(const mat_t* A, const mat_t* B) {
    ASSERT(A->n == B->m);

    mat_t C = mat_create(A->m, B->n);

    dgemm("n", "n", &A->m, &B->n, &A->n, &fone, A->p, &A->m, B->p, &B->n, &fzero, C.p, &C.m);

    return C;
}

// A <- diag([s_1, ..., s_m]) * A
INLINE void mat_scale_diag(const mat_t* S, mat_t* A) {
    ASSERT(S->m == A->m);
    ASSERT(S->n == 1);

    for (lapackint_t i = 0; i < A->m; ++i) {
        dscal(&A->n, S->p + i, A->p + i, &A->ld);
    }

    // There is a lapack function DLASCL2 to perform this task, but it is not available in mwlapack.
    // dlascl2(&A.m, &A.n, S.p, A.p, &A.m);
}

// A <- diag([1 / s_1, ..., 1 / s_m]) * A
INLINE void mat_scale_recip_diag(const mat_t* S, mat_t* A) {
    ASSERT(S->m == A->m);
    ASSERT(S->n == 1);

    for (lapackint_t i = 0; i < A->m; ++i) {
        double s = 1 / (S->p[i]);
        dscal(&A->n, &s, A->p + i, &A->ld);
    }

    // There is a lapack function DLARSCL2 to perform this task, but it is not available in
    // mwlapack.
    // dlarscl2(&A.m, &A.n, S.p, A.p, &A.m);
}

// b <- a + b
INLINE void vec_add(const mat_t* a, mat_t* b) {
    ASSERT(a->m == b->m);
    ASSERT(a->n == 1);
    ASSERT(b->n == 1);

    daxpy(&a->m, &fone, a->p, &ione, b->p, &ione);
}

// b <- -a + b
INLINE void vec_minus(const mat_t* a, mat_t* b) {
    ASSERT(a->m == b->m);
    ASSERT(a->n == 1);
    ASSERT(b->n == 1);

    daxpy(&a->m, &fmone, a->p, &ione, b->p, &ione);
}

INLINE void vec_minus_2(const mat_t* a, mat_t b) { vec_minus(a, &b); }

// b <- s * a + b
INLINE void vec_add_factor(double s, const mat_t* a, mat_t* b) {
    ASSERT(a->m == b->m);
    ASSERT(a->n == 1);
    ASSERT(b->n == 1);

    daxpy(&a->m, &s, a->p, &ione, b->p, &ione);
}

// a <- s * a
INLINE void vec_scalar_mult(double s, mat_t* a) {
    ASSERT(a->n == 1);

    dscal(&a->m, &s, a->p, &ione);
}

// B <- A + B, for contiguous matrices (i.e. matrices with .ld == .m)
INLINE void mat_add_cont(const mat_t* A, mat_t* B) {
    ASSERT(A->ld == A->m);
    ASSERT(B->ld == B->m);
    ASSERT(A->m == B->m);
    ASSERT(A->n == B->n);

    const lapackint_t n = A->m * A->n;
    daxpy(&n, &fone, A->p, &ione, B->p, &ione);
}

// ret <- max_i |v_i|
// (There is the BLAS function IDAMAX that returns the index of the maximum element.
// It should be tested if this is faster when proper benchmarking is set up.)
INLINE double vec_inf_norm(const mat_t* v) {
    ASSERT(v->n == 1);

    double gam = fabs(v->p[0]);

    for (lapackint_t i = 1; i < v->m; ++i) {
        double c = fabs(v->p[i]);
        gam = (gam > c) ? gam : c;
    }

    return gam;
}

// ret <- a^T * b
INLINE double vec_dot(const mat_t* a, const mat_t* b) {
    ASSERT(a->m == b->m);
    ASSERT(a->n == 1);
    ASSERT(b->n == 1);

    return ddot(&a->m, a->p, &ione, b->p, &ione);
}

// upper(A) <- U where A = U^T U
// ret <- 0 if succesful
INLINE lapackint_t chol(mat_t* A) {
    lapackint_t info;

    dpotrf("U", &A->m, A->p, &A->ld, &info);

    return info;
}

// b_i <- a_i * b_i
INLINE void vec_mult(const mat_t* a, mat_t* b) {
    ASSERT(a->m == b->m);
    ASSERT(a->n == 1);
    ASSERT(b->n == 1);

    for (lapackint_t i = 0; i < a->m; ++i) {
        b->p[i] *= a->p[i];
    }
}

typedef struct {
    double c;
    double s;
} givens_t;

// Returns the Given rotation for the two scalars a and b, such that
// %R^\top [a;b] = [c; 0], with R = [c s; -s c]. It follows [Golub, Alg. 5.1.3]
INLINE givens_t givens(double a, double b) {
    givens_t res;

    if (b == 0) {
        res.c = 1.0;
        res.s = 0.0;
    } else {
        if (fabs(b) > fabs(a)) {
            const double t = -a / b;
            res.s = 1 / sqrt(1 + t * t);
            res.c = res.s * t;
        } else {
            const double t = -b / a;
            res.c = 1 / sqrt(1 + t * t);
            res.s = res.c * t;
        }
    }

    return res;
}
