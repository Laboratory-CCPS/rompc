#pragma once

#include "mex.h"

#include "algorithm/linalg.h"

static mat_t mview_matrix_from_mx(mxArray* mxa) {
    return mat_from_raw_buffer(mxGetPr(mxa), mxGetM(mxa), mxGetN(mxa));
}

static const mat_t cview_matrix_from_mx(const mxArray* mxa) {
    return mat_from_raw_buffer(mxGetPr(mxa), mxGetM(mxa), mxGetN(mxa));
}

static const mat_t cview_matrix_from_mx_field(const mxArray* mxs, const char* fieldname) {
    const mxArray* mxa = mxGetField(mxs, 0, fieldname);
    return cview_matrix_from_mx(mxa);
}

static unsigned long long get_assert_positive_integer(const mxArray* mx, int arg_id,
                                                      const char* const mex_file_id) {

    if (!mxIsNumeric(mx) || (mxGetClassID(mx) != mxDOUBLE_CLASS) || mxIsComplex(mx) ||
        !mxIsScalar(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar real value expected",
                          arg_id + 1);
    }

    const double x = *mxGetPr(mx);
    const unsigned long long n = (unsigned long long)x;

    if ((x <= 0) || mxIsNaN(x) || mxIsInf(x) || (x != n)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: postive integer value expected",
                          arg_id + 1);
    }

    return n;
}

static const mat_t cview_assert_double_matrix(const mxArray* mx, int arg_id,
                                              const char* const mex_file_id) {

    if (!mxIsNumeric(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: numeric array expected", arg_id + 1);
    }

    if (mxGetClassID(mx) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: double array expected", arg_id + 1);
    }

    if (mxIsComplex(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no complex data expected", arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mx) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no more than two dimensions expected",
                          arg_id + 1);
    }

    return cview_matrix_from_mx(mx);
}

static const mat_t cview_assert_nonempty_double_matrix(const mxArray* mx, int arg_id,
                                                       const char* const mex_file_id) {

    if (!mxIsNumeric(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: numeric array expected", arg_id + 1);
    }

    if (mxGetClassID(mx) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: double array expected", arg_id + 1);
    }

    if (mxIsComplex(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no complex data expected", arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mx) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no more than two dimensions expected",
                          arg_id + 1);
    }

    if (mxGetM(mx) == 0 || mxGetN(mx) == 0) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: non-empty matrix expected",
                          arg_id + 1);
    }

    return cview_matrix_from_mx(mx);
}

static const mat_t cview_assert_nonempty_square_double_matrix(const mxArray* mx, int arg_id,
                                                              const char* const mex_file_id) {

    if (!mxIsNumeric(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: numeric array expected", arg_id + 1);
    }

    if (mxGetClassID(mx) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: double array expected", arg_id + 1);
    }

    if (mxIsComplex(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no complex data expected", arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mx) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no more than two dimensions expected",
                          arg_id + 1);
    }

    if (mxGetM(mx) == 0 || mxGetM(mx) != mxGetN(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: non-empty square matrix expected",
                          arg_id + 1);
    }

    return cview_matrix_from_mx(mx);
}

static const mat_t cview_assert_nonempty_double_vector(const mxArray* mx, int arg_id,
                                                       const char* const mex_file_id) {

    if (!mxIsNumeric(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: numeric vector expected", arg_id + 1);
    }

    if (mxGetClassID(mx) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: double vector expected", arg_id + 1);
    }

    if (mxIsComplex(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: no complex data expected", arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mx) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: vector expected", arg_id + 1);
    }

    const size_t m = mxGetM(mx);
    const size_t n = mxGetN(mx);

    if (m != 1 && n != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: vector expected", arg_id + 1);
    }

    if (m == 0 || n == 0) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: non-empty vector expected",
                          arg_id + 1);
    }

    const size_t mv = (m > n) ? m : n;
    const mat_t v = { .p = mxGetPr(mx), .m = mv, .n = 1, .ld = 1 };

    return v;
}

static double get_assert_positive_double_from_struct_field(const mxArray* mx, const char* fieldname,
                                                           int arg_id,
                                                           const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa) || (mxGetClassID(mxa) != mxDOUBLE_CLASS) || mxIsComplex(mxa) ||
        !mxIsScalar(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: scalar real value expected",
                          fieldname, arg_id + 1);
    }

    const double x = *mxGetPr(mxa);

    if (x <= 0 || mxIsNaN(x) || mxIsInf(x)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: postive value expected",
                          fieldname, arg_id + 1);
    }

    return x;
}

static unsigned long long
get_assert_positive_integer_from_struct_field(const mxArray* mx, const char* fieldname, int arg_id,
                                              const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa) || (mxGetClassID(mxa) != mxDOUBLE_CLASS) || mxIsComplex(mxa) ||
        !mxIsScalar(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: scalar real value expected",
                          fieldname, arg_id + 1);
    }

    const double x = *mxGetPr(mxa);
    const unsigned long long n = (unsigned long long)x;

    if ((x <= 0) || mxIsNaN(x) || mxIsInf(x) || (x != n)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: postive integer value expected",
                          fieldname, arg_id + 1);
    }

    return n;
}

static const mat_t cview_assert_double_matrix_from_struct_field(const mxArray* mx,
                                                                const char* fieldname, int arg_id,
                                                                const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: numeric array expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetClassID(mxa) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: double array expected",
                          fieldname, arg_id + 1);
    }

    if (mxIsComplex(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no complex data expected", fieldname,
                          arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mxa) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no more than two dimensions expected",
                          fieldname, arg_id + 1);
    }

    return cview_matrix_from_mx(mxa);
}

static const mat_t
cview_assert_nonempty_double_matrix_from_struct_field(const mxArray* mx, const char* fieldname,
                                                      int arg_id, const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: numeric array expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetClassID(mxa) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: double array expected",
                          fieldname, arg_id + 1);
    }

    if (mxIsComplex(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no complex data expected", fieldname,
                          arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mxa) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no more than two dimensions expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetM(mxa) == 0 || mxGetN(mxa) == 0) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: non-empty matrix expected", fieldname,
                          arg_id + 1);
    }

    return cview_matrix_from_mx(mxa);
}

static const mat_t cview_assert_nonempty_square_double_matrix_from_struct_field(
    const mxArray* mx, const char* fieldname, int arg_id, const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: numeric array expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetClassID(mxa) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: double array expected",
                          fieldname, arg_id + 1);
    }

    if (mxIsComplex(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no complex data expected", fieldname,
                          arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mxa) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no more than two dimensions expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetM(mx) == 0 || mxGetM(mx) != mxGetN(mx)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: non-empty square matrix expected",
                          fieldname, arg_id + 1);
    }

    return cview_matrix_from_mx(mxa);
}

static const mat_t
cview_assert_nonempty_double_vector_from_struct_field(const mxArray* mx, const char* fieldname,
                                                      int arg_id, const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mx, 0, fieldname);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          fieldname);
    }

    if (!mxIsNumeric(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: numeric vector expected",
                          fieldname, arg_id + 1);
    }

    if (mxGetClassID(mxa) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: double vector expected",
                          fieldname, arg_id + 1);
    }

    if (mxIsComplex(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: no complex data expected", fieldname,
                          arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mxa) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: vector expected",
                          fieldname, arg_id + 1);
    }

    const size_t m = mxGetM(mxa);
    const size_t n = mxGetN(mxa);

    if (m != 1 && n != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: vector expected",
                          fieldname, arg_id + 1);
    }

    if (m == 0 || n == 0) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s' of argument %d invalid: non-empty vector expected", fieldname,
                          arg_id + 1);
    }

    const size_t mv = (m > n) ? m : n;
    const mat_t v = { .p = mxGetPr(mxa), .m = mv, .n = 1, .ld = 1 };

    return v;
}

static const mat_t
cview_assert_double_matrix_from_struct_struct_field(const mxArray* mx, const char* field_outer,
                                                    const char* field_inner, int arg_id,
                                                    const char* const mex_file_id) {

    if (!mxIsStruct(mx)) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: struct expected", arg_id + 1);
    }

    if (mxGetNumberOfElements(mx) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: scalar struct expected", arg_id + 1);
    }

    const mxArray* mxo = mxGetField(mx, 0, field_outer);

    if (mxo == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "argument %d invalid: field '%s' not found", arg_id + 1,
                          field_outer);
    }

    if (!mxIsStruct(mxo)) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: struct expected",
                          field_outer, arg_id + 1);
    }

    if (mxGetNumberOfElements(mxo) != 1) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: scalar struct expected",
                          field_outer, arg_id + 1);
    }

    const mxArray* mxa = mxGetField(mxo, 0, field_inner);

    if (mxa == NULL) {
        mexErrMsgIdAndTxt(mex_file_id, "field '%s' of argument %d invalid: field '%d' not found",
                          field_outer, arg_id + 1, field_inner);
    }

    if (!mxIsNumeric(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s'.'%s' of argument %d invalid: numeric array expected",
                          field_outer, field_inner, arg_id + 1);
    }

    if (mxGetClassID(mxa) != mxDOUBLE_CLASS) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s'.'%s' of argument %d invalid: double array expected",
                          field_outer, field_inner, arg_id + 1);
    }

    if (mxIsComplex(mxa)) {
        mexErrMsgIdAndTxt(mex_file_id,
                          "field '%s'.'%s' of argument %d invalid: no complex data expected",
                          field_outer, field_inner, arg_id + 1);
    }

    if (mxGetNumberOfDimensions(mxa) != 2) {
        // number of dimension is at least 2, even for scalar values
        mexErrMsgIdAndTxt(
            mex_file_id,
            "field '%s'.'%s' of argument %d invalid: no more than two dimensions expected",
            field_outer, field_inner, arg_id + 1);
    }

    return cview_matrix_from_mx(mxa);
}

static void assert_mex(bool cond, const char* const msg, const char* const mex_file_id) {
    if (!cond) {
        mexErrMsgIdAndTxt(mex_file_id, msg);
    }
}