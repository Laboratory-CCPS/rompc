#pragma once

#ifdef COMPILE_MEX_FUNCTION
	#include "mex.h"

	#define ASSERT(expr) mxAssert(expr, "")
	#define ASSERT_MSG(expr, msg) mxAssert(expr, msg)
#else
	#include <assert.h>

	#define ASSERT(expr) assert(expr)
	#define ASSERT_MSG(expr, msg) assert((expr) && (msg))
#endif
