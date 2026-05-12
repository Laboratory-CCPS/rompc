cmake_minimum_required(VERSION 3.26)
project(condensed_solver_mex)

set(CMAKE_CXX_STANDARD 14)

# Headers are listed here as well to tell the IDE that
# they are part of the project.
add_executable(condensed_solver_mex
        src/algorithm/alg.c
        src/algorithm/alg.h
        src/algorithm/assertions.h
        src/algorithm/linalg.c
        src/algorithm/linalg.h
        src/algorithm/utils.h
        src/mexinterface/mx_helper.h
        src/mexinterface/rompc_mex_helper.h
        src/mexinterface/condensed_predictor_corrector.c
        src/mexinterface/condensed_rompc_step.c
        src/mexinterface/forward_backward_solve.c
        src/mexinterface/get_condensed_adjoint.c
        src/mexinterface/get_condensed_gradient.c
        src/mexinterface/get_condensed_hessian.c
        src/mexinterface/get_condensed_primal_feasibility.c
        src/mexinterface/get_condensed_quadratic_constraints.c
        src/mexinterface/multiplier_line_search.c
        src/mexinterface/update_qr.c
)

target_include_directories(condensed_solver_mex
        PRIVATE
        src
)

target_include_directories(condensed_solver_mex
        PRIVATE
        "C:/Program Files/MATLAB/R2025a/extern/include"
)
