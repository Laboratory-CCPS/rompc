# A Condensed and Efficient Interior Point Solver for Reduced-Order Model Predictive Control

This code accompanies the paper *A Condensed and Efficient Interior Point Solver for Reduced-Order Model Predictive Control* [1] presented at IFAC World Congress 2026.
It provides a C implementation of an interior point solver tailored to reduced-order model predictive control (MPC) as well as the code to reproduce the results presented in the aforementioned paper.

The following optimisation problem is solved:

$$
\mathscr{P}_N(x,\ \tilde z): \qquad \min_{\alpha,\ \tau} x_N^T P_f x_N + \sum_{k=0}^{N-1} x_k^T Q_c x_k + u_k^T R_c u_k
$$

$$
\mathrm{s.t.}\hspace{1em}
\begin{array}{rll}
    x_0 & = x \\
    x_{k+1} & = A x_k + B u_k & k = 0,\ \ldots,\ N-1\\
    u_k & = K x_k + z_k & k = 0,\ \ldots,\ N-1\\
    G_x x_k & \leq f_x & k = 0,\ \ldots,\ N-1\\
    G_f x_N & \leq f_f \\
    G_u u_k & \leq f_u & k = 0,\ \ldots,\ N-1\\
\end{array}
$$

$$
z = U \alpha + \tau \xi + (1 - \tau) \tilde z \qquad \mathrm{with} \quad  z = (z_0^T \ \ldots \ z_{N-1}^T)^T
$$

- Non-fixed parameters (MPC application):
    - The current (measured) state $x$ (that must also obey the state constraints described below).
    - An initial admissible guess $\tilde z$ for the non-reduced optimisation variable (see below).
- Fixed parameters:
    - $N$ is the prediction horizon.
    - $P_f$, $Q_c$, and $R_c$ are positive definite weighting matrices.
    - $A$ and $B$ describe a linear, time-discrete, time-invariant dynamic system with $n$ states and $m$ inputs.
    - $K$ is a stabilising, fixed state feedback gain.
    - The states and inputs must lie within compact polyhedral sets defined by $(G_f,\ f_f)$, $(G_x,\ f_x)$, and $(G_u,\ f_u)$, each of which must contain the origin in its interior.
    - The definition of the subspace $(U,\ \xi)$ (see below) used for order-reduction.
- Optimisation variables:
    - $z \in \mathbb{R}^{Nm}$ is the non-reduced optimisation variable consisting of all input values (before adding the output of the stabilising feedback $K$) for all time steps of the prediction horizon.
    - The order-reduction is achieved by constraining $z$ to an $r$-dimensional affine set described by $U$ and $\xi$, $\\{U p + \xi\ | \ p  \in \mathbb{R}^r \\}$
    - Furthermore, to provide feasibility and stability guarantees for the reduced-order problem, an initial admissible guess $\tilde z$ for the non-reduced optimisation variable is incorporated with an additional degree of freedom [2].
        In this setup the non-reduced optimisation variable is given by

        $z = U \alpha + \tau \xi + (1 - \tau) \tilde z$

        with $\alpha \in \mathbb{R}^{r}$ and $\tau \in \mathbb{R}$ being the optimisation variables for the reduced-order problem $\mathscr{P}_N(x,\ \tilde z)$.


## Usage

### Prerequisites
- Matlab must be provided with a C-compiler to build mex-files.
    To check if a compiler is provided, type
    ```matlab
    mex -setup
    ```
    into the Matlab prompt.

    If no compiler is selected and no compiler is offered for selection also,
    the (free) Matlab-Add-On `MATLAB Support for MinGW-w64 C/C++/Fortran Compiler` can be installed to provide one.

- [Yalmip](https://yalmip.github.io/) including the [MOSEK solver](https://www.mosek.com/) must be installed.
    (This is not necessary for compiling and using our solver, but to run the example that reproduces the results presented in [1].)

- To generate the problem data, using the method presented in [2], the Matlab toolbox [MPT3](https://www.mpt3.org) must be installed.
    Alternatively the problem data is also provided in `six_masses_condensed.mat`.


### Run the example

1. Building the mex files
    (Must be done once.)

    - Open the folder `./c` in Matlab.
    - Run the script `build_mex.m`

        This will create the mex file for the complete algorithm as well as separate mex files
        for individual substeps of the algorithm.
        The mex files are stored in the `./mex` folder, that initially only contains the m-files that
        provide the documentation.

2. (optional) Generating the problem data (Must be done once.)
    - Open the root of this repository in Matlab.
    - Run the script `six_masses_condensed.m`.
        (This script requires MPT3, see above.)

        This will overwrite the provided file `six_masses_condensed.mat` with an equivalent (but not necessarily identical) representation of the problem data.

3. Run the example
    - Open the root of this repository in Matlab.
    - Add the path `./mex` to the Matlab search path
        ```matlab
        addpath("./mex");
        ```
    - Run the script `paper_ifac_2026.m`

        There are deliberately no outputs to the console to reduce influences on the time measurement.
        On our test machine it runs for about 30 seconds.

## Remarks for development

- The project is set up in a way that the algorithm itself does not depend on Matlab.
    However in the current state of this repository, only the setup for a compilation as mex file is provided.

- The `CMakeLists.txt` defines a "sham" project that is not intended to be actually compiled, but serves only to setup the paths for an IDE correctly.

    - Within `CMakeLists.txt` the path to the local Matlab installation must be adapted accordingly.

- A Matlab-implementation of the algorithm is provided within the `./c/tests_mex` folder that is used to run tests.
    To run the tests (or benchmarks, which compares the execution times of the different implementations), change into the folder `./c` and run `add_dev_paths`.
    Then (after building the mex-files) the test and benchmarks can be executed by running `run_tests` and `run_benchmarks`, resp.

- When compiling the code as mex function (in debug mode), the preprocessor constant `COMPILE_FOR_MEX` has to be defined (`-DCOMPILE_FOR_MEX`).
    This makes the assertions within the program use the mex-versions.
    If the default C `assert`s are used instead, triggering an assertion will end the Matlab process.
    (This is done automatically if the provided `build_mex.m` is used.)


## References

[1] R. Schurig, E. Lenz, and R. Findeisen (2026). *A Condensed and Efficient Interior Point Solver for Reduced-Order Model Predictive Control*. In 23rd IFAC World Congress 2026, *accepted*

[2] R. Schurig, A. Himmel, and R. Findeisen (2024). *Geometric data-driven dimensionality reduction in MPC with guarantees.* In 2024 European Control Conference (ECC), 3164–3170

## License

This software is licensed under the MIT license.
