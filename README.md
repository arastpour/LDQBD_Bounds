The program in this repository is an implementation of Algorithm 5 reported in the following paper: Rastpour, A. A. Ingolfsson, B. Sandıkcı. (2021) Algorithms for queueing systems with reneging and priorities modeled as quasi-birth-death processes. *INFORMS Journal on Computing*; Under Minor Revision. This algorithm is developed for calculateing bounds on performance measures of an Erlang A system with two priority classes.

Software requirement: this program has to be run on MATLAB R2019, or newer versions of MATLAB.

This program includes 8 MATLAB `.m` files as described below:
1. `Main.m` is the main file, which is used to take user inputs and calls all other functions to calculate the bounds.
2. `function_BlockTridiagonalInverse2.m` calculates the inverse of a block-tridiagonal matrix using its properties.
3. `function_class2waiting.m` calculates bounds on the Class-2 Wait Probability. If a user is interested in another performance measure, this function should be adjusted accordingly.
4. `function_ErlangA_truncation_level.m` calculates the truncation phase, $p^\ast$.
5. `function_Q.m` forms block matrices of the infiniteseimal matrix at any desired level ell, for ell = 0, 1, ... .
6. `function_R_low_and_up.m` calculates bounds on the rate matrices of level ell, for ell = 0, 1, ... .
7. `function_tau_cond.m` calcualtes the threshold level $\ell^\tau$, above which rate matrices become non-increasing.
8. `function_x0.m` calculates the unnormalized stationary distribution of Level 0.
