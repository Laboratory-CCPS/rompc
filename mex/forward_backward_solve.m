function x = forward_backward_solve(G, b)
%FORWARD_BACKWARD_SOLVE Solves the linear system of equation Ax = b, where
%A = G' * G with G upper triangular
%   G:  Defines matrix A
%   b:  rhs of the linear system of equations
%   x:  Vector that solves the system of equations