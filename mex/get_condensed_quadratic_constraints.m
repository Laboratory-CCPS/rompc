function D = get_condensed_quadratic_constraints(W, Gx, Gu, Gf, y, lambda, N)
%GET_CONDENSED_QUADRATIC_CONSTRAINTS Returns the matrix G^\top S G
%   W:      First N*(n+m) rows of Q_2
%           (N*(n+m) by s) matrix
%   Gx:     lhs of polyhedral state constraints
%           (qx by n) matrix
%   Gu:     lhs of polyhedral control constraints
%           (qu by m) matrix
%   Gf:     lhs of polyhedral state terminal constraints
%           (qf by n) matrix
%   y:      Slack variable
%           (q by 1) vector
%   lambda: Lagrange multiplier
%           (q by 1) vector
%   N:      Prediction horizon
%           Natural number