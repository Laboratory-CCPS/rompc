function theta = get_condensed_adjoint(lambda, W, Gx, Gu, Gf, N)
%GET_CONDENSED_ADJOINT Evaluates the "adjoint" equation G^\top \lambda
%   lambda: Lagrange multiplier
%           (q by 1) vector
%   W:  First N*(n+m) rows of Q_2
%       (N*(n+m) by s) matrix
%   Gx:     lhs of polyhedral state constraints
%           (qx by n) matrix
%   Gu:     lhs of polyhedral control constraints
%           (qu by m) matrix
%   Gf:     lhs of polyhedral state terminal constraints
%           (qf by n) matrix
%   N:      Prediction horizon
%           Natural number