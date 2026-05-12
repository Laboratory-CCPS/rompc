function delta = get_condensed_primal_feasibility(w, W, Gx, Gu, Gf, fx, fu, ff, N)
%GET_CONDENSED_PRIMAL_FEASIBILITY Summary of this function goes here
%   w:  Control-state sequence that corresponds to p, i.e. w = Zg + Wp
%       (N*(n+m) by 1) vector
%   W:  First N*(n+m) rows of Q_2
%       (N*(n+m) by s) matrix
%   Gx:     lhs of polyhedral state constraints
%           (qx by n) matrix
%   Gu:     lhs of polyhedral control constraints
%           (qu by m) matrix
%   Gf:     lhs of polyhedral state terminal constraints
%           (qf by n) matrix
%   qx:     rhs of polyhedral state constraints
%           (qx by 1) vector
%   qu:     rhs of polyhedral control constraints
%           (qu by 1) vector
%   qf:     rhs of polyhedral state terminal constraints
%           (qf by 1) vector
%   N:      Prediction horizon
%           Natural number