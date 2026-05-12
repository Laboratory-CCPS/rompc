function eta = get_condensed_gradient(w, W, Q, R, P, N)
%GET_CONDENSED_GRADIENT Returns the gradient eta = Hp + F at the current point
%   w:  Control-state sequence that corresponds to p, i.e. w = Zg + Wp
%       (N*(n+m) by 1) vector
%   W:  First N*(n+m) rows of Q_2
%       (N*(n+m) by s) matrix
%   Q:  State stage cost
%       (n by n) positive definite matrix
%   R:  Control stage cost
%       (m by m) positive definite matrix
%   P:  State terminal cost
%       (n by n) positive definite matrix
%   N:  Prediction horizon
%       Natural number
