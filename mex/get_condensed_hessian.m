function H = get_condensed_hessian(W, Q, R, P, N)
%GET_CONDENSED_HESSIAN Returns the Hessian that result from condensing
%via QR decomposition
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