function theta = get_condensed_adjoint(lambda,W,Gx,Gu,Gf,N) %# codegen
%GET_ADJOINT Evaluates the "adjoint" equation G^\top \lambda
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

% dimensions
[qx,n] = size(Gx); [qu,m] = size(Gu); qf = size(Gf,1);
s = size(W,2);

% initialise
theta = zeros(s,1);

% loop over horizon
for k=1:N
    % indices for weighted control & state variables
    Iu = (k-1)*(m+n)+1:k*m+(k-1)*n;
    Ix = k*m+(k-1)*n+1:k*(m+n);
    % indices for Lagrange multipliers related to constraints
    Icu = (k-1)*(qu+qx)+1:k*qu+(k-1)*qx;
    if k==N
        Icx = k*qu+(k-1)*qx+1:k*qu+(k-1)*qx+qf;
        Gk = Gf;
    else
        Icx = k*qu+(k-1)*qx+1:k*(qu+qx);
        Gk = Gx;
    end
    % add control component
    theta = theta + W(Iu,:)' * (Gu' * lambda(Icu));
    % add state component
    theta = theta + W(Ix,:)' * (Gk' * lambda(Icx));
end
end