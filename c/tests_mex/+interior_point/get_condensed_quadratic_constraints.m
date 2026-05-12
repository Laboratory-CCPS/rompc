function D = get_condensed_quadratic_constraints(W,Gx,Gu,Gf,y,lambda,N) %# codegen
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

% dimensions
[qx,n] = size(Gx); [qu,m] = size(Gu); qf = size(Gf,1);
s = size(W,2);

% initialise
D = zeros(s);

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
    % add
    D = D + (y(Icu) .\ (Gu * W(Iu,:)))' * (lambda(Icu) .* (Gu * W(Iu,:))) ...
        + (y(Icx) .\ (Gk * W(Ix,:)))' * (lambda(Icx) .* (Gk * W(Ix,:)));
end
end