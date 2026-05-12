function delta = get_condensed_primal_feasibility(w,W,Gx,Gu,Gf,fx,fu,ff,N) %# codegen
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

% dimensions
[qx,n] = size(Gx); [qu,m] = size(Gu); qf = size(Gf,1);
q = N*qu + (N-1)*qx + qf;
s = size(W,2);

% initialise
delta = NaN(q,1);

% loop over horizon
for k=1:N
    % indices for weighted control & state variables
    Iu = (k-1)*(m+n)+1:k*m+(k-1)*n;
    Ix = k*m+(k-1)*n+1:k*(m+n);
    % indices for Lagrange multipliers related to constraints
    Icu = (k-1)*(qu+qx)+1:k*qu+(k-1)*qx;
    if k==N
        Icx = k*qu+(k-1)*qx+1:k*qu+(k-1)*qx+qf;
        Gk = Gf; fk = ff;
    else
        Icx = k*qu+(k-1)*qx+1:k*(qu+qx);
        Gk = Gx; fk = fx;
    end
    % evaluate current control constraint
    delta(Icu) = Gu * w(Iu) - fu;

    % evaluate current state constraint
    delta(Icx) = Gk * w(Ix) - fk;
end
end