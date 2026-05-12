function eta = get_condensed_gradient(w,W,Q,R,P,N) %# codegen
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

% dimensions
n = length(Q); m = length(R);
s = size(W,2);

% initialise
eta = zeros(s,1);

% loop over horizon
for k=1:N
    % indices for weighted control & state variables
    Iu = (k-1)*(m+n)+1:k*m+(k-1)*n;
    Ix = k*m+(k-1)*n+1:k*(m+n);
    % terminal cost?
    if k==N
        D = P;
    else
        D = Q;
    end

    % add control cost
    eta = eta + W(Iu,:)' * (R * w(Iu));
    % add state cost
    eta = eta + W(Ix,:)' * (D * w(Ix));
end
end