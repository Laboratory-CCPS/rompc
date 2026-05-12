function H = get_condensed_hessian(W,Q,R,P,N) %# codegen
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

% dimensions
n = length(Q); m = length(R);
s = size(W,2);
% initialise 
H = zeros(s);

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
    % add
    H = H + W(Ix,:)' * (D * W(Ix,:)) + W(Iu,:)' * (R * W(Iu,:));
end
end