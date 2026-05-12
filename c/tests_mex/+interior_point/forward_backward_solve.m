function x = forward_backward_solve(G,b) %# codegen
%FORWARD_BACKWARD_SOLVE Solves the linear system of equation Ax = b, where
%A = G' * G with G upper triangular
%   G:  Defines matrix A
%   b:  rhs of the linear system of equations
%   x:  Vector that solves the system of equations

% dimension
n = length(b);

% first step: forward solve
y = NaN(n,1); % initialize
for i=1:n
    % evaluate rhs
    y(i) = b(i);
    % add "past" information
    for j=1:i-1
        y(i) = y(i) - G(j,i) * y(j);
    end
    % final step
    y(i) = y(i) / G(i,i);
end

% second step: backward solve
x = NaN(n,1); % initialize
for i=1:n
    % backward time
    k = n - i + 1; % we start at k = n
    % evaluate rhs
    x(k) = y(k);
    % add "past" information
    for j=k+1:n
        x(k) = x(k) - G(k,j) * x(j);
    end
    % final step
    x(k) = x(k) / G(k,k);
end
end