function alpha = multiplier_line_search(nu,Delta_nu) %# codegen
%MULTIPLIER_LINE_SEARCH Solves the problem
%   max nu + alpha .* Delta_nu  s.t. alpha \in (0,1]
% under the assumption that all entries of nu are positive.
%   Detailed explanation goes here

% negative update directions
ind = (Delta_nu < 0);

if all(~ind) % all update directions are nonnegative
    alpha = 1;
else
    kappa = nu(ind); Delta_kappa = Delta_nu(ind);
    mu = kappa ./ (-Delta_kappa);
    alpha = min(mu);
end
end

