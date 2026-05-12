function alpha = multiplier_line_search(nu, Delta_nu)
%MULTIPLIER_LINE_SEARCH Solves the problem
%   max nu + alpha .* Delta_nu  s.t. alpha \in (0,1]
% under the assumption that all entries of nu are positive.
%   Detailed explanation goes here