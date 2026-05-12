function [w_opt,beta_opt,flag] = condensed_predictor_corrector(x0,zt,sys,rompc,p_init) %# codegen
%CONDENSED_PREDICTOR_CORRECTOR Implements the Predictor-Algorithm as in
%[Alg. 16.4, Nocedal/Wright] for the condensed roMPC-QP
%   x0:         Initial state
%               (n by 1) vector
%   zt:         Admissible guess
%               (N*m by 1) vector
%
%   sys:        System description, as a structure. With fields:
%
%               A:  system matrix
%                   (n by n) matrix
%               B:  input matrix
%                   (n by m) matrix
%               Q:  Weight matrix Q
%                   (n by n) matrix
%               R:  Weight matrix R
%                   (m by m) matrix
%               Gx: lhs of plyhedral state constraints
%                   (qx by n) matrix
%               fx: rhs of polyhedral state constraints
%                   (qx by 1 vector)
%               Gu: lhs of polyhedral control constraints
%                   (qu by m) matrix
%               fu: rhs of polyhedral control constraints
%                   (qu by 1) vector
%
%   rompc:       Reduced-order MPC problem, as a structure. With fields:
%
%               N:  Prediction horizon
%                   Natural number
%               P:  Terminal state cost
%                   (n by n) matrix
%               K:  pre-stabilising feedback gain
%                   (m by n) matrix
%               Gf: lhs of polyhedral terminal state coonstraints
%                   (qf by n) matrix
%               ff: rhs of polyhedral terminal state coonstraints
%                   (qf by 1) vector
%               eq: Represents the equality constraints, as a structure.
%                   With fields:
%
%                   Q1, (N(n+m)+s by N(n+m)) matrix
%                   Q2, (N(n+m)+s by s) matrix
%                   R,  (N(n+m) by N*(n+m)) matrix

TOL = 1e-6;
MAX_ITER = 50;

% dimensions
n = length(sys.Q); m = length(sys.R);
N = rompc.N;
q = N*size(sys.Gu,1) + (N-1)*size(sys.Gx,1) + size(rompc.Gf,1);

% partition
Z = rompc.eq.Q1(1:N*(n+m),:);
W = rompc.eq.Q2(1:N*(n+m),:);

%%% FIND PARTICULAR SOLUTION TO EQUALITY CONSTRAINT %%% [Golub, Alg. 3.1.1]
% initialise
g = [sys.A*x0; zeros((N-1)*n,1); ...
    [-rompc.K*x0; zeros((N-1)*m,1)] - zt];
g(1) = g(1) / rompc.eq.R(1,1);
for i=2:N*(n+m)
    g(i) = (g(i) - rompc.eq.R(1:i-1,i)'*g(1:i-1)) / rompc.eq.R(i,i);
end

%%% Build condensed Hessian %%%
H = interior_point.get_condensed_hessian(W,sys.Q,sys.R,rompc.P,N);

% initialise parameters
k = 0;
p = p_init; % = beta_k
y = ones(q,1); % = y_k
l = ones(q,1); % = lambda_k

w_opt = zeros(size(W, 1), 1);
beta_opt = zeros(size(rompc.eq.Q1, 1) - N * (n + m), 1);

% solve KKT system iteratively
while true
    %%% Current control-state sequence %%%
    w = Z*g + W*p;

    %%% Residuals %%%
    eta = interior_point.get_condensed_gradient(w,W,sys.Q,sys.R,rompc.P,rompc.N);
    delta = interior_point.get_condensed_primal_feasibility(w,W,sys.Gx,sys.Gu,rompc.Gf, ...
        sys.fx, sys.fu, rompc.ff, N);
    rd = eta + interior_point.get_condensed_adjoint(l,W,sys.Gx,sys.Gu,rompc.Gf,N);
    rp = -(delta+y);

    %%% Stopping criterion %%%
    if (norm(rp,'Inf') < TOL) && (norm(rd,'Inf') < TOL)
        flag = 1;
        break;
    elseif k == MAX_ITER
        flag = -1;
        break;
    end

    %%% Build normal equaltion %%%
    % update L
    L = H ...
        + interior_point.get_condensed_quadratic_constraints(W,sys.Gx,sys.Gu,rompc.Gf,y,l,N);
    % decompose
    [Lc, chol_flag] = chol(L);

    if chol_flag ~= 0
        flag = -2;
        return;
    end

    % build rhs
    nu = (y .\ (l .* delta)) + l;
    theta = interior_point.get_condensed_adjoint(nu,W,sys.Gx,sys.Gu,rompc.Gf,N);
    % final rhs
    phat = -(eta + theta);

    %%% Find affine update directions %%%
    Dp_aff = interior_point.forward_backward_solve(Lc,phat);
    % initialise
    waff = w + W*Dp_aff;
    Dl_aff = interior_point.get_condensed_primal_feasibility(waff,W,sys.Gx,sys.Gu,rompc.Gf, ...
        sys.fx, sys.fu, rompc.ff, N);
    % affine slack step
    Dy_aff = -(Dl_aff + y);
    % affine dual step
    Dl_aff = (y .\ (l .* Dl_aff));

    %%% Centering parameter %%%
    mu = (y'*l) / q;
    alpha_aff = interior_point.multiplier_line_search([y;l],[Dy_aff;Dl_aff]);
    mu_aff = ((y + alpha_aff.*Dy_aff)'*(l + alpha_aff.*Dl_aff)) / q;
    sigma = (mu_aff / mu)^3;

    %%% Total step %%%
    % correction term
    kappa = y .\ (Dl_aff .* Dy_aff - sigma*mu);

    % total primal step
    phat = phat + ...
        interior_point.get_condensed_adjoint(kappa,W,sys.Gx,sys.Gu,rompc.Gf,N);
    Dp = interior_point.forward_backward_solve(Lc,phat);
    % initialise
    what = w + W*Dp;
    Dl = interior_point.get_condensed_primal_feasibility(what,W,sys.Gx,sys.Gu,rompc.Gf, ...
        sys.fx, sys.fu, rompc.ff, N);
    % total slack step
    Dy = -(Dl + y);
    % total dual step
    Dl = (y .\ (l .* Dl)) - kappa;

    %%% STEP LENGTH %%%
    tau = 0.5; % this can be improved on!
    alpha_primal = interior_point.multiplier_line_search(tau .* y, Dy);
    alpha_dual = interior_point.multiplier_line_search(tau .* l, Dl);
    alpha = min(alpha_dual,alpha_primal); %max(alpha_dual,alpha_primal);

    %%% UPDATE %%%
    p = p + alpha .* Dp;
    y = y + alpha .* Dy;
    l = l + alpha .* Dl;
    % increase step counter
    k = k+1;
end
% optimiser
w_opt = w;
beta_opt = rompc.eq.Q1(N*(n+m)+1:end,:)*g + rompc.eq.Q2(N*(n+m)+1:end,:)*p;
end