clc; clearvars; close all;

% This script creates the problem data used in paper_ifac_2026.m
% When completed successfully, this script overwrote the provided
% six_masses_condensed.mat with an equivalent representation of the
% example problem.

% Prerequisities: MPT3 toolbox (https://www.mpt3.org)

%% define system as in [Wang/Boyd 2010]
k = 1;          % spring constant
lam = 0;        % damping constant
a = -2 * k;
b = -2 * lam;
c = k;
d = lam;

n = 12; % state dimension
m = 3; % input dimension

Acts = [zeros(n / 2), eye(n / 2); ...
        [a, c, 0, 0, 0, 0, b, d, 0, 0, 0, 0; ...
         c, a, c, 0, 0, 0, d, b, d, 0, 0, 0; ...
         0, c, a, c, 0, 0, 0, d, b, d, 0, 0; ...
         0, 0, c, a, c, 0, 0, 0, d, b, d, 0; ...
         0, 0, 0, c, a, c, 0, 0, 0, d, b, d; ...
         0, 0, 0, 0, c, a, 0, 0, 0, 0, d, b]];

Bcts = [zeros(n / 2, m); ...
        [1,  0,  0; ...
        -1,  0,  0; ...
         0,  1,  0; ...
         0,  0,  1; ...
         0, -1,  0; ...
         0,  0, -1]];

% convert to discrete-time system
Ts = 0.5;       % sampling time
A = expm(Ts * Acts);
B = (Acts \ (A - eye(n))) * Bcts;
[n, m] = size(B);

% Quadratic state penalty with weighting matrix Q
Qc = eye(n);
% Quadratic input penalty with weighting matrix R
Rc = eye(m);

% box constraint on the state
xbar = repmat(4, n, 1);
Gx = [eye(n); -eye(n)];
fx = repmat(xbar, 2, 1);

% box constraint on the control value
ubar = repmat(0.5, m, 1);
Gu = [eye(m); -eye(m)];
fu = repmat(ubar, 2, 1);

%% terminal ingredients based on LQR

% termianl cost and controller
[Pf, Kf] = idare(A, B, Qc, Rc);
Kf = -Kf; % positive feedback sign

% set up LQR closed-loop system in MPT3
model = LTISystem('A', A + B * Kf);
% extend constraints for the invariant set to capture the control value
% constraints
Ge = [Gx; Gu * Kf];
fe = [fx; fu];

% compute max. positive invariant set
Xf = model.invariantSet('X', Polyhedron('A', Ge, 'b', fe));
Xf = Xf.minHRep();
% extract matrix represnetation
Gf = Xf.A;
ff = Xf.b;

%% Define MPC problem
% box initial set
X_0 = [0, -2, 0, 0, 0, 0, zeros(1, n / 2); ...
       0, 2, 0, 0, 0, 0, zeros(1, n / 2); ...
       0, 0, 0, -2, 0, 0, zeros(1, n / 2); ...
       0, 0, 0, 2, 0, 0, zeros(1, n / 2)]';

% prediction horizon
N = 30;

% pre-stabilise w/ terminal controller
K = Kf;
A_cl = A + B * K;

%% Define MPC problem in Yalmip

% states over the horizon
xi = sdpvar(n,  N + 1);

% controls over the horizon
omega = sdpvar(m, N);
zeta = sdpvar(m, N);

J = 0;
Con = [];
for t = 0:N-1
    % from time to index
    k = t + 1;
    % get control & state
    ut = omega(:, k); % = u_t
    zt = zeta(:, k); % = z_t
    xt = xi(:, k); % = x_t
    xt1 = xi(:, k + 1); % = x_{t+1}

    % add cost
    J = J + xt' * (Qc * xt) + ut' * (Rc * ut);

    % add continuity constraint
    Con = [Con,  xt1 == A * xt + B * ut]; %#ok<AGROW>
    % add prestabilising constraint
    Con = [Con,  ut == Kf * xt + zt]; %#ok<AGROW>
    % add state constraints
    Con = [Con,  Gx * xt <= fx]; %#ok<AGROW>
    % add control constrainits
    Con = [Con,  Gu * ut <= fu]; %#ok<AGROW>
end
% add terminal constraints
xN = xi(:, N + 1);
J = J + xN' * (Pf * xN);
Con = [Con, Gf * xN <= ff];

ops = sdpsettings('verbose',  0);
sol_yal = optimizer(Con, J, ops, xi(:, 1), zeta);

%% Data collection for PCA
M = 100;
X_cl = NaN(n, 0);
Z_cl = NaN(N * m, 0);
% sample initial states
for ii=1:M
    x24init = 1 + rand(2, 1);
    % construct initial state
    x_init = [0, x24init(1), 0, x24init(2), zeros(1, n / 2 + 2)]';
    % simulate
    x = x_init;
    while max(Gf * x-ff) > 1e-2 % while we're not in the terminal set...
        % solve ocp
        zeta_opt = sol_yal(x);
        z_opt = reshape(zeta_opt,  N * m,  1);
        % store
        X_cl(:, end + 1) = x; %#ok<SAGROW>
        Z_cl(:, end + 1) = z_opt; %#ok<SAGROW>
        % update state
        u = K * x + zeta_opt(:, 1);
        x = A * x + B * u;
    end
end

%% PCA subspace
r = 5;
b0 = mean(Z_cl, 2);
Del = Z_cl - b0;
[U_svd, si, ~] = svd(Del, 0, 'vector');
U = U_svd(:, 1:r);
% project offset to orthogonal complement of ran U
xi = b0 - U * (U' * b0);

%% condense the equality contraints offline

% initialise
C = zeros(N * n, N * (n + m));
D = zeros(N * m, N * (n + m));

% loop over block rows
rowIndex = 1;
colIndex = 1;
row1Index = 1;
col1Index = 1;

for k = 1:N
    % matrices to insert
    if k == 1
        F1 = [-B, eye(n)];
        F2 = eye(m);
    else
        F1 = [-A, -B, eye(n)];
        F2 = [-K, eye(m)];
    end
    q = size(F1, 2);
    l = size(F2, 2);
    % insert
    C(rowIndex : rowIndex + n - 1, colIndex : colIndex + q - 1) = F1;
    D(row1Index : row1Index + m - 1, col1Index : col1Index + l - 1) = F2;

    % update pointer
    rowIndex = rowIndex + n;
    colIndex = colIndex + q - n;
    row1Index = row1Index + m;
    col1Index = col1Index + l;
end

% build offline equality constraints
E = [C, zeros(N * n, r); ...
     -D, U];

% factor
[Qcheck, Rcheck] = qr(E');
Rcheck = Rcheck(1 : N * (n + m), :); % delete zero rows

%% store as structures

% data on system
sys = struct();
sys.A = A;
sys.B = B;
sys.Gx = Gx;
sys.fx = fx;
sys.Gu = Gu;
sys.fu = fu;
sys.Q = Qc;
sys.R = Rc;

% data on MPC
rompc = struct();
rompc.N = N;
rompc.K = Kf;
rompc.P = Pf;
rompc.Gf = Gf;
rompc.ff = ff;

% data on subspace
subspace = struct();
subspace.U = U;
subspace.xi = xi;

save('six_masses_condensed.mat', ...
     'sys', 'rompc', 'subspace', 'Qcheck', 'Rcheck');
