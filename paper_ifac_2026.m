clc; clearvars;
%% load the offline problem data

load('six_masses_condensed.mat');

%% dimensions

n = length(sys.Q);
m = length(sys.R);
N = rompc.N;
r = size(subspace.U, 2);
s = r + 1;
qu = size(sys.Gu, 1);
qx = size(sys.Gx, 1);
qf = size(rompc.Gf, 1);
q = N * qu + (N - 1) * qx + qf;


%% original problem in Yalmip

ops = sdpsettings('verbose', 0);
K_f = -dlqr(sys.A, sys.B, sys.Q, sys.R);

% states over the horizon
xi = sdpvar(n, N + 1);

% controls over the horizon
omega = sdpvar(m, N);
zeta = sdpvar(m, N);

J = 0;
Con = [];
for t = 0:rompc.N-1
    % from time to index
    k = t + 1;
    % get control & state
    ut = omega(:, k); % = u_t
    zt = zeta(:, k); % = z_t
    xt = xi(:, k); % = x_t
    xt1 = xi(:, k + 1); % = x_{t+1}

    % add cost
    J = J + xt' * (sys.Q * xt) + ut' * (sys.R * ut);

    % add continuity constraint
    Con = [Con, xt1 == sys.A * xt + sys.B * ut]; %#ok<AGROW>
    % add prestabilising constraint
    Con = [Con, ut == K_f * xt + zt]; %#ok<AGROW>
    % add state constrainits
    Con = [Con, sys.Gx * xt <= sys.fx]; %#ok<AGROW>
    % add control constrainits
    Con = [Con, sys.Gu * ut <= sys.fu]; %#ok<AGROW>
end

% add terminal constraints
xN = xi(:, N + 1);
J = J + xN' * (rompc.P * xN);
Con = [Con, rompc.Gf * xN <= rompc.ff];

sol_yal = optimizer(Con, J, ops, xi(:, 1), omega(:, 1));

%% MPC loop

% number of mpc steps (time steps) simulated for each run
nsteps = 50;

n_warmup_runs = 2;
n_mpc_runs = 20;

n_max_run_attemps = ceil(1.5 * (n_warmup_runs + n_mpc_runs));

all_solver_times = nan(2, nsteps, n_mpc_runs);

i_successfull_runs = 0;

for i_run = 1:n_max_run_attemps

    % define initial state & admissible guess
    % (masses two and four are displaced)
    x_init = [ 0; 1.5; 0; 1; zeros(n - 4, 1) ];
    
    x_init = x_init + (0.2 * rand(12, 1) - 0.1);

    zt = zeros(N * m, 1);

    xt = x_init;

    t = 0;
    Xcl = NaN(length(xt), nsteps + 1);
    Xcl(:, 1) = xt;

    Ucl = NaN(m, nsteps);
    t_mpc = -1;
    costs = NaN(1, nsteps);
    err = NaN(1, nsteps);
    solver_times = zeros(2, nsteps);

    beta_init = [zeros(r, 1); 1];
    w_init = zeros(N * (n + m), 1);
    y_init = ones(q, 1);
    l_init = ones(q, 1);

    success = true;
    for i_step = 1:nsteps

        tic();
        [ut, zt, flag, stats] = condensed_rompc_step(xt, zt, Rcheck, Qcheck, sys, rompc, subspace);
        t_solve_rompc = toc();
        solver_times(1, i_step) = t_solve_rompc;

        if flag < 0
            fprintf('Interior point solver did not converge at time t = %d: %s\n', t, stats.msg);

            if i_step > 1
                error('Interior point solver couldn''t find a solution for t > 0');
            end

            success = false;
            break;
        end

        % % solve!
        tic();
        u_yal = sol_yal{xt};
        t_solve_yalmip = toc();
        solver_times(2, i_step) = t_solve_yalmip;

        % store cost
        costs(i_step) = ut' * (sys.R * ut) + xt' * (sys.Q * xt);

        % apply control
        xt = sys.A * xt + sys.B * ut;

        % store
        Xcl(:, i_step + 1) = xt;
        Ucl(:, i_step) = ut;
        t = t + 1;
    end

    if ~success
        continue;
    end

    i_successfull_runs = i_successfull_runs + 1;

    if i_successfull_runs > n_warmup_runs
        all_solver_times(:, :, i_successfull_runs - n_warmup_runs) = solver_times;
    end

    if i_successfull_runs == n_warmup_runs + n_mpc_runs
        break;
    end
end

if i_successfull_runs < n_warmup_runs + n_mpc_runs
    error('Too many failed runs');
end

%%

% get actual solver used by yalmip
solver_tag = regexp(evalc('sol_yal'), 'Solver:\s*(\S+)', 'tokens');
solver_tag = solver_tag{1}{1};

solver_times = mean(all_solver_times, 3);
solver_times_std = std(all_solver_times, [], 3);

tx = (0:t-1);

our_mean_ms = solver_times(1, :) .* 1e3;
yalmip_mean_ms = solver_times(2, :) .* 1e3;

our_lower_ms = (solver_times(1, :) - solver_times_std(1, :)) * 1e3;
our_upper_ms = (solver_times(1, :) + solver_times_std(1, :)) * 1e3;

yalmip_lower_ms = (solver_times(2, :) - solver_times_std(2, :)) * 1e3;
yalmip_upper_ms = (solver_times(2, :) + solver_times_std(2, :)) * 1e3;

% Create vectors for manual "stairs" plot to use them with the patch
% function. (These vectors end in "_s".)
tx_s = reshape([tx(1:end-1); tx(2:end)], 1, []);

our_lower_ms_s = kron(our_lower_ms(1:end-1), [1, 1]);
our_upper_ms_s = kron(our_upper_ms(1:end-1), [1, 1]);

yalmip_lower_ms_s = kron(yalmip_lower_ms(1:end-1), [1, 1]);
yalmip_upper_ms_s = kron(yalmip_upper_ms(1:end-1), [1, 1]);

figure();
hold('on');
stairs(tx, our_mean_ms);
p = patch(...
    'XData', [tx_s, flip(tx_s)], ...
    'YData', [our_lower_ms_s, flip(our_upper_ms_s)], ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'SeriesIndex', 1, ...
    'FaceAlpha', 0.2);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';

stairs(tx, yalmip_mean_ms);
p = patch(...
    'XData', [tx_s, flip(tx_s)], ...
    'YData', [yalmip_lower_ms_s, flip(yalmip_upper_ms_s)], ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'none', ...
    'SeriesIndex', 2, ...
    'FaceAlpha', 0.2);
p.Annotation.LegendInformation.IconDisplayStyle = 'off';

xlabel('t');
ylabel('solvertime (ms)');
legend({'our', solver_tag});
ylims = ylim();
ylim([0, ylims(2)]);
