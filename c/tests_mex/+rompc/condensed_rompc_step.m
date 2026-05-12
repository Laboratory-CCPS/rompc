function [ut, zt, flag, net_duration] = condensed_rompc_step(xt, zt, Rcheck, Qcheck, sys, rompc, subspace) %# codegen

    n = length(sys.Q);
    m = length(sys.R);
    N = rompc.N;
    r = size(subspace.U,2);
    s = r+1;

    ut = zeros(3, 1);
    net_duration = [nan, nan, nan];

    % update equality constraints
    z0 = subspace.xi - zt;
    % upper Hessenberg matrix
    Hcheck = [z0'; Rcheck(N*n+1:end,N*n+1:end)];

    % update factorisation
    tic();
    [Rtilde,Qhat2] = condensing.update_qr(Hcheck,Qcheck(:,N*n+1:N*(n+m)));
    net_duration(1) = toc();

    Qhat1 = [Qcheck(:,1:N*n); zeros(1,N*n)];
    Qhat3 = [Qcheck(:,N*(n+m)+1:end); zeros(1,r)];

    % form full matrices
    % eq = struct();
    Q = [Qhat1, Qhat2, Qhat3];
    % rompc.eq.R = [Rcheck(1:N*n,:); ...
    %     zeros(size(Rtilde,1),N*n), Rtilde];
    rompc.eq.R(1:N*n,:) = Rcheck(1:N*n,:);
    rompc.eq.R(N*n+1:end,1:N*n) = 0;
    rompc.eq.R(N*n+1:end,N*n+1:end) = Rtilde;

    % partition
    rompc.eq.Q1(:,1:N*(n+m)) = Q(:,1:N*(n+m));
    rompc.eq.Q2(:,1:s) = Q(:,N*(n+m)+1:end);
    net_duration(2) = toc() - net_duration(1);

    % store
    %rompc.eq = eq;

    tic();
    [w_opt,beta_opt,flag] = interior_point.condensed_predictor_corrector(xt,zt,sys,rompc,zeros(s,1));
    net_duration(3) = toc() - net_duration(2);

    if flag ~= 1
        return;
    end

    % extract control
    ut = w_opt(1:m);

    % recover optimal pre-stabilised trajectory
    z_opt = [subspace.U, z0]*beta_opt + zt;
    % update admissible guess
    zt(1:end-m) = z_opt(m+1:end);
    zt(end-m+1:end) = zeros(m,1); % works iff K = K_f
end
