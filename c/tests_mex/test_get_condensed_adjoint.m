function result = test_get_condensed_adjoint(run)

    if run == "test"
        %%
        q = 958;
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        qf = 82;
        
        lambda = ones(q, 1);
        W = randn(N * (n + m), s);
        Gx = [eye(n); -eye(n)];
        Gu = [eye(m); -eye(m)];
        Gf = randn(qf, n);
        
        theta_r = interior_point.get_condensed_adjoint(lambda, W, Gx, Gu, Gf, N);
        theta = get_condensed_adjoint(lambda, W, Gx, Gu, Gf, N);
        
        errors = norm(theta_r - theta) / norm(theta_r);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        q = 958;
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        qf = 82;
        
        lambda = ones(q, 1);
        W = randn(N * (n + m), s);
        Gx = [eye(n); -eye(n)];
        Gu = [eye(m); -eye(m)];
        Gf = randn(qf, n);

        Nruns = 10000;
        tic()
        for i = 1:Nruns
            theta_r = interior_point.get_condensed_adjoint(lambda, W, Gx, Gu, Gf, N);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            theta = get_condensed_adjoint(lambda, W, Gx, Gu, Gf, N);
        end
        t = toc();

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(theta_r - theta) / norm(theta_r);
        result.passed = all(result.errors < 1e-10);
    end

end
