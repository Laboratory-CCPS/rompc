function result = test_get_condensed_quadratic_contraints(run)

    if run == "test"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        qx = 24;
        qu = 6;
        qf = 82;
        
        q = qu * N + qx * (N - 1) + qf;
        
        W = randn(N * (n + m), s);
        y = rand(q, 1) + 0.5;
        lambda = rand(q, 1) + 0.5;
        Gx = randn(qx, n);
        Gu = randn(qu, m);
        Gf = randn(qf, n);
        
        D_r = interior_point.get_condensed_quadratic_constraints(W, Gx, Gu, Gf, y, lambda, N);
        D = get_condensed_quadratic_constraints(W, Gx, Gu, Gf, y, lambda, N);
        
        errors = norm(D_r - D) / norm(D_r);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        qx = 24;
        qu = 6;
        qf = 82;
        
        q = qu * N + qx * (N - 1) + qf;
        
        W = randn(N * (n + m), s);
        y = rand(q, 1) + 0.5;
        lambda = rand(q, 1) + 0.5;
        Gx = randn(qx, n);
        Gu = randn(qu, m);
        Gf = randn(qf, n);

        Nruns = 1000;
        tic()
        for i = 1:Nruns
            D_r = interior_point.get_condensed_quadratic_constraints(W, Gx, Gu, Gf, y, lambda, N);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            D = get_condensed_quadratic_constraints(W, Gx, Gu, Gf, y, lambda, N);
        end
        t = toc();

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(D_r - D) / norm(D_r);
        result.passed = all(result.errors < 1e-10);
    end

end
