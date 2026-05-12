function result = test_get_condensed_primal_feasibility(run)

    if run == "test"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        qx = 24;
        qu = 6;
        qf = 82;
        
        w = randn(N * (n + m), 1);
        W = randn(N * (n + m), s);
        Gx = randn(qx, n);
        Gu = randn(qu, m);
        Gf = randn(qf, n);
        fx = randn(qx, 1);
        fu = randn(qu, 1);
        ff = randn(qf, 1);        

        delta_r = interior_point.get_condensed_primal_feasibility(w, W, Gx, Gu, Gf, fx, fu, ff, N);
        delta = get_condensed_primal_feasibility(w, W, Gx, Gu, Gf, fx, fu, ff, N);
        
        errors = norm(delta_r - delta) / norm(delta_r);
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
        
        w = randn(N * (n + m), 1);
        W = randn(N * (n + m), s);
        Gx = randn(qx, n);
        Gu = randn(qu, m);
        Gf = randn(qf, n);
        fx = randn(qx, 1);
        fu = randn(qu, 1);
        ff = randn(qf, 1);

        Nruns = 10000;
        tic()
        for i = 1:Nruns
            delta_r = interior_point.get_condensed_primal_feasibility(w, W, Gx, Gu, Gf, fx, fu, ff, N);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            delta = get_condensed_primal_feasibility(w, W, Gx, Gu, Gf, fx, fu, ff, N);
        end
        t = toc();
        
        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(delta_r - delta) / norm(delta_r);
        result.passed = all(result.errors < 1e-10);
    end

end
