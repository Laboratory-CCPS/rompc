function result = test_get_condensed_gradient(run)

    if run == "test"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        
        w = randn(N * (n + m), 1);
        W = randn(N * (n + m), s);
        Q = eye(n);
        R = eye(m);
        P = eye(n);
        
        eta_r = interior_point.get_condensed_gradient(w, W, Q, R, P, N);
        eta = get_condensed_gradient(w, W, Q, R, P, N);            
        
        errors = norm(eta_r - eta) / norm(eta_r);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        
        w = randn(N * (n + m), 1);
        W = randn(N * (n + m), s);
        Q = eye(n);
        R = eye(m);
        P = eye(n);

        Nruns = 10000;
        tic()
        for i = 1:Nruns
            eta_r = interior_point.get_condensed_gradient(w, W, Q, R, P, N);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            eta = get_condensed_gradient(w, W, Q, R, P, N);
        end
        t = toc();

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(eta_r - eta) / norm(eta_r);
        result.passed = all(result.errors < 1e-10);
    end

end
