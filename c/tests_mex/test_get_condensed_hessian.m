function result = test_get_condensed_hessian(run)

    if run == "test"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        
        W = randn(N * (n + m), s);
        Q = randn(n, n);
        Q = Q * Q.';
        R = randn(m, m);
        R = R * R.';
        P = randn(n, n);
        P = P * P.';
        
        H_r = interior_point.get_condensed_hessian(W, Q, R, P, N);
        H = get_condensed_hessian(W, Q, R, P, N);        
        
        errors = norm(H_r - H) / norm(H_r);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        N = 30;
        n = 12;
        m = 3;
        s = 6;
        
        W = randn(N * (n + m), s);
        Q = randn(n, n);
        Q = Q * Q.';
        R = randn(m, m);
        R = R * R.';
        P = randn(n, n);
        P = P * P.';

        Nruns = 10000;
        tic()
        for i = 1:Nruns
            H_r = interior_point.get_condensed_hessian(W, Q, R, P, N);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            H = get_condensed_hessian(W, Q, R, P, N);
        end
        t = toc();
        
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(H_r - H) / norm(H_r);
        result.passed = all(result.errors < 1e-10);
    end

end
