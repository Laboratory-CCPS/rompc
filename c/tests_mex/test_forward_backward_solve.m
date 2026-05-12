function result = test_forward_backward_solve(run)

    if run == "test"
        %%
        G = [2, 4; 0, 3];
        b = [4; 12];
        
        xr = interior_point.forward_backward_solve(G, b);
        x = forward_backward_solve(G, b);

        errors = norm(xr - x) / norm(xr);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        n = 6;
        G = triu(randn(n, n));
        b = randn(n, 1);
        
        Nruns = 10000;
        tic()
        for i = 1:Nruns
            xr = interior_point.forward_backward_solve(G, b);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            x = forward_backward_solve(G, b);
        end
        t = toc();

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;

        result.errors = norm(xr - x) / norm(xr);
        result.passed = all(result.errors < 1e-10);
    end

end
