function result = test_multiplier_line_search(run)

    if run == "test"
        %%
        nu = [1, 2, 0.2].';
        Delta_nu = [1, 1, 1].';
        
        alphar = interior_point.multiplier_line_search(nu, Delta_nu);
        alpha = multiplier_line_search(nu, Delta_nu);
        
        errors = norm(alphar - alpha) / norm(alphar);
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        n = 1916;
        nu = rand(n, 1) + 0.5;
        Delta_nu = randn(n, 1);
        
        
        Nruns = 10000;
        tic()
        for i = 1:Nruns
            alphar = interior_point.multiplier_line_search(nu, Delta_nu);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            alpha = multiplier_line_search(nu, Delta_nu);
        end
        t = toc();
        
        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = norm(alphar - alpha) / norm(alphar);
        result.passed = all(result.errors < 1e-10);
    end

end
