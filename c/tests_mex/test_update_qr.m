function result = test_update_qr(run)

    if run == "test"
        %%
        
        [ar, br] = condensing.update_qr([1, 1, 1; 2, 2, 2; 0, 3, 3; 0, 0, 4], [10, 10, 10; 20, 20, 20]);
        [a, b] = update_qr([1, 1, 1; 2, 2, 2; 0, 3, 3; 0, 0, 4], [10, 10, 10; 20, 20, 20]);
        
        errors = [norm(ar - a) / norm(ar), norm(br - b) / norm(br)];
        result.passed = all(errors < 1e-10);
        
    elseif run == "benchmark"
        %%
        cols = 90;
        rows = 455;
        
        Hcheck = hess(randn(cols, cols));
        Hcheck(end + 1, :) = 0;
        Hcheck(end, end) = randn();
        Qcheck2 = randn(rows, cols);
        
        Nruns = 10000;
        tic()
        for i = 1:Nruns
            [ar, br] = condensing.update_qr(Hcheck, Qcheck2);
        end
        tr = toc();
        
        tic();
        for i = 1:Nruns
            [a, b] = update_qr(Hcheck, Qcheck2);
        end
        t = toc();
        
        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;
        
        result.errors = [norm(ar - a) / norm(ar), norm(br - b) / norm(br)];
        result.passed = all(result.errors < 1e-10);
    end

end
