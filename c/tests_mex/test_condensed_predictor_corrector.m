function result = test_condensed_predictor_corrector(run)

    if run == "test"
        %%
        data = load("condensed_predictor_corrector_input.mat");
        s = 6;
        p_init = zeros(s, 1);

        [w_opt_r, beta_opt_r, flag_r] = interior_point.condensed_predictor_corrector(data.xt, data.zt, data.sys, data.rompc, p_init);
        [w_opt, beta_opt, flag, ~] = condensed_predictor_corrector(data.xt, data.zt, data.sys, data.rompc, p_init);

        errors = [norm(w_opt_r - w_opt) / norm(w_opt_r), norm(beta_opt_r - beta_opt) / norm(beta_opt_r), abs(flag_r - flag)];
        result.passed = all(errors < 1e-10);

    elseif run == "benchmark"
        %%
        data = load("condensed_predictor_corrector_input.mat");
        s = 6;
        p_init = zeros(s, 1);

        Nruns = 10;
        tr = 0;
        for i = 1:Nruns
            tic();
            [w_opt_r, beta_opt_r, flag_r] = interior_point.condensed_predictor_corrector(data.xt, data.zt, data.sys, data.rompc, p_init);
            tr = tr + toc();
        end

        t = 0;
        for i = 1:Nruns
            tic();
            [w_opt, beta_opt, flag] = condensed_predictor_corrector(data.xt, data.zt, data.sys, data.rompc, p_init);
            t = t + toc();
        end

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;

        result.errors = [norm(w_opt_r - w_opt) / norm(w_opt_r), norm(beta_opt_r - beta_opt) / norm(beta_opt_r), abs(flag_r - flag)];
        result.passed = all(result.errors < 1e-10);
    end
end
