function result = test_condensed_rompc_step(run)

    if run == "test"
        %%
        data = load("condensed_rompc_step_input.mat");

        [ut_r, zt_r, flag_r, ~] = rompc.condensed_rompc_step(data.xt, data.zt, data.Rcheck, data.Qcheck, data.sys, data.rompc, data.subspace);
        [ut, zt, flag, ~] = condensed_rompc_step(data.xt, data.zt, data.Rcheck, data.Qcheck, data.sys, data.rompc, data.subspace);


        errors = [norm(ut_r - ut) / norm(ut_r), norm(zt_r - zt) / norm(zt_r), norm(flag_r - flag)];
        result.passed = all(errors < 1e-10);

    elseif run == "benchmark"
        %%
        data = load("condensed_rompc_step_input.mat");

        Nruns = 10;
        tr = 0;
        for i = 1:Nruns
            tic();
            [ut_r, zt_r, flag_r, ~] = rompc.condensed_rompc_step(data.xt, data.zt, data.Rcheck, data.Qcheck, data.sys, data.rompc, data.subspace);
            tr = tr + toc();
        end

        t = 0;
        for i = 1:Nruns
            tic();
            [ut, zt, flag, ~] = condensed_rompc_step(data.xt, data.zt, data.Rcheck, data.Qcheck, data.sys, data.rompc, data.subspace);
            t = t + toc();
        end

        result = [];
        result.t_ml = tr / Nruns;
        result.t_mex = t / Nruns;

        result.errors = [norm(ut_r - ut) / norm(ut_r), norm(zt_r - zt) / norm(zt_r), norm(flag_r - flag)];
        result.passed = all(result.errors < 1e-10);
    end
end
