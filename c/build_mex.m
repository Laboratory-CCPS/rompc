function build_mex()
    clc();

    % specify mex flags
    flags = {};
    % flags{end + 1} = '-g';     % debug mode (activates assertions)
    % flags{end + 1} = '-v';     % verbose

    % specify compiler flags (without the leading "-" or "/")
    add_comp_flags = {};
    % add_comp_flags{end + 1} = 'Wall';

    % always use this flag!
    add_comp_flags{end + 1} = 'DCOMPILE_MEX_FUNCTION';

    outputdir = "../mex";

    mex_files = { ...
        "./src/mexinterface/update_qr.c", ...
        "./src/mexinterface/forward_backward_solve.c", ...
        "./src/mexinterface/multiplier_line_search.c", ...
        "./src/mexinterface/get_condensed_adjoint.c", ...
        "./src/mexinterface/get_condensed_gradient.c", ...
        "./src/mexinterface/get_condensed_hessian.c", ...
        "./src/mexinterface/get_condensed_primal_feasibility.c", ...
        "./src/mexinterface/get_condensed_quadratic_constraints.c", ...
        "./src/mexinterface/condensed_predictor_corrector.c", ...
        "./src/mexinterface/condensed_rompc_step.c", ...
    };

    add_source_files = {"./src/algorithm/alg.c", "./src/algorithm/linalg.c"};
    add_include_directories = {"./src", "C:/Program Files/MATLAB/R2025a/extern/include"};
    add_libraries = {"mwblas", "mwlapack"};

    cc = mex.getCompilerConfigurations("c");

    args = {};

    for i = 1:length(flags)
        args{end + 1} = flags{i}; %#ok<AGROW>
    end

    if ~isempty(outputdir)
        args{end + 1} = "-outdir";
        args{end + 1} = sprintf('"%s"', outputdir);
    end

    for i = 1:length(add_source_files)
        args{end + 1} = add_source_files{i}; %#ok<AGROW>
    end

    for i = 1:length(add_include_directories)
        args{end + 1} = sprintf('-I"%s"', add_include_directories{i}); %#ok<AGROW>
    end

    for i = 1:length(add_libraries)
        args{end + 1} = sprintf('-l"%s"', add_libraries{i}); %#ok<AGROW>
    end

    if ~isempty(add_comp_flags)
        if strcmp(cc.Manufacturer, 'GNU')
            cflags = [' -', strjoin(add_comp_flags, ' -')];
            args{end + 1} = sprintf('CFLAGS="$CFLAGS %s"', cflags);
        
        elseif strcmp(cc.Manufacturer, 'Microsoft')
            cflags = [' /', strjoin(add_comp_flags, ' /')];
            args{end + 1} = sprintf('COMPFLAGS="$COMPFLAGS %s"', cflags);
        
        else
            error("unknown compiler manufacturer: %s", cc.Manufacturer);
        end
    end


    for i = 1:length(mex_files)
        fprintf('(%d/%d) "%s"\n', i, length(mex_files), mex_files{i});
        mex(mex_files{i}, args{:});
    end
end
