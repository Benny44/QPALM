%A test to to compare the qpalm mex and matlab versions 
warning ('off','all');
% warning('off','MATLAB:mex:GccVersion_link')

qpalm_make;

for k=1:100
    %% Generate data
    m = 300;n = 50;
    A = sprandn(m, n, 1e-1, 1e-4); 

    lb = -2*ones(m,1);
    ub =  2*ones(m,1);
    Q = sprandsym(n, 1e-1, 1e-4, 1);
    q = 10*randn(n,1);

    % load('qpalm_diff9')
    %% Solve with qpalm mex
    solver = qpalm;
    settings = solver.default_settings();
    settings.delta = 10;
    settings.proximal = true;
    settings.scaling = 2;
    settings.max_iter = 300;
    settings.eps_abs = 1e-4;
    settings.eps_rel = 1e-4;
    settings.gamma_max = 1e7;

    % settings.print_iter = 10;
    % settings.reset_newton_iter = 1000;
    settings.verbose = false;
    % settings.sigma_max = 1e0;

    solver.setup(Q, q, A, lb, ub, settings); 
    res = solver.solve();

    %% Solve with qpalm matlab
    opts.Delta   = 10;
    opts.eps_abs = 1e-4;
    opts.eps_rel = 1e-4;
    opts.maxiter = 300;
    opts.scaling = 'simple';
    opts.scaling_iter = 2;
    opts.gammaMax = 1e7;

    % opts.reset_newton_iter = 1000;
    % opts.sigma_max = 1e9;

    opts.solver  = 'newton';
    opts.scalar_sig = false;
    opts.proximal = true;
    tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);qpalm_time = toc;

    %% Compare
%     format longe
%     fprintf('           |      QPALM mex       |     QPALM Matlab \n') 
%     fprintf('Iterations |  %10d   \t  |  %10d \n', res.info.iter, stats_qpalm.iter-1);
%     fprintf('Pri_res    | %.13e  | %.13e \n', res.info.pri_res_norm, stats_qpalm.nrm_rp(end));
%     fprintf('Dua_res    | %.13e  | %.13e \n', res.info.dua_res_norm, stats_qpalm.nrm_rd(end));
%     fprintf('Time (ms)  |      %.7f       |      %.7f \n', res.info.run_time*1000, qpalm_time*1000);
%     fprintf('\n Difference in x: %.16e\n', norm(res.x-x_qpalm,inf));

    solver.delete();

    if(norm(res.x-x_qpalm,inf)/norm(res.x, inf) > 1e-10)
        fprintf('\n\n ERROR: Mex and Matlab versions might have diverged\n\n ');
        error('');
    end

end

fprintf('\n\n SUCCESS: Mex and Matlab do the same thing!\n\n');
