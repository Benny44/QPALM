%A demo to show how the matlab interface for qpalm is called for a
%nonconvex QP

%% Generate data
m = 300;n = 50;
% m = 1; n = 2000;
A = sprandn(m, n, 1e-1, 1e-4); 

lb = -2*ones(m,1);
ub =  2*ones(m,1);
% ub = lb;
Q = sprandsym(n, 1e-1, 1e-4);
q = 10*randn(n,1);


%% Solve with qpalm mex
solver = qpalm;
settings = solver.default_settings();
%IMPORTANT: set nonconvex to true for nonconvex QPs
settings.nonconvex = true;

settings.delta = 100;
settings.proximal = true;
settings.scaling = 2;
settings.max_iter = 5000;
settings.eps_abs = 1e-4;
settings.eps_rel = 1e-10;
settings.gamma_init = 1e1;
settings.gamma_max = 1e7;
settings.verbose = false;
% settings.reset_newton_iter = 1;

solver.setup(Q, q, A, lb, ub, settings); 
res = solver.solve();

%% Solve with qpalm matlab
opts.nonconvex = true;
opts.nonconvex_approx = false;

opts.Delta   = 100;
opts.eps_abs = 1e-4;
opts.eps_rel = 1e-10;
opts.maxiter = 5000;
opts.scaling = 'simple';
opts.scaling_iter = 2;
opts.gamma = 1e1;
opts.gammaMax = 1e7;
opts.verbose = false;
% opts.reset_newton_iter = 1;

opts.solver  = 'newton';
opts.scalar_sig = false;
opts.proximal = true;
tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);qpalm_time = toc;

%% Compare
format longe
fprintf('           |      QPALM mex       |     QPALM Matlab \n') 
fprintf('Iterations |  %10d   \t  |  %10d \n', res.info.iter, stats_qpalm.iter-1);
fprintf('Pri_res    | %.13e  | %.13e \n', res.info.pri_res_norm, stats_qpalm.nrm_rp(end));
fprintf('Dua_res    | %.13e  | %.13e \n', res.info.dua_res_norm, stats_qpalm.nrm_rd(end));
fprintf('Time (ms)  |      %.7f       |      %.7f \n', res.info.run_time*1000, qpalm_time*1000);
fprintf('Objective  | %.13e  | %.13e \n', 0.5*res.x'*Q*res.x + q'*res.x, 0.5*x_qpalm'*Q*x_qpalm + q'*x_qpalm )
fprintf('Status     | %20s | %20s \n', res.info.status, stats_qpalm.status);
fprintf('\n Difference in x: %.16e\n', norm(res.x-x_qpalm,inf)/norm(res.x,inf));

solver.delete();





