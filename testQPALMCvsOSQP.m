% clear;close all;clc
m = 3000;n = 500;
% rng(6)
A = sprandn(m, n, 1e-1, 1e-4); 

lb = -2*ones(m,1);
ub =  2*ones(m,1);
Q = sprandsym(n, 1e-1, 1e-4, 1); %Q=sparse(n,n);
q = 10*randn(n,1);

clear;close all;clc
m = 500;n = 100;
A = sprandn(m, n, 1e-1,1e-8); 
lb = -1e3*rand(m,1);
ub =  rand(m,1);
Q = sprandsym(n, 9e-1, 1e-8, 1); 
% Q=sparse(n,n);
q = 100*randn(n,1);

load('dual-bug.mat')

fprintf('nnz A: %d, nnz Q: %d\n', nnz(A), nnz(Q));

%% QPALM C
% 
solver = qpalm;
settings = solver.default_settings();
settings.verbose = false;
settings.scaling = 10;
settings.max_iter = 10000;
settings.eps_abs = 1e-4;
settings.eps_rel = 1e-4;
settings.delta   = 10;
% settings.memory  = 10;
settings.proximal = false;
solver.setup(Q, q, A, lb, ub, settings);
tic
% profile on;
res_qpalm = solver.solve();
% profile viewer;
QPALMCtime = toc;
%% Quadprog

% tic
% [xs,fs,es] = quadprog(Q,q,[A;-A],[ub;-lb]);
% QPtime = toc;
% 
% xs-res.x


%% OSQP
% 
solver = osqp;
osqp_settings = solver.default_settings();
% settings.verbose = true;
% osqp_settings.scaling = 10;
osqp_settings.scaling = settings.scaling;
osqp_settings.max_iter = settings.max_iter;
osqp_settings.eps_abs = settings.eps_abs;
osqp_settings.eps_rel = settings.eps_rel;
osqp_settings.verbose = settings.verbose; osqp_settings.verbose = false;
% osqp_settings.adaptive_rho = 1;
% osqp_settings.adaptive_rho_tolerance = 5;
% osqp_settings.linsys_solver ='mkl pardiso';
% osqp_settings.check_termination = 25;
solver.setup(Q, q, A, lb, ub, osqp_settings);
tic
res_osqp = solver.solve();
OSQPtime = toc;

% %% QPALM MATLAB
%Copy the settings
opts.Delta   = 10;
opts.eps_abs = settings.eps_abs;
opts.eps_rel = settings.eps_rel;
opts.eps_abs_in = settings.eps_abs_in;
opts.eps_rel_in = settings.eps_rel_in;
% opts.memory  = settings.memory;
opts.maxiter = settings.max_iter;
opts.rho     = settings.rho;
opts.theta   = settings.theta;
opts.scaling = 'simple';
opts.scaling_iter = settings.scaling; 
% opts.scaling = 2;

opts.solver  = 'newton';
% opts.solver = 'newton';
opts.scalar_sig = false;
opts.lbfgs_precon = false;
opts.proximal = settings.proximal;
% opts.scalar_sig = true;
tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);qpalm_time = toc
display(stats_qpalm.status)
% 
% opts.proximal = true;
% opts.gamma    = 1e4;
% opts.gammaUpd = 10;
% opts.gammaMax = 1e8;
% % opts.scalar_sig = true;
% tic;[x_qpalm2,y_qpalm2,stats_qpalm2] = qpalm(Q,q,A,lb,ub,[],[],opts);toc
% % stats_qpalm2.status
% %% OUTPUT
% 
% OSQPfeas = norm([min(A*res.x-lb,0);min(ub-A*res.x,0)],inf);
% QPALMfeas = norm([min(A*x_qpalm-lb,0);min(ub-A*x_qpalm,0)],inf);
% QPALMfeas2 = norm([min(A*x_qpalm2-lb,0);min(ub-A*x_qpalm2,0)],inf);
% 
% OSQPobj = 1/2*res.x'*Q*res.x + q'*res.x;
% QPALMobj = 1/2*x_qpalm'*Q*x_qpalm + q'*x_qpalm;
% QPALMobj2 = 1/2*x_qpalm2'*Q*x_qpalm2 + q'*x_qpalm2;
% 
QPALMCfeas = norm([min(A*res_qpalm.x-lb,0);min(ub-A*res_qpalm.x,0)],inf);
QPALMfeas = norm([min(A*res_qpalm.x-lb,0);min(ub-A*res_qpalm.x,0)],inf);
OSQPfeas = norm([min(A*x_qpalm-lb,0);min(ub-A*x_qpalm,0)],inf);
QPALMCobj = 1/2*res_qpalm.x'*Q*res_qpalm.x + q'*res_qpalm.x;
QPALMobj = 1/2*x_qpalm'*Q*x_qpalm + q'*x_qpalm;
OSQPobj = 1/2*res_osqp.x'*Q*res_osqp.x + q'*res_osqp.x;

fprintf('           |    QPALM-C   |    OSQP         |   QPALM  \n')
fprintf('Iterations |    %3d       |    %3d          |    %3d   \n',...
    res_qpalm.info.iter,...
    res_osqp.info.iter,...
    stats_qpalm.iter)
fprintf('Runtime    |   %3.2e   |    %3.2e     |    %3.2e   \n',...
    res_qpalm.info.run_time,...
    res_osqp.info.run_time,...
    qpalm_time)
fprintf('Setup time |   %3.2e   |    %3.2e     |    %3.2e   \n',...
    res_qpalm.info.setup_time,...
    res_osqp.info.setup_time,...
    0)
fprintf('Solve time |   %3.2e   |    %3.2e     |    %3.2e   \n',...
    res_qpalm.info.solve_time,...
    res_osqp.info.solve_time,...
    0)
fprintf('Violation  |   %3.2e   |    %3.2e     |    %3.2e   \n', ...
    QPALMCfeas, ...
    OSQPfeas,...
    0)
fprintf('Objective  | %3.2e    |    %3.2e    |  %3.2e \n', ...
    QPALMCobj, ...
    OSQPobj,...
    QPALMobj)


