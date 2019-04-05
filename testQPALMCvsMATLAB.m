clear;close all;clc
m = 50;n = 30;
rng(1)
A = sprandn(m, n, 5e-1,1e-1);
nnz(A)
% A = sparse(ones(m,n));
lb = -2*ones(m,1);
ub =  2*ones(m,1);
% Q = sparse(ones(n,n));
Q = sprandsym(n, 5e-1, 1e-1, 1); %Q=sparse(n,n);
q = 10*randn(n,1);

% full(A)
%% QPALM C
% 
solver = qpalm;
settings = solver.default_settings();
% settings.verbose = true;
settings.proximal = true;
settings.scaling = 10;
settings.max_iter = 1000;
settings.eps_abs = 1e-4;
settings.eps_rel = 1e-4;
settings.tau_init = 1.5;
solver.setup(Q, q, A, lb, ub, settings);
tic
res = solver.solve();
QPALMtime = toc
%% Quadprog

% tic
% [xs,fs,es] = quadprog(Q,q,[A;-A],[ub;-lb]);
% QPtime = toc;
% 
% xs-res.x

%% QPALM MATLAB
%Copy the settings
opts.Delta   = settings.delta;
opts.eps_abs = settings.eps_abs;
opts.eps_rel = settings.eps_rel;
opts.eps_abs_in = settings.eps_abs_in;
opts.eps_rel_in = settings.eps_rel_in;
opts.memory  = settings.memory;
opts.maxiter = settings.max_iter;
opts.rho     = settings.rho;
opts.theta   = settings.theta;
opts.scaling = 'simple';
opts.scaling_iter = settings.scaling;

% opts.solver  = 'lbfgs';
opts.solver = 'newton';
opts.scalar_sig = false;
opts.lbfgs_precon = false;
opts.proximal = settings.proximal;
% opts.scalar_sig = true;
tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);toc
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
QPALMCfeas = norm([min(A*res.x-lb,0);min(ub-A*res.x,0)],inf);
QPALMfeas = norm([min(A*x_qpalm-lb,0);min(ub-A*x_qpalm,0)],inf);
QPALMCobj = 1/2*res.x'*Q*res.x + q'*res.x;
QPALMobj = 1/2*x_qpalm'*Q*x_qpalm + q'*x_qpalm;
% 
fprintf('           |   QPALM (C)   |   QPALM  \n')
fprintf('Iterations |   %3d      |    %3d   \n',...
    res.info.iter,...
    stats_qpalm.iter...
    )
fprintf('Violation  | %3.2e |  %3.2e \n', ...
    QPALMCfeas, ...
    QPALMfeas...
    )
fprintf('Objective  | %3.2e |  %3.2e \n', ...
    QPALMCobj, ...
    QPALMobj...    
    )


