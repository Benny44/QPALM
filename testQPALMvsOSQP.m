clear;close all;clc
m = 1000;n = 100;
A = sprandn(m, n, 1e-1,1e-8); 
lb = -1e3*rand(m,1);
ub =  rand(m,1);
Q = sprandsym(n, 9e-1, 1e-8, 1); 
% Q=sparse(n,n);
q = 100*randn(n,1);

%% OSQP
% 
solver = osqp;
settings = solver.default_settings();
settings.verbose = true;
settings.scaling = 10;
settings.max_iter = 10000;
settings.eps_abs = 1e-6;
settings.eps_rel = 1e-6;
% settings.adaptive_rho = false;
% settings.sigma = 0;
% settings.alpha = 1;
solver.setup(Q, q, A, lb, ub, settings);
tic
res = solver.solve();
OSQPtime = toc;

%% Quadprog

% tic
% [xs,fs,es] = quadprog(Q,q,[A;-A],[ub;-lb]);
% QPtime = toc;

%% QPALM
opts.Delta   = 10;
opts.eps_abs = 1e-6;
opts.eps_rel = 1e-6;
% opts.eps_abs_in = 1e-1;
% opts.eps_rel_in = 1e-1;
opts.solver  = 'newton';
opts.memory  = 10;
opts.maxiter = 10000;
opts.scaling = 'simple';
opts.scaling_iter = 0;
opts.scalar_sig = false;
opts.lbfgs_precon = false;
opts.proximal = true;
opts.gamma    = 1e4;
opts.gammaUpd = 10;
opts.gammaMax = 1e8;
tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);QPALMtime=toc;
display(stats_qpalm.status)

opts.Delta   = 10;
opts.scaling_iter = 1;
opts.solver = 'newton';
tic;[x_qpalm2,y_qpalm2,stats_qpalm2] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);QPALM2time=toc;
% stats_qpalm2.status
%% OUTPUT

OSQPfeas = norm([min(A*res.x-lb,0);min(ub-A*res.x,0)],inf);
QPALMfeas = norm([min(A*x_qpalm-lb,0);min(ub-A*x_qpalm,0)],inf);
QPALMfeas2 = norm([min(A*x_qpalm2-lb,0);min(ub-A*x_qpalm2,0)],inf);

OSQPobj = 1/2*res.x'*Q*res.x + q'*res.x;
QPALMobj = 1/2*x_qpalm'*Q*x_qpalm + q'*x_qpalm;
QPALMobj2 = 1/2*x_qpalm2'*Q*x_qpalm2 + q'*x_qpalm2;

fprintf('           |   OSQP       |   QPALM         |   QPALM2\n')
fprintf('Iterations |   %3d        |    %3d          |    %3d\n',...
    res.info.iter,...
    stats_qpalm.iter,...
    stats_qpalm2.iter...
    )
fprintf('Time       |  %3.2e    |   %3.2e      |  %3.2e\n',...
    OSQPtime,...
    QPALMtime,...
    QPALM2time...
    )

fprintf('Violation  |  %3.2e    |   %3.2e      |  %3.2e\n', ...
    OSQPfeas, ...
    QPALMfeas,...
    QPALMfeas2...    
    )
fprintf('Objective  | %3.2e    |  %3.2e      |  %3.2e\n', ...
    OSQPobj, ...
    QPALMobj,...
    QPALMobj2...    
    )


