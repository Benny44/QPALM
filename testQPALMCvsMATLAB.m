clear;close all;clc
m = 500;n = 300;
% rng(1)
A = sprandn(m, n, 5e-1,1e-1);
% A = sparse(ones(m,n));
lb = -2*ones(m,1);
ub =  2*ones(m,1);
% Q = sparse(ones(n,n));
Q = sprandsym(n, 5e-1, 1e-1, 1); %Q=sparse(n,n);
q = 10*randn(n,1);

% load('dual-bug-nl.mat')


x0 = [18.141378485712202
  36.896595037645973
   2.582191656442478
   9.253341674991914
   1.149763491253560
 -10.369845078121662
  -1.775841906763338
 -18.950809638303912
   9.727965340368442
  -4.711294877824136
  13.423053204405463
   4.372518471362451
   8.638069164603243
 -16.912721174120403
  37.696702236698634
 -15.203181801217713
  72.141375485189641
   9.271831651995498
 -10.120014445106698
   3.280167865124274
   5.166872269095653
   0.171050963269970
  -6.256644774019049
  25.007290602139534
 -29.394588506027084
 -36.780109978661081
 -37.902675481941095
  22.342733346164763
   6.079305025380402
   1.639091818564381];

y0 = [0
  -0.281967068395563
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
  -0.088170297273864
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
   0.007317885802707
                   0
   2.188830991476911
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
   0.152287493281371
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
   0.026971053167083
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
   0.073767167982456
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
   0.007317885802707
                   0
  -0.088170297273864
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
  -0.088170297273864
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.264283251391454
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   2.012369632548829
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
   0.435142150036130
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
   0.033242498093382
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -3.201221626658378
                   0
                   0
                   0
                   0
                   0
  -0.088170297273864
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
  -0.000261008765564
                   0
                   0
                   0
                   0
                   0
   0.007317885802707
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0
                   0]*1e2;

% full(A)
%% QPALM C
% 
solver = qpalm;
settings = solver.default_settings();
% settings.verbose = true;
settings.proximal = false;
settings.scaling = 2;
settings.max_iter = 100;
settings.eps_abs = 1e-4;
settings.eps_rel = 1e-4;
settings.tau_init = 1.5;

% tic
solver.setup(Q, q, A, lb, ub, settings); 
res = solver.solve();
fprintf('QPALM C \n');
fprintf('Elapsed time is %f seconds\n', res.info.run_time);
% QPALMtime = toc
% tic
% res = solver.optimize(Q, q, A, lb, ub, settings);
% QPALMtime = toc
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
fprintf('QPALM MATLAB \n');

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


