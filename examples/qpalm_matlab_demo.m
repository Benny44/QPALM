m = 3000;n = 500;
A = sprandn(m, n, 1e-1, 1e-4); 

lb = -2*ones(m,1);
ub =  2*ones(m,1);
Q = sprandsym(n, 1e-1, 1e-4, 1);
q = 10*randn(n,1);


%% QPALM MATLAB
opts.Delta   = 10;
opts.eps_abs = 1e-4;
opts.eps_rel = 1e-4;
opts.maxiter = 1000;
opts.scaling = 'simple';
opts.scaling_iter = 2;

opts.solver  = 'newton';
opts.scalar_sig = false;
opts.proximal = true;
tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(Q,q,A,lb,ub,[],[],opts);qpalm_time = toc
display(stats_qpalm.status)




