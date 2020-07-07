%Portfolio optimization
%Sparse regressor selection
clear; close all;

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;
options.gurobi = false;
options.VERBOSE = false;
options.SCALING_ITER = 10;
options.EPS_ABS = 1e-6;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];
% rng(1)

nb_gamma = 2;
n_values = 100:20:1000;
nb_n = length(n_values);

for i = 1:nb_n
    n = n_values(i);
    m = n+1;
    
    k = ceil(n/10);
    F = sprandn(n,k,5e-1,1e-4);
    D = diag(rand(n,1)*sqrt(k));
    mu = randn(n,1);
%     mu(1:2:end) = mu(1:2:end)/1000;
    Q = blkdiag(D, eye(k));
    Q = sparse(Q);
    
%     m = m-1;
%     A = [F' -eye(k);
%         eye(n) zeros(n,k)];
%     A = sparse(A);
%     lb = [zeros(k+n,1)];
%     ub = [zeros(k,1); ones(n,1)];
    A = [ones(1,n) zeros(1,k); 
        F' -eye(k);
        eye(n) zeros(n,k)];
    A = sparse(A);
    lb = [1; zeros(k+n,1)];
    ub = [1; zeros(k,1); 1e20*ones(n,1)];

    
    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
    
    for gamma = logspace(-2,2,nb_gamma);
        q = [-1/2/gamma*mu; zeros(k,1)];
        prob.q = q;
        [ X, timings, iter, status, options, stats ] = compare_QP_solvers(prob, options);
%         [X, timings, options] = compare_QP_solvers(prob, options);
        if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
        if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
        if options.osqp, osqp_time = osqp_time + timings.osqp; end
        if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
        if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time/nb_gamma; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time/nb_gamma; end
    if options.osqp, Tosqp(i) = osqp_time/nb_gamma; end
    if options.qpoases, Tqpoases(i) = qpoases_time/nb_gamma; end
    if options.gurobi, Tgurobi(i) = gurobi_time/nb_gamma; end
    
    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
end

% save('output/Portfolio', 'n_values','Tqpalm_matlab','Tqpalm_c','Tosqp','Tqpoases','Tgurobi');

%% Plot results

% plot_QP_comparison('output/Portfolio')
    
