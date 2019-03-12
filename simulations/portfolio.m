%Portfolio optimization
%Sparse regressor selection
clear; close all;

options.qpalm = true;
options.osqp = true;
options.qpoases = true;

nb_gamma = 11;
n_values = 10:400;
nb_n = length(n_values);

Tqpalm = [];
Tosqp = [];
Tqpoases = [];

for i = 1:nb_n
    n = n_values(i);
    m = n+1;
    
    k = ceil(n/10);
    F = sprandn(n,k,5e-1);
    D = diag(rand(n,1)*sqrt(k));
    mu = randn(n,1);
    Q = blkdiag(D, eye(k));
    Q = sparse(Q);
    A = [ones(1,n) zeros(1,k); 
        F' -eye(k);
        eye(n) zeros(n,k)];
    A = sparse(A);
    lb = [1; zeros(k+n,1)];
    ub = [1; zeros(k,1); ones(n,1)];
    
    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub;
    
    qpalm_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    
    for gamma = logspace(-2,2,nb_gamma);
        q = [-1/2/gamma*mu; zeros(k,1)];
        prob.q = q;
        [X, timings, options] = compare_QP_solvers(prob, options);
        if options.qpalm, qpalm_time = qpalm_time + timings.qpalm; end
        if options.osqp, osqp_time = osqp_time + timings.osqp; end
        if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    end
    
    if options.qpalm, Tqpalm(i) = qpalm_time/nb_gamma; end
    if options.osqp, Tosqp(i) = osqp_time/nb_gamma; end
    if options.qpoases, Tqpoases(i) = qpoases_time/nb_gamma; end
    
end

save('Portfolio', 'n_values','Tqpalm','Tosqp','Tqpoases');

%% Plot results

plot_portfolio
    
