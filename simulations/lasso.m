%Sparse regressor selection
clear; close all;

options.qpalm = true;
options.osqp = true;
options.qpoases = true;

nb_gamma = 21;
n_values = 10:200;
nb_n = length(n_values);


for i = 1:nb_n
    n = n_values(i);
    m = 2*n;
    C = sprandn(m,n,4e-4);
    x_hat = sprandn(n,1,5e-1)/n;
    eps = randn(m,1)/4;
    d = C*x_hat+eps;
    
    Q = blkdiag(zeros(n,n), eye(m), zeros(n,n));
    Q = sparse(Q);
    A = [C -eye(m) zeros(m,n); 
        eye(n), zeros(n,m), eye(n);
        -eye(n), zeros(n,m), eye(n)];
    A = sparse(A);
    lb = [d; zeros(n+n,1)];
%     lb = sparse(lb);
    ub = [d; inf*ones(n+n,1)];
%     ub = sparse(ub);
    
    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub;
    
    qpalm_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    
    for gamma = logspace(-2,2,nb_gamma);
        q = [zeros(m+n,1); gamma*ones(n,1)];
        prob.q = q;
        [X, timings, options] = compare_QP_solvers(prob, options);
        qpalm_time = qpalm_time + timings.qpalm;
        osqp_time = osqp_time + timings.osqp;
        qpoases_time = qpoases_time + timings.qpoases;
    end
    
    if options.qpalm, Tqpalm(i) = qpalm_time/nb_gamma; end
    if options.osqp, Tosqp(i) = osqp_time/nb_gamma; end
    if options.qpoases, Tqpoases(i) = qpoases_time/nb_gamma; end
    
end

save('Lasso', 'n_values','Tqpalm','Tosqp','Tqpoases');

%% Plot results

plot_lasso
    