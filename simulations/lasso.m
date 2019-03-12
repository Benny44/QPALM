%Sparse regressor selection
clear; close all;

options.qpalm_matlab = true;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;

nb_gamma = 21;
n_values = 50:51;
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
    ub = [d; inf*ones(n+n,1)];

    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    
    for gamma = logspace(-2,2,nb_gamma);
        q = [zeros(m+n,1); gamma*ones(n,1)];
        prob.q = q;
        [X, timings, options] = compare_QP_solvers(prob, options);
        if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
        if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
        if options.osqp, osqp_time = osqp_time + timings.osqp; end
        if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time/nb_gamma; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time/nb_gamma; end
    if options.osqp, Tosqp(i) = osqp_time/nb_gamma; end
    if options.qpoases, Tqpoases(i) = qpoases_time/nb_gamma; end
    
end

save('output/Lasso', 'n_values','Tqpalm_matlab','Tqpalm_c','Tosqp','Tqpoases');

%% Plot results

plot_QP_comparison('output/Lasso')
    