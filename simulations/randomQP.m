%Sparse regressor selection
clear; close all;

options.qpalm_matlab = true;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;

n_values = 200:220;
nb_n = length(n_values);

for i = 1:nb_n
    n = n_values(i);
    m = 6*n;
    
    M = sprandn(n, n, 5e-1);
    Q = M*M';
    
    A = sprandn(m,n,5e-1);
    q = randn(n,1);
    lb = -rand(m,1);
    ub = rand(m,1);
    
    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub; prob.q = q;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
      
    [X, timings, options] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    
end

save('output/randomQP', 'n_values','Tqpalm_matlab','Tqpalm_c','Tosqp','Tqpoases');

%% Plot results

plot_QP_comparison('output/randomQP')
    