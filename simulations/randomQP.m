%Random QP
clear; close all;

current = fileparts(mfilename('fullpath'));
cd(current);

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;
options.gurobi = true;
options.qrqp = true;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tqrqp = [];
Tgurobi = [];

% n_values = 3:4;
n_values = 20:20:100;
nb_n = length(n_values);

options.EPS_ABS = 1e-6;
options.SCALING_ITER = 2;

for i = 1:nb_n
    n = n_values(i);
    m = 1*n;
    
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
    qrqp_time = 0;
    gurobi_time = 0;
      
    [X, timings, iter, status, options] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.qrqp, qrqp_time = qrqp_time + timings.qrqp; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    if options.qrqp, Tqrqp(i) = qrqp_time; end
    if options.gurobi, Tgurobi(i) = gurobi_time; end

    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.qrqp, Iter_qrqp(i) = iter.qrqp; end %iteration number not available
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
    if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.qrqp, Status_qrqp{i} = status.qrqp; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
    if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.qrqp, X_qrqp{i} = X.qrqp; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    
end

save('output/randomQP');

%% Plot results

plot_QP_comparison('output/randomQP')
    