close all;
clear;

load('convexQPs/AUG2D')

addpath('/home/ben/Documents/Projects/QPALM/simulations')
addpath('/home/ben/Documents/Projects/QPALM/interfaces/mex')

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = false;
options.qpoases = false;
options.gurobi = false;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];

Iter_qpalm_matlab = [];
Iter_qpalm_c = [];
Iter_osqp = [];
Iter_qpoases = [];
Iter_gurobi = [];

maros_files = {};
Stats_qpalm_matlab = {};

options.SCALING_ITER=10;
options.EPS_ABS=1e-6;
options.VERBOSE = false;
options.TIME_LIMIT = 3600; 

for i = 1

    matData{i}   = Data;
    n            = matData{i}.n;
    m            = matData{i}.m;
    Q            = matData{i}.Q;
    q            = matData{i}.q;
    l            = matData{i}.bl;
    u            = matData{i}.bu;
    lb           = matData{i}.cl;
    ub           = matData{i}.cu;
    A            = matData{i}.A;
    
    prob.Q = Q; prob.q = q; prob.l = l; prob.u = u; prob.lb = lb; prob.ub = ub; prob.A = A;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
      
    [X, timings, iter, status, options, stats] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    if options.gurobi, Tgurobi(i) = gurobi_time; end

    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_matlab, Iter_out_qpalm_matlab(i) = stats.qpalm_matlab.iter_out; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
    if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
    if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    
    if options.qpalm_matlab, Stats_qpalm_matlab{i} = stats.qpalm_matlab; end
end

