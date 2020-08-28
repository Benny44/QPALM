close all;
clear;

current_path = fileparts(mfilename('fullpath'));
cd(current_path);

myFolder = './';

if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = false;
options.qpoases = false;
options.gurobi = false;
options.NewtonKKTqp = false;
options.ipopt = false;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];
T_NewtonKKTqp = [];
T_ipopt = [];

Iter_qpalm_matlab = [];
Iter_qpalm_c = [];
Iter_osqp = [];
Iter_qpoases = [];
Iter_gurobi = [];
Iter_NewtonKKTqp = [];
Iter_ipopt = [];

files = {};
Stats_qpalm_matlab = {};

options.SCALING_ITER=10;
options.EPS_ABS=1e-6;
options.VERBOSE = true;
options.TIME_LIMIT = 3600;
options.MAXITER = 1e9;
options.NONCONVEX = true;

Data = {};
ll = 1;
load('QPNBAND');
for i = 1:ll    
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
    
    %delete empty rows of A (assuming lb <= 0 <= ub for these rows)
    b = A*rand(n,1);
    I = find(abs(b) == 0);
    A(I,:) = [];
    lb(I,:) = [];
    ub(I,:) = [];
    
    prob.Q = Q; prob.q = q; prob.l = l; prob.u = u; prob.lb = lb; prob.ub = ub; prob.A = A;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
    NewtonKKTqp_time = 0;
    ipopt_time = 0;
    
    
    [X, timings, iter, status, options, stats] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    if options.NewtonKKTqp, NewtonKKTqp_time = NewtonKKTqp_time + timings.NewtonKKTqp; end
    if options.ipopt, ipopt_time = ipopt_time + timings.ipopt; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    if options.gurobi, Tgurobi(i) = gurobi_time; end
    if options.NewtonKKTqp, T_NewtonKKTqp(i) = NewtonKKTqp_time; end
    if options.ipopt, Tipopt(i) = ipopt_time; end

    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_matlab, Iter_out_qpalm_matlab(i) = stats.qpalm_matlab.iter_out; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    if options.NewtonKKTqp, Iter_NewtonKKTqp(i) = iter.NewtonKKTqp; end
    if options.ipopt, Iter_ipopt(i) = iter.ipopt; end

    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
    if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    if options.NewtonKKTqp, Status_NewtonKKTqp{i} = status.NewtonKKTqp; end
    if options.ipopt, Status_ipopt{i} = status.ipopt; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
    if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    if options.NewtonKKTqp, X_NewtonKKTqp{i} = X.NewtonKKTqp; end
    if options.ipopt, X_ipopt{i} = X.ipopt; end

    if options.qpalm_matlab, Stats_qpalm_matlab{i} = stats.qpalm_matlab; end
    if options.NewtonKKTqp, Stats_NewtonKKTqp{i} = stats.NewtonKKTqp; end
 
end

n_values = 1:ll;

clear matData;
% save('../output/cutest');
