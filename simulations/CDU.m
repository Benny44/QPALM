%Sparse regressor selection
clear; close all;

current = fileparts(mfilename('fullpath'));
cd(current);

loadBenchmarkCDU

nQP = 1;

options.qpalm_matlab = false;
options.qpalm_c = false;
options.osqp = false;
options.qpoases = false;
options.gurobi = true;

options.VERBOSE = false;
options.SCALING_ITER = 10;
options.MAXITER = 10000;
options.EPS_ABS = 1e-6;
options.EPS_REL = 1e-6;
options.MAX_TIME = 100;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];

osqp_solver = osqp;
osqp_settings = osqp_solver.default_settings();

osqp_settings.scaling = options.SCALING_ITER;
osqp_settings.max_iter = options.MAXITER;
osqp_settings.eps_abs = options.EPS_ABS;
osqp_settings.eps_rel = options.EPS_REL;
osqp_settings.verbose = options.VERBOSE;

%% qpoases (can deal directly with a sequence of QPs)
if options.qpoases
    qpoases_options = qpOASES_options('default', 'printLevel', 0, 'terminationTolerance', 1e-6);
    Iter_qpoases = [];

%     [X_qpoases,fval,Status_qpoases,Iter_qpoases,lambda,auxOutput] = qpOASES(H,g',A,lb',ub',lbA',ubA',qpoases_options);
%     [X_qpoases,fval,Status_qpoases,Iter_qpoases,lambda,auxOutput] = qpOASES(H,g',A,lb',ub',lbA',ubA',qpoases_options);
    [QP,x,fval,status, iter, ~, auxOutput] = qpOASES_sequence( 'i',sparse(H),g(1,:)',sparse(A),lb(1,:)',ub(1,:)',lbA(1,:)',ubA(1,:)',qpoases_options );
    Status_qpoases{1} = status;
    X_qpoases{1} = x;
    Iter_qpoases = [Iter_qpoases iter];
    Tqpoases = [Tqpoases auxOutput.cpuTime];
    for i = 1:nQP
        [x,fval,status, iter, ~, auxOutput] = qpOASES_sequence( 'h',QP,g(i,:)',lb(i,:)',ub(i,:)',lbA(i,:)',ubA(i,:)',qpoases_options );
        Status_qpoases{i} = status;
        X_qpoases{i} = x;
        Iter_qpoases = [Iter_qpoases iter];
        Tqpoases = [Tqpoases auxOutput.cpuTime];
    end
    qpOASES_sequence( 'c',QP );    
%     Tqpoases = auxOutput.cpuTime;

end

%% other solvers
for i = 1:nQP
    options.i = i;
    
    prob.Q = sparse(H);
    prob.A = sparse(A);
    prob.A_combined = sparse([A; eye(nV)]); 
    prob.lbA_combined = [lbA(i,:)'; lb(i,:)'];
    prob.ubA_combined = [ubA(i,:)'; ub(i,:)'];
    prob.lbA = lbA(i,:)';
    prob.ubA = ubA(i,:)';
    prob.lb = lb(i,:)';
    prob.ub = ub(i,:)';
    prob.q = g(i,:)';
    
    if i==1
        osqp_solver.setup(prob.Q, prob.q, prob.A_combined,prob.lbA_combined,prob.ubA_combined, osqp_settings);
    end 
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
      
    
    [X, timings, iter, status, options] = compare_QP_solvers_sequential(prob, options, osqp_solver);
    
    %Shift solution for warm_starting
%     options.x = [options.x(33:end); zeros(32,1)];
%     options.y = [options.y(33:end); zeros(32,1)]; %equally many constraints as vars
%     options.sig = [options.sig(33:end); 2e1*ones(32,1)];
    
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
%     if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
%     if options.qpoases, Tqpoases(i) = qpoases_time; end
    if options.gurobi, Tgurobi(i) = gurobi_time; end

    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
%     if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
%     if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
%     if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    
end

osqp_solver.delete();

save('output/CDU');

%% Plot results

% Iter_qpalm_matlab
    