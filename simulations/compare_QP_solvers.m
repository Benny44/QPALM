function [ x, timings, iter, status, options, stats ] = compare_QP_solvers( prob, options )
%Run QPALM (Matlab), QPALM (C), OSQP, qpoases, and GUROBI on the given problem 
%n times and return the solution and timings

n = 1; %to get an average timing
t = zeros(n,1);

x = {};
stats = {};
timings = {};
iter = {};
status = {};


VERBOSE = false;
SCALING_ITER = 10;
MAXITER = 5000;
MAX_TIME = 10000;

if ~isfield(options, 'EPS_ABS')
    EPS_ABS = 1e-8;
    EPS_REL = EPS_ABS;
else
    EPS_ABS = options.EPS_ABS;
    EPS_REL = EPS_ABS;
end
%% QPALM Matlab

if isfield(options, 'x')
    x_warm_start = x;
else 
    x_warm_start = [];
end
if isfield(options, 'y')
    y_warm_start = y;
else 
    y_warm_start = [];
end

%Prep A and lbA and ubA for QPALM and OSQP
if isfield(prob, 'l') || isfield(prob, 'u')
    A = [prob.A; eye(size(prob.A,2))];
    A_combined = true;
else
    A = prob.A;
    lbA = prob.lb;
    ubA = prob.ub;
end
if isfield(prob, 'l') && isfield(prob, 'u')
    lbA = [prob.lb; prob.l];
    ubA = [prob.ub; prob.u];
elseif isfield(prob, 'l') && ~isfield(prob, 'u')
    lbA = [prob.lb; prob.l];
    ubA = [prob.ub; inf*ones(length(prob.ub),1)];
elseif ~isfield(prob, 'l') && isfield(prob, 'u')
    lbA = [prob.lb; -inf*ones(length(prob.ub),1)];
    ubA = [prob.ub; prob.u];
end

if options.qpalm_matlab       
    
    for k = 1:n
        opts.solver = 'newton';
        opts.scalar_sig = false;
        opts.maxiter = 50000;
        opts.eps_abs = EPS_ABS;
        opts.eps_rel = EPS_REL;
        opts.eps_abs_in = min(EPS_ABS*1e6, 1);
        opts.eps_rel_in = min(EPS_REL*1e6, 1);
        opts.eps_pinf = EPS_ABS;
        opts.eps_dinf = EPS_ABS;
        opts.proximal = true;
        opts.gamma    = 1e1;
        opts.gammaUpd = 10;
        opts.gammaMax = 1e6;
%         opts.sig = 5;
        opts.Delta   = 100;
        opts.scaling = 'simple';
        opts.scaling_iter = 10;
        tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(prob.Q,prob.q,A,lbA,ubA,x_warm_start,y_warm_start,opts);
        qpalm_time = toc;
        t(k) = qpalm_time;
        
    end
    
    status.qpalm_matlab = stats_qpalm.status;
    iter.qpalm_matlab = stats_qpalm.iter;
    timings.qpalm_matlab = sum(t)/n;
    x.qpalm_matlab = x_qpalm;
    stats.qpalm_matlab = stats_qpalm;
    
    
%     if timings.qpalm_matlab > MAX_TIME
%         options.qpalm_matlab = false;
%     end
    
end

%% QPALM C

if options.qpalm_c
    
    for k = 1:n
        solver = qpalm;
        settings = solver.default_settings();
        
        settings.verbose = VERBOSE;
        settings.scaling = 10;
        settings.max_iter = 30000;
        settings.eps_abs_in = min(EPS_ABS*1e6, 1);
        settings.eps_rel_in = min(EPS_REL*1e6, 1);
        settings.eps_abs = EPS_ABS;
        settings.eps_rel = EPS_REL;
        settings.eps_prim_inf = EPS_ABS;
        settings.eps_dual_inf = EPS_ABS;
        settings.delta   = 100;
%         settings.memory  = 20;
        
        solver.setup(prob.Q, prob.q, A,lbA,ubA, settings);
        res_qpalm = solver.solve();
        t(k) = res_qpalm.info.run_time;
    end

    status.qpalm_c = res_qpalm.info.status;
    iter.qpalm_c = res_qpalm.info.iter;
    timings.qpalm_c = sum(t)/n;
    x.qpalm_c = res_qpalm.x;
    
    if timings.qpalm_c > MAX_TIME
        options.qpalm_c = false;
    end
    
end
%% OSQP
% 
if options.osqp
    for k = 1:n
        solver = osqp;
        osqp_settings = solver.default_settings();

        osqp_settings.scaling = SCALING_ITER;
        osqp_settings.max_iter = MAXITER;
        osqp_settings.eps_abs = EPS_ABS;
        osqp_settings.eps_rel = EPS_REL;
        osqp_settings.verbose = VERBOSE;
        osqp_settings.polish = true;
        solver.setup(prob.Q, prob.q, A,lbA,ubA, osqp_settings);
        res_osqp = solver.solve();
        solver.delete();
        t(k) = res_osqp.info.run_time;
        
    end

    status.osqp = res_osqp.info.status;
    iter.osqp = res_osqp.info.iter;
    timings.osqp = sum(t)/n;
    if timings.osqp > MAX_TIME
        options.osqp = false;
    end
    x.osqp = res_osqp.x;
end
%% qpoases
if isfield(prob, 'l') && isfield(prob, 'u')
    l = prob.l;
    u = prob.u;
elseif isfield(prob, 'l') && ~isfield(prob, 'u')
    l = prob.l;
    u = [];
elseif ~isfield(prob, 'l') && isfield(prob, 'u')
    l = [];
    u = prob.u;
else
    l = [];
    u = [];
end

if options.qpoases
    qpoases_options = qpOASES_options('default', 'printLevel', 0, 'terminationTolerance', 1e-6);

    for k = 1:n
        [x.qpoases,fval,status.qpoases,iter.qpoases,lambda,auxOutput] = qpOASES(prob.Q,prob.q,prob.A,l,u,prob.lb,prob.ub,qpoases_options);
        t(k) = auxOutput.cpuTime;
    end
    timings.qpoases = sum(t)/n;
    
    if timings.qpoases > MAX_TIME
        options.qpoases = false;
    end

end

if options.gurobi
    if isfield(options, 'lp') && options.lp
       %do not define model.Q 
    else
        model.Q = 0.5*prob.Q;
    end
    model.obj = prob.q;
    model.A = [prob.A;-prob.A];
    model.rhs = [prob.ub;-prob.lb];
    model.lb = -inf*ones(size(prob.q)); %gurobi uses default lb = 0 on the vars
    if isfield(prob, 'l') && isfield(prob, 'u')
        model.lb = prob.l;
        model.ub = prob.u;
    elseif isfield(prob, 'l') && ~isfield(prob, 'u')
        model.lb = prob.l;
    elseif ~isfield(prob, 'l') && isfield(prob, 'u')
        model.ub = prob.u;
    end
    
    model.sense = '<';
    
    params.outputflag = 0;
    params.OptimalityTol = EPS_ABS;
    params.FeasibilityTol = EPS_ABS;
    
    for k = 1:n
        results = gurobi(model,params);
        t(k) = results.runtime;
    end
    timings.gurobi = sum(t)/n;
    status.gurobi = results.status;
    iter.gurobi = results.baritercount;
    if isfield(results,'x')
        x.gurobi = results.x;
    else
        x.gurobi = nan;
    end
    if timings.gurobi > MAX_TIME
        options.gurobi = false;
    end
    

end

