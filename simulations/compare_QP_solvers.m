function [ x, timings, iter, status, options ] = compare_QP_solvers( prob, options )
%Run QPALM (Matlab), QPALM (C), OSQP, qpoases, and GUROBI on the given problem 
%n times and return the solution and timings

n = 2; %to get an average timing
t = zeros(n-1,1);

VERBOSE = false;
SCALING_ITER = 10;
MAXITER = 10000;
EPS_ABS = 1e-6;
EPS_REL = 1e-6;
MAX_TIME = 10000;
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

A_combined = false; %if A is already combined in qpalm_matlab, don't do it again for OSQP

if options.qpalm_matlab

    
    for k = 1:n
        opts.solver = 'newton';
        opts.scalar_sig = false;
        opts.maxiter = 500;
        opts.eps_abs = EPS_ABS;
        opts.eps_rel = EPS_REL;
        opts.proximal = true;
        opts.gamma    = 1e4;
        opts.gammaUpd = 10;
        opts.gammaMax = 1e6;
        opts.Delta   = 10;
        opts.scaling = 'simple';
        opts.scaling_iter = SCALING_ITER; opts.scaling_iter = 1;
        tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(prob.Q,prob.q,prob.A,prob.lb,prob.ub,[],[],opts);
        qpalm_time = toc;
        if k > 1
            t(k-1) = qpalm_time;
        end
    end
    status.qpalm_matlab = stats_qpalm.status;
    iter.qpalm_matlab = stats_qpalm.iter;
    timings.qpalm_matlab = sum(t)/(n-1);
    x.qpalm_matlab = x_qpalm;
    
    if timings.qpalm_matlab > MAX_TIME
        options.qpalm_matlab = false;
    end
    
end

%% QPALM C

if options.qpalm_c
    for k = 1:n
        solver = qpalm;
        settings = solver.default_settings();
        
        settings.verbose = VERBOSE;
        settings.scaling = SCALING_ITER;
        settings.max_iter = MAXITER;
        settings.eps_abs = EPS_ABS;
        settings.eps_rel = EPS_REL;
        settings.delta   = 1.2;
        settings.memory  = 20;
        
        solver.setup(prob.Q, prob.q, prob.A, prob.lb, prob.ub, settings);
        res_qpalm = solver.solve();
        if k > 1
            t(k-1) = res_qpalm.info.run_time;
        end
    end

    status.qpalm_c = res_qpalm.info.status;
    iter.qpalm_c = res_qpalm.info.iter;
    timings.qpalm_c = sum(t)/(n-1);
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
        solver.setup(prob.Q, prob.q, prob.A, prob.lb, prob.ub, osqp_settings);
        res_osqp = solver.solve();
        solver.delete();
        if k > 1
            t(k-1) = res_osqp.info.run_time;
        end
    end

    status.osqp = res_osqp.info.status;
    iter.osqp = res_osqp.info.iter;
    timings.osqp = sum(t)/(n-1);
    if timings.osqp > MAX_TIME
        options.osqp = false;
    end
    x.osqp = res_osqp.x;
end
%% qpoases
if options.qpoases
    qpoases_options = qpOASES_options('default', 'printLevel', 0);

    for k = 1:n
        [x.qpoases,fval,status.qpoases,iter.qpoases,lambda,auxOutput] = qpOASES(prob.Q,prob.q,prob.A,[],[],prob.lb,prob.ub,qpoases_options);
        if k > 1
            t(k-1) = auxOutput.cpuTime;
        end
    end
    timings.qpoases = sum(t)/(n-1);
    
    if timings.qpoases > MAX_TIME
        options.qpoases = false;
    end

end

if options.gurobi
    model.Q = 0.5*prob.Q;
    model.obj = prob.q;
    model.A = [prob.A;-prob.A];
    model.rhs = [prob.ub;-prob.lb];
    model.lb = -inf*ones(size(prob.q)); %gurobi uses default lb = 0 on the vars
    model.sense = '<';
    
    params.outputflag = 0;
    params.OptimalityTol = EPS_ABS;
    params.FeasibilityTol = EPS_ABS;
    
    for k = 1:n
        results = gurobi(model,params);
        if k > 1
            t(k-1) = results.runtime;
        end
    end
    timings.gurobi = sum(t)/(n-1);
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

