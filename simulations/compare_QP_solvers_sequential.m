function [ x, timings, iter, status, options, osqp_solver] = compare_QP_solvers_sequential( prob, options, qpalm_solver, osqp_solver )
%Run QPALM (Matlab), QPALM (C), OSQP, qpoases, and GUROBI on the given problem 
% as one problem in a sequence of problems

% n = 2; %to get an average timing
% t = zeros(n-1,1);


%% QPALM Matlab

if isfield(options, 'x')
    x_warm_start = options.x;
else 
    x_warm_start = [];
end
if isfield(options, 'y')
    y_warm_start = options.y;
else 
    y_warm_start = [];
end

% %Prep A and lbA and ubA for QPALM and OSQP
% if isfield(prob, 'l') || isfield(prob, 'u')
%     A = [prob.A; eye(size(prob.A,2))];
%     A_combined = true;
% else
%     A = prob.A;
%     lbA = prob.lb;
%     ubA = prob.ub;
% end
% if isfield(prob, 'l') && isfield(prob, 'u')
%     lbA = [prob.lb; prob.l];
%     ubA = [prob.ub; prob.u];
% elseif isfield(prob, 'l') && ~isfield(prob, 'u')
%     lbA = [prob.lb; prob.l];
%     ubA = [prob.ub; inf*ones(length(prob.ub),1)];
% elseif ~isfield(prob, 'l') && isfield(prob, 'u')
%     lbA = [prob.lb; -inf*ones(length(prob.ub),1)];
%     ubA = [prob.ub; prob.u];
% end

if options.qpalm_matlab     
    opts.solver = 'newton';
    opts.scalar_sig = false;
    opts.maxiter = 500;
    opts.eps_abs = options.EPS_ABS;
    opts.eps_rel = options.EPS_REL;
    opts.proximal = true;
    opts.gamma    = 1e4;
    opts.gammaUpd = 10;
    opts.gammaMax = 1e6;
    opts.Delta   = 100;
%     opts.scaling = 'simple';
    if isfield(options, 'sig')
%         opts.sig = options.sig;
%         opts.gamma = options.gamma;
%         opts.LD = options.LD;
%         opts.active_cnstrs = options.active_cnstrs;
    end

    opts.scaling = '';
    opts.scaling_iter = options.SCALING_ITER; opts.scaling_iter = 1;
    tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(prob.Q,prob.q,prob.A_combined,prob.lbA_combined,prob.ubA_combined,x_warm_start,y_warm_start,opts);
    timings.qpalm_matlab = toc;
    
    options.LD = stats_qpalm.LD;
    options.active_cnstrs = stats_qpalm.active_cnstrs;
    options.sig = stats_qpalm.sig;
    options.gamma = stats_qpalm.gamma(end);
    
    status.qpalm_matlab = stats_qpalm.status;
    iter.qpalm_matlab = stats_qpalm.iter;
    x.qpalm_matlab = x_qpalm;
    
    if timings.qpalm_matlab > options.MAX_TIME
        options.qpalm_matlab = false;
    end
    
    %Set x and y for warm starting next time
    options.x = x_qpalm;
    options.y = y_qpalm;
    
end

% %% QPALM C
% 
% if options.qpalm_c && false
%     
%     for k = 1:n
%         solver = qpalm;
%         settings = solver.default_settings();
%         
%         settings.verbose = VERBOSE;
%         settings.scaling = SCALING_ITER;
%         settings.max_iter = MAXITER;
%         settings.eps_abs = EPS_ABS;
%         settings.eps_rel = EPS_REL;
%         settings.delta   = 1.2;
%         settings.memory  = 20;
%         
%         solver.setup(prob.Q, prob.q, A,lbA,ubA, settings);
%         
%         res_qpalm = solver.solve();
%         if k > 1
%             t(k-1) = res_qpalm.info.run_time;
%         end
%     end
% 
%     status.qpalm_c = res_qpalm.info.status;
%     iter.qpalm_c = res_qpalm.info.iter;
%     timings.qpalm_c = sum(t)/(n-1);
%     x.qpalm_c = res_qpalm.x;
%     
%     if timings.qpalm_c > MAX_TIME
%         options.qpalm_c = false;
%     end
%     
% end
%% QPALM C
if options.qpalm_c
    if options.i ~= 1
%         qpalm_solver.update('q', prob.q, 'l', prob.lbA_combined, 'u', prob.ubA_combined);
        qpalm_solver.update_q(prob.q);
        qpalm_solver.update_bounds(prob.lbA_combined, prob.ubA_combined);
        qpalm_solver.warm_start(x_warm_start, y_warm_start);
    end
    
    res_qpalm = qpalm_solver.solve();
    status.qpalm_c = res_qpalm.info.status;
    iter.qpalm_c = res_qpalm.info.iter;
    timings.qpalm_c = res_qpalm.info.run_time;
    x.qpalm_c = res_qpalm.x;
    
    options.x = res_qpalm.x;
    options.y = res_qpalm.y;
    
%     res_qpalm.info;
end

%% OSQP
% 
if options.osqp
    if options.i ~= 1
        osqp_solver.update('q', prob.q, 'l', prob.lbA_combined, 'u', prob.ubA_combined);
    end
    res_osqp = osqp_solver.solve(); %osqp should warm start automatically

    status.osqp = res_osqp.info.status;
    iter.osqp = res_osqp.info.iter;
    timings.osqp = res_osqp.info.run_time;
    if timings.osqp > options.MAX_TIME
        options.osqp = false;
    end
    x.osqp = res_osqp.x;
end


%% gurobi

if options.gurobi
    model.Q = 0.5*prob.Q;
    model.obj = prob.q;
    model.A = [prob.A;-prob.A];
    model.rhs = [prob.ubA;-prob.lbA];
% %     model.lb = -inf*ones(size(prob.q)); %gurobi uses default lb = 0 on the vars
% %     if isfield(prob, 'l') && isfield(prob, 'u')
%         model.lb = prob.l;
%         model.ub = prob.u;
%     elseif isfield(prob, 'l') && ~isfield(prob, 'u')
%         model.lb = prob.l;
%     elseif ~isfield(prob, 'l') && isfield(prob, 'u')
%         model.ub = prob.u;
%     end
    model.lb = prob.lb;
    model.ub = prob.ub;
    
    model.sense = '<';
    
    params.outputflag = 0;
    params.OptimalityTol = options.EPS_ABS;
    params.FeasibilityTol = options.EPS_ABS;
    
    results = gurobi(model,params);
   
    timings.gurobi = results.runtime;
    status.gurobi = results.status;
    iter.gurobi = results.baritercount;
    if isfield(results,'x')
        x.gurobi = results.x;
    else
        x.gurobi = nan;
    end
    if timings.gurobi > options.MAX_TIME
        options.gurobi = false;
    end
    

end

