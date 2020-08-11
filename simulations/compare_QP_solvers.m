function [ x, timings, iter, status, options, stats ] = compare_QP_solvers( prob, options )
%Run QPALM (Matlab), QPALM (C), OSQP, qpoases, and GUROBI on the given problem 
%n times and return the solution and timings

% if isfield(options, 'qrqp') && options.qrqp
%     import casadi.*
% end

n = 1; %to get an average timing
t = zeros(n,1);

x = {};
stats = {};
timings = {};
iter = {};
status = {};

if ~isfield(options, 'VERBOSE')
    VERBOSE = false;
else
    VERBOSE = options.VERBOSE;
end
if ~isfield(options, 'SCALING_ITER')
    SCALING_ITER = 10;
else
    SCALING_ITER = options.SCALING_ITER;
end
if ~isfield(options, 'MAXITER')
    MAXITER = 5000;
else
    MAXITER = options.MAXITER;
end

if ~isfield(options, 'TIME_LIMIT')
    TIME_LIMIT = 3600;
else
    TIME_LIMIT = options.TIME_LIMIT;
end

if ~isfield(options, 'EPS_ABS')
    EPS_ABS = 1e-6;
    EPS_REL = EPS_ABS;
else
    EPS_ABS = options.EPS_ABS;
    EPS_REL = EPS_ABS;
end

if ~isfield(options, 'NONCONVEX')
    NONCONVEX = false;
else
    NONCONVEX = options.NONCONVEX;
end

if ~isfield(options, 'update')
    update = false;
else
    update = options.update;
end

if ~isfield(options, 'return_solvers')
    return_solvers = false;
else
    return_solvers = options.return_solvers;
end

if ~isfield(options, 'solvers_setup')
    solvers_setup = false;
else
    solvers_setup = options.solvers_setup;
end

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

%Prep A and lbA and ubA for QPALM and OSQP
if isfield(prob, 'l') || isfield(prob, 'u')
    A = [prob.A; speye(size(prob.A,2))];
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

if isfield(options, 'qpalm_matlab') && options.qpalm_matlab     
    
    for k = 1:n
        opts.solver = 'newton';
        opts.scalar_sig = false;
        opts.maxiter = MAXITER;
        opts.eps_abs = EPS_ABS;
        opts.eps_rel = EPS_REL;
        opts.eps_abs_in = min(EPS_ABS*1e6, 1);
        opts.eps_rel_in = min(EPS_REL*1e6, 1);
        opts.eps_pinf = EPS_ABS;
        opts.eps_dinf = EPS_ABS;
        opts.proximal = true;
        opts.gamma    = 1e7;
        opts.gammaUpd = 10;
        opts.gammaMax = 1e7;
        opts.verbose = VERBOSE;
        opts.Delta   = 100;
        opts.scaling = 'simple';
        opts.scaling_iter = SCALING_ITER;
        opts.linsys = 10;
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

if isfield(options, 'qpalm_c') && options.qpalm_c
    
    for k = 1:n
        if ~solvers_setup
            solver = qpalm;
        else
            solver = options.qpalm_solver;
        end
        
        settings = solver.default_settings();
        settings.verbose = VERBOSE;
        settings.print_iter = 100;
        settings.scaling = 10;
        settings.max_iter = MAXITER;
        settings.eps_abs_in = min(EPS_ABS*1e6, 1);
        settings.eps_rel_in = min(EPS_REL*1e6, 1);
        settings.rho = 0.1;
        settings.eps_abs = EPS_ABS;
        settings.eps_rel = EPS_REL;
        settings.eps_prim_inf = EPS_ABS;
        settings.eps_dual_inf = EPS_ABS;
        settings.time_limit = TIME_LIMIT;
        settings.proximal = true;
        settings.gamma_init = 1e7;
        settings.gamma_max = 1e7;
        settings.sigma_init = 2e1;
        settings.delta = 100;
        settings.factorization_method = 2; %0: KKT, 1: SCHUR
        settings.nonconvex = NONCONVEX;
            
        if ~update
            solver.setup(prob.Q, prob.q, A,lbA,ubA, settings);
        else
%             solver.update_q(prob.q);
            solver.update_bounds(lbA, ubA);
        end
        
        if (~isempty(x_warm_start) || ~isempty(y_warm_start))
            solver.warm_start(x_warm_start, y_warm_start);
        end
        try
            res_qpalm = solver.solve();
            t(k) = res_qpalm.info.run_time;
        catch
            t(k) = TIME_LIMIT;
            res_qpalm.info.status = 'failed';
            res_qpalm.info.iter = inf;
            res_qpalm.x = nan;
        end
        
        if ~return_solvers
            solver.delete();
        elseif ~solvers_setup
            options.qpalm_solver = solver;
            options.solvers_setup = true;
        end
        
%         options.y = res_qpalm.y;
            
    end

    status.qpalm_c = res_qpalm.info.status;
    iter.qpalm_c = res_qpalm.info.iter;
    timings.qpalm_c = sum(t)/n;
    x.qpalm_c = res_qpalm.x;
    
%     if timings.qpalm_c > MAX_TIME
%         options.qpalm_c = false;
%     end
    
end
%% OSQP
% 
if isfield(options, 'osqp') && options.osqp
    for k = 1:n
        
        if ~solvers_setup
            solver = osqp;
        else
            solver = options.osqp_solver;
        end

        osqp_settings = solver.default_settings();
        osqp_settings.scaling = SCALING_ITER;
        osqp_settings.max_iter = 10000000000;
        osqp_settings.time_limit = TIME_LIMIT;
        osqp_settings.eps_abs = EPS_ABS;
        osqp_settings.eps_rel = EPS_REL;
        osqp_settings.eps_prim_inf = EPS_ABS;
        osqp_settings.eps_dual_inf = EPS_ABS;
        osqp_settings.verbose = VERBOSE;
        
        if ~update
            solver.setup(prob.Q, prob.q, A,lbA,ubA, osqp_settings);
        else
            new.l = lbA;
            new.u = ubA;
            solver.update(new);
        end
        
        if (~isempty(x_warm_start) || ~isempty(y_warm_start))
            solver.warm_start('x', x_warm_start, 'y', y_warm_start);
        end
        try
            res_osqp = solver.solve();
            t(k) = res_osqp.info.run_time;
        catch
            t(k) = TIME_LIMIT;
            res_osqp.info.status = 'failed';
            res_osqp.info.iter = inf;
            res_osqp.x = nan;
        end

        if ~return_solvers
            solver.delete();
        elseif ~solvers_setup
            options.osqp_solver = solver;
            options.solvers_setup = true;
        end 
        
    end

    status.osqp = res_osqp.info.status;
    iter.osqp = res_osqp.info.iter;
    timings.osqp = sum(t)/n;
%     if timings.osqp > MAX_TIME
%         options.osqp = false;
%     end
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

if isfield(options, 'qpoases') && options.qpoases
    qpoases_options = qpOASES_options('default', 'printLevel', 0, 'terminationTolerance', EPS_ABS, 'maxCpuTime', TIME_LIMIT, 'maxIter', MAXITER);

    for k = 1:n
        if (~isempty(x_warm_start))
            auxInput = qpOASES_auxInput('x0', x_warm_start);
            [x.qpoases,fval,status.qpoases,iter.qpoases,lambda,auxOutput] = qpOASES(prob.Q,prob.q,prob.A,l,u,prob.lb,prob.ub,qpoases_options, auxInput);
        else
            [x.qpoases,fval,status.qpoases,iter.qpoases,lambda,auxOutput] = qpOASES(prob.Q,prob.q,prob.A,l,u,prob.lb,prob.ub,qpoases_options);
        end
        t(k) = auxOutput.cpuTime;
    end
    timings.qpoases = sum(t)/n;
    
%     if timings.qpoases > MAX_TIME
%         options.qpoases = false;
%     end

end

if isfield(options, 'qrqp') && options.qrqp
    qrqp_options = struct;
    qrqp_options.print_iter = false;
    qrqp_options.record_time = true;
    qrqp_options.print_header = false;
%     qrqp_options.max_iter = 1000;
    qp_struct = struct;
    qp_struct.h = sparsity(casadi.DM(prob.Q));
    qp_struct.a = sparsity(casadi.DM(prob.A));
    solver_qrqp = casadi.conic('solver','qrqp',qp_struct,qrqp_options);
    prob_qrqp.g = prob.q;
    prob_qrqp.lbx = -inf*ones(size(prob.q));
    prob_qrqp.ubx = inf*ones(size(prob.q));
    prob_qrqp.lba = prob.lb;
    prob_qrqp.uba = prob.ub;
    prob_qrqp.h = prob.Q;
    prob_qrqp.a = prob.A;
    
    
    for k = 1:n
        try
            res_qrqp = solver_qrqp.call(prob_qrqp);
        catch
            fprintf('QRQP failed: %s\n', solver_qrqp.stats.return_status);
            res_qrqp.x = 0;
        end
        t(k) = solver_qrqp.stats.t_wall_total;
    end
    x.qrqp = full(res_qrqp.x);
    status.qrqp = solver_qrqp.stats.return_status;
    iter.qrqp = solver_qrqp.stats.n_call_total; %not available
    timings.qrqp = sum(t)/n;
end
    
if isfield(options, 'gurobi') && options.gurobi
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
    
    params.OutputFlag = 0;
    params.OptimalityTol = EPS_ABS;
    params.FeasibilityTol = EPS_ABS;
    params.TimeLimit = TIME_LIMIT;
    
    for k = 1:n
        try
            results = gurobi(model,params);
            t(k) = results.runtime;
        catch
            results.status = 'FAILED';
            results.baritercount = inf;
            t(k) = TIME_LIMIT;
        end
        
    end
    timings.gurobi = sum(t)/n;
    status.gurobi = results.status;
    iter.gurobi = results.baritercount;
    if isfield(results,'x')
        x.gurobi = results.x;
    else
        x.gurobi = nan;
    end
%     if timings.gurobi > MAX_TIME
%         options.gurobi = false;
%     end
    
end

cu = prob.ub;
cl = prob.lb;
if isfield(prob, 'u')
    bu = prob.u;
else
    bu = [];
end
if isfield(prob, 'l')
    bl = prob.l;
else
    bl = [];
end


if isfield(options, 'NewtonKKTqp') && options.NewtonKKTqp
    act_cu = cu ~= 1e20;
    act_cl = cl ~= -1e20;
    act_bu = bu ~= 1e20;
    act_bl = bl ~= -1e20;
    
    Id = speye(size(prob.Q));
    newtonA = [prob.A(act_cu,:); -prob.A(act_cl,:); Id(act_bu,:); -Id(act_bl,:)];
    newtonA = full(newtonA);
    newtonB = [cu(act_cu)+1e-10;       -cl(act_cl)-1e-10;       bu(act_bu)+1e-10;   -bl(act_bl)-1e-10]; %very bad with eq. constr
    
    opts.maxit = MAXITER;
    opts.verbose = VERBOSE;
    opts.barrier = 0;
    opts.tol_cost_fn = EPS_ABS;
    opts.tol = EPS_ABS;
    
    try
        Q_full = full(prob.Q);
        x0 = zeros(size(prob.q));
        tic; [x_final,logs] = NewtonKKTqp(Q_full,prob.q,newtonA,newtonB,x0,opts); newton_time = toc;
    
        x.NewtonKKTqp = x_final;
        stats.NewtonKKTqp = logs;
        iter.NewtonKKTqp = logs.nb_iter;
        if (newton_time > 10*logs.cpu_total)
            timings.NewtonKKTqp = newton_time;
        else
            timings.NewtonKKTqp = logs.cpu_total;
        end
        status.NewtonKKTqp = 'unknown';
    catch
        disp('NewtonKKTqp failed for some reason');
        status.NewtonKKTqp = 'failed';
        x.NewtonKKTqp = [];
        iter.NewtonKKTqp = inf;
        timings.NewtonKKTqp = TIME_LIMIT;
        stats.NewtonKKTqp = [];
    end
end

if isfield(options, 'ipopt') && options.ipopt
    ubg = prob.ub;
    lbg = prob.lb;
    if isfield(prob, 'u')
        ubx = prob.u;
    else
        ubx = [];
    end
    if isfield(prob, 'l')
        lbx = prob.l;
    else
        lbx = [];
    end
    
    import casadi.*
    n = length(prob.q);
    xi = MX.sym('x', n);
    f = 0.5*xi'*prob.Q*xi+prob.q'*xi;
    g = prob.A*xi;
    prob_struct = struct('x',xi,'f',f,'g',g); %
    % Options for ipopt
    opt_ipopt.jit = false;
    opt_ipopt.jit_options.flags = {'-O3'};
    opt_ipopt.jit_options.verbose = true;
    opt_ipopt.compiler = 'shell';
    %opt_ipopt.jit_options.compiler = 'ccache gcc';
%     opt_ipopt.ipopt.hessian_approximation = 'limited-memory';
%     opt_ipopt.ipopt.limited_memory_max_history = 30;
    opt_ipopt.ipopt.tol = EPS_ABS;
    opt_ipopt.ipopt.constr_viol_tol = EPS_ABS;
    opt_ipopt.ipopt.print_level = 0;
    opt_ipopt.ipopt.max_cpu_time = TIME_LIMIT;
    opt_ipopt.ipopt.max_iter = 1000000000;
    opt_ipopt.ipopt.hessian_constant = 'yes';
    opt_ipopt.print_time = false;

    S = nlpsol('S', 'ipopt', prob_struct, opt_ipopt);
    r = S('x0',x_warm_start,'lbx',lbx,'ubx',...
        ubx,'lbg',lbg,'ubg',ubg); %Solve  
    timings.ipopt = S.stats.t_wall_S;
    status.ipopt = S.stats.return_status;
    x.ipopt = r.x;
    iter.ipopt = S.stats.iter_count;
    
%     results.ts_casadi(k) = S.stats.t_wall_S;
%     results.ts_casadi(k) = toc(start);
%     results.ts_casadi(k) = S.stats.t_wall_mainloop;
end
