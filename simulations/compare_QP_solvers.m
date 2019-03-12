function [ x, timings, options ] = compare_QP_solvers( prob, options )
%Run QPALM, OSQP, qpoases, and GUROBI on the given problem n times and return the
%solution and timings

n = 2; %to get an average timing
t = zeros(n-1,1);
%% QPALM C
%
if options.qpalm
    for k = 1:n
        solver = qpalm;
        settings = solver.default_settings();
        settings.verbose = false;
        settings.scaling = 10;
        settings.max_iter = 10000;
        settings.eps_abs = 1e-4;
        settings.eps_rel = 1e-4;
        settings.delta   = 1.2;
        settings.memory  = 20;
        solver.setup(prob.Q, prob.q, prob.A, prob.lb, prob.ub, settings);
        res_qpalm = solver.solve();
        if k > 1
            t(k-1) = res_qpalm.info.run_time;
        end
    end

    timings.qpalm = sum(t)/(n-1);
    x.qpalm = res_qpalm.x;
    
    if timings.qpalm > 10
        options.qpalm = false;
    end
    
end
%% OSQP
% 
if options.osqp
    for k = 1:n
        solver = osqp;
        osqp_settings = solver.default_settings();
        % settings.verbose = true;
        % osqp_settings.scaling = 10;
        osqp_settings.scaling = settings.scaling;
        osqp_settings.max_iter = settings.max_iter;
        osqp_settings.eps_abs = settings.eps_abs;
        osqp_settings.eps_rel = settings.eps_rel;
        osqp_settings.verbose = settings.verbose; %osqp_settings.verbose = false;
    %     osqp_settings.adaptive_rho = 1;
    %     osqp_settings.adaptive_rho_tolerance = 5;
        % osqp_settings.linsys_solver ='mkl pardiso';
    %     osqp_settings.check_termination = 25;
        solver.setup(prob.Q, prob.q, prob.A, prob.lb, prob.ub, osqp_settings);
        res_osqp = solver.solve();
        solver.delete();
        if k > 1
            t(k-1) = res_osqp.info.run_time;
        end
    end

    timings.osqp = sum(t)/(n-1);
    if timings.osqp > 10
        options.osqp = false;
    end
    x.osqp = res_osqp.x;
end
%% qpoases
if options.qpoases
    qpoases_options = qpOASES_options('default', 'printLevel', 0);

    for k = 1:n
        [x.qpoases,fval,exitflag,iter,lambda,auxOutput] = qpOASES(prob.Q,prob.q,prob.A,[],[],prob.lb,prob.ub,qpoases_options);
        if k > 1
            t(k-1) = auxOutput.cpuTime;
        end
    end
    timings.qpoases = sum(t)/(n-1);
    
    if timings.qpoases > 10
        options.qpoases = false;
    end

end

end

