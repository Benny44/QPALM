close all;
clear;

myFolder = './';

if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
end

[Filename, pathname] = uigetfile('*.mat',myFolder,'MultiSelect','on');
if ~iscell(Filename), Filenamec{1} = Filename;Filename = Filenamec;end
matFiles             = fullfile(pathname,Filename);
ll                   = length(matFiles);
out_maha             = cell(1, ll);
out_osqp             = cell(1, ll);
new                  = {};


for k = 1:ll
    baseFileName = Filename{k};
    new{k}       = char(baseFileName(1:end-4));
    fullFileName = fullfile(baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    matData{k}   = load(fullFileName);
    n            = matData{k}.n;
    m            = matData{k}.m;
    P            = matData{k}.P + 1e-3*speye(n);
    q            = matData{k}.q;
    lb           = matData{k}.l;
    ub           = matData{k}.u;
    A            = matData{k}.A;
    
        %% QPALM C
    % 
    solver = qpalm;
    settings = solver.default_settings();
    settings.verbose = false;
    settings.scaling = 10;
    settings.max_iter = 10000;
    settings.eps_abs = 1e-6;
    settings.eps_rel = 1e-6;
    settings.delta   = 1.2;
    settings.memory  = 20;
    solver.setup(P, q, A, lb, ub, settings);
    res_qpalm = solver.solve();
    results(k).res_qpalm = res_qpalm;
    
    %% OSQP

    solver = osqp;
    osqp_settings = solver.default_settings();
    % settings.verbose = true;
    % osqp_settings.scaling = 10;
    osqp_settings.scaling = settings.scaling;
    osqp_settings.max_iter = settings.max_iter;
    osqp_settings.eps_abs = settings.eps_abs;
    osqp_settings.eps_rel = settings.eps_rel;
    osqp_settings.verbose = settings.verbose;
    osqp_settings.check_termination = 25;
    solver.setup(P, q, A, lb, ub, osqp_settings);
    res_osqp = solver.solve();
    results(k).res_osqp = res_osqp;

    % %% QPALM MATLAB
    %Copy the settings
    opts.Delta   = settings.delta;
    opts.eps_abs = settings.eps_abs;
    opts.eps_rel = settings.eps_rel;
    opts.eps_abs_in = settings.eps_abs_in;
    opts.eps_rel_in = settings.eps_rel_in;
    opts.memory  = settings.memory;
    opts.maxiter = settings.max_iter;
    opts.rho     = settings.rho;
    opts.theta   = settings.theta;
%     opts.scaling = 'simple';
    opts.scaling_iter = settings.scaling;
    opts.solver  = 'lbfgs';
    % opts.solver = 'newton';
    opts.proximal = true;
    tic;[x_qpalm,y_qpalm,stats_qpalm] = qpalm_matlab(P,q,A,lb,ub,[],[],opts);qpalm_time = toc
    display(stats_qpalm.status)

    
    QPALMfeas = norm([min(A*x_qpalm-lb,0);min(ub-A*x_qpalm,0)],inf);    
    QPALMCfeas = norm([min(A*res_qpalm.x-lb,0);min(ub-A*res_qpalm.x,0)],inf);
    OSQPfeas = norm([min(A*res_osqp.x-lb,0);min(ub-A*res_osqp.x,0)],inf);
    QPALMCobj = 1/2*res_qpalm.x'*P*res_qpalm.x + q'*res_qpalm.x;
    QPALMobj = 1/2*x_qpalm'*P*x_qpalm + q'*x_qpalm;
    OSQPobj = 1/2*res_osqp.x'*P*res_osqp.x + q'*res_osqp.x;

    fprintf('           |   QPALM    |   OSQP     |   QPALM-N \n')
    fprintf('Iterations |   %3d      |   %3d      |    %3d   \n',...
        res_qpalm.info.iter,...
        res_osqp.info.iter,...
        stats_qpalm.iter...
        )
    fprintf('Runtime    |   %3.2e  |   %3.2e  |    %3.2e   \n',...
        res_qpalm.info.run_time,...
        res_osqp.info.run_time,...
        qpalm_time...
        )
    fprintf('Setup time |   %3.2e  |   %3.2e  |    %3.2e   \n',...
        res_qpalm.info.setup_time,...
        res_osqp.info.setup_time,...
        qpalm_time...
        )
    fprintf('Solve time |   %3.2e  |   %3.2e  |    %3.2e   \n',...
        res_qpalm.info.solve_time,...
        res_osqp.info.solve_time,...
        qpalm_time...
        )
    fprintf('Violation  | %3.2e | %3.2e |  %3.2e \n', ...
        QPALMCfeas,...
        OSQPfeas,...
        QPALMfeas...
        )
    fprintf('Objective  | %3.2e | %3.2e |  %3.2e \n', ...
        QPALMCobj, ...
        OSQPobj,...   
        QPALMCobj...
        )
    
end
