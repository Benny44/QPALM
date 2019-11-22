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

[Filename, pathname] = uigetfile('*.mat',myFolder,'MultiSelect','on');
if ~iscell(Filename), Filenamec{1} = Filename;Filename = Filenamec;end
matFiles             = fullfile(pathname,Filename);
ll                   = length(matFiles);
out_maha             = cell(1, ll);
out_osqp             = cell(1, ll);
new                  = {};

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

options.SCALING_ITER=2;
options.EPS_ABS=1e-6;
options.VERBOSE=true;

for i = 1:ll
    baseFileName = Filename{i};
    new{i}       = char(baseFileName(1:end-4));
    fullFileName = fullfile(baseFileName);
    maros_files{i} = fullFileName;
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    matData{i}   = load(fullFileName);
    n            = matData{i}.n;
    m            = matData{i}.m;
    P            = matData{i}.P;% + 1e-3*speye(n);
    q            = matData{i}.q;
    lb           = matData{i}.l;
    ub           = matData{i}.u;
    A            = matData{i}.A;
    
    ind = ~(lb==-1e20 & ub==1e20);
    m = sum(ind);
    lb = lb(ind);
    ub = ub(ind);
    A = A(ind, :);
    
    prob.Q = P; prob.q = q; prob.lb = lb; prob.ub = ub; prob.A = A;
    
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

n_values = 1:ll;

save('../output/maros_meszaros');

