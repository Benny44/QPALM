TIME_LIMIT = 3600;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_OSQP_BOYD.mat')
Status_osqp_boyd = Status_osqp;
Tosqp_boyd = Tosqp;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_OSQP.mat')
Tosqp = [Tosqp, Tosqp_boyd];
Status_osqp = [Status_osqp, Status_osqp_boyd];

[gs_osqp, fail_rate_osqp] = compute_geometric_mean(Tosqp, Status_osqp, 'solved', TIME_LIMIT);



load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_GUROBI_BOYD.mat')
Status_gurobi_boyd = Status_gurobi;
Tgurobi_boyd = Tgurobi;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_GUROBI.mat')
Tgurobi = [Tgurobi, Tgurobi_boyd];
Status_gurobi = [Status_gurobi, Status_gurobi_boyd];

[gs_gurobi, fail_rate_gurobi] = compute_geometric_mean(Tgurobi, Status_gurobi, 'OPTIMAL', TIME_LIMIT);

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_QPALM_BOYD.mat')

Status_qpalm_c_boyd = Status_qpalm_c;
Tqpalm_c_boyd = Tqpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_QPALM.mat')
Tqpalm_c = [Tqpalm_c, Tqpalm_c_boyd];
Status_qpalm_c = [Status_qpalm_c, Status_qpalm_c_boyd];

[gs_qpalm, fail_rate_qpalm] = compute_geometric_mean(Tqpalm_c, Status_qpalm_c, 'solved', TIME_LIMIT);

gs_min = min([gs_gurobi, gs_qpalm, gs_osqp]);
gs_qpalm = gs_qpalm/gs_min;
gs_osqp = gs_osqp/gs_min;
gs_gurobi = gs_gurobi/gs_min;

% fprintf('gs_qpalm: %.4f, gs_osqp: %.4f, gs_gurobi: %.4f\n', gs_qpalm, gs_osqp, gs_gurobi);
fprintf('Shifted geometric means & %.4f & %.4f & %.4f\\\\\n', gs_qpalm, gs_osqp, gs_gurobi);

% fprintf('fail_qpalm: %.4f, fail_osqp: %.4f, fail_gurobi: %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);
fprintf('Failure rate [\\%%] & %.4f & %.4f & %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);