close all

TIME_LIMIT = 3600;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_OSQP_BOYD.mat')
Status_osqp_boyd = Status_osqp;
Tosqp_boyd = Tosqp;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_OSQP.mat')
Tosqp = [Tosqp, Tosqp_boyd];
Status_osqp = [Status_osqp, Status_osqp_boyd];

[gs_osqp, fail_rate_osqp, Tosqp] = compute_geometric_mean(Tosqp, Status_osqp, 'solved', TIME_LIMIT);
Tosqp_hat = Tosqp;


load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_GUROBI_BOYD.mat')
Status_gurobi_boyd = Status_gurobi;
Tgurobi_boyd = Tgurobi;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_GUROBI.mat')
Tgurobi = [Tgurobi, Tgurobi_boyd]; 
Status_gurobi = [Status_gurobi, Status_gurobi_boyd];

[gs_gurobi, fail_rate_gurobi, Tgurobi] = compute_geometric_mean(Tgurobi, Status_gurobi, 'OPTIMAL', TIME_LIMIT);
Tgurobi_hat = Tgurobi;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_QPALM_BOYD.mat')

Status_qpalm_c_boyd = Status_qpalm_c;
Tqpalm_c_boyd = Tqpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/MM_QPALM.mat')
Tqpalm_c = [Tqpalm_c, Tqpalm_c_boyd];
Status_qpalm_c = [Status_qpalm_c, Status_qpalm_c_boyd];

[gs_qpalm, fail_rate_qpalm, Tqpalm_c] = compute_geometric_mean(Tqpalm_c, Status_qpalm_c, 'solved', TIME_LIMIT);

gs_min = min([gs_gurobi, gs_qpalm, gs_osqp]);
gs_qpalm = gs_qpalm/gs_min;
gs_osqp = gs_osqp/gs_min;
gs_gurobi = gs_gurobi/gs_min;

% fprintf('gs_qpalm: %.4f, gs_osqp: %.4f, gs_gurobi: %.4f\n', gs_qpalm, gs_osqp, gs_gurobi);
fprintf('Shifted geometric means & %.4f & %.4f & %.4f\\\\\n', gs_qpalm, gs_osqp, gs_gurobi);

% fprintf('fail_qpalm: %.4f, fail_osqp: %.4f, fail_gurobi: %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);
fprintf('Failure rate [\\%%] & %.4f & %.4f & %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);


%% Performance profiles
Tgurobi = Tgurobi_hat;
Tosqp = Tosqp_hat;
Tosqp(Tosqp==TIME_LIMIT) = inf;
Tqpalm_c(Tqpalm_c==TIME_LIMIT) = inf;
Tgurobi(Tgurobi==TIME_LIMIT) = inf;

%Compare QPALM and Gurobi
Tmin = min(Tqpalm_c, Tgurobi);
r_qpalm_c = Tqpalm_c./Tmin;
r_gurobi = Tgurobi./Tmin;

rmax = max(max(r_gurobi(r_gurobi~=inf & ~isnan(r_gurobi))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
end


[xgu, ygu, fgu] = make_performance_profile(r_gurobi);
if fgu
    xgu = [xgu rmax];
    ygu = [ygu ygu(end)];
end

figure
plot(log10(xqp), yqp, 'b', log10(xgu) ,ygu, 'k')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'Gurobi','Location', 'SouthEast')

%Compare QPALM and OSQP
Tmin = min(Tqpalm_c, Tosqp);
r_qpalm_c = Tqpalm_c./Tmin;
r_osqp = Tosqp./Tmin;

rmax = max(max(r_osqp(r_osqp~=inf & ~isnan(r_osqp))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
end


[xos, yos, fos] = make_performance_profile(r_osqp);
if fos
    xos = [xos rmax];
    yos = [yos yos(end)];
end

figure
plot(log10(xqp), yqp, 'b', log10(xos) ,yos, 'r')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'OSQP','Location', 'SouthEast')
