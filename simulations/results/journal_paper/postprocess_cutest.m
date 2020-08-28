close all

TIME_LIMIT = 3600;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final.mat')
Tqpalm_c_A = Tqpalm_c;
Status_qpalm_c_A = Status_qpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final.mat')
i = [25, 42, 60]; %LINCONT, NASH and STATIC3 (PI, PI and DI)
Tqpalm_c(i) = [];
Status_qpalm_c(i) = [];

Tqpalm = [Tqpalm_c_A, Tqpalm_c];
Status_qpalm_c = [Status_qpalm_c_A, Status_qpalm_c];

[gs_qpalm, fail_rate_qpalm, Tqpalm] = compute_geometric_mean(Tqpalm, Status_qpalm_c, 'solved', TIME_LIMIT);

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final_ipopt.mat')
Tipopt_A = Tipopt;
Status_ipopt_A = Status_ipopt;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final_ipopt.mat')
i = [25, 42, 60]; %LINCONT, NASH and STATIC3 (PI, PI and DI)
Tipopt(i) = [];
Status_ipopt(i) = [];

Tipopt = [Tipopt_A, Tipopt];
Status_ipopt = [Status_ipopt_A, Status_ipopt];

[gs_ipopt, fail_rate_ipopt, Tipopt] = compute_geometric_mean(Tipopt, Status_ipopt, 'Solve_Succeeded', TIME_LIMIT);

gs_min = min([gs_ipopt, gs_qpalm]);
gs_qpalm = gs_qpalm/gs_min;
gs_ipopt = gs_ipopt/gs_min;
% gs_gurobi = gs_gurobi/gs_min;

% fprintf('gs_qpalm: %.4f, gs_osqp: %.4f, gs_gurobi: %.4f\n', gs_qpalm, gs_osqp, gs_gurobi);
fprintf('Shifted geometric means & %.4f & %.4f\\\\\n', gs_qpalm, gs_ipopt);

% fprintf('fail_qpalm: %.4f, fail_osqp: %.4f, fail_gurobi: %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);
fprintf('Failure rate [\\%%] & %.4f & %.4f\n', fail_rate_qpalm, fail_rate_ipopt);


%% Performance profiles
% Tgurobi = Tgurobi_hat;
% Tosqp = Tosqp_hat;
% Tosqp(Tosqp==TIME_LIMIT) = inf;
Tqpalm(Tqpalm==TIME_LIMIT) = inf;
Tipopt(Tipopt==TIME_LIMIT) = inf;

%Compare QPALM and Gurobi
Tmin = min(Tqpalm, Tipopt);
r_qpalm_c = Tqpalm./Tmin;
r_ipopt = Tipopt./Tmin;

rmax = max(max(r_ipopt(r_ipopt~=inf & ~isnan(r_ipopt))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
end


[xgu, ygu, fgu] = make_performance_profile(r_ipopt);
if fgu
    xgu = [xgu rmax];
    ygu = [ygu ygu(end)];
end

figure
plot(log10(xqp), yqp, 'b', log10(xgu) ,ygu, 'k')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'IPOPT','Location', 'SouthEast')

% %Compare QPALM and OSQP
% Tmin = min(Tqpalm_c, Tosqp);
% r_qpalm_c = Tqpalm_c./Tmin;
% r_osqp = Tosqp./Tmin;
% 
% rmax = max(max(r_osqp(r_osqp~=inf & ~isnan(r_osqp))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 
% 
% [xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
% if fqp
%     xqp = [xqp rmax];
%     yqp = [yqp yqp(end)];
% end
% 
% 
% [xos, yos, fos] = make_performance_profile(r_osqp);
% if fos
%     xos = [xos rmax];
%     yos = [yos yos(end)];
% end
% 
% figure
% plot(log10(xqp), yqp, 'b', log10(xos) ,yos, 'r')
% set(gca,'fontsize',14)
% xlabel('log_{10}(f)')
% ylabel('fraction of solver within f of best')
% legend('QPALM', 'OSQP','Location', 'SouthEast')
