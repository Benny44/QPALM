close all

TIME_LIMIT = 3600;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final.mat')
Tqpalm_c_A = Tqpalm_c;
Status_qpalm_c_A = Status_qpalm_c;
X_qpalm_c_A = X_qpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final.mat')
i = [25, 42, 60]; %LINCONT, NASH and STATIC3 (PI, PI and DI)
Tqpalm_c(i) = [];
Status_qpalm_c(i) = [];
X_qpalm_c(i) = [];

Tqpalm = [Tqpalm_c_A, Tqpalm_c];
Status_qpalm_c = [Status_qpalm_c_A, Status_qpalm_c];
X_qpalm_c = [X_qpalm_c_A, X_qpalm_c];

X_qpalm_c_bar = X_qpalm_c;
Tqpalm_c_bar = Tqpalm;

[gs_qpalm, fail_rate_qpalm, Tqpalm] = compute_geometric_mean(Tqpalm, Status_qpalm_c, 'solved', TIME_LIMIT);

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final_ipopt.mat')
Tipopt_A = Tipopt;
Status_ipopt_A = Status_ipopt;
X_ipopt_A = X_ipopt;
files_A = files;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final_ipopt.mat')
i = [25, 42, 60]; %LINCONT, NASH and STATIC3 (PI, PI and DI)
Tipopt(i) = [];
Status_ipopt(i) = [];
X_ipopt(i) = [];
files(i) = [];

Tipopt = [Tipopt_A, Tipopt];
Status_ipopt = [Status_ipopt_A, Status_ipopt];
X_ipopt = [X_ipopt_A, X_ipopt];
files = [files_A, files];

% [gs_ipopt, fail_rate_ipopt, Tipopt] = compute_geometric_mean(Tipopt, Status_ipopt, 'Solve_Succeeded', TIME_LIMIT);
% 
% gs_min = min([gs_ipopt, gs_qpalm]);
% gs_qpalm = gs_qpalm/gs_min;
% gs_ipopt = gs_ipopt/gs_min;
% gs_gurobi = gs_gurobi/gs_min;

% % fprintf('gs_qpalm: %.4f, gs_osqp: %.4f, gs_gurobi: %.4f\n', gs_qpalm, gs_osqp, gs_gurobi);
% fprintf('Shifted geometric means & %.4f & %.4f\\\\\n', gs_qpalm, gs_ipopt);
% 
% % fprintf('fail_qpalm: %.4f, fail_osqp: %.4f, fail_gurobi: %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);
% fprintf('Failure rate [\\%%] & %.4f & %.4f\n', fail_rate_qpalm, fail_rate_ipopt);


%% Performance profile (timing)
% Tgurobi = Tgurobi_hat;
% Tosqp = Tosqp_hat;
% Tosqp(Tosqp==TIME_LIMIT) = inf;
% Tqpalm(Tqpalm==TIME_LIMIT) = inf;
% Tipopt(Tipopt==TIME_LIMIT) = inf;
Tqpalm(~strcmp(Status_qpalm_c, 'solved')) = inf; Tqpalm_c = Tqpalm;
Tipopt(~strcmp(Status_ipopt, 'Solve_Succeeded')) = inf;

%% Performance profile (objective)
% Tgurobi = Tgurobi_hat;
% Tosqp = Tosqp_hat;
% Tosqp(Tosqp==TIME_LIMIT) = inf;
X_qpalm_c = X_qpalm_c_bar;
% Tqpalm_c = Tqpalm_c_bar;

Obj_qpalm = [];
Obj_ipopt = [];
Tqp = [];
Tip = [];

x_tol = 1e-5;
nb_same = 0;
nb_dif = 0;
nb_both_fail = 0;
nb_qpalm_wins = 0;
nb_ipopt_wins = 0;

for i = 1:length(Tqpalm_c)
    xqp = X_qpalm_c{i};
    xip = full(X_ipopt{i});
    load(files{i});
    obj_qp = 0.5*xqp'*Data.Q*xqp + Data.q'*xqp;
    obj_ip = 0.5*xip'*Data.Q*xip + Data.q'*xip;
    feas_qp = norm([max(Data.A*xqp - Data.cu, 0)+min(Data.A*xqp - Data.cl, 0); max(xqp - Data.bu, 0)+min(xqp - Data.bl, 0)]);
    feas_ip = norm([max(Data.A*xip - Data.cu, 0)+min(Data.A*xip - Data.cl, 0); max(xip - Data.bu, 0)+min(xip - Data.bl, 0)]);
    diff_x = norm(xip-xqp)/norm(xqp);    
    if (Tqpalm_c(i) == inf && Tipopt(i) == inf)
        if (feas_qp < feas_ip && obj_qp < obj_ip)
            nb_qpalm_wins = nb_qpalm_wins + 1;
        elseif (feas_qp < feas_ip && obj_qp < obj_ip)
            nb_ipopt_wins = nb_ipopt_wins + 1;
        else
%             nb_same = nb_same + 1;
%             Tqp(nb_same) = Tqpalm_c(i);
%             Tip(nb_same) = Tipopt(i);
        end
    elseif diff_x < 1e-6
        nb_same = nb_same + 1;
        Tqp(nb_same) = Tqpalm_c(i);
        Tip(nb_same) = Tipopt(i);
    else
        nb_dif = nb_dif + 1;
        Obj_qpalm(nb_dif) = obj_qp;
        if (Tqpalm_c(i) ~= inf && Tipopt(i) == inf)
            Obj_ipopt(nb_dif) = inf;
        else
            Obj_ipopt(nb_dif) = obj_ip;
        end
    end
end

%Compare objectives for different solutions
nb_qpalm_wins = sum(Obj_qpalm < Obj_ipopt);
nb_ipopt_wins = sum(Obj_qpalm > Obj_ipopt);

%Compare timings for same solution of QPALM and Ipopt
Tmin = min(Tqp, Tip);
r_qpalm_c = Tqp./Tmin;
r_ipopt = Tip./Tmin;

rmax = max(max(r_ipopt(r_ipopt~=inf & ~isnan(r_ipopt))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xgu, ygu, fgu] = make_performance_profile(r_ipopt);
% if fgu
    xgu = [xgu rmax];
    ygu = [ygu ygu(end)];
% end


[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
% if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
% end

figure
plot(log10(xqp), yqp, 'b', log10(xgu) ,ygu, 'k')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'IPOPT','Location', 'SouthEast')

Tip(Tip==inf) = TIME_LIMIT;
Tqp(Tqp==inf) = TIME_LIMIT;
n = length(Tip);
all_success = {};
for k = 1:n
    all_success{k} = 'success';
end
[gs_ipopt, ~, ~] = compute_geometric_mean(Tip, all_success, 'success', TIME_LIMIT);

[gs_qpalm, ~, ~] = compute_geometric_mean(Tqp, all_success, 'success', TIME_LIMIT);

% gs_min = min([gs_ipopt, gs_qpalm]);
% gs_qpalm = gs_qpalm/gs_min;
% gs_ipopt = gs_ipopt/gs_min;
fprintf('Runtime (sgm) & %.4f & %.4f\\\\\n', gs_qpalm, gs_ipopt);
fprintf('Failure rate [\\%%] & %.4f & %.4f\\\\\n', fail_rate_qpalm, fail_rate_ipopt);
fprintf('Optimal & %4d & %4d \n', nb_qpalm_wins, nb_ipopt_wins);
