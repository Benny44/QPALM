clc


load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final.mat')
% Iter_qpalm = Iter_qpalm_c;
Tqpalm = Tqpalm_c;
Status_qpalm = Status_qpalm_c;
X_qpalm = X_qpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_A_final_ipopt.mat')
i = find(Tipopt > TIME_LIMIT);
if (~isempty(i))
    for j = i
        Status_ipopt{j} = 'Maximum_CpuTime_Exceeded';
    end
end

[~, nb_files] = size(files);

for i = 1:nb_files
    load(matFiles{i});
    
    fprintf('%s & %5ld & %5ld & ', Data.name, Data.n, Data.m);
    
    if (strcmp(Status_qpalm{i}, 'solved'))
        qpalm_status = 'S';
    elseif (strcmp(Status_qpalm{i}, 'primal infeasible'))
        qpalm_status = 'PI';
    elseif (strcmp(Status_qpalm{i}, 'dual infeasible'))
        qpalm_status = 'DI';
    else
        qpalm_status = 'F';
    end
    
%     fprintf('%2s & ', qpalm_status);
    
    if (~strcmp(qpalm_status, 'S'))
        qpalm_time = qpalm_status;
        qpalm_obj = '/';
        fprintf('%s & ', qpalm_time);
    else
        qpalm_time = Tqpalm(i);
        qpalm_obj = 0.5*X_qpalm{i}'*Data.Q*X_qpalm{i} + Data.q'*X_qpalm{i};
        fprintf('%.3e & ', qpalm_time);
    end
    
    if (strcmp(Status_ipopt{i}, 'Solve_Succeeded'))
        ipopt_status = 'S';
    elseif (strcmp(Status_ipopt{i}, 'Infeasible_Problem_Detected'))
        ipopt_status = 'PI';
    elseif (strcmp(Status_ipopt{i}, 'Diverging_Iterates'))
        ipopt_status = 'DI';
    else
        ipopt_status = 'F';
    end
    
%     fprintf('%2s & ', ipopt_status);
    
    if (~strcmp(ipopt_status, 'S'))
        ipopt_time = ipopt_status;
        ipopt_obj = '/';
        fprintf('%s & ', ipopt_time);
    else
        ipopt_time = Tipopt(i);
        ipopt_obj = full(0.5*X_ipopt{i}'*Data.Q*X_ipopt{i} + Data.q'*X_ipopt{i});
        fprintf('%.3e & ', ipopt_time);
    end
    
    if (strcmp(qpalm_status, 'S'))
        if (qpalm_obj >= 0)
            fprintf('\\phantom{-}');
        end
        fprintf('%.3e & ', qpalm_obj);
    else
        fprintf('%s & ', qpalm_obj);
    end
    
    if (strcmp(ipopt_status, 'S'))
        if (ipopt_obj >= 0)
            fprintf('\\phantom{-}');
        end
        fprintf('%.3e', ipopt_obj);
    else
        fprintf('%s', ipopt_obj);
    end
    
    fprintf('\\\\ \n');     
end

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final.mat')
% Iter_qpalm = Iter_qpalm_c;
Tqpalm = Tqpalm_c;
Status_qpalm = Status_qpalm_c;
X_qpalm = X_qpalm_c;

load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/full_noA_final_ipopt.mat')
i = find(Tipopt > TIME_LIMIT);
if (~isempty(i))
    for j = i
        Status_ipopt{j} = 'Maximum_CpuTime_Exceeded';
    end
end

[~, nb_files] = size(files);

for i = 1:nb_files
    load(matFiles{i});
    
    fprintf('%s & %5ld & %5ld & ', Data.name, Data.n, Data.m);
    
    if (strcmp(Status_qpalm{i}, 'solved'))
        qpalm_status = 'S';
    elseif (strcmp(Status_qpalm{i}, 'primal infeasible'))
        qpalm_status = 'PI';
    elseif (strcmp(Status_qpalm{i}, 'dual infeasible'))
        qpalm_status = 'DI';
    else
        qpalm_status = 'F';
    end
    
%     fprintf('%2s & ', qpalm_status);
    
    if (~strcmp(qpalm_status, 'S'))
        qpalm_time = qpalm_status;
        qpalm_obj = '/';
        fprintf('%s & ', qpalm_time);
    else
        qpalm_time = Tqpalm(i);
        qpalm_obj = 0.5*X_qpalm{i}'*Data.Q*X_qpalm{i} + Data.q'*X_qpalm{i};
        fprintf('%.3e & ', qpalm_time);
    end
    
    if (strcmp(Status_ipopt{i}, 'Solve_Succeeded'))
        ipopt_status = 'S';
    elseif (strcmp(Status_ipopt{i}, 'Infeasible_Problem_Detected'))
        ipopt_status = 'PI';
    elseif (strcmp(Status_ipopt{i}, 'Diverging_Iterates'))
        ipopt_status = 'DI';
    else
        ipopt_status = 'F';
    end
    
%     fprintf('%2s & ', ipopt_status);
    
    if (~strcmp(ipopt_status, 'S'))
        ipopt_time = ipopt_status;
        ipopt_obj = '/';
        fprintf('%s & ', ipopt_time);
    else
        ipopt_time = Tipopt(i);
        ipopt_obj = full(0.5*X_ipopt{i}'*Data.Q*X_ipopt{i} + Data.q'*X_ipopt{i});
        fprintf('%.3e & ', ipopt_time);
    end
    
    if (strcmp(qpalm_status, 'S'))
        if (qpalm_obj >= 0)
            fprintf('\\phantom{-}');
        end
        fprintf('%.3e & ', qpalm_obj);
    else
        fprintf('%s & ', qpalm_obj);
    end
    
    if (strcmp(ipopt_status, 'S'))
        if (ipopt_obj >= 0)
            fprintf('\\phantom{-}');
        end
        fprintf('%.3e', ipopt_obj);
    else
        fprintf('%s', ipopt_obj);
    end
    
    fprintf('\\\\ \n');    
end