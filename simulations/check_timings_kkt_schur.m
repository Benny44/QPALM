close all;

schur_approx = true;
% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);

load('mm_ladel_kkt.mat')
T_qpalm_c_kkt = Tqpalm_c;
Iter_qpalm_c_kkt = Iter_qpalm_c;

load('mm_ladel_schur_limit_2.mat')
T_qpalm_c_schur = Tqpalm_c;
Iter_qpalm_c_schur = Iter_qpalm_c;

load('mm_cholmod.mat')
T_qpalm_c_cholmod = Tqpalm_c;

figure 
loglog([1e-3 1e3],[1 1],'b--');
hold on
for i = 1:136
    if ismember(i, [13,14,15,21])
        continue;
    else
        prob = matData{i};
        nnz_kkt = nnz(triu(prob.P)) + nnz(prob.A) + prob.m;
        if schur_approx
            Annz = zeros(prob.m,1);
            At = prob.A';
            for j = 1:prob.m
                Annz(j) = nnz(At(:,j));
            end
            nnz_schur_approx = nnz(triu(prob.P)) + 0.5*sum(Annz.*(Annz-1));
            loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'b*');
%             loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_kkt(i)/T_qpalm_c_cholmod(i), 'r*');
%             loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_schur(i)/T_qpalm_c(i), 'gsq');

        else
            nnz_schur = nnz(triu(prob.P+prob.A'*prob.A));
            loglog(nnz_kkt/nnz_schur, T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'b*');
%             loglog(nnz_kkt/nnz_schur, T_qpalm_c_kkt(i)/T_qpalm_c_cholmod(i), 'r*');
%             loglog(nnz_kkt/nnz_schur, T_qpalm_c_schur(i)/T_qpalm_c_cholmod(i), 'gsq');
        end
    end
end

title('Comparison KKT and SCHUR methods');
ylabel('T_{KKT} / T_{SCHUR}');
% title('Comparison LADEL and CHOLMOD');
% ylabel('T_{SCHUR_{LADEL}} / T_{SCHUR_{CHOLMOD}}');

if schur_approx
    xlabel('nnz_{KKT} / nnz_{SCHUR_{approx}}');
else
    xlabel('nnz_{KKT} / nnz_{SCHUR}');
end