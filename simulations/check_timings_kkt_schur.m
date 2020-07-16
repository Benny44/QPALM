close all;

schur_approx = true;
% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);

% load('mm_ladel_kkt.mat')
load('full_newest_kkt.mat')
% load('kkt_reference')
T_qpalm_c_kkt = Tqpalm_c;
Iter_qpalm_c_kkt = Iter_qpalm_c;

% load('mm_ladel_schur_limit_2.mat')
load('full_newest_schur.mat');
% load('schur_reference')

T_qpalm_c_schur = Tqpalm_c;
Iter_qpalm_c_schur = Iter_qpalm_c;

% load('mm_cholmod.mat')
% T_qpalm_c_cholmod = Tqpalm_c;

figure 
loglog([1e-4 1e3],[1 1],'k--');
hold on
% loglog([2 2], [1e-5, 1e1],'g--')  

NNZ_KKT = zeros(size(T_qpalm_c_kkt));
NNZ_SCHUR = zeros(size(T_qpalm_c_kkt));
for i = 1:length(T_qpalm_c_kkt)
%     if ismember(i, [13,14,15,21])
%         continue;
%     else
        prob = matData{i};
        nnz_kkt = nnz(triu(prob.P,1)) + prob.n + nnz(prob.A) + prob.m;
        NNZ_KKT(i) = nnz_kkt;
        if schur_approx
            nnz_schur_approx = nnz(triu(prob.P,1)) + prob.n;
            
            At = prob.A';
            max_Annz = 0;
            n = matData{i}.n;
            for j = 1:prob.m
                Annz = nnz(At(:,j));
                if Annz + max_Annz <= n
                    nnz_schur_approx = nnz_schur_approx + 0.5*Annz*(Annz-1);
                else
                    nnz_schur_approx = nnz_schur_approx + (n-max_Annz)*(Annz-(n-max_Annz+1)/2);
                end
                max_Annz = max(max_Annz, Annz);
            end
            
            
            NNZ_SCHUR(i) = nnz_schur_approx;
            loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'bx');
%             loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_kkt(i)/T_qpalm_c_cholmod(i), 'r*');
%             loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_schur(i)/T_qpalm_c(i), 'gsq');

        else
            try
                nnz_schur = nnz(triu(prob.P+prob.A'*prob.A, 1)) + prob.n;
            catch
                nnz_schur = matData{i}.n^2;
            end
            loglog(nnz_kkt/nnz_schur, T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'b*');
%             loglog(nnz_kkt/nnz_schur, T_qpalm_c_kkt(i)/T_qpalm_c_cholmod(i), 'r*');
%             loglog(nnz_kkt/nnz_schur, T_qpalm_c_schur(i)/T_qpalm_c_cholmod(i), 'gsq');
        end
%     end
end

set(gca,'fontsize',12)
title('Comparison KKT and SCHUR methods');
ylabel('T_{KKT} / T_{SCHUR}');
% title('Comparison LADEL and CHOLMOD');
% ylabel('T_{SCHUR_{LADEL}} / T_{SCHUR_{CHOLMOD}}');

if schur_approx
    xlabel('nnz_{KKT} / nnz_{SCHUR_{approx}}');
else
    xlabel('nnz_{KKT} / nnz_{SCHUR}');
end