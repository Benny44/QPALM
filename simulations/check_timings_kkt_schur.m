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
loglog([1e-8 1e3],[1 1],'k--');
hold on
loglog([2 2], [1e-5, 1e1],'r--')  

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
                max_Annz = max(max_Annz, Annz);
            end
            
            for j = 1:prob.m
                Annz = nnz(At(:,j));
                if Annz + max_Annz <= n
                    nnz_schur_approx = nnz_schur_approx + 0.5*Annz*(Annz-1);
                else
                    nnz_schur_approx = nnz_schur_approx + (n-max_Annz)*(Annz-(n-max_Annz+1)/2);
                end
                max_Annz = max(max_Annz, Annz);
            end
            if (2*max_Annz > n)
                nnz_schur_approx = nnz_schur_approx + 0.5*max_Annz*(max_Annz-1) - (n-max_Annz)*(max_Annz-(n-max_Annz+1)/2);
            end
            nnz_schur_approx = min(n*(n-1)/2, nnz_schur_approx);
            
            NNZ_SCHUR(i) = nnz_schur_approx;
%             if (nnz_schur_approx > 10000) 
%                 loglog(nnz_kkt/nnz_schur_approx, T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'rx');
%             else
                loglog((nnz_kkt/nnz_schur_approx)^2*prob.n/(prob.n+prob.m), T_qpalm_c_kkt(i)/T_qpalm_c_schur(i), 'bx', 'MarkerSize', 12);
%             end
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

set(gca,'fontsize',16)
title('Comparison KKT and Schur methods');
ylabel('T_{KKT} / T_{Schur}');
% title('Comparison LADEL and CHOLMOD');
% ylabel('T_{SCHUR_{LADEL}} / T_{SCHUR_{CHOLMOD}}');

if schur_approx
    xlabel('$\displaystyle\frac{n}{n+m} \displaystyle\frac{|{\mathcal K}|^2}{|\widetilde{H}|^2} $','interpreter','latex')
%     xlabel('n/(n+m) nnz_{KKT} / nnz_{SCHUR_{approx}} ');
else
    xlabel('nnz_{KKT} / nnz_{SCHUR}');
end
set(gca,'xtick',logspace(-8,3,12))
axis tight