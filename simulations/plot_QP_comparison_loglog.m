function plot_QP_comparison_loglog( file )
%Helper function to plot results from simulation

load(file)
close all
figure

loglog(condition_number(1:length(Tqpalm_matlab)), Tqpalm_matlab, 'b',...
    condition_number(1:length(Tosqp)), Tosqp, 'r',...
    condition_number(1:length(Tqpoases)), Tqpoases, 'g',...
    condition_number(1:length(Tgurobi)), Tgurobi, 'k');

grid on
set(gca,'fontsize',14)
xlabel('Condition number')
ylabel('Runtime (s)')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','northeast')

figure

loglog(condition_number(1:length(Tqpalm_matlab)), Iter_qpalm_matlab, 'b',...
    condition_number(1:length(Tosqp)), Iter_osqp, 'r',...
    condition_number(1:length(Tqpoases)), Iter_qpoases, 'g',...
    condition_number(1:length(Tgurobi)), Iter_gurobi, 'k');

grid on
set(gca,'fontsize',14)
xlabel('Condition number')
ylabel('Number of iterations')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','northeast')

end

