function plot_QP_comparison_loglog( file )
%Helper function to plot results from simulation

load(file)
close all
figure

loglog(rc(1:length(Tqpalm_matlab)), Tqpalm_matlab, 'b',...
    rc(1:length(Tosqp)), Tosqp, 'r',...
    rc(1:length(Tqpoases)), Tqpoases, 'g',...
    rc(1:length(Tgurobi)), Tgurobi, 'k');

grid on
set(gca,'fontsize',14)
xlabel('Condition number')
ylabel('Runtime (s)')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','northeast')

figure

loglog(rc(1:length(Tqpalm_matlab)), Iter_qpalm_matlab, 'b',...
    rc(1:length(Tosqp)), Iter_osqp, 'r',...
    rc(1:length(Tqpoases)), Iter_qpoases, 'g',...
    rc(1:length(Tgurobi)), Iter_gurobi, 'k');

grid on
set(gca,'fontsize',14)
xlabel('Condition number')
ylabel('Number of iterations')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','northeast')

end

