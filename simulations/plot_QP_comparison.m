function plot_QP_comparison( file )
%Helper function to plot results from simulation

load(file)
close all
figure

semilogy(n_values(1:length(Tqpalm_matlab)), Tqpalm_matlab, 'm',...
    n_values(1:length(Tqpalm_c)), Tqpalm_c, 'b',...
    n_values(1:length(Tosqp)), Tosqp, 'r',...
    n_values(1:length(Tqpoases)), Tqpoases, 'g',...
    n_values(1:length(Tgurobi)), Tgurobi, 'k')%,...
%     n_values(1:length(Tqrqp)), Tqrqp, 'c');
%      
grid on
set(gca,'fontsize',14)
xlabel('Number of primal variables')
ylabel('Runtime (s)')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','southeast')

end

