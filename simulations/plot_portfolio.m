load('output/Portfolio')
close all
figure

semilogy(n_values(1:length(Tqpalm_matlab)), Tqpalm_matlab, 'b',...
    n_values(1:length(Tqpalm_c)), Tqpalm_c, 'm',...
    n_values(1:length(Tosqp)), Tosqp, 'r',...
    n_values(1:length(Tqpoases)), Tqpoases, 'g');
grid on
legend('QPALM (Matlab)','QPALM (C)', 'OSQP', 'qpOASES','Location','southeast')