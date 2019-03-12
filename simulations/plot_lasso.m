load('Lasso')
close all
figure

% hold on
semilogy(n_values(1:length(Tqpalm)), Tqpalm, 'b',...
    n_values(1:length(Tosqp)), Tosqp, 'r',...
    n_values(1:length(Tqpoases)), Tqpoases, 'g');
% semilogy(n_values(1:length(Tosqp)), Tosqp, 'r');
% semilogy(n_values(1:length(Tqpoases)), Tqpoases, 'g');
grid on
legend('QPALM', 'OSQP', 'qpOASES','Location','southeast')