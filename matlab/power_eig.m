function [ lambda ] = power_eig(A)
%Power method to compute largest (in absolute value) eigenvalue of A
TOL = 1e-5;

[m,n] = size(A);

assert(m==n, 'A must be square');

x = rand(n,1);
x = x/norm(x);

lambda = 0;
% lambda_prev = -inf;

iter = 0;
while norm(A*x - lambda*x, inf) > TOL
    iter = iter+1;
%     lambda_prev = lambda;
    
    x_next = A*x;
    
    x_next_norm = norm(x_next,inf);
    if (x_next_norm == 0) %A=0, or Ax and x are orthogonal. TODO, catch if orthogonal
        lambda = 0;
        return;
    end
    
    lambda = (x_next'*x)/(x'*x);
    x = x_next/x_next_norm;
end

fprintf('Power iterations: %d\n', iter);

end

