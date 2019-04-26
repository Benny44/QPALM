function [ lambda ] = power_eig(A)
%Power method to compute largest (in absolute value) eigenvalue of A
TOL = 0.005;

[m,n] = size(A);

assert(m==n, 'A must be square');

x = ones(n,1);

lambda = 0;
lambda_prev = -inf;

while abs(lambda-lambda_prev) > TOL
    lambda_prev = lambda;
    
    x_next = A*x;
    
    x_next_norm = norm(x_next,inf);
    if (x_next_norm == 0) %A=0, or Ax and x are orthogonal. TODO, catch if orthogonal
        lambda = 0;
        return;
    end
    
    lambda = (x_next'*x)/(x'*x);
    x = x_next/x_next_norm;
end

end

