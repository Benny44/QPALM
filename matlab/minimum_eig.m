function [ lambda ] = minimum_eig( A )
%Calculate the minimum eigenvalue of A. This can be done in two steps:
%1) Calculate the maximum (absolute value) eigenvalue of A, lambda_max. 
%If this is negative, we are done. Else
%2) Calculate the maximum (absolute value) eigenvalue of A-lambda_max*eye,
%lambda_min. The result is then lambda = lambda_min + lambda_max

[m,n] = size(A);
assert(m==n, 'A must be square');

lambda_max = power_eig(A);
if lambda_max <= 0
    lambda = lambda_max;
else
    lambda_min = power_eig(A-lambda_max*eye(n));
    lambda = lambda_min + lambda_max; 
end

end

