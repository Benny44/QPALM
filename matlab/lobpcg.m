function lambda_min = lobpcg( A )
%Calculate the minimum eigenvalue of the matrix A using the LOBPCG
%algorithm

TOL = 1e-5;

[m,n] = size(A);

assert(m==n, 'A must be square');

x = rand(n,1);
x = x/norm(x);

max_iter = 1000;
for i = 1:max_iter
    w = A*x - rayleigh(A,x)*x;
    if i == 1
        V = orthonormal_basis([x, w]);
    else
        V = orthonormal_basis([x, w, x_prev]);
    end
    x_prev = x;
    
    B = V'*A*V;
    [eig_vec_B, lambda_B] = eig(B);
    [lambda_min, index] = min(diag(lambda_B));
    x = V*eig_vec_B(:,index);
    if norm(A*x - lambda_min*x, inf) < TOL
        fprintf('Lambda: %.4e, Iter: %d \n', lambda_min, i);
        break;
    end 
   
end

end

function lambda = rayleigh(A,x)
%Helper function to compute the rayleigh quotient
lambda = (x'*A*x)/(x'*x);
end

function V = orthonormal_basis(vectors)
%Assume vectors are LI
n = size(vectors,2);
V=zeros(size(vectors));
for k = 1:n
    v = vectors(:,k);
    for j = k-1:-1:1
        v = v - (V(:,j)'*v)*V(:,j);
    end
    V(:,k) = v/norm(v);
end

end

