function lambda_min = lobpcg( A, x )
%Calculate the minimum eigenvalue of the matrix A using the LOBPCG
%algorithm

TOL = 1e-5;

[m,n] = size(A);

assert(m==n, 'A must be square');

% rng(1);

if nargin < 2
    x = rand(n,1);
    % x = ones(n,1);
    x = x/norm(x);
end

max_iter = 10000;
Ax = A*x;

lambda_min = Ax'*x;
w = Ax - lambda_min*x;

%get residual orthonormal wrt x
w = w - x*(x'*w);
w = w/norm(w);
Aw = A*w;
xAw = Ax'*w;
wAw = w'*Aw;

[V, lambda] = eig([lambda_min xAw; xAw wAw]);
lambda_min = lambda(1,1); %The first eigenvalue is the smallest
eig_vec = V(:,1);

p = w*eig_vec(2);
Ap = Aw*eig_vec(2);
x = x*eig_vec(1)+p;
Ax = Ax*eig_vec(1)+Ap;

for i = 1:max_iter
    w = Ax-lambda_min*x;
    if norm(w, inf) < TOL 
        lambda_min = lambda_min - sqrt(2)*norm(w); %theoretical bound on lambda
        if n <= 3 %This means we have found the exact eigenvalue
            lambda_min = lambda_min - 1e-6;
        end
        break;
    end
    w = w-x*(x'*w);
    w = w/norm(w);
    Aw = A*w;
    
    xAw = Ax'*w;
    wAw = w'*Aw;
    
    norm_p = norm(p);
    p = p/norm_p;
    Ap = Ap/norm_p;
    
    xAp = Ax'*p;
    wAp = Aw'*p;
    pAp = Ap'*p;
    xp = x'*p;
    wp = w'*p;
    
    B = [lambda_min xAw xAp; xAw wAw wAp; xAp wAp pAp];
    C = [1 0 xp; 0 1 wp; xp wp 1];
 
    [eig_vec_B, lambda_B] = eig(B, C);
    lambda_min = lambda_B(1,1);
    eig_vec = eig_vec_B(:,1);
    
    p = w*eig_vec(2) + p*eig_vec(3);
    Ap = Aw*eig_vec(2) + Ap*eig_vec(3);
    x = x*eig_vec(1) + p;
    Ax = Ax*eig_vec(1) + Ap;
       
end

end
