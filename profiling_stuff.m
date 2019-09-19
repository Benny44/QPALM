load('profiling_qp4')

y = rand(m,1);
Aty = zeros(n,1);


tic;
for i=1:200
    Aty = A'*y;
end
toc


x = rand(n,1);
Qx = zeros(n,1);

tic;
for i=1:200
    Qx = Q*x;
end
toc