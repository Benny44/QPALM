Q = Data.Q;
A = [Data.A; speye(size(Data.A,2))];
At = A';
n = Data.n;
m = Data.m + n;
max_Annz = 0;

nnz_schur_approx = nnz(triu(Q,1)) + n;
for j = 1:m
    Annz = nnz(At(:,j));
    if Annz + max_Annz <= n
        nnz_schur_approx = nnz_schur_approx + 0.5*Annz*(Annz-1);
    else
        nnz_schur_approx = nnz_schur_approx + (n-max_Annz)*(Annz-(n-max_Annz+1)/2);
    end
    max_Annz = max(max_Annz, Annz);
end

nnz_kkt = nnz(triu(Q,1)) + n + nnz(At) + m;