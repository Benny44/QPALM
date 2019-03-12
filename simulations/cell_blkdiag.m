function X = cell_blkdiag( Xk, n, Xf )
%Helper function to create a block diagonal matrix with n times Xk and one 
%Xf on the diagonal. With only two arguments, Xf = [] is assumed
if (nargin == 2)    Xf = []; end

X_diag = cell(1,n+1);
for k = 1:n     
    X_diag{k} = Xk; 
end
X_diag{n+1} = Xf;
X = blkdiag(X_diag{:});


end

