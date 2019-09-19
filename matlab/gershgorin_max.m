function ub_eig = gershgorin_max(M)
% Use Gershgorin to upper bound the eigenvalues of the matrix M.

center = diag(M);
radius = sum(abs(M), 2) - abs(center);
ub_eig = max(center+radius);


