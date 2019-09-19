function lb_eig = gershgorin_min(M)
% Use Gershgorin to lower bound the eigenvalues of the matrix M.

center = diag(M);
radius = sum(abs(M), 2) - abs(center);
lb_eig = min(center-radius);
