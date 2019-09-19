function lb_eig = gershgorin(M)
% Use Gershgorin to lower bound the eigenvalues of the matrix M.

for row = 1:size(M,1)
    center = M(row,row);
    radius = sum(abs(M(row,:))) - abs(center);
    if row==1
        lb_eig = center - radius;
    else
        lb_eig = min(center-radius, lb_eig);
    end
end
