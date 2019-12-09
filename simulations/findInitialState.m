function x0 = findInitialState(A, B, T, x_upper, u_upper, F, f )
%Find a tight initial state for the given problem by solving an LP

nx = size(B,1);
nu = size(B,2);

c = zeros(nx*(T+1)+nu*T,1);
c(1:nx) = rand(nx,1);

nF = size(F,1);
Aineq = zeros(nF, nx*(T+1)+nu*T);
Aineq(:, nx*T+nu*T+1:end) = F;
bineq = f;



M = zeros(nx*T, nx*(T+1) + nu*T);
for k = 0:T-1
    M(k*nx+1:(k+1)*nx,k*(nx+nu)+1:k*(nx+nu)+nx) = A;
    M(k*nx+1:(k+1)*nx,(k+1)*nx+k*nu+1:(k+1)*(nx+nu)) = B;
    M(k*nx+1:(k+1)*nx,(k+1)*(nx+nu)+1:(k+1)*(nx+nu)+nx) = -eye(nx);
end 

Aeq = M;
beq = zeros(nx*T,1);

lb = zeros(nx*(T+1)+nu*T,1);
ub = zeros(nx*(T+1)+nu*T,1);

e = 0.000;

for k=0:T-1
   lb(k*(nx+nu)+1:k*(nx+nu)+nx) = -x_upper*(1-e);
   ub(k*(nx+nu)+1:k*(nx+nu)+nx) = x_upper*(1-e);
   lb(k*(nx+nu)+nx+1:k*(nx+nu)+nx+nu) = -u_upper*(1-e);
   ub(k*(nx+nu)+nx+1:k*(nx+nu)+nx+nu) = u_upper*(1-e);
end

[X,~,EXITFLAG] = linprog(c, Aineq, bineq, Aeq, beq, lb, ub);
fprintf('Exitflag: %d\n', EXITFLAG);

x0 = X(1:nx);


% M = zeros(nx*T+nx+nu*T+nx*T+nf, nx*(T+1)+nu*T);
%     
% %     M = zeros(nx*(T+1) + nx*(T+1) + nu*T + nxF, );
%     lb = zeros(nx*T+nx+nu*T+nx*T+nf, 1);
%     ub = zeros(nx*T+nx+nu*T+nx*T+nf, 1);
%     
%     for k = 0:T-1
%         M(k*nx+1:(k+1)*nx,k*(nx+nu)+1:k*(nx+nu)+nx) = A;
%         M(k*nx+1:(k+1)*nx,(k+1)*nx+k*nu+1:(k+1)*(nx+nu)) = B;
%         M(k*nx+1:(k+1)*nx,(k+1)*(nx+nu)+1:(k+1)*(nx+nu)+nx) = -eye(nx);
%     end
%     M(T*nx+1:T*nx+nx*(T+1)+nu*T, 1:nx*(T+1)+nu*T) = eye(nx*(T+1)+nu*T); %state and input constraints
%     M(T*nx+nx*(T+1)+nu*T+1:end, nx*T+nu*T+1:end) = F; %terminal constraints
%     
%     
%     lb(T*nx+1:(T+1)*nx) = x_init;
%     ub(T*nx+1:(T+1)*nx) = x_init;
%     
%     for k = 0:T-1
%         lb((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = -u_upper;
%         ub((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = u_upper;
%         lb((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = -x_upper;
%         ub((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = x_upper;
%     end
%     
%     %Terminal constraints
%     ub((T+1)*nx+T*(nx+nu)+1:end) = f;
%     lb((T+1)*nx+T*(nx+nu)+1:end) = -inf;



end

