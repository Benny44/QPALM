%Sparse regressor selection
clear; close all;

options.qpalm = true;
options.osqp = true;
options.qpoases = true;

% nb_gamma = 21;
T_values = 10:20;
nx = 10;
nu = 5;

n_values = (T_values+1)*nx+T_values*nu;

nb_n = length(n_values);


for i = 1:nb_n
    T = T_values(i);
    
    delta = randn(nx, 1)*0.01;
    A = eye(nx)+diag(delta);
    %stabilize A
    lambda_max = max(eig(A))-1;
    if lambda_max > 0
        A = A - diag(ones(nx,1)*lambda_max + 1e-6);
    end
    B = randn(nx, nu);
    Q = diag(sprand(nx,1,7e-1)*10);
    R = 0.1*eye(nu);
    QT = dare(A,B,full(Q),R);
    x_upper = rand(nx,1)+1;
    u_upper = rand(nu,1)*0.1;
    
    x_init = rand(nx,1).*x_upper*0.5 - x_upper;
    
    M = zeros(nx*(T+1) + nx*T + nu*T, nx*(T+1)+nu*T);
    lb = zeros(nx*(T+1) + nx*T + nu*T, 1);
    ub = zeros(nx*(T+1) + nx*T + nu*T, 1);
    
    for k = 0:T-1
        M(k*nx+1:(k+1)*nx,k*nx+1:(k+1)*nx) = A;
        M(k*nx+1:(k+1)*nx,(k+1)*nx+1:(k+1)*nx+nu) = B;
        M(k*nx+1:(k+1)*nx,(k+1)*nx+nu+1:(k+1)*nx+nu+nx) = -eye(nx);
    end
    M(T*nx+1:end, :) = eye(nx*(T+1)+nu*T);
    
    lb(T*nx+1:(T+1)*nx) = x_init;
    ub(T*nx+1:(T+1)*nx) = x_init;
    
    for k = 0:T-1
        lb((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = -u_upper;
        ub((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = u_upper;
        lb((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = -x_upper;
        ub((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = x_upper;
    end
           
    q = zeros(nx*(T+1)+nu*T,1);
    Q = cell_blkdiag(blkdiag(Q, R), T, QT);
    
    Q = sparse(Q);
    M = sparse(M);
    prob.Q = Q; prob.A = M; prob.lb = lb; prob.ub = ub; prob.q = q;
    
    qpalm_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
      
    [X, timings, options] = compare_QP_solvers(prob, options);
    if options.qpalm, qpalm_time = timings.qpalm; end
    if options.osqp, osqp_time = timings.osqp; end
    if options.qpoases, qpoases_time = timings.qpoases; end
    
    if options.qpalm, Tqpalm(i) = qpalm_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    
end

save('MPC', 'n_values','Tqpalm','Tosqp','Tqpoases');

%% Plot results

plot_mpc
    