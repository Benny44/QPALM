%Sparse regressor selection
clear; close all;

current = fileparts(mfilename('fullpath'));
cd(current);

options.qpalm_matlab = true;
options.qpalm_c = false;
options.osqp = true;
options.qpoases = true;
options.gurobi = true;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];

% nb_gamma = 21;
T_values = 10:60;
nx = 10;
nu = 5;

n_values = (T_values+1)*nx+T_values*nu;

nb_n = length(n_values);

rng(1)

for i = 1:nb_n
    T = T_values(i);
    
    delta = randn(nx, 1)*0.01;
%     A = eye(nx)+diag(delta);
    A = randn(nx,nx);
    %stabilize A
    A = A/(max(abs(eig(A)))*2);
%     lambda_max = max(eig(A))-1;
%     if lambda_max > 0
%         A = A - diag(ones(nx,1)*lambda_max + 1e-6);
%     end
    B = randn(nx, nu);
%     Q = diag(sprand(nx,1,7e-1)*10);
    M = sprandn(nx, nx, 5e-1);
    Q = M*M';
%     Q = diag(rand(nx,1)*10);

    R = 0.1*eye(nu);
    QT = dare(A,B,full(Q),R);
    x_upper = rand(1,1)+1; x_upper = x_upper*ones(nx,1);
    u_upper = rand(1,1)*0.1; u_upper = u_upper*ones(nu,1);
    
    x_init = rand(nx,1).*x_upper*0.5 - x_upper;
    
    M = zeros(nx*(T+1) + nx*T + nu*T, nx*(T+1)+nu*T);
    lb = zeros(nx*(T+1) + nx*T + nu*T, 1);
    ub = zeros(nx*(T+1) + nx*T + nu*T, 1);
    
    for k = 0:T-1
        M(k*nx+1:(k+1)*nx,k*(nx+nu)+1:k*(nx+nu)+nx) = A;
        M(k*nx+1:(k+1)*nx,(k+1)*nx+k*nu+1:(k+1)*(nx+nu)) = B;
        M(k*nx+1:(k+1)*nx,(k+1)*(nx+nu)+1:(k+1)*(nx+nu)+nx) = -eye(nx);
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
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
      
    [X, timings, iter, status, options] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time; end
    if options.osqp, Tosqp(i) = osqp_time; end
    if options.qpoases, Tqpoases(i) = qpoases_time; end
    if options.gurobi, Tgurobi(i) = gurobi_time; end

    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
    if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
    if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    
end

save('output/MPC');

%% Plot results

plot_QP_comparison('output/MPC')
    