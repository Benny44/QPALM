%Random MPC (for one problem, solve the OCPs sequentially)
clear; close all;

current = fileparts(mfilename('fullpath'));
cd(current);

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;
options.gurobi = true;
options.VERBOSE = false;
options.MAXITER = 10000;
options.TIME_LIMIT = 100;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];

% nb_gamma = 21;
% rng(1)
T = 30;
nx = 10;
nu = 5;

n_values = 1:T;
nb_n = T;

% rng(1)
A = 2*sprandsym(nx,1,1e-2);
B = 1*randn(nx, nu);
M = 5*sprandn(nx, nx, 5e-1, 1e-4);
Q = M*M';

R = 0.01*eye(nu);
QT = dare(full(A),B,full(Q),R);
K = dlqr(full(A),B,full(Q),R);
K = -K; % Note that matlab defines K as -K
Ak = A+B*K; % Closed-loop dynamics
x_upper = 2*rand(1,1)+10; x_upper = x_upper*ones(nx,1);
u_upper = 2*rand(1,1)+10; u_upper = u_upper*ones(nu,1);

H = [eye(nx); -eye(nx)];
h = [x_upper; x_upper];
X = polytope(H,h);
Hu = [eye(nu); -eye(nu)];
hu = [u_upper; u_upper];
U = polytope(Hu,hu);

HH = [H;Hu*K]; hh = [h;hu]; % State and input constraints

% Compute the maximal invariant set
O = polytope(HH,hh);
while 1
    Oprev = O;
    [F,f] = double(O);  
    % Compute the pre-set
    O = polytope([F;F*Ak],[f;f]);
    if O == Oprev, break; end
end
%     termSet = O;
[F,f] = double(O);
[nf, nx] = size(F); 

%     x_init = rand(nx,1).*x_upper*0.5 - x_upper;
x_init = findInitialState(A, B, T, x_upper, u_upper, F, f);

M = zeros(nx*T+nx+nu*T+nx*T+nf, nx*(T+1)+nu*T);

%     M = zeros(nx*(T+1) + nx*(T+1) + nu*T + nxF, );
lb = zeros(nx*T+nx+nu*T+nx*T+nf, 1);
ub = zeros(nx*T+nx+nu*T+nx*T+nf, 1);

for k = 0:T-1
    M(k*nx+1:(k+1)*nx,k*(nx+nu)+1:k*(nx+nu)+nx) = A;
    M(k*nx+1:(k+1)*nx,(k+1)*nx+k*nu+1:(k+1)*(nx+nu)) = B;
    M(k*nx+1:(k+1)*nx,(k+1)*(nx+nu)+1:(k+1)*(nx+nu)+nx) = -eye(nx);
end
M(T*nx+1:T*nx+nx*(T+1)+nu*T, 1:nx*(T+1)+nu*T) = eye(nx*(T+1)+nu*T); %state and input constraints
M(T*nx+nx*(T+1)+nu*T+1:end, nx*T+nu*T+1:end) = F; %terminal constraints


lb(T*nx+1:(T+1)*nx) = x_init;
ub(T*nx+1:(T+1)*nx) = x_init;

for k = 0:T-1
    lb((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = -u_upper;
    ub((T+1)*nx+k*(nx+nu)+1: (T+1)*nx+k*(nx+nu)+nu) = u_upper;
    lb((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = -x_upper;
    ub((T+1)*nx+k*(nx+nu)+nu+1: (T+1)*nx+(k+1)*(nx+nu)) = x_upper;
end

%Terminal constraints
ub((T+1)*nx+T*(nx+nu)+1:end) = f;
lb((T+1)*nx+T*(nx+nu)+1:end) = -inf;

q = zeros(nx*(T+1)+nu*T,1);
Q = cell_blkdiag(blkdiag(Q, R), T, QT);

Q = sparse(Q);
M = sparse(M);
prob.Q = Q; prob.A = M; prob.lb = lb; prob.ub = ub; prob.q = q;

% Initialize with the autonomous trajectory
x_warm = zeros(nx*(T+1)+nu*T,1);
x_warm(1:nx) = x_init;
x_k = x_init;
for k = 1:T
    x_k = A*x_k;
    x_warm(k*(nx+nu)+1:k*nu+(k+1)*nx) = x_k;
end
options.x = x_warm;

options.return_solvers = true;

for i = 1:T
       
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
    
    %% Apply the first input (of the qpalm solution), shift the vectors for the next OCP and fix the bounds
    if strcmp(status.qpalm_c, 'solved')
        
        %Apply the first input (and add a small disturbance)
        w = normrnd(0, 1e-2, nx, 1);
        x_warm = X.qpalm_c;
        u = x_warm(nx+1:nx+nu);
        x_init = A*x_init + B*u + w;
        
        % Shift the solution to warm start. Add zero input and autonomous
        % state evolution for the last input/state
        x_warm = [x_warm(nx+nu+1:end); zeros(nu,1); A*x_warm(end-nx+1:end)];
        
        % Update bounds for initial state
        prob.lb(T*nx+1:(T+1)*nx) = x_init;
        prob.ub(T*nx+1:(T+1)*nx) = x_init;
        
        options.x = x_warm;
        options.update = true;
        
    else
        error('Failed to solve for some reason');
    end
    
end

if options.osqp, Tosqp(strcmp(Status_osqp, 'run time limit reached')) = options.TIME_LIMIT; end

save('output/MPC_sequential');

%% Plot results

plot_QP_comparison('output/MPC_sequential')
    