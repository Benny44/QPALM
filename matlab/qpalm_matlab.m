function [x,yh,stats] = qpalm_matlab(Q,q,A,bmin,bmax,x,y,opts)
[m,n] = size(A);
if nargin<7 || isempty(y)
    y = zeros(m,1);
end

if nargin<6 || isempty(x)
    x = zeros(n,1);
end

if nargin<8 || ~isfield(opts,'verbose')
    verbose = false; 
else
    verbose = opts.verbose;
end

if nargin<8 || ~isfield(opts,'scaling')
    scaling = ''; %'ruiz', 'simple', 'outer_iter'
else
    scaling = opts.scaling;
end

if nargin<8 || ~isfield(opts,'scaling_iter')
    scaling_iter = 2; %'ruiz'
else
    scaling_iter = opts.scaling_iter;
end

D_scale = 1;
E_scale = 1;
c_scale = 1;
if strcmp(scaling, 'ruiz')
    [x,y,Q,q,A,bmin,bmax,D_scale,E_scale,c_scale] = modified_ruiz_equilibration(x,y,Q,q,A,bmin,bmax,scaling_iter);
elseif strcmp(scaling, 'simple') %unit rows A, and unit diagonal elements in Q
    [x,y,Q,q,A,bmin,bmax,D_scale,E_scale,c_scale] = simple_equilibration(x,y,Q,q,A,bmin,bmax,scaling_iter);
end 
[m,n] = size(A);

x_prev = x; x0 = x;
Ax = A*x;
Qx = Q*x;
Aty = A'*y;

if nargin<8 || ~isfield(opts,'Delta')
    Delta = 10;
else
    Delta = opts.Delta;
end

if nargin<8 || ~isfield(opts,'sig')
    f = 0.5*(x'*Qx)+q'*x;
    dist = Ax-min(max(Ax,bmin),bmax);
%     nrm_rp_unscaled = norm(dist,inf);
    dist2 = dist'*dist;
    sig = max(1e-4,min(2e1*max(1,abs(f))/max(1,0.5*dist2),1e4))*ones(m,1);
%     sig = max(1, Delta*abs(dist)/(nrm_rp_unscaled+1e-6)).*sig;
%     sig = 5*ones(m,1);
else
    sig = opts.sig.*ones(m,1);
end

if nargin<8 || ~isfield(opts,'sigma_max')
    sigma_max = 1e9;
else
    sigma_max = opts.sigma_max;
end

if nargin<8 || ~isfield(opts,'eps_abs')
    eps_abs = 1e-4;
else
    eps_abs = opts.eps_abs;
end

if nargin<8 || ~isfield(opts,'eps_rel')
    eps_rel = 1e-4;
else
    eps_rel = opts.eps_rel;
end

if nargin<8 || ~isfield(opts,'eps_pinf')
    eps_pinf = 1e-5;
else
    eps_pinf = opts.eps_pinf;
end

if nargin<8 || ~isfield(opts,'eps_dinf')
    eps_dinf = 1e-5;
else
    eps_dinf = opts.eps_dinf;
end

if nargin<8 || ~isfield(opts,'rho')
    rho = 1e-1;
else
    rho = opts.rho;
end

if nargin<8 || ~isfield(opts,'theta')
    theta = 0.25;
else
    theta = opts.theta;
end

if nargin<8 || ~isfield(opts,'eps_abs_in')
    eps_abs_in = 1;
else
    eps_abs_in = opts.eps_abs_in;
end

if nargin<8 || ~isfield(opts,'eps_rel_in')
    eps_rel_in = 1;
else
    eps_rel_in = opts.eps_rel_in;
end

if nargin<8 || ~isfield(opts,'scalar_sig')
    scalar_sig = false;
else
    scalar_sig = opts.scalar_sig;
end

if nargin<8 || ~isfield(opts,'maxiter')
    maxiter = 1000;
else
    maxiter = opts.maxiter;
end

if nargin<8 || ~isfield(opts,'inner_maxiter')
    inner_maxiter = 100;
else
    inner_maxiter = opts.inner_maxiter;
end

if nargin<8 || ~isfield(opts,'reset_newton_iter')
    reset_newton_iter = 100;
else
    reset_newton_iter = opts.reset_newton_iter;
end

if nargin<8 || ~isfield(opts,'print_iter')
    print_iter = 1;
else
    print_iter = opts.print_iter;
end


if nargin<8 || ~isfield(opts,'solver')
    solver = 'newton';
else
    solver = opts.solver;
end

if nargin<8 || ~isfield(opts,'memory')
    memory = 10;
else
    memory = opts.memory;
end

if nargin<8 || ~isfield(opts,'lbfgs_precon') %Birgin/Martinez 8.2.3
    lbfgs_precon = false; 
else
    lbfgs_precon = opts.lbfgs_precon;
end

if nargin<8 || ~isfield(opts,'proximal')
    proximal = false;
else
    proximal = opts.proximal;
end

if nargin<8 || ~isfield(opts,'gamma')
    gamma = 1e1;
else
    gamma = opts.gamma;
end

if nargin<8 || ~isfield(opts,'gammaUpd')
    gammaUpd = 1e1;
else
    gammaUpd = opts.gammaUpd;
end

if nargin<8 || ~isfield(opts,'gammaMax')
    gammaMax = 1e7;
else
    gammaMax = opts.gammaMax;
end

%Factorization caching
if nargin<8 || ~isfield(opts,'active_cnstrs') 
    active_cnstrs_old = false(m,1);
    active_cnstrs = false(m,1);
    LD = [];
else
    active_cnstrs_old = opts.active_cnstrs;
    active_cnstrs = opts.active_cnstrs;
    LD = opts.LD;
end
    

if nargin<8 || ~isfield(opts,'nonconvex')
    nonconvex = false;
else
    nonconvex = opts.nonconvex;
end

if nargin<8 || ~isfield(opts,'nonconvex_approx')
    nonconvex_approx = false;
else
    nonconvex_approx = opts.nonconvex_approx;
end

gamma_maxed = false;

if nonconvex 
    if nonconvex_approx
        %do nothing for now
        gamma_maxed = true;
    else
        lambda = lobpcg(Q);
    %     lambda_adj = lambda - 1e-3; %adjust for tolerance
        if lambda < 0
            proximal = true;
            gamma = 1/abs(lambda);
            gammaMax = gamma;
            gamma_maxed = true;
        end
    end
end

% Q = (Q+Q')/2;
% lmin = min(eig(Q));
% fprintf('Lowest eigenvalue of Q is: %.4f\n', lmin);
% if lmin < 0
%     fprintf('The maximal gamma should therefore be: %.4f\n', 1/abs(lmin));
% end

if proximal
    Q = Q+1/gamma*speye(n);
    Qx = Q*x;
else
    gamma = inf;
end

if strcmp(scaling,'outer_iter') %perform this here when sig is known
    [x,Q,q,A,D_scale] = outer_iter_scaling(x,Q,q,A,sig,D_scale);
    nrmq = norm(D_scale\q,inf);
end

K = 1; 
sig_updated = true; %Reset lbfgs initially and perform gradient descent step
reset_newton = true;

%Initialization for Qdx and Adx used in is_dual_infeasible;
tau = 0; Qd = zeros(n,1); Ad = zeros(m,1); d = zeros(n,1);
Qdx = 0; Adx = 0;

%Precompute for ldlupdate
Asqrtsigt = (sparse(1:m,1:m,sqrt(sig),m,m)*A)';
Asig  = (sparse(1:m,1:m,sig,m,m)*A);

y_next = zeros(m,1);
k_prev = 0; 

na = 0;

stats.nact(1) = 0;
stats.nact_changed(1) = 0;



for k = 1:maxiter
    
    if verbose && k > 1 && mod(k,print_iter)==0
        fprintf('iter: %5d,\t nrm_rp: %e,\t nrm_rd: %e,\t nrm_rd2: %e,\t tau: %e \n', k-2, nrm_rp, nrm_rd, nrm_rd2, tau);
    end
    
    stats.gamma(k) = gamma;
%     stats.sigma(:,k) = sig;
    
   Axys = Ax+y./sig;
   z    = min(max(Axys,bmin),bmax);                         % z-update 
   rp   = Ax-z;                                             % primal residual
   yh   = y+sig.*rp;                                        % candidate dual update
   df   = Qx+q-x0/gamma;                                    % cost gradient
   Atyh = A'*yh;
   dphi = df+Atyh;                                          % Augmented Lagrangian gradient

   if proximal
        nrm_rd = 1/c_scale*norm((dphi-1/gamma*(x-x0))./D_scale,inf); %Residual of original problem
        nrm_rd2 = 1/c_scale*norm(dphi./D_scale,inf); %Residual of the proximal problem
   else
       nrm_rd = norm(dphi./D_scale,inf)/c_scale;
       nrm_rd2 = nrm_rd;
   end
   nrm_rp = norm(rp./E_scale,inf);
   nrm_rp_unscaled = norm(rp,inf);
   stats.nrm_rd(k) = nrm_rd;
   stats.nrm_rp(k) = nrm_rp;

   eps_primal   = eps_abs + eps_rel*max(norm(Ax./E_scale,inf),norm(z./E_scale,inf));            % primal eps
   rel_d        = norm([Qx./D_scale;Atyh./D_scale;q./D_scale],inf)/c_scale;
   eps_dual     = eps_abs + eps_rel*rel_d;    % dual eps
   eps_dual_in  = eps_abs_in + eps_rel_in*rel_d;  % inner dual eps
   
   dy = yh-y; Atdy = Atyh - Aty;
   dx = x-x_prev; 
   if proximal 
       Qdx2 = Qdx - tau/gamma*d;
   else
       Qdx2 = Qdx;
   end
          
   if nrm_rd <= eps_dual && nrm_rp <= eps_primal
       stats.status = 'solved';
       if verbose && k > 1 && mod(k,print_iter)==0
        fprintf('iter: %5d,\t nrm_rp: %e,\t nrm_rd: %e,\t nrm_rd2: %e,\t tau: %e \n', k-1, nrm_rp, nrm_rd, nrm_rd2, tau);
       end
       break
   elseif is_primal_infeasible(dy, Atdy, bmin, bmax, D_scale, E_scale, eps_pinf)
       stats.status = 'primal_infeasible';
       stats.pinf_certificate = 1/c_scale*(E_scale.*dy);
       break
   elseif is_dual_infeasible(dx, Qdx2, q, Adx, bmin, bmax, D_scale, E_scale, c_scale, eps_dinf)
       stats.status = 'dual_infeasible';
       stats.dinf_certificate = D_scale.*dx;
       break
   elseif k == k_prev + inner_maxiter %inner problem maxiter termination
       %do outer update except dual and tolerance updates
       k_prev = k;
       gamma_changed = proximal && gamma < gammaMax;
           
       if K > 1 && nrm_rp > eps_primal
           if scalar_sig
               adj_sig  = norm(rp,inf)>theta*norm(rpK,inf) &active_cnstrs;
           else
               adj_sig  = abs(rp)>theta*abs(rpK) &active_cnstrs;
           end
    %            sig  = min((1-(1-Delta).*adj_sig).*sig,1e8);
            prev_sig = sig;
            sig(adj_sig) = min(sigma_max, max(1, Delta*abs(rp(adj_sig))/(nrm_rp_unscaled+1e-6)).*sig(adj_sig));
            sig_changed = sig ~= prev_sig;
            nb_sig_changed = sum(sig_changed);

            if gamma_changed 
                reset_newton = true; 
            elseif nb_sig_changed == 0
                %do nothing
            elseif nb_sig_changed <= 40 
                LD = ldlupdate(LD, (sparse(1:nb_sig_changed, 1:nb_sig_changed, ...
                    sqrt(sig(sig_changed)-prev_sig(sig_changed)), nb_sig_changed, nb_sig_changed)...
                    *A(sig_changed,:))','+');
            else
                reset_newton = true; 
            end            

           Asqrtsigt = (sparse(1:m,1:m,sqrt(sig),m,m)*A)';
           Asig      = (sparse(1:m,1:m,sig,m,m)*A);
       end
       
       rpK = rp;
       if K>1
           stats.iter_in(K) = k-(sum(stats.iter_in));
       else
           stats.iter_in(K) = k;
       end
       K   = K+1;
       if proximal
           x0=x;
           if gamma < gammaMax
               Q=Q-1/gamma*speye(n);
               Qx = Qx-1/gamma*x; 
               gamma=min(gamma*gammaUpd, gammaMax);
               Q=Q+1/gamma*speye(n); %Q = original Q + 1/gamma*eye
               Qx = Qx+1/gamma*x; 
               reset_newton = true; 
           end
       end
       
   elseif nrm_rd2 <= eps_dual_in
       if verbose; fprintf('-------------------------------------------------------------\n'); end;
       k_prev = k;
%        if K == 1 || eps_abs_in ~= eps_abs
           y = yh; Aty = Atyh;
%        else %nesterov acceleration
%            y = yh + (K-2)/(K+1)*(yh-y);
%            Aty = Atyh + (K-2)/(K+1)*(Atyh-Aty);
%        end
       eps_abs_in = max(rho*eps_abs_in,eps_abs);
       eps_rel_in = max(rho*eps_rel_in,eps_rel);
       if K > 1 && nrm_rp > eps_primal

           gamma_changed = proximal && gamma < gammaMax;
           
           if scalar_sig
               adj_sig  = norm(rp,inf)>theta*norm(rpK,inf)&active_cnstrs;
           else
               adj_sig  = abs(rp)>theta*abs(rpK)&active_cnstrs;
           end
           prev_sig = sig;
%            sig  = min((1-(1-Delta).*adj_sig).*sig,1e8);
            sig(adj_sig) = min(sigma_max, max(1, Delta*abs(rp(adj_sig))/(nrm_rp_unscaled+1e-6)).*sig(adj_sig));
            sig_changed = sig ~= prev_sig;
            nb_sig_changed = sum(sig_changed);
            
           if gamma_changed 
            reset_newton = true; 
           elseif nb_sig_changed == 0
                %do nothing
           elseif nb_sig_changed <= 40
               LD = ldlupdate(LD, (sparse(1:nb_sig_changed, 1:nb_sig_changed, ...
                    sqrt(sig(sig_changed)-prev_sig(sig_changed)), nb_sig_changed, nb_sig_changed)...
                    *A(sig_changed,:))','+');
           else
                reset_newton = true; 
           end            
           
           Asqrtsigt = (sparse(1:m,1:m,sqrt(sig),m,m)*A)';
           Asig      = (sparse(1:m,1:m,sig,m,m)*A);
  
       end
       rpK = rp;
       if K>1
           stats.iter_in(K) = k-(sum(stats.iter_in));
       else
           stats.iter_in(K) = k;
       end
       K   = K+1;
       if proximal
           x0=x;
           prev_gamma = gamma;
           active_cnstrs_old = active_cnstrs;
           Axys = Ax+y./sig;
           active_cnstrs = (Axys<=bmin | Axys>=bmax);

           stats.nact_changed(k+1) = sum(abs(active_cnstrs_old - active_cnstrs));
           if K > 2 && ~gamma_maxed && stats.nact_changed(k) == 0 && stats.nact_changed(k+1) == 0 && nrm_rp < eps_primal 
               fprintf('Gamma boosted on iter %d\n', k);
               if na == 0
                   gamma=1e12;
               else
                   gamma=max(1e14/gershgorin_max(A(active_cnstrs,:)'*Asig(active_cnstrs,:)), gammaMax);
                   gamma_maxed = true;
               end
%                
           elseif gamma < gammaMax
               gamma=min(gamma*gammaUpd, gammaMax);
           end
           
           if prev_gamma ~= gamma       
               Q = Q - 1/prev_gamma*speye(n);
               Q = Q + 1/gamma*speye(n);
               Qx = Qx - 1/prev_gamma*x;
               Qx = Qx + 1/gamma*x;
               reset_newton = true;
           end
           
       end
        
       stats.nact(k+1) = na;
       
   else
      
      if strcmp(solver, 'newton') 
          % Newton direction
          active_cnstrs = (Axys<=bmin | Axys>=bmax);
          na = nnz(active_cnstrs);
          stats.nact(k+1) = na;
          if (isempty(active_cnstrs_old))
              stats.nact_changed(k+1) = na;
          else
              stats.nact_changed(k+1) = sum(abs(active_cnstrs-active_cnstrs_old));
              
              if reset_newton
                  stats.nact_changed(k+1) = -1*stats.nact_changed(k+1);
              end
          end
          
          if mod(k,reset_newton_iter)==0
              reset_newton = true;
          end
          prev_gamma = gamma;
          [d,LD,gamma, gammaMax] = computedir(LD,Q,A,Asqrtsigt,Asig,-dphi,active_cnstrs,active_cnstrs_old, reset_newton, na, gamma, gammaMax, nonconvex_approx); %TODO: Refactor
          reset_newton = false;
          if prev_gamma ~= gamma       
               Q = Q - 1/prev_gamma*speye(n);
               Q = Q + 1/gamma*speye(n);
               Qx = Qx - 1/prev_gamma*x;
               Qx = Qx + 1/gamma*x;
%                reset_newton = true;
           end
              
          active_cnstrs_old = active_cnstrs;

      end
      
      Qd = Q*d;
      Ad = A*d;
      % Exact line search
      eta = d'*Qd;
      beta = d'*df; 
      
      if nonconvex_approx 
          if eta < 0 %a direction of negative curvature (after regularization)
              Q = Q - 1/gamma*speye(n);
              Qx = Qx - 1/gamma*x;
              gamma = gamma/10; %Alternatively, use eta -1/gamma*d'*d / norm(d)^2 as an estimate for lambda_min
              gammaMax = gamma;
              reset_newton = true;
              fprintf('Gamma updated (neg curve) to: %.4f\n', gamma);
              Q = Q + 1/gamma*speye(n);
              Qx = Qx + 1/gamma*x;
              continue;
          end
      end
      
      sig_sqr = sqrt(sig);
      delta = -sig_sqr.*Ad;
      delta = [delta;-delta];
      alpha = [(y+sig.*(Ax-bmin))./sig_sqr;(sig.*(bmax-Ax)-y)./sig_sqr];
%       tau = PWALineSearch(eta,beta,delta,alpha);
%       if (stats.nact_changed(k) > 40)
          tau = PWAlinesearch_mex(eta,beta,delta,alpha,int64(2*m));
%       else
%           tau = LineSearchArmijo(eta,beta,Ad,yh,Axys,z,sig,bmax,bmin);
%       end
%       tau = NewtonLS(eta,beta,delta,alpha);
%       tau = BPLS(eta,beta,delta,alpha);
      stats.tau(k) = tau;

      % Store previous values (lbfgs)
      x_prev  = x;
      % Update
      x     = x  + tau*d;
      Adx   = tau*Ad;
      Ax    = Ax + Adx;
      Qdx   = tau*Qd;
      Qx    = Qx + Qdx;
      
   end
end
if k == maxiter
    stats.status = 'maximum number of iterations';
end
stats.iter = k;
stats.iter_out = K;
stats.sig = sig;
stats.LD = LD;
stats.active_cnstrs = active_cnstrs;
if isfield(stats, 'nact_changed')
    stats.sum_nact_changed = sum(abs(stats.nact_changed));
end
if K>1
    stats.iter_in(K) = k-(sum(stats.iter_in));
else
    stats.iter_in(K) = k;
end

%objective
stats.obj = (0.5*x'*(Q-1/gamma*speye(n))*x + q'*x)/c_scale;

%unscale the solution
x = D_scale.*x;
yh = (E_scale.*yh)/c_scale;



end

%% ========================================================================

function is_pinf = is_primal_infeasible(dy, Atdy, bmin, bmax, D_scale, E_scale, eps_pinf) %OSQP
    is_pinf = false;
    eps_pinf_norm_Edy = eps_pinf*norm(E_scale.*dy,inf);
    if eps_pinf_norm_Edy>0
        if norm((Atdy)./D_scale,inf) <= eps_pinf_norm_Edy
            if (bmax'*max(dy,0) + bmin'*min(dy,0)) <= -eps_pinf_norm_Edy;
               is_pinf = true;
            end
        end
    end
end

%% ========================================================================

function is_dinf = is_dual_infeasible(dx, Qdx, q, Adx, bmin, bmax, D_scale, E_scale, c_scale, eps_dinf) %OSQP

eps_dinf_norm_Ddx = eps_dinf*norm(D_scale.*dx,inf);
if eps_dinf_norm_Ddx == 0
    is_dinf = false;return
end
Adx = Adx./E_scale;
if  any((bmax < inf & Adx >= eps_dinf_norm_Ddx) | (bmin > -inf & Adx <= -eps_dinf_norm_Ddx))
    is_dinf = false;return
end

dxQdx = dx'*Qdx;

is_dinf = dxQdx <= -c_scale*eps_dinf_norm_Ddx^2 || ...
    (dxQdx <= c_scale*eps_dinf_norm_Ddx^2 && q'*dx <= -c_scale*eps_dinf_norm_Ddx);

end
%% ========================================================================

function tf = PWALineSearch(eta,beta,delta,gamma)

% inz    = abs(delta)>1e-12 & ~isinf(gamma); 
% if sum(inz) ~= size(delta,1)
%     fprintf('inz \n');
% end
% t = 0:0.001:1;
% psi = zeros(length(t),1);
% for k = 1:length(t)
%     tau = t(k);
%     psi(k) = eta*tau + beta + delta'*max(delta*tau-gamma,0);
% end
% if any(diff(psi) < 0)
%     fprintf('nonconvexity in the linesearch\n');
% end
inz = 1:length(gamma);
gamma  = gamma(inz);
delta  = delta(inz);
s      = gamma./delta;
ns     = length(s);
indVec = 1:ns;
P      = delta>0;
L = s>0;
J = (P&~L)|(~P&L);
a      = eta+delta(J)'*delta(J);
b      = beta-delta(J)'*gamma(J);
indL   = indVec(L);
s = s(indL);
nL = length(indL);
[s,is] = sort(s);
if isempty(is) || a*s(1)+b>0
    tf = -b/a;
    return;
end
i=1;
while i<nL
    iz = indL(is(i));
    if  P(iz)
        a = a+delta(iz)^2;
        b = b-delta(iz)*gamma(iz);
    else
        a = a-delta(iz)^2;
        b = b+delta(iz)*gamma(iz);
    end
    i=i+1;
    if a*s(i)+b>0
        tf=-b/a;
        return;
    end
end
iz = indL(is(i));
if  P(iz)
    a = a+delta(iz)^2;
    b = b-delta(iz)*gamma(iz);
else
    a = a-delta(iz)^2;
    b = b+delta(iz)*gamma(iz);
end
tf=-b/a;

end
%% ========================================================================
function tau = LineSearchArmijo(eta,beta,Ad,yh,Axys,z,sig,bmax,bmin)
    rhs = 1e-4*(beta+Ad'*yh); %c1 = 1e-4
    dist2 = 0.5*sig'*(Axys-z).^2;
    tau = 1; %tau_init = 1
    while true
        lhs = 0.5*tau^2*eta + tau*beta - dist2;
        temp = Axys+tau*Ad-bmax;
        lhs = lhs+0.5*temp(temp>0)'*temp(temp>0);
        temp = temp+bmax-bmin;
        lhs = lhs+0.5*temp(temp<0)'*temp(temp<0);
        if (lhs <= tau*rhs) 
            break; 
        end
        tau = 0.5*tau;
    end  
end

%% ========================================================================
function [x, y, Q, q, A, bmin, bmax, D, E, c] = modified_ruiz_equilibration(x, y, Q,q,A,bmin,bmax,scaling_iter)
    Anzr = max(abs(A),[],2)>0;
    A = A(Anzr,:);
    bmin = bmin(Anzr,1);
    bmax = bmax(Anzr,1);
    [m,n] = size(A);
    c = 1;
    S = ones(n+m,1); 


    for k = 1:scaling_iter
        delta   = 1./sqrt([max(max(abs(Q),[],1),max(abs(A),[],1))';max(abs(A),[],2)]);
        ddelta1 = sparse(1:n,1:n,delta(1:n),n,n);
        ddelta2 = sparse(1:m,1:m,delta(n+1:n+m),m,m);
        Q = ddelta1*(Q*ddelta1); q = ddelta1*q; 
        A = ddelta2*(A*ddelta1); 
%         gamma = 1/max(mean(max(abs(Q),[],2)), norm(q,Inf));
%         Q = gamma*Q; q = gamma*q;c = gamma*c;
        S = delta.*S; 
    end
    
    D = S(1:n);
    E = S(n+1:n+m);
    bmin = E.*bmin; bmax = E.*bmax;

    x = x./D;
    y = c*(y./E);

end

%% ========================================================================
function [x,y,Q,q,A,bmin,bmax,D,E,c] = simple_equilibration(x,y,Q,q,A,bmin,bmax,scaling_iter)
    Anzr = max(abs(A),[],2)>0;
    A = A(Anzr,:);
    bmin = bmin(Anzr,1);
    bmax = bmax(Anzr,1);
    y = y(Anzr,1);
    [m,n] = size(A);
    
    %Make Q have unit diagonal elements
%     D = sparse(diag(1./sqrt(diag(Q)))); 
%     Q = D*Q*D; q = D*q; A = A*D;
    
%     Ruiz scaling on A (works well)
    E = ones(m,1);
    D = ones(n,1);
    colinds = 1:m;
    rowinds = 1:n;
    for k=1:scaling_iter
        Aabs = abs(A);
        r = 1./sqrt(max(Aabs,[],2));
        cint = sqrt(max(Aabs,[],1));
        cint(cint < 1e-12) = 1;
        c = 1./cint';
        R = sparse(colinds,colinds,r,m,m);
        C = sparse(rowinds,rowinds,c,n,n);
        A = R*A*C;
        E = E.*r;
        D = D.*c;
    end
    Dm = sparse(rowinds,rowinds,D,n,n);
    Q = Dm*(Q*Dm);
    q = D.*q;
    bmin = E.*bmin; bmax = E.*bmax;

    x = x./D; y = (y./E);

    c = 1;
    if (scaling_iter)
        c = 1/max(1, norm(Q*x+q,inf)); %Add cost scaling 10.2.2 Birgin/Martinez
    end
    Q = c*Q; q=c*q;
   
    y = c*y;

    
end

%% ========================================================================
function [x,Q,q,A,D] = outer_iter_scaling(x,Q,q,A,sig,D_scale)

    D = diag(1./sqrt(diag(Q+A'*diag(sig)*A))); D = sparse(D);
    x = D\x; Q = D*Q*D; q = D*q; A = A*D;
    
    D = D_scale*D; %Keep track of total scaling

end