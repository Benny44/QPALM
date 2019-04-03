function [x,yh,stats] = qpalm_matlab(Q,q,A,bmin,bmax,x,y,opts)
[m,n] = size(A);
if nargin<7 || isempty(y)
    y = zeros(m,1);
end

if nargin<6 || isempty(x)
    x = zeros(n,1);
end

if nargin<8 || ~isfield(opts,'scaling')
    scaling = ''; %'ruiz', 'simple', 'outer_iter'
else
    scaling = opts.scaling;
end

if nargin<8 || ~isfield(opts,'scaling_iter')
    scaling_iter = 10; %'ruiz'
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

if nargin<8 || ~isfield(opts,'sig')
    f = 0.5*(x'*Qx)+q'*x;
    dist = Ax-min(max(Ax,bmin),bmax);
    dist2 = dist'*dist;
    sig = max(1e-8,min(2e1*max(1,abs(f))/max(1,0.5*dist2),1e8))*ones(m,1);
else
    sig = opts.sig.*ones(m,1);
end

if nargin<8 || ~isfield(opts,'linsys')
    linsys = 2;
else
    linsys = opts.linsys;
end

if linsys == 1
    linopts.POSDEF = true;
    linopts.SYM = true;
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
    eps_pinf = 1e-4;
else
    eps_pinf = opts.eps_pinf;
end

if nargin<8 || ~isfield(opts,'eps_dinf')
    eps_dinf = 1e-4;
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

if nargin<8 || ~isfield(opts,'Delta')
    Delta = 10;
else
    Delta = opts.Delta;
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
    gamma = 1e6;
else
    gamma = opts.gamma;
end

if nargin<8 || ~isfield(opts,'gammaUpd')
    gammaUpd = 1e1;
else
    gammaUpd = opts.gammaUpd;
end

if nargin<8 || ~isfield(opts,'gammaMax')
    gammaMax = 1e8;
else
    gammaMax = opts.gammaMax;
end

%Factorization caching
if nargin<8 || ~isfield(opts,'active_cnstrs') 
    active_cnstrs_old = [];
    active_cnstrs = [];
    LD = [];
else
    active_cnstrs_old = opts.active_cnstrs;
    active_cnstrs = opts.active_cnstrs;
    LD = opts.LD;
end
    

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


%Initialization for Qdx and Adx used in is_dual_infeasible;
tau = 0; Qd = zeros(n,1); Ad = zeros(m,1);

%Precompute for ldlupdate
Asqrtsigt = (sparse(1:m,1:m,sqrt(sig),m,m)*A)';
Asig  = (sparse(1:m,1:m,sig,m,m)*A);

for k = 1:maxiter
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
   stats.nrm_rd(k) = nrm_rd;
   stats.nrm_rp(k) = nrm_rp;

   eps_primal   = eps_abs + eps_rel*max(norm(Ax./E_scale,inf),norm(z./E_scale,inf));            % primal eps
   rel_d        = norm([Qx./D_scale;Atyh./D_scale;q./D_scale],inf)/c_scale;
   eps_dual     = eps_abs + eps_rel*rel_d;    % dual eps
   eps_dual_in  = eps_abs_in + eps_rel_in*rel_d;  % inner dual eps
   
   dy = yh-y; Atdy = Atyh - Aty;
   dx = x-x_prev; Qdx = Qd*tau; Adx = Ad*tau;
   
   if nrm_rd <= eps_dual && nrm_rp <= eps_primal
       stats.status = 'solved';
       break
   elseif is_primal_infeasible(dy, Atdy, bmin, bmax, D_scale, E_scale, eps_pinf)
       stats.status = 'primal_infeasible';
       stats.pinf_certificate = 1/c_scale*(E_scale.*dy);
       break
   elseif is_dual_infeasible(dx, Qdx, q, Adx, bmin, bmax, D_scale, E_scale, c_scale, eps_dinf)
       stats.status = 'dual_infeasible';
       stats.dinf_certificate = D_scale.*dx;
       break
   elseif nrm_rd2 <= eps_dual_in
       y = yh; Aty = Atyh;
       eps_abs_in = max(rho*eps_abs_in,eps_abs);
       eps_rel_in = max(rho*eps_rel_in,eps_rel);
       if K > 1 && nrm_rp > eps_primal

           if scalar_sig
               adj_sig  = norm(rp,inf)>theta*norm(rpK,inf);
           else
               adj_sig  = abs(rp)>theta*abs(rpK);
           end
           sig  = min((1-(1-Delta).*adj_sig).*sig,1e8);
           sig_updated = true;
           
           Asqrtsigt = (sparse(1:m,1:m,sqrt(sig),m,m)*A)';
           Asig      = (sparse(1:m,1:m,sig,m,m)*A);
           
%            if strcmp(scaling,'outer_iter') %perform this here when sig is known
%                [x,Q,q,A,D_scale] = outer_iter_scaling(x,Q,q,A,sig,D_scale);
%                Qx = Q*x; 
%            end
  
       end
       rpK = rp;
       if K>1
           stats.iter_in(K) = k-(sum(stats.iter_in));
       else
           stats.iter_in(K) = k;
       end
       K   = K+1;
       active_cnstrs_old = [];
       if proximal
           x0=x;
           if gamma ~= gammaMax
               Q=Q-1/gamma*speye(n);
               gamma=min(gamma*gammaUpd, gammaMax);
               Q=Q+1/gamma*speye(n); %Q = original Q + 1/gamma*eye
               Qx = Q*x; 
           end
       end
   else
      if strcmp(solver, 'newton') 
          % Newton direction
          active_cnstrs = (Axys<bmin | Axys>bmax);
          na = nnz(active_cnstrs);
          stats.nact(k) = na;
          if na
              switch linsys
                  case 0 % sparse backslash
                      d = -(Q+A(active_cnstrs,:)'*sparse(1:na,1:na,sig(active_cnstrs),na,na)*A(active_cnstrs,:))\dphi;
                  case 1% Cholesky via backslash 1
                      if na>1
                          d = -linsolve((Q+A(active_cnstrs,:)'*diag(sig(active_cnstrs))*A(active_cnstrs,:)),dphi,linopts);
                      else
                          d = -(Q+A(active_cnstrs,:)'*diag(sig(active_cnstrs))*A(active_cnstrs,:))\dphi;
                      end
                  case 2% ldlchol 2
                      [d,LD] = computedir(LD,Q,A,Asqrtsigt,Asig,-dphi,active_cnstrs,active_cnstrs_old);
    %                   LD = ldlchol(Q+A(active_cnstrs,:)'*spdiags(sig(active_cnstrs),0,na,na)*A(active_cnstrs,:));
    %                   d  = -ldlsolve (LD,dphi);
                  case 3% lchol  3
                      [L,~,p] = lchol(Q+A(active_cnstrs,:)'*sparse(1:na,1:na,sig(active_cnstrs),na,na)*A(active_cnstrs,:));
                      d(p,:)  = -L'\(L\(dphi(p,:)));
                  case 4
                      % chol 4
                      [L,~,p]   = chol(Q+A(active_cnstrs,:)'*sparse(1:na,1:na,sig(active_cnstrs),na,na)*A(active_cnstrs,:),'lower','vector');
                      d(p,:)    = -L'\(L\(dphi(p,:)));
                  case 5% Backslash 5
                      dlam  = -[Q A(active_cnstrs,:)';A(active_cnstrs,:) -sparse(1:na,1:na,1./sig(active_cnstrs),na,na)]\[dphi;zeros(na,1)];
                      d     = dlam(1:n,1);
                  case 6% Indefinite LDL 6
                      [L, D, p] = ldl([Q A(active_cnstrs,:)';A(active_cnstrs,:) -sparse(1:na,1:na,1./sig(active_cnstrs),na,na)], 'vector');
                      dlam      = -[dphi;zeros(na,1)];
                      dlam(p,:) = L'\(D\(L\(dlam(p,:))));
                      d         = dlam(1:n,1);
                  case 7% ldlsparse
                      dlam  = -ldlsparse([Q A(active_cnstrs,:)';A(active_cnstrs,:) -sparse(1:na,1:na,1./sig(active_cnstrs),na,na)],[],[dphi;zeros(na,1)]);
                      d     = dlam(1:n,1);
              end
              active_cnstrs_old = active_cnstrs;
          else
              LD = ldlchol(Q);
              d = -ldlsolve (LD,dphi);
          end
          
      elseif strcmp(solver, 'lbfgs')
          % lbfgs direction
          if sig_updated || memory==0 %do gradient descent and start over in lbfgs after sigma is updated
              d = -dphi;
              sig_updated = false;
              mem = 0; 
              curridx = 0;
              Sbuffer = [];Ybuffer = []; YSbuffer = [];
          else
              lbfgs_s = x - x_prev;
              lbfgs_y = dphi-dphi_prev;

              lbfgs_ys  = lbfgs_y'*lbfgs_s;
              
              if abs(lbfgs_ys) >= 1e-8*norm(lbfgs_s,inf)*norm(lbfgs_y,inf)
                  
                  if mem == memory
                      curridx = 1;
                      Sbuffer   = [Sbuffer(:,2:end), lbfgs_s];
                      Ybuffer   = [Ybuffer(:,2:end), lbfgs_y];
                      YSbuffer  = [YSbuffer(2:end), lbfgs_ys];
                  else
                      mem = mem + 1;
                      curridx = curridx+1;
                      Sbuffer(:,mem) = lbfgs_s  ;
                      Ybuffer(:,mem) = lbfgs_y  ;
                      YSbuffer(mem)  = lbfgs_ys;
                  end
              end
              
              
              if lbfgs_precon
                  %Birgin/Martinez 8.2.3
                  active_cnstrs = (Axys<bmin | Axys>bmax);
                  Ac = A(active_cnstrs,:)'*diag(sig(active_cnstrs))*A(active_cnstrs,:);
                  Dc = diag(Ac);
                  sigma_spec = (lbfgs_y-Dc.*lbfgs_s)'*lbfgs_s/(lbfgs_s'*lbfgs_s);
                  sigma = max(1e-12, min(sigma_spec, 1e12));
                  Dc = Dc + sigma;
                  Dcinv = diag(1./Dc);
                  if lbfgs_ys <= 1e-8*norm(lbfgs_s)*norm(lbfgs_y)
                      Hp = Dcinv;
                  else
                      Hp = Dcinv + ((lbfgs_s+Dcinv*lbfgs_y)*lbfgs_s' + ...
                          lbfgs_s*(lbfgs_s-Dcinv*lbfgs_y)')/lbfgs_ys - ...
                          (lbfgs_s-Dcinv*lbfgs_y)'*lbfgs_y*(lbfgs_s*lbfgs_s')/lbfgs_ys^2; 
                  end
              else
                  Hp = lbfgs_ys/(lbfgs_y'*lbfgs_y);
              end

              alpha = zeros(1, mem);
              qq = -dphi;
              for i = mem:-1:1
                  alpha(i)  = 1/YSbuffer(i) * (Sbuffer(:,i)'*qq);
                  qq = qq - alpha(i)*Ybuffer(:,i);
              end

              d = Hp*qq;
              for i = 1:mem
                  bet = 1/YSbuffer(i) * Ybuffer(:,i)'*d;
                  d   = d + (alpha(i)-bet)*Sbuffer(:,i);
              end
            
%               d = lbfgs(Sbuffer, Ybuffer, YSbuffer, Hp,-dphi, int32(mem), int32(memory));
          end
      end
      
      Qd = Q*d;
      Ad = A*d;
      % Exact line search
      eta = d'*Qd;
      beta = d'*df;
      sig_sqr = sqrt(sig);
      delta = -sig_sqr.*Ad;
      delta = [delta;-delta];
      alpha = [(y+sig.*(Ax-bmin))./sig_sqr;(sig.*(bmax-Ax)-y)./sig_sqr];
%       tau = PWALineSearch(eta,beta,delta,alpha);
      tau = PWAlinesearch_mex(eta,beta,delta,alpha,int64(2*m));
%       tau = NewtonLS(eta,beta,delta,alpha);
%       tau = BPLS(eta,beta,delta,alpha);
      stats.tau(k) = tau;

      % Store previous values (lbfgs)
      x_prev  = x;
      dphi_prev = dphi;
      % Update
      x     = x  + tau*d;
      Ax    = Ax + tau*Ad;
      Qx    = Qx + tau*Qd;
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
if K>1
    stats.iter_in(K) = k-(sum(stats.iter_in));
else
    stats.iter_in(K) = k;
end

%unscale
x = D_scale.*x;
y = (E_scale.*yh)/c_scale;

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
%     is_pinf = eps_pinf_norm_Edy > 0 ... %dy must be nonzero
%         && norm((Atdy)./D_scale,inf) <= eps_pinf_norm_Edy ...
%         && (bmax'*max(dy,0) + bmin'*min(dy,0)) <= -eps_pinf_norm_Edy;

end

%% ========================================================================

function is_dinf = is_dual_infeasible(dx, Qdx, q, Adx, bmin, bmax, D_scale, E_scale, c_scale, eps_dinf) %OSQP

eps_dinf_norm_Ddx = eps_dinf*norm(D_scale.*dx,inf);
if eps_dinf_norm_Ddx == 0
    is_dinf = false;return
end
Adx = Adx./E_scale;
if  any(bmax < inf & Adx >= eps_dinf_norm_Ddx) | (bmin > -inf & Adx <= -eps_dinf_norm_Ddx)
    is_dinf = false;return
end
% for k = 1:length(bmax),if (bmax(k) < 1e20 && Adx(k) >= eps_dinf_norm_Ddx) || (bmin(k) > -1e20 && Adx(k) <= -eps_dinf_norm_Ddx),is_dinf = false; return;end,end

is_dinf = norm(Qdx./D_scale,inf) <= c_scale*eps_dinf_norm_Ddx ...
    && q'*dx <= -c_scale*eps_dinf_norm_Ddx;

end
%% ========================================================================

function tf = PWALineSearch(eta,beta,delta,gamma)

inz    = abs(delta)>1e-12 & ~isinf(gamma); 
% if sum(inz) ~= size(delta,1)
%     fprintf('inz \n');
% end
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
end
tf=-b/a;

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
        c = 1./sqrt(max(Aabs,[],1))';
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
% Make A have unit row norms
% D = speye(n);
% E = sparse(1:m,1:m,1./vecnorm(A,2,2),m,m);
% A = E*A; bmin = E*bmin; bmax = E*bmax;
    c = 1;
    c = 1/max(1, norm(Q*x+q,inf)); %Add cost scaling 10.2.2 Birgin/Martinez
    Q = c*Q; q=c*q;
   

    x = x./D; y = c*(y./E);
end

%% ========================================================================
function [x,Q,q,A,D] = outer_iter_scaling(x,Q,q,A,sig,D_scale)

    D = diag(1./sqrt(diag(Q+A'*diag(sig)*A))); D = sparse(D);
    x = D\x; Q = D*Q*D; q = D*q; A = A*D;
    
    D = D_scale*D; %Keep track of total scaling

end