function tau = NewtonLS(eta,beta,delta,alpha)
% Newton's method starting from right-most breakpoint
% Adopted from Nesterov, "Lectures on Convex Optimization", Appendix A

zd = (delta==0);
alpha(zd) = [];
delta(zd) = [];
s = alpha./delta;
Pd = delta>0;
Nd = ~Pd;
tau = max(s);
for k=1:100
   A = (tau>s & Pd)|(tau<s & Nd);
   intercept = beta - delta(A)'*alpha(A);
   slope     = eta  + delta(A)'*delta(A);
   if abs(slope*tau+intercept)<1e-12
       break
   end
   tau = -intercept/slope;
end

