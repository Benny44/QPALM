function tau = BPLS(eta,beta,delta,alpha)
tL = -inf;tU = inf;
s = alpha./delta;
T = s;
ip = delta>0;
np = ~ip;

for k=1:1000
%     tau = median(T);
    tau = fast_median_ip(T);
%     tau = fast_median(T);
%     Ip = (tau>s & ip)|(tau<s & np);
%     slope = eta+delta(Ip)'*delta(Ip);
%     intercept = beta-delta(Ip)'*alpha(Ip);
%     psi = slope*tau+intercept;
    psi = eta*tau+beta+delta'*max(delta*tau-alpha,0);
    if psi<0
        tL = tau;
        T = T(tau<T);
    else
        tU = tau;
        T = T(tau>T);
    end
    if isempty(T)
        psiL = eta*tL+beta+delta'*max(delta*tL-alpha,0);
        psiU = eta*tU+beta+delta'*max(delta*tU-alpha,0);
        tau = tL-psiL*(tU-tL)/(psiU-psiL);
        break
    end
end
