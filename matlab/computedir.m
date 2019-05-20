function [d,LD] = computedir(LD,Q,A,Asqrtsigt,Asig,b,active_cnstrs,active_cnstrs_old, reset_newton)
enter = active_cnstrs & (~active_cnstrs_old);
leave = active_cnstrs_old & (~active_cnstrs);
ne = sum(enter);
nl = sum(leave);
if ~reset_newton && (ne+nl)<=40
    if ne>0
        Ae = Asqrtsigt(:,enter);
        LD = ldlupdate(LD,Ae,'+');
    end
    if nl>0
        Al = Asqrtsigt(:,leave);
        LD = ldlupdate(LD,Al,'-');
    end
else
    Asig_active_cnstrs = Asig(active_cnstrs,:);
    LD = ldlchol(Q+A(active_cnstrs,:)'*Asig_active_cnstrs);
end
d = ldlsolve (LD,b);
