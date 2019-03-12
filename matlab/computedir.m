function [d,LD] = computedir(LD,Q,A,Asqrtsigt,Asig,b,active_cnstrs,active_cnstrs_old)
na = sum(active_cnstrs);
if ~isempty(active_cnstrs_old)
    enter = active_cnstrs & (~active_cnstrs_old);
    leave = active_cnstrs_old & (~active_cnstrs);
    ne = nnz(enter);
    nl = sum(leave);
    if ne>0
%         Ap = (sparse(1:ne,1:ne,sqrt(sig(enter)),ne,ne)*A(enter,:));
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
