function [d,LD] = computedir(LD,Q,A,sig,b,active_cnstrs,active_cnstrs_old)
na = nnz(active_cnstrs);
if ~isempty(active_cnstrs_old)
    enter = active_cnstrs & (~active_cnstrs_old);
    leave = active_cnstrs_old & (~active_cnstrs);
    ne = nnz(enter);
    nl = nnz(leave);
    if ne>0
        LD = ldlupdate(LD,(sparse(1:ne,1:ne,sqrt(sig(enter)),ne,ne)*A(enter,:))','+');
    end
    if nl>0
        LD = ldlupdate(LD,(sparse(1:nl,1:nl,sqrt(sig(leave)),nl,nl)*A(leave,:))','-');
    end
else
    LD = ldlchol(Q+A(active_cnstrs,:)'*(sparse(1:na,1:na,sig(active_cnstrs),na,na)*A(active_cnstrs,:)));
end
d = ldlsolve (LD,b);
