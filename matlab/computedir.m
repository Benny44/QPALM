function [d,LD] = computedir(LD,Q,A,Asqrtsigt,Asig,b,active_cnstrs,active_cnstrs_old, reset_newton, na)
if na
    enter = active_cnstrs & (~active_cnstrs_old);
    leave = active_cnstrs_old & (~active_cnstrs);
    ne = sum(enter);
    nl = sum(leave);
    if ~reset_newton && (ne+nl)<=160
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
        AsigA = A(active_cnstrs,:)'*Asig_active_cnstrs;
        LD = ldlchol(Q+AsigA);
    end
else
    LD = ldlchol(Q);
end
d = ldlsolve (LD,b);
