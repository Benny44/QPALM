function [d,LD,gamma, gammaMax] = computedir(LD,Q,A,Asqrtsigt,Asig,b,active_cnstrs,active_cnstrs_old, reset_newton, na, gamma, gammaMax, nonconvex_approx)
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
        if nonconvex_approx
            [LD,p] = ldlchol(Q+AsigA);
            while p || (condest(LD)*eps > .1)
            Q = Q - 1/gamma*speye(size(Q,1));
            gamma = gamma/10; %regularize more heavily. TODO: find a smart way to determine how heavy to regularize
            gammaMax = gamma;
            fprintf('Gamma updated to: %.4f\n', gamma);
            pause(0.001)
            Q = Q + 1/gamma*speye(size(Q,1)); 
            [LD, p] = ldlchol(Q+AsigA);
            end
        else
            LD = ldlchol(Q+AsigA);
        end
    end
else
    if nonconvex_approx
        [LD, p] = ldlchol(Q);
        while p || (condest(LD)*eps > 1)
            Q = Q - 1/gamma*speye(size(Q,1));
            gamma = gamma/10; %regularize more heavily. TODO: find a smart way to determine how heavy to regularize
            gammaMax = gamma;
            Q = Q + 1/gamma*speye(size(Q,1)); 
            [LD, p] = ldlchol(Q);
        end
    else
        LD = ldlchol(Q);
    end
end
d = ldlsolve (LD,b);

