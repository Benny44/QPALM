function [gs, fail_rate] = compute_geometric_mean(T, Status, success, max_time)
gs = 0;
fail_rate = 0;

shift = 1;
n = length(T);
for k = 1:n
    if ~strcmp(Status{k}, success)
        T(k) = max_time;
        fail_rate = fail_rate + 1;
    elseif T(k) > max_time
        T(k) = max_time;
        fail_rate = fail_rate + 1;
    end
    gs = gs + log(T(k)+shift)/n;

end

gs = exp(gs) - shift;

fail_rate = 100*fail_rate/n; %in %