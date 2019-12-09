function [ x, y, failure_detected] = make_performance_profile( r )
%make the performance profile (stepwise graph)
x = [1];
y = [0];

r = sort(r);

failure_detected = false;

y_temp = 0;

T = length(r);

for k = 1:T
    if (r(k) == inf || isnan(r(k)))
        failure_detected = true;
    elseif k == 1
        x_temp = r(k);
        y_temp = 1/T;
        if r(k) ~= 1
            x = [x x_temp];
            y = [y 0];
        end
            
    else
        if r(k) ~= r(k-1)
            x = [x x_temp r(k)];
            y = [y y_temp y_temp];
            
            y_temp = y_temp + 1/T;
            x_temp = r(k);
        else
            y_temp = y_temp + 1/T;
        end
    end
end


end

