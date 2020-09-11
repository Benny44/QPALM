function [ x, y, failure_detected] = make_performance_profile( r )
%make the performance profile (stepwise graph)
x = [1];
y = [0];

r = sort(r);

failure_detected = false;

T = length(r);

for k = 1:T
    if (r(k) == inf || isnan(r(k)))
        failure_detected = true;
        if (r(k) == inf) 
            %deal with the plotting extension outside of this function
            %that is, extend x = [x rmax]; y = [y y(end)];
            break; 
        end
    elseif k == 1
        if r(k) ~= 1
            x = [x r(k) r(k)];
            y = [y 0 k/T];
        else
            x = [x r(k)];
            y = [y k/T];
        end
            
    else
        if r(k) ~= r(k-1)
            x = [x r(k) r(k)];
            y = [y (k-1)/T k/T];
        else
            x = [x r(k)];
            y = [y k/T];
        end
    end
end

