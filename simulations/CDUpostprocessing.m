%Calculate mean and max runtime and iterations

load('results/CDUqpalm')
meanIterQpalm = mean(Iter_qpalm_matlab);
maxIterQpalm = max(Iter_qpalm_matlab);

meanTimingQpalm = mean(Tqpalm_matlab)*1000; %in ms
maxTimingQpalm = max(Tqpalm_matlab)*1000; %in ms

fprintf('QPALM: %d %d %f %f\n', meanIterQpalm, maxIterQpalm, meanTimingQpalm, maxTimingQpalm);

load('results/CDU');
meanIterQpoases = mean(Iter_qpoases);
maxIterQpoases = max(Iter_qpoases);

meanTimingQpoases = mean(Tqpoases)*1000; %in ms
maxTimingQpoases = max(Tqpoases)*1000; %in ms

fprintf('qpOASES: %d %d %f %f\n', meanIterQpoases, maxIterQpoases, meanTimingQpoases, maxTimingQpoases);

meanIterOsqp = mean(Iter_osqp);
maxIterOsqp = max(Iter_osqp);

meanTimingOsqp = mean(Tosqp)*1000; %in ms
maxTimingOsqp = max(Tosqp)*1000; %in ms

fprintf('OSQP: %d %d %f %f\n', meanIterOsqp, maxIterOsqp, meanTimingOsqp, maxTimingOsqp);
