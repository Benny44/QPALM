% Load the QP via the cutest interface, save all the relevant parameters 
% and store the result as a mat file for later
% This script is supposed to be run in load_cutest_and_save.sh

Data = cutest_setup();
Data.name = Data.name(~isspace(Data.name));
folder = '/home/ben/Documents/Projects/QPALM/simulations/cutest_set';

Data.Q = cutest_isphess(zeros(Data.n,1), 0);
Data.q = cutest_grad(zeros(Data.n,1));

if Data.m > 0
    [~, Data.A] = cutest_scons(zeros(Data.n,1)); %Does not include bounds on variables
    c = cutest_cons(zeros(Data.n,1));
    up = Data.cu ~= 1e20;
    lo = Data.cl ~= -1e20;
    Data.cu(up) = Data.cu(up) - c(up);
    Data.cl(lo) = Data.cl(lo) - c(lo);
else
    Data.A = [];
end

fudge = 0;
if (nnz(Data.Q) > 10000)
    addpath('/home/ben/Documents/Projects/QPALMmatlab/matlab');
    try
        lam_min = lobpcg(Data.Q);
        fudge = -1e-5;
    catch
        lam_min = min(eig(Data.Q));
    end
else
    lam_min = min(eig(Data.Q));
end
if lam_min < fudge
    output_file = [folder '/' 'nonconvexQPs' '/' Data.name];
else
    output_file = [folder '/' 'convexQPs' '/' Data.name];
end
save(output_file, 'Data');

quit;
