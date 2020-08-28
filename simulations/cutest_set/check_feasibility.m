% load('/home/ben/Documents/Projects/QPALM/simulations/cutest_set/nonconvexQPs/A/A2NNDNIL.mat')
load('/home/ben/Documents/Projects/QPALM/simulations/cutest_set/nonconvexQPs/A/A5NNDNIL.mat')
% load('/home/ben/Documents/Projects/QPALM/simulations/cutest_set/nonconvexQPs/A/A2ENSNDL.mat')

A = [Data.A; speye(Data.n)];
A = [A; -A];
b = [Data.cu; Data.bu; -Data.cl; -Data.bl];

[x,fval,exitflag,OUTPUT] = linprog(zeros(Data.n,1),A,b);