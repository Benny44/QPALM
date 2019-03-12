%Run setup for the varied solvers

current = fileparts(mfilename('fullpath'));
cd(current);
%% qpOASES
cd('../external/qpOASES-3.2.1/interfaces/matlab')
make

%% OSQP
cd(current);
cd('../external/osqp-matlab')
make_osqp

%% QPALM
cd(current);
cd('..')
qpalm_setup

cd(current);
