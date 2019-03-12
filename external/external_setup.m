%Setup LBFGS and OSQP
clear all
close all

% %% LBFGS setup
% 
% current_path = fileparts(mfilename('fullpath'));
% LBFGS_path = fullfile(current_path, 'external','lbfgs');
% 
% error_msg = 'The C compiler could not succesfully compile ';
% if mex('-outdir', LBFGS_path, fullfile(LBFGS_path,'lbfgs.c'), '-compatibleArrayDims'), error([error_msg, LBFGS_path]); end

%% OSQP setup

external_path = fullfile(current_path, 'external');
addpath(external_path);
cd(external_path);
websave('install_osqp.m', 'https://dl.bintray.com/bstellato/generic/OSQP/0.5.0/install_osqp.m');
install_osqp
!git clone --recurse-submodules https://github.com/oxfordcontrol/osqp-matlab
cd osqp-matlab
make_osqp

%% Savepath
addpath(genpath(current_path));
savepath;

cd(current_path);

clear current_path LBFGS_path external_path error_msg
clc