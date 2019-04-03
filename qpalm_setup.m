% clear
close all

%% Change to correct directory
current_path = fileparts(mfilename('fullpath'));
cd(current_path);

%% Create the output directories if they don't exist
if ~exist('./obj','dir')
    mkdir('obj');
end

if ~exist('./lib','dir')
    mkdir('lib');
end
%% Build QPALM library

system('make clean');
system('make lib');

%% Mex interface
mex_path = fullfile(current_path, 'mexInterface');
addpath(mex_path);
savepath;

remove_existing_mex = sprintf('rm -f %s/*.mex*', mex_path);
system(remove_existing_mex);

error_msg = 'The C compiler could not succesfully compile ';
if mex('-outdir', mex_path, fullfile(mex_path,'qpalm_mex.c'),...
        '-Iinclude', '-lqpalm','-Llib', '-lm', '-largeArrayDims',...
        '-ISuiteSparse/include', '-LSuiteSparse/lib',...
        '-lcholmod', '-lamd', '-lcolamd', '-lsuitesparseconfig',...
        'CFLAGS="\$CFLAGS -std=c99 -fPIC -DMATLAB -O3 -DPROFILING -DPRINTING"')
    error([error_msg, mex_path]);
end
   
 
 
 
 
