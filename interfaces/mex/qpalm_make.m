function qpalm_make (solver)
%QPALM_MAKE compiles the QPALM mexFunction
%
% Example:
%   qpalm_make
%   qpalm_make('ladel')
%   qpalm_make('cholmod')
%
% QPALM relies either on LADEL or on CHOLMOD, AMD, COLAMD and optionally CCOLAMD, CAMD, and METIS.
warning('off','MATLAB:mex:GccVersion_link')

close all

% store the path from which this function is called
current_path = pwd;

% Change to directory containing this function
this_path = fileparts(mfilename('fullpath'));
cd(this_path);

!git submodule update --init
!git submodule update --recursive --remote


%% Determine the linear systems solver and specify the include paths
fprintf ('Compiling QPALM mex function\n')

details = 0 ;	    % 1 if details of each command are to be printed

v = version ;
try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
    mac = ismac ;
catch                                                                       %#ok
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
    mac = 0 ;
end

flags = '' ;
is64 = ~isempty (strfind (computer, '64')) ;
if (is64)
    % 64-bit MATLAB
    flags = '-largeArrayDims' ;
end

% MATLAB 8.3.0 now has a -silent option to keep 'mex' from burbling too much
if (~verLessThan ('matlab', '8.3.0'))
    flags = ['-silent ' flags] ;
end

if nargin < 1
    solver = 'ladel';
    fprintf('No linear systems solver selected. Using LADEL as the default.\n');
    fprintf('Use qpalm_make(''cholmod'') if you wish to compile with CHOLMOD instead.\n');
elseif strcmp(solver, 'ladel')
    fprintf('Using LADEL as the linear systems solver.\n');
elseif strcmp(solver, 'cholmod');
    fprintf('Using CHOLMOD as the linear systems solver.\n');
else
    error('Unrecognised linear systems solver. Use ''ladel'' or ''cholmod''.\n');
end

if strcmp(solver, 'ladel')
    include = '-I. -I../../include -I../../LADEL/include -I../../LADEL/amd/Include -I../../LADEL/interfaces/mex/include';
    flags = [flags ' -DUSE_LADEL'];
elseif strcmp(solver, 'cholmod')
    include = '-I. -I../../include -I./ -I../../suitesparse/AMD/Include -I../../suitesparse/COLAMD/Include -I../../suitesparse/CCOLAMD/Include -I../../suitesparse/CAMD/Include -I../../suitesparse/CHOLMOD/Include -I../../suitesparse/CHOLMOD/MATLAB -I../../suitesparse/SuiteSparse_config' ;
    flags = [flags ' -DUSE_CHOLMOD'];
    if (verLessThan ('matlab', '7.0'))
        % do not attempt to compile CHOLMOD with large file support
        include = [include ' -DNLARGEFILE'] ;
    elseif (~pc)
        % Linux/Unix require these flags for large file support
        include = [include ' -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE'] ;
    end
    
    if (verLessThan ('matlab', '6.5'))
        % logical class does not exist in MATLAB 6.1 or earlie
        include = [include ' -DMATLAB6p1_OR_EARLIER'] ;
    end

     % Determine if METIS is available
    if strcmp(solver, 'cholmod')
        metis_path = '../../suitesparse/metis-5.1.0' ;
    end
    have_metis = exist (metis_path, 'dir') ;
    
    if (have_metis)
    %     fprintf ('Compiling CHOLMOD with METIS 5.1.0 for MATLAB Version %s\n', v) ;
        include = [include ' -I' metis_path '/include'] ;
        include = [include ' -I' metis_path '/GKlib'] ;
        include = [include ' -I' metis_path '/libmetis'] ;
    else
    %     fprintf ('Compiling CHOLMOD without METIS for MATLAB Version %s\n', v) ;
        include = ['-DNPARTITION ' include] ;
    end

    flags = [flags ' -DMATLAB_MEX_FILE -DNMALLOC'];
end

 %---------------------------------------------------------------------------
 % BLAS option
 %---------------------------------------------------------------------------

 % This is exceedingly ugly.  The MATLAB mex command needs to be told where to
 % find the LAPACK and BLAS libraries, which is a real portability nightmare.

if (pc)
    if (verLessThan ('matlab', '7.5'))
        lapack = 'libmwlapack.lib' ;
    else
        lapack = 'libmwlapack.lib libmwblas.lib' ;
    end
else
    if (verLessThan ('matlab', '7.5'))
        lapack = '-lmwlapack' ;
    else
        lapack = '-lmwlapack -lmwblas' ;
    end
end

if (is64 && ~verLessThan ('matlab', '7.8'))
    % versions 7.8 and later on 64-bit platforms use a 64-bit BLAS
%     fprintf ('with 64-bit BLAS\n') ;
    flags = [flags ' -DBLAS64'] ;
end

if (~(pc || mac))
    % for POSIX timing routine
    lapack = [lapack ' -lrt'] ;
end

%% List all the source files
qpalm_src = { ...
    'qpalm',...
    'solver_interface',...
    'lin_alg',...
    'linesearch',...
    'newton',...
    'scaling',...
    'termination',...
    'util',...
    'validate',...
    'nonconvex',...
    'iteration',...
    };
    
qpalm_src_path = '../../src';
for i = 1:length (qpalm_src)
    qpalm_src {i} = [qpalm_src_path '/' qpalm_src{i}] ;
end

if strcmp(solver, 'ladel')
    
    ladel_src = { ...
        'ladel_col_counts', ...
        'ladel_copy', ...
        'ladel_debug_print', ...
        'ladel_etree', ...
        'ladel_global', ...
        'ladel', ...
        'ladel_ldl_numeric', ...
        'ladel_ldl_symbolic', ...
        'ladel_matvec', ...
        'ladel_matmat', ...
        'ladel_submatrix', ...
        'ladel_add', ...
        'ladel_pattern', ...
        'ladel_permutation', ...
        'ladel_postorder', ...
        'ladel_rank1_mod', ...
        'ladel_row_mod', ...
        'ladel_scale', ...
        'ladel_transpose', ...
        'ladel_upper_diag', ...
        };

    ladel_src_path = '../../LADEL/src';
    for i = 1:length(ladel_src)
        ladel_src{i} = [ladel_src_path '/' ladel_src{i}] ;
    end

    ladel_mex_util_src = ['../../LADEL/interfaces/mex/ladel_mex_util']; 

    amd_src = { ...
        'amd_1', ...
        'amd_2', ...
        'amd_aat', ...
        'amd_control', ...
        'amd_defaults', ...
        'amd_dump', ...
        'amd_global', ...
        'amd_info', ...
        'amd_order', ...
        'amd_postorder', ...
        'amd_post_tree', ...
        'amd_preprocess', ...
        'amd_valid', ...
        'SuiteSparse_config', ...
        };

    amd_src_path = '../../LADEL/amd/Source';
    for i = 1:length(amd_src)
        amd_src{i} = [amd_src_path '/' amd_src{i}];
    end
    
    source = [amd_src ladel_src ladel_mex_util_src qpalm_src] ;
    
elseif strcmp(solver, 'cholmod')
    config_src = { '../../suitesparse/SuiteSparse_config/SuiteSparse_config' } ;

    ordering_src = { ...
        '../../suitesparse/AMD/Source/amd_1', ...
        '../../suitesparse/AMD/Source/amd_2', ...
        '../../suitesparse/AMD/Source/amd_aat', ...
        '../../suitesparse/AMD/Source/amd_control', ...
        '../../suitesparse/AMD/Source/amd_defaults', ...
        '../../suitesparse/AMD/Source/amd_dump', ...
        '../../suitesparse/AMD/Source/amd_global', ...
        '../../suitesparse/AMD/Source/amd_info', ...
        '../../suitesparse/AMD/Source/amd_order', ...
        '../../suitesparse/AMD/Source/amd_postorder', ...
        '../../suitesparse/AMD/Source/amd_post_tree', ...
        '../../suitesparse/AMD/Source/amd_preprocess', ...
        '../../suitesparse/AMD/Source/amd_valid', ...
        '../../suitesparse/CAMD/Source/camd_1', ...
        '../../suitesparse/CAMD/Source/camd_2', ...
        '../../suitesparse/CAMD/Source/camd_aat', ...
        '../../suitesparse/CAMD/Source/camd_control', ...
        '../../suitesparse/CAMD/Source/camd_defaults', ...
        '../../suitesparse/CAMD/Source/camd_dump', ...
        '../../suitesparse/CAMD/Source/camd_global', ...
        '../../suitesparse/CAMD/Source/camd_info', ...
        '../../suitesparse/CAMD/Source/camd_order', ...
        '../../suitesparse/CAMD/Source/camd_postorder', ...
        '../../suitesparse/CAMD/Source/camd_preprocess', ...
        '../../suitesparse/CAMD/Source/camd_valid', ...
        '../../suitesparse/COLAMD/Source/colamd', ...
        '../../suitesparse/CCOLAMD/Source/ccolamd' } ;

    if (have_metis)

        metis_src = {
            'GKlib/b64', ...
            'GKlib/blas', ...
            'GKlib/csr', ...
            'GKlib/error', ...
            'GKlib/evaluate', ...
            'GKlib/fkvkselect', ...
            'GKlib/fs', ...
            'GKlib/getopt', ...
            'GKlib/gkregex', ...
            'GKlib/graph', ...
            'GKlib/htable', ...
            'GKlib/io', ...
            'GKlib/itemsets', ...
            'GKlib/mcore', ...
            'GKlib/memory', ...
            'GKlib/omp', ...
            'GKlib/pdb', ...
            'GKlib/pqueue', ...
            'GKlib/random', ...
            'GKlib/rw', ...
            'GKlib/seq', ...
            'GKlib/sort', ...
            'GKlib/string', ...
            'GKlib/timers', ...
            'GKlib/tokenizer', ...
            'GKlib/util', ...
            'libmetis/auxapi', ...
            'libmetis/balance', ...
            'libmetis/bucketsort', ...
            'libmetis/checkgraph', ...
            'libmetis/coarsen', ...
            'libmetis/compress', ...
            'libmetis/contig', ...
            'libmetis/debug', ...
            'libmetis/fm', ...
            'libmetis/fortran', ...
            'libmetis/frename', ...
            'libmetis/gklib', ...
            'libmetis/graph', ...
            'libmetis/initpart', ...
            'libmetis/kmetis', ...
            'libmetis/kwayfm', ...
            'libmetis/kwayrefine', ...
            'libmetis/mcutil', ...
            'libmetis/mesh', ...
            'libmetis/meshpart', ...
            'libmetis/minconn', ...
            'libmetis/mincover', ...
            'libmetis/mmd', ...
            'libmetis/ometis', ...
            'libmetis/options', ...
            'libmetis/parmetis', ...
            'libmetis/pmetis', ...
            'libmetis/refine', ...
            'libmetis/separator', ...
            'libmetis/sfm', ...
            'libmetis/srefine', ...
            'libmetis/stat', ...
            'libmetis/timing', ...
            'libmetis/util', ...
            'libmetis/wspace', ...
        } ;

        for i = 1:length (metis_src)
            metis_src {i} = [metis_path '/' metis_src{i}] ;
        end
    end

    cholmod_matlab = { '../../suitesparse/CHOLMOD/MATLAB/cholmod_matlab' } ;

    cholmod_src = {
        '../../suitesparse/CHOLMOD/Core/cholmod_aat', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_add', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_band', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_change_factor', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_common', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_complex', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_copy', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_dense', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_error', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_factor', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_memory', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_sparse', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_transpose', ...
        '../../suitesparse/CHOLMOD/Core/cholmod_triplet', ...
        '../../suitesparse/CHOLMOD/Check/cholmod_check', ...
        '../../suitesparse/CHOLMOD/Check/cholmod_read', ...
        '../../suitesparse/CHOLMOD/Check/cholmod_write', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_amd', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_analyze', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_colamd', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_etree', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_factorize', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_postorder', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_rcond', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_resymbol', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_rowcolcounts', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_rowfac', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_solve', ...
        '../../suitesparse/CHOLMOD/Cholesky/cholmod_spsolve', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_drop', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_horzcat', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_norm', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_scale', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_sdmult', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_ssmult', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_submatrix', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_vertcat', ...
        '../../suitesparse/CHOLMOD/MatrixOps/cholmod_symmetry', ...
        '../../suitesparse/CHOLMOD/Modify/cholmod_rowadd', ...
        '../../suitesparse/CHOLMOD/Modify/cholmod_rowdel', ...
        '../../suitesparse/CHOLMOD/Modify/cholmod_updown', ...
        '../../suitesparse/CHOLMOD/Supernodal/cholmod_super_numeric', ...
        '../../suitesparse/CHOLMOD/Supernodal/cholmod_super_solve', ...
        '../../suitesparse/CHOLMOD/Supernodal/cholmod_super_symbolic', ...
        '../../suitesparse/CHOLMOD/Partition/cholmod_ccolamd', ...
        '../../suitesparse/CHOLMOD/Partition/cholmod_csymamd', ...
        '../../suitesparse/CHOLMOD/Partition/cholmod_camd', ...
        '../../suitesparse/CHOLMOD/Partition/cholmod_metis', ...
        '../../suitesparse/CHOLMOD/Partition/cholmod_nesdis' } ;
    
    if (pc)
        % Windows does not have drand48 and srand48, required by METIS.  Use
        % drand48 and srand48 in CHOLMOD/MATLAB/Windows/rand48.c instead.
        % Also provide Windows with an empty <strings.h> include file.
        cholmod_matlab = [cholmod_matlab {'Windows/rand48'}] ;
        include = [include ' -IWindows'] ;
    end
    
    source = [ordering_src config_src cholmod_src cholmod_matlab qpalm_src] ;
    if (have_metis)
        source = [metis_src source] ;
    end
end

%% Compile everything
if (pc)
    obj_extension = '.obj' ;
else
    obj_extension = '.o' ;
end
obj = '' ;

kk = 0 ;
cflags = 'CFLAGS="\$CFLAGS -std=c99 -fPIC -DMATLAB -DDAMD -O3 -DPROFILING -DPRINTING"';
flags = [cflags ' ' flags];

for f = source
    ff = f {1} ;
    if (strcmp(solver, 'cholmod') && isequal (ff, [metis_path '/GKlib/util']))
        % special case, since a file with the same name also exists in libmetis
        copyfile ([ff '.c'], 'GKlib_util.c', 'f') ;
        ff = 'GKlib_util' ;
        o = 'GKlib_util' ;
    elseif (strcmp(solver, 'cholmod') && isequal (ff, [metis_path '/GKlib/graph']))
        % special case, since a file with the same name also exist in libmetis
        copyfile ([ff '.c'], 'GKlib_graph.c', 'f') ;
        ff = 'GKlib_graph' ;
        o = 'GKlib_graph' ;
    elseif (isequal (ff, [qpalm_src_path '/util']))
        % special case, since a file with the same name also exists in
        % qpalm
        copyfile ([ff '.c'], 'qpalm_util.c', 'f') ;
        ff = 'qpalm_util' ;
        o = 'qpalm_util' ;
    else
        slash = strfind (ff, '/') ;
        if (isempty (slash))
            slash = 1 ;
        else
            slash = slash (end) + 1 ;
        end
        o = ff (slash:end) ;
    end
    % fprintf ('%s\n', o) ;
    o = [o obj_extension] ;
    obj = [obj  ' ' o] ;					            
    s = sprintf ('mex %s -DDLONG -O %s -c %s.c', flags, include, ff) ;
    s = [s obj ' ' lapack] ;
    kk = do_cmd (s, kk, details) ;
end

s = sprintf ('mex %s -DDLONG -O %s %s.c', flags, include, 'qpalm_mex') ;
s = [s obj ' ' lapack] ;	
kk = do_cmd (s, kk, details) ;

 % clean up
s = ['delete ' obj] ;
do_cmd (s, kk, details) ;
fprintf ('\nQPALM successfully compiled\n') ;

% remove the renamed METIS files, if they exist
if (exist ('GKlib_util.c', 'file'))
    delete ('GKlib_util.c') ;
end
if (exist ('GKlib_graph.c', 'file'))
    delete ('GKlib_graph.c') ;
end
if (exist ('qpalm_util.c', 'file'))
    delete ('qpalm_util.c') ;
end

%change to qpalm main path
% cd(current_path);
qpalm_path = fullfile(this_path, '..');
cd(qpalm_path);

function kk = do_cmd (s, kk, details)
 %DO_CMD: evaluate a command, and either print it or print a "."
if (details)
    fprintf ('%s\n', s) ;
else
    if (mod (kk, 60) == 0)
	fprintf ('\n') ;
    end
    kk = kk + 1 ;
    fprintf ('.') ;
end
eval (s) ;
