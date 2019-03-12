current_path = fileparts(mfilename('fullpath'));

error_msg = 'The C compiler could not succesfully compile ';
if mex('-outdir', current_path, fullfile(current_path,'PWAlinesearch_mex.c'),...
         '-lm', '-largeArrayDims',...
        'CFLAGS="\$CFLAGS -std=c99 -O3"')
    error([error_msg, mex_path]);
end