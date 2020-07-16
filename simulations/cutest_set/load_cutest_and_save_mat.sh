#!/bin/bash

# curdir=`pwd`
# echo $curdir
data_folder=/home/ben/Downloads/cutest/mastsif/qpsif_crashing
# data_folder=/home/ben/Downloads/cutest/mastsif/qpsif
# data_folder=/home/ben/Downloads/cutest/mastsif/qpsif/test

for f in `ls $data_folder | sort -V`
do
    cutest2matlab $data_folder/$f    
    LD_PRELOAD="/usr/local/MATLAB/R2020a/bin/glnxa64/mkl.so:/usr/local/MATLAB/R2020a/bin/glnxa64/mklcompat.so:${LD_PRELOAD}" ${MYMATLAB}/bin/matlab -nosplash -nojvm -D"valgrind --error-limit=no --tool=memcheck -v --leak-check=yes" -r "run('load_cutest_and_save_mat.m');"
    # LD_PRELOAD="/usr/local/MATLAB/R2020a/bin/glnxa64/mkl.so:/usr/local/MATLAB/R2020a/bin/glnxa64/mklcompat.so:${LD_PRELOAD}" ${MYMATLAB}/bin/matlab -nosplash -nojvm -r "run('load_cutest_and_save_mat.m');"
done
