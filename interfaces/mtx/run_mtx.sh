#!/bin/bash
cd $PWD
cd ..
make clean
make mtxformat
cd mtx
./qpalm_mtx mtx_files/a.mtx mtx_files/h.mtx mtx_files/g.mtx mtx_files/lba.mtx mtx_files/uba.mtx
