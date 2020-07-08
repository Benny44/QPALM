#!/bin/bash

#data_folder
folder=/home/ben/Documents/Projects/QPALM/interfaces/qps/all_small

#settings_file (file containing qpalm settings)
settings_file=/home/ben/Documents/Projects/QPALM/interfaces/qps/sample_settings.txt

#output folder for logs
output_folder=/home/ben/Documents/Projects/QPALM/logs

#folder with qpalm_qps executable
builddir=/home/ben/Documents/Projects/QPALM/build/debug/bin
cd $builddir

#for f in /home/ben/Documents/Projects/QPALM/simulations/maros_meszaros/qps/qpdata/
#copies=find /home/ben/Documents/Projects/QPALM/qps/all/ -iwholename "*_copy.qps"
#if [ -z "$copies" ]
#then
#else
  for f in $folder/*_copy.qps
  do
    if [ -z "$f" ]
    then
      rm $f
    fi
  done
#fi


for f in `ls $folder | sort -V`
#for f in `ls /home/ben/Documents/Projects/QPALM/qps/all/hs21.qps | sort -V`
#for f in /home/ben/Documents/Projects/QPALM/qps/problem/*.qps
#for f in /home/ben/Documents/Projects/QPALM/simulations/maros_meszaros/qps/qpdata/cute_big_qps/liswet*.qps
do
  if [[ $f != *"_copy.qps" ]]
  then
    #valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./qpalm_qps $f
    #./qpalm_qps $f 
    # { # try
        length=${#f}
        f2=${f:0:length-4} 
        touch $output_folder/$f2.log
        #sudo nice -n -20 ./qpalm_qps /home/ben/Documents/Projects/QPALM/qps/all/$f > $curdir/logs/$f2.log
        #sudo nice -n -20 ./qpalm_qps /home/ben/Documents/Projects/QPALM/qps/all/$f
        ./qpalm_qps $folder/$f $settings_file > $output_folder/$f2.log
        # ./qpalm_qps $folder/$f $settings_file
        # valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --quiet --suppressions=/home/ben/Documents/Projects/QPALM/valgrind/dl_open.supp ./qpalm_qps $folder/$f $settings_file 
        #./qpalm_qps $f
        #save your output
    #} || { # catch
    #    echo "$f2 &  & &  &  &  & f & \\" >> $curdir/build/debug/bin/out.tex
        # save log for exception 
    #}
    #./qpalm_qps /home/ben/Documents/Projects/QPALM/qps/all/$f > logs/$f.log
  fi
done
#./qpalm_qps /home/ben/Documents/Projects/QPALM/simulations/maros_meszaros/qps/hs21.qps


#simulations/maros_meszaros/qps/qpdata/brunel





