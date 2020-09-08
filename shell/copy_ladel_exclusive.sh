#!/bin/bash

# Copy everything to a new directory
if [ -d "../QPALM_ladel" ]; then
    rm -r ../QPALM_ladel
fi

mkdir ../QPALM_ladel

rsync -a . ../QPALM_ladel \
    --exclude .git/ \
    --exclude '*.mat' \
    --exclude '*.qps' \
    --exclude '*.SIF' \
    --exclude '*.oqp' \
    --exclude profiling \
    --exclude 'jit_tmp.c' \
    --exclude build/ \
    --exclude suitesparse/ \
    --exclude cmakeSuitesparse/ \
    --exclude docs/ \

# Remove all the cholmod code in C
cd ../QPALM_ladel
re_ifdef="#ifdef.*"
re_ifndef="#ifndef.*"
re_endif="#endif.*"
deleting=1
for f in `rgrep -l "#elif defined USE_CHOLMOD" --include="*.c" --include="*.h"` 
do
    i=0
    ladel_i=0
    declare -i lineno=1
    while read -r line
    do
        if [ "$line" == "#ifdef USE_LADEL" ]; then
            ladel_i=1
            sed -i "$lineno d" "$f"
        else
            if [ "$line" == "#elif defined USE_CHOLMOD" ]; then
                ((i=i+1))
                ladel_i=0
                echo "Text read at $lineno from $f: $line"
            fi

            if [ "$ladel_i" -ge "$deleting" ]; then
                if [[ $line =~ $re_ifdef || $line =~ $re_ifndef ]]; then
                    ((ladel_i=ladel_i+1))
                fi
            fi 

            if [ "$i" -ge "$deleting" ]; then
                if [[ $line =~ $re_ifdef || $line =~ $re_ifndef ]]; then
                    ((i=i+1))
                fi
                echo "Deleting line $lineno"
                sed -i "$lineno d" "$f"
            else
                let lineno++
            fi

            if [[ $line =~ $re_endif ]]; then
                if [ "$i" -ge "$deleting" ]; then 
                    ((i=i-1))
                fi
                if [ "$ladel_i" -ge "$deleting" ]; then 
                    ((ladel_i=ladel_i-1))
                    if [ "$ladel_i" == "0" ]; then
                        let lineno--
                        sed -i "$lineno d" "$f"
                    fi
                fi
            fi
        fi
        
        # sed -i "1 d" "$f"
    done < "$f"
    if [ "$i" -ge "$deleting" ]; then
        sed -i "$lineno d" "$f"
    fi
    # sed -i '/./{H;$!d} ; x ; s/#elif defined USE_CHOLMOD.*#endif//' $f
    # sed -i '/./{H;$!d} ; x ; s/#ifdef USE_LADEL\n#include "ladel.h"\n#endif/#include "ladel.h"/' $f
    # sed -i 's/#ifdef USE_LADEL//' $f
done


# Clean up any externals
