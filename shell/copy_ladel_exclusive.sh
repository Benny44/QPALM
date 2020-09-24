#!/bin/bash

# Copy everything to a new directory
# if [ -d "../QPALM_ladel" ]; then
#     rm -r ../QPALM_ladel
# fi

# mkdir ../QPALM_ladel

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
    --exclude LADEL/ \
    --exclude docs/ \

if [ -d "../QPALM_ladel/docs" ]; then 
    rsync -a docs/Doxyfile ../QPALM_ladel/docs
    rsync -a docs/ref.bib ../QPALM_ladel/docs
else
    mkdir ../QPALM_ladel/docs
    rsync -a docs/ref.bib ../QPALM_ladel/docs
fi

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
                # echo "Text read at $lineno from $f: $line"
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
                # echo "Deleting line $lineno"
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

sed -i "s/%   qpalm_make('cholmod')//" interfaces/mex/qpalm_make.m
sed -i "s/    fprintf(\'No linear systems solver selected. Using LADEL as the default.\n\');//" interfaces/mex/qpalm_make.m
sed -i "s/    fprintf('Use qpalm_make(''cholmod'') if you wish to compile with CHOLMOD instead.\n');//" interfaces/mex/qpalm_make.m
sed -i "s/    fprintf('Using CHOLMOD as the linear systems solver.\n');/error('This version of QPALM only works with LADEL.\n');/" interfaces/mex/qpalm_make.m

## Clean up any externals

#License
cp private/LICENSE ../QPALM_ladel/

#Matlab
sed -i 's/No linear systems solver selected. Using LADEL as the default.\\n//' interfaces/mex/qpalm_make.m
sed -i 's/if you wish to compile with CHOLMOD instead.\\n//' interfaces/mex/qpalm_make.m

sed -i "s/    fprintf('');//" interfaces/mex/qpalm_make.m
sed -i "s/    fprintf('Use qpalm_make(''cholmod'') ');//" interfaces/mex/qpalm_make.m

sed -i "s/fprintf('Using CHOLMOD as the linear systems solver./error('CHOLMOD is not available in this version of QPALM/" interfaces/mex/qpalm_make.m
sed -i "s/or ''cholmod''.//" interfaces/mex/qpalm_make.m

ifs=0
declare -i lineno=1
f="interfaces/mex/qpalm_make.m"
re_if="if.*"
re_for="for i.*"
re_modify="Modify.*"
re_elseif="elseif.*"
re_ifmetis="if METIS.*"
re_end="end.*"
while read -r line
do
    if [ "$line" == "elseif strcmp(solver, 'cholmod')" ]; then
        ((ifs=ifs+1))
        sed -i "$lineno d" "$f"
    else
        if [ "$ifs" -ge "$deleting" ]; then
        # echo "Text read at $lineno from $f ($ifs ifs): $line"
            if [[ $line =~ $re_if || $line =~ $re_for ]]; then
                if ! [[ $line =~ $re_elseif || $line =~ $re_ifmetis || $line =~ $re_modify ]]; then
                    ((ifs=ifs+1))
                fi
            fi
            if [[ $line =~ $re_end ]]; then
                ((ifs=ifs-1))
            fi
            if [ "$ifs" -ge "$deleting" ]; then
                sed -i "$lineno d" "$f"
            else
                let lineno++
            fi
        else
            let lineno++
        fi
    fi
done < "$f"

#Readme
sed -i 's/## Check out QPALM_vLADEL//' README.md
sed -i 's#Check out \[this\](https://github.com/Benny44/QPALM_vLADEL) for the main LGPL-v3 licensed version of QPALM based on LADEL. This repo is only maintained because it provides an interface also to CHOLMOD, which might be more useful than LADEL for dense problems.##' README.md 

sed -i 's#* Suitesparse: authored by Tim Davis. Each of its modules is licensed separately, see \[suitesparse/LICENSE.txt\](https://github.com/jluttine/suitesparse/blob/e409f9fb39181ea86718dbf91ce39c2c7e6c3dcd/LICENSE.txt). The main module used in QPALM is CHOLMOD.#\* LADEL: authored by Ben Hermans and licensed under \[LGPL-v3\](https://github.com/Benny44/LADEL/blob/master/LICENSE).#' README.md
sed -i 's#* Intel MKL: authored by the Intel Corporation and licensed under the Intel Simplified Software License.##' README.md

# sed -i 's#QPALM.svg#QPALM_vLADEL.svg#' README.md
sed -i 's#Benny44/QPALM#Benny44/QPALM_vLADEL#' README.md
sed -i 's#io/QPALM#io/QPALM_vLADEL#' README.md

#Documentation
sed -i 's#github.com/Benny44/QPALM#github.com/Benny44/QPALM_vLADEL#' doxypages/mainpage.dox
sed -i 's# * GPL 3.0# * LGPL-3.0#' doxypages/mainpage.dox
sed -i 's#ladel/cholmod#ladel#' include/solver_interface.h
sed -i 's#cholmod/ladel#ladel#' include/solver_interface.h
sed -i 's#Finally, all the settings relevant to cholmod (and suitesparse) are included##' include/solver_interface.h
sed -i 's#in this file as well.##'  include/solver_interface.h
sed -i 's#ladel/cholmod#ladel#' src/solver_interface.c
sed -i 's#cholmod/ladel#ladel#' src/solver_interface.c
sed -i 's#Finally, all the settings relevant to cholmod (and suitesparse) are included##' src/solver_interface.c
sed -i 's#in this file as well.##'  src/solver_interface.c
sed -i 's#@note The function in this file makes use of the cholmod scale routines.##' include/scaling.h
sed -i 's#@note This function makes use of the cholmod scale routines.##' include/scaling.h
sed -i 's#@note The function in this file makes use of the cholmod scale routines.##' src/scaling.c
sed -i 's#cholmod_factor#factor#' src/iteration.c
sed -i 's#cholmod_factor#factor#' include/iteration.h


cd docs
doxygen Doxyfile
cd ..