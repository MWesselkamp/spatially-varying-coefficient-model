#!/bin/bash

cd ./svcm

input=("1" "2" "3" "4")

for i in ${input[*]}
do
    echo "submitting fgamSVC${i}.sh"
    msub compute_svcm_hpc${i}.sh
done
