#!/usr/bin/env bash

for vcf in $(ls -1 ${HOME}/scratch/pop-gen-02-results/16-freebayes-theta-*/*filt*.vcf); do
    # echo "${vcf}"
    ./21-pairwise-pi-from-vcf.R "${vcf}" &
done

