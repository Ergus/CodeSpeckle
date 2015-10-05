#!/bin/bash

for a in 1 2 4 8; do
    qsub -N plasma_$a -l select=$a:ppn=20 -q regular -v solver=plasma ./job.sh
    qsub -N magma_$a -l select=$a:ppn=20 -q gpu -v solver=magma ./job.sh
done
