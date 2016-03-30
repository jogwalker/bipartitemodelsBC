#!/bin/bash
#
#PBS -N jgwZINB
#PBS -o output
#PBS -e error
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

outdir="~/bipartitemodelsBC/results/finaltrunc"
runlength=300000
burn=75000
thin=50

R --no-save --args ${outdir} ${runlength} ${burn} ${thin} < "$HOME/bipartitemodelsBC/finalmodels/finaltrunc-march.R"
