#!/bin/bash
#
#PBS -N jgwUpdate
#PBS -o output
#PBS -e error
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

outdir="~/bipartitemodelsBC/results/finaltrunc"
runlength=150000
thin=150

R --no-save --args ${outdir} ${runlength} ${thin} < "$HOME/bipartitemodelsBC/finalmodels/update-march.R"
