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

R --no-save --args ${outdir} ${runlength} < "$HOME/bipartitemodelsBC/finalmodels/finaltrunc.R"
