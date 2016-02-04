#!/bin/bash
#
#PBS -N jgwZINB
#PBS -o output
#PBS -e error
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash

outdir="~/bipartitemodelsBC/results/"
runlength=50000

R --no-save --args ${outdir} ${runlength} < "$HOME/bipartitemodelsBC/finalmodels/cnj-mat.R"
