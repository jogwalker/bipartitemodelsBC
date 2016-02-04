#!/bin/bash
#
#PBS -N jgwZINBup
#PBS -o output/newwild
#PBS -e error/newwild
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 0-3

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

declare -a files=("newwild/cnj_output_mat.RData" "newwild/cnj_output.RData" "newprior/cnj_output.RData" "newprior/sd/cnj_output50000.RData")

outdir="~/bipartitemodelsBC/results"
runlength=1000
reps=10
model=${files[${PBS_ARRAYID}]}


R --no-save --args ${outdir} ${runlength} ${reps} ${model} < "$HOME/bipartitemodelsBC/updatemodel.R"