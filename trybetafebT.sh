#!/bin/bash
#
#PBS -N jgwZINBfeb
#PBS -o output/beta
#PBS -e error/beta
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 0-16

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/betaT"
outdir="~/bipartitemodelsBC/results/betaT"
iter=105000
burnin=5000
thin=100

ls -1 ${filepath}

i=0
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})

echo ${filepath}/${array[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength[${PBS_ARRAYID}]} ${iter} ${burnin} ${thin} < ${filepath}/${array[${PBS_ARRAYID}]}
