#!/bin/bash
#
#PBS -N jgwZINBfeb
#PBS -o output
#PBS -e error
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 0-3

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/alreadyfeb"
outdir="~/bipartitemodelsBC/results/alreadyfeb"
iter=105000

ls -1 ${filepath}

i=0
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})

echo ${filepath}/${array[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength[${PBS_ARRAYID}]} ${iter} < ${filepath}/${array[${PBS_ARRAYID}]}
