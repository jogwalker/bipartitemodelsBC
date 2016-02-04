#!/bin/bash
#
#PBS -N jgwZINB
#PBS -o output/newwild
#PBS -e error/newwild
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 1-2

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/newprior"
outdir="~/bipartitemodelsBC/results/newprior/long"
runlength=100000

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})
echo ${filepath}/${array[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength} < ${filepath}/${array[${PBS_ARRAYID}]}
