#!/bin/bash
#
#PBS -N jgwZINBbeta
#PBS -o output/newwild
#PBS -e error/newwild
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 4-6

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/newspec-pete"
outdir="~/bipartitemodelsBC/results/newspec"
runlength=3000

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})
echo ${filepath}/${array[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength} < ${filepath}/${array[${PBS_ARRAYID}]}
