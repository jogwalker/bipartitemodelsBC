#!/bin/bash
#
#PBS -N jgwZIPphase1
#PBS -o output/ZIPphase1
#PBS -e error/ZIPphase1
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 9

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/ZIPphase1"
outdir="~/bipartitemodelsBC/results/ZIPphase1_19Nov"
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
