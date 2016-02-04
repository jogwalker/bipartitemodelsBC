#!/bin/bash
#
#PBS -N jgwZINBtrunc
#PBS -o output/trunc
#PBS -e error/trunc
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 0-11

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/newspec19Jan"
outdir="~/bipartitemodelsBC/results/newspec19Jan"

declare -a runlength=(3000 3000 3000 3000 50000 50000 50000 50000 100000 100000 100000 100000)

ls -1 ${filepath}

i=0
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})

declare -a file=(${array[@]} ${array[@]} ${array[@]})

echo ${filepath}/${file[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength[${PBS_ARRAYID}]} < ${filepath}/${file[${PBS_ARRAYID}]}
