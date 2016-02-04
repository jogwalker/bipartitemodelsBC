#/bin/bash
#
#PBS -N jgwZINB
#PBS -o output/ZINB5
#PBS -e error/ZINB5
#PBS -l walltime=144:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash
#PBS -t 6-8

if [ -n "${1}" ]; then
  echo "${1}"
  PBS_ARRAYID=${1}
fi

filepath="$HOME/bipartitemodelsBC/finalmodels/ZINBphase5"
outdir="~/bipartitemodelsBC/results/ZINB5_27Nov"
runlength=200000

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})
echo ${filepath}/${array[${PBS_ARRAYID}]}

R --no-save --args ${outdir} ${runlength} < ${filepath}/${array[${PBS_ARRAYID}]}
