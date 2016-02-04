#!/bin/bash
#
#PBS -N jgwZINBconvert
#PBS -o output
#PBS -e error
#PBS -l walltime=12:00:00
#PBS -l nodes=1:ppn=1
#PBS -S /bin/bash


filepath="$HOME/bipartitemodelsBC/finalmodels/ZIPphase4forNB"
outdir="$HOME/bipartitemodelsBC/finalmodels/ZINBphase4"

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})


for j in "${array[@]}"
do
	sed -e 's/file=(\(.*\)))/endhere\n\n\1/' -e '/cat(/,/endhere/D' -e 's_/jags/_/jagsNB/_' < ${filepath}/${j} > ${outdir}/${j}
done

