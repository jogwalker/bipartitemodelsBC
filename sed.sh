#!/bin/bash


filepath="${HOME}/Documents/bipartitemodelsBC/finalmodels/beta"
outdir="${HOME}/Documents/bipartitemodelsBC/finalmodels/betaT"

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})


for j in "${array[@]}"
do
	echo $j
	sed -e 's/\/beta\//\/betaT\//' < ${filepath}/${j} > ${outdir}/${j}
done

