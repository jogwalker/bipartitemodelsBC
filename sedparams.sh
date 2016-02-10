#!/bin/bash


filepath="${HOME}/bipartitemodelsBC/finalmodels/jagsNB/beta"
outdir="${HOME}/bipartitemodelsBC/finalmodels/jagsNB/betaFull"

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
	sed -e 's/D_byhost/PD_host/' < ${filepath}/${j} > ${outdir}/${j}
done

