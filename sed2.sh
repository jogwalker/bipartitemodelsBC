#!/bin/bash


filepath="${HOME}/Documents/bipartitemodelsBC/finalmodels/jagsNB/beta"
outdir="${HOME}/Documents/bipartitemodelsBC/finalmodels/jagsNB/betaT"

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})


for j in "${array[@]}"
do
	sed -e 's/T(0.00001,0.99999)/T(0.0001,0.9999)/' < ${filepath}/${j} > ${outdir}/${j}
done

