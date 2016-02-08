#!/bin/bash


filepath="${HOME}/Documents/bipartitemodelsBC/finalmodels/beta"
outdir="${HOME}/Documents/bipartitemodelsBC/finalmodels/betaFull"

ls -1 ${filepath}

i=1
while read line
do
    array[ $i ]="$line"        
    (( i++ ))
done < <(ls -1 ${filepath})


for j in "${array[@]}"
do
	sed -e 's/'mu'/'beta\',\'alpha\',\'r\',\'HB_invert\',\'PD_host\',\'hosts'/' < ${filepath}/${j} > ${outdir}/${j}
done

